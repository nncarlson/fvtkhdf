"""Independent checks of HTG numeric overloads and writer integration."""
import sys
from pathlib import Path

try:
    import h5py
    import numpy as np
except ImportError as exc:
    print(exc)
    raise SystemExit(77)

operation, directory = sys.argv[1:]
directory = Path(directory)
if operation == "reader":
    try:
        from vtkmodules.vtkIOHDF import vtkHDFReader
        from vtkmodules.vtkCommonCore import vtkVersion
        from vtkmodules.vtkCommonDataModel import vtkHyperTreeGrid, vtkHyperTreeGridNonOrientedGeometryCursor
        from vtkmodules.util.numpy_support import vtk_to_numpy
    except ImportError as exc:
        print(exc)
        raise SystemExit(77)
    if (vtkVersion.GetVTKMajorVersion(), vtkVersion.GetVTKMinorVersion()) < (9, 7):
        raise SystemExit(77)

KINDS = dict(int8=np.dtype("u1"), int32=np.dtype("i4"), int64=np.dtype("i8"),
             real32=np.dtype("f4"), real64=np.dtype("f8"))


def ids_by_part(nproc, phase=0):
    return [np.array([i - 1 for i in range(1, 4) if (i + phase) % nproc == rank], dtype=int)
            for rank in range(nproc)]


def numeric_values(kind, components, ids, phase):
    scalar = 10 * phase + 3 * ids + (.25 if kind.startswith("real") else 0)
    values = scalar if components == 1 else np.stack((scalar + 1, scalar + 2), axis=-1)
    return np.asarray(values, dtype=KINDS[kind])


def reader_for(path):
    reader = vtkHDFReader()
    errors = []
    reader.AddObserver("ErrorEvent", lambda *args: errors.append(args))
    reader.SetFileName(str(path))
    reader.UpdateInformation()
    return reader, errors


def pieces(output):
    if isinstance(output, vtkHyperTreeGrid):
        return [output]
    return [output.GetPartitionAsDataObject(i) for i in range(output.GetNumberOfPartitions())]


def cell_indices(htg, ids):
    cursor = vtkHyperTreeGridNonOrientedGeometryCursor()
    indices = []
    for tree in ids:
        htg.InitializeNonOrientedGeometryCursor(cursor, int(tree), False)
        assert cursor.HasTree() and cursor.IsLeaf()
        indices.append(cursor.GetGlobalNodeIndex())
    assert htg.GetNumberOfCells() == len(ids)
    return np.array(indices, dtype=int)


def check_numeric(mode, suffix, chunk_bytes):
    path = directory / f"htg_numeric_{suffix}.vtkhdf"
    with h5py.File(path) as file:
        root = file["VTKHDF"]
        nproc = int(root["Steps/NumberOfParts"][0]) if mode else len(root["NumberOfTrees"])
        phases = [0, 0, 2, 2] if mode else [0]
        names = {}
        for kind in KINDS:
            for components, shape in ((1, "scalar"), (2, "vector")):
                if mode != 3:
                    names[f"static_{kind}_{shape}"] = (kind, components, False)
                if mode:
                    names[f"{kind}_{shape}"] = (kind, components, True)
        assert set(root["CellData"]) == set(names)
        for name, (kind, components, temporal) in names.items():
            dataset = root[f"CellData/{name}"]
            assert dataset.dtype == KINDS[kind], (path, name, dataset.dtype)
            assert dataset.attrs["Name"] == name.encode()
            assert dataset.shape == ((6 if temporal else 3),) + (() if components == 1 else (2,))
            if temporal:
                count = max(1, chunk_bytes // (KINDS[kind].itemsize * components))
                assert dataset.chunks == (count,) + (() if components == 1 else (2,)), (path, name, dataset.chunks)
                np.testing.assert_array_equal(root[f"Steps/CellDataOffsets/{name}"], [0, 0, 3, 3])
            else:
                assert dataset.chunks is None
            for step, phase in enumerate(phases):
                ownership = 1 if mode == 3 and step >= 2 else 0
                ids = np.concatenate(ids_by_part(nproc, ownership if temporal else 0))
                offset = (3 if step >= 2 else 0) if temporal else 0
                np.testing.assert_array_equal(dataset[offset:offset+3],
                                              numeric_values(kind, components, ids, phase if temporal else 0),
                                              err_msg=f"{path}:{name} step {step}")
        if mode:
            np.testing.assert_array_equal(root["Steps/Values"], [0, 1, 2, 3])
            assert int(root["Steps"].attrs["NSteps"]) == 4
            for name, value_bytes in (("XCoordinates", 8), ("TreeIds", 8), ("Descriptors", 1)):
                assert root[name].chunks == (max(1, chunk_bytes // value_bytes),), (path, name)
            assert root["Steps/CellDataOffsets/real64_scalar"].chunks == (512,)
            assert root["Steps/PartOffsets"].chunks == (512,)
        for step in range(len(phases)):
            ownership = 1 if mode == 3 and step >= 2 else 0
            want = ids_by_part(nproc, ownership)
            start = 3 if mode == 3 and step >= 2 else 0
            np.testing.assert_array_equal(root["TreeIds"][start:start+3], np.concatenate(want))

    if operation == "reader":
        reader, errors = reader_for(path)
        assert reader.GetNumberOfSteps() == len(phases)
        for step in list(range(len(phases))) + list(range(len(phases) - 1, -1, -1)):
            reader.SetStep(step)
            reader.Update()
            assert not errors, (path, step, errors)
            output = pieces(reader.GetOutputDataObject(0))
            ownership = 1 if mode == 3 and step >= 2 else 0
            assert len(output) == nproc
            for htg, ids in zip(output, ids_by_part(nproc, ownership)):
                indices = cell_indices(htg, ids)
                shift = int(mode >= 2 and step >= 2)
                np.testing.assert_array_equal(vtk_to_numpy(htg.GetXCoordinates()), np.arange(4) + shift)
                cell_data = htg.GetCellData()
                assert {cell_data.GetArrayName(i) for i in range(cell_data.GetNumberOfArrays())} == set(names)
                for name, (kind, components, temporal) in names.items():
                    array = cell_data.GetArray(name)
                    assert array.GetNumberOfComponents() == components
                    actual = vtk_to_numpy(array)
                    assert actual.dtype == KINDS[kind]
                    np.testing.assert_array_equal(actual[indices],
                                                  numeric_values(kind, components, ids, phases[step] if temporal else 0),
                                                  err_msg=f"{path}:{name} step {step}")


def check_names():
    path = directory / "htg_names.vtkhdf"
    # Internal HDF5 name: original label, base value, temporal.
    names = {
        "a_b": ("a/b", 10, False),
        "a_b_1": ("a.b", 20, True),
        "a_b_1_1": ("a_b_1", 30, False),
        "Unnamed_CellData": ("Unnamed_CellData", 40, False),
        "label": ("label", 50, False),
        "a_b_2": ("a/b", 60, True),
        "reverse_one": ("reverse.one", 80, True),
        "reverse_one_1": ("reverse one", 70, False),
    }
    with h5py.File(path) as file:
        root = file["VTKHDF"]
        nproc = int(root["Steps/NumberOfParts"][0])
        ids = np.concatenate(ids_by_part(nproc))
        assert set(root["CellData"]) == set(names)
        for name, (label, base, temporal) in names.items():
            dataset = root[f"CellData/{name}"]
            assert dataset.attrs["Name"] == label.encode()
            expected = np.concatenate((ids + base, ids + base + 100)) if temporal else ids + base
            np.testing.assert_array_equal(dataset, expected)
            if temporal:
                np.testing.assert_array_equal(root[f"Steps/CellDataOffsets/{name}"], [0, 3, 3])
    if operation == "reader":
        reader, errors = reader_for(path)
        for step in (0, 1, 2, 1, 0):
            reader.SetStep(step)
            reader.Update()
            assert not errors, errors
            output = pieces(reader.GetOutputDataObject(0))
            assert len(output) == nproc
            for htg, ids in zip(output, ids_by_part(nproc)):
                indices = cell_indices(htg, ids)
                for name, (_, base, temporal) in names.items():
                    np.testing.assert_array_equal(vtk_to_numpy(htg.GetCellData().GetArray(name))[indices],
                                                  ids + base + (100 * min(step, 1) if temporal else 0))


def check_field_lifecycle():
    path = directory / "htg_field_lifecycle.vtkhdf"
    constants = dict(before=10, during=20, between=30, after=40, last_open=50)
    with h5py.File(path) as file:
        root = file["VTKHDF"]
        assert set(root["FieldData"]) == set(constants) | {"time"}
        np.testing.assert_array_equal(root["Steps/Values"], [0, 1, 2, 3])
        assert int(root["Steps"].attrs["NSteps"]) == 4
        for name, value in constants.items():
            np.testing.assert_array_equal(root[f"FieldData/{name}"], [[value]])
            np.testing.assert_array_equal(root[f"Steps/FieldDataOffsets/{name}"], [0, 0, 0, 0])
            np.testing.assert_array_equal(root[f"Steps/FieldDataSizes/{name}"], [[1, 1]] * 4)
        np.testing.assert_array_equal(root["FieldData/time"], [[11], [33]])
        np.testing.assert_array_equal(root["Steps/FieldDataOffsets/time"], [0, 0, 0, 1])
        np.testing.assert_array_equal(root["Steps/FieldDataSizes/time"], [[1, 0], [1, 1], [1, 1], [1, 1]])
        assert root["FieldData/time"].chunks == (8, 1)  # forwarded 64-byte chunk policy
    if operation == "reader":
        # Check each stored state independently of the installed reader's
        # cache: VTK 9.7.2 can drop constant field arrays on repeated HTG reads.
        for step in range(4):
            reader, errors = reader_for(path)
            reader.SetStep(step)
            reader.Update()
            assert not errors, errors
            fields = reader.GetOutputDataObject(0).GetFieldData()
            for name, value in constants.items():
                array = fields.GetArray(name)
                assert array is not None, (path, step, name)
                np.testing.assert_array_equal(vtk_to_numpy(array), [value])
            expected = [] if step == 0 else [33 if step == 3 else 11]
            np.testing.assert_array_equal(vtk_to_numpy(fields.GetArray("time")), expected)


for mode in range(4):
    check_numeric(mode, str(mode), 64)
check_numeric(1, "default", 128 * 1024)
check_names()
check_field_lifecycle()
print(f"HTG API {operation} checks passed")
