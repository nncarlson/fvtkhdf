"""Check HTG block metadata, shared time stepping, and VTK composite readback."""
import sys
from pathlib import Path

try:
    import h5py
    import numpy as np
except ImportError as exc:
    print(exc)
    raise SystemExit(77)

operation, directory = sys.argv[1:3]
prefix = sys.argv[3] if len(sys.argv) > 3 else "mb"
path = Path(directory) / f"{prefix}_hypertree.vtkhdf"
static_path = Path(directory) / f"{prefix}_hypertree_static.vtkhdf"
composite_type = b"PartitionedDataSetCollection" if prefix == "pdc" else b"MultiBlockDataSet"
interface_names = ("vtkInterfaceNormals", "vtkInterfaceIntercepts")
names = ("fixed", "moving", "dynamic")


def tree_ids(nproc, phase):
    active = max(nproc - 1, 1) if phase else nproc
    return [[g for g in range(2) if (g + phase) % active == rank] for rank in range(nproc)]


def cell_values(mode, step, ids):
    phase = int(mode == 3 and step > 0)
    return [1000 * mode + 100 * g + 10 * min(step, 1) + j
            for g in ids for j in range(3 if phase and g == 0 else 1)]


def interface_values(values):
    values = np.asarray(values, dtype=np.float64)[:, None]
    return (values + .125*np.arange(1, 4), -values - .25*np.arange(1, 4))


def check_interface_layout(block, values, offset=0):
    for name, expected in zip(interface_names, interface_values(values)):
        assert block.attrs[name[3:]+"Name"] == name.encode()
        dataset = block[f"CellData/{name}"]
        assert dataset.attrs["Name"] == name.encode()
        assert dataset.dtype == np.dtype("f8")
        assert dataset.shape == block["CellData/value"].shape + (3,)
        np.testing.assert_array_equal(dataset[offset:offset+len(values)], expected)
        if "Steps" in block:
            np.testing.assert_array_equal(block[f"Steps/CellDataOffsets/{name}"], [0, 2, 2])


def check_parent(root):
    assert root.attrs["Type"] == composite_type
    assert "Steps" not in root
    for name in ("InterfaceNormalsName", "InterfaceInterceptsName"):
        assert name not in root.attrs  # the attributes belong to each HTG block


with h5py.File(path) as file:
    root = file["VTKHDF"]
    check_parent(root)
    assert set(root["Assembly"]) == set(names)
    nproc = int(root["Group_0/Steps/NumberOfParts"][0])
    for mode, name in enumerate(names, 1):
        block = root[f"Group_{mode - 1}"]
        assert block.attrs["Type"] == b"HyperTreeGrid"
        assert block.attrs["Name"] == name.encode()
        assert int(np.asarray(block.attrs["Index"]).flat[0]) == mode - 1
        link = root["Assembly"].get(name, getlink=True)
        assert isinstance(link, h5py.SoftLink) and link.path == f"/VTKHDF/Group_{mode - 1}"
        assert "PointData" not in block
        np.testing.assert_array_equal(block.attrs["Dimensions"], [3, 1, 1])
        np.testing.assert_array_equal(block["FieldData/id"], [[mode]])
        steps = block["Steps"]
        assert int(steps.attrs["NSteps"]) == 3  # includes the step finalized by the owner
        np.testing.assert_array_equal(steps["Values"], [0, 1, 2])
        np.testing.assert_array_equal(steps["NumberOfParts"], [nproc] * 3)
        np.testing.assert_array_equal(steps["PartOffsets"], [0, nproc, nproc] if mode == 3 else [0] * 3)
        np.testing.assert_array_equal(steps["TreeIdsOffsets"], [0, 2, 2] if mode == 3 else [0] * 3)
        np.testing.assert_array_equal(steps["CellDataOffsets/value"], [0, 2, 2])
        np.testing.assert_array_equal(steps["FieldDataOffsets/id"], [0, 0, 0])
        np.testing.assert_array_equal(steps["FieldDataSizes/id"], [[1, 1]] * 3)
        np.testing.assert_array_equal(steps["XCoordinatesOffsets"], [0, 3*nproc, 3*nproc] if mode > 1 else [0] * 3)
        for dataset in ("XCoordinates", "TreeIds", "CellData/value"):
            assert block[dataset].chunks == (8,), (name, dataset, block[dataset].chunks)
        assert steps["CellDataOffsets/value"].chunks == (512,)
        for step in range(3):
            phase = int(mode == 3 and step > 0)
            ids = tree_ids(nproc, phase)
            part0 = int(steps["PartOffsets"][step])
            np.testing.assert_array_equal(block["NumberOfTrees"][part0:part0+nproc], list(map(len, ids)))
            counts = [sum(3 if phase and g == 0 else 1 for g in part) for part in ids]
            np.testing.assert_array_equal(block["NumberOfCells"][part0:part0+nproc], counts)
            tree0 = int(steps["TreeIdsOffsets"][step])
            np.testing.assert_array_equal(block["TreeIds"][tree0:tree0+2], [g for part in ids for g in part])
            offset = int(steps["CellDataOffsets/value"][step])
            np.testing.assert_array_equal(block["CellData/value"][offset:offset+sum(counts)],
                                          cell_values(mode, step, [g for part in ids for g in part]))
            check_interface_layout(block, cell_values(mode, step, [g for part in ids for g in part]), offset)
            offset = int(steps["XCoordinatesOffsets"][step])
            shift = int(mode > 1 and step > 0)
            np.testing.assert_array_equal(block["XCoordinates"][offset:offset+3*nproc],
                                          np.tile(np.arange(3) + shift, nproc))

with h5py.File(static_path) as file:
    root = file["VTKHDF"]
    check_parent(root)
    assert list(root["Assembly"]) == ["static"]
    block = root["Group_0"]
    assert block.attrs["Type"] == b"HyperTreeGrid"
    assert block.attrs["Name"] == b"static"
    assert "Steps" not in block
    assert len(block["NumberOfTrees"]) == nproc
    ids = [g for part in tree_ids(nproc, 0) for g in part]
    np.testing.assert_array_equal(block["CellData/value"], cell_values(0, 0, ids))
    check_interface_layout(block, cell_values(0, 0, ids))

if operation == "reader":
    try:
        from vtkmodules.vtkIOHDF import vtkHDFReader
        from vtkmodules.vtkCommonCore import vtkVersion
        from vtkmodules.vtkCommonDataModel import vtkHyperTreeGridNonOrientedGeometryCursor
        from vtkmodules.util.numpy_support import vtk_to_numpy
    except ImportError as exc:
        print(exc)
        raise SystemExit(77)
    if (vtkVersion.GetVTKMajorVersion(), vtkVersion.GetVTKMinorVersion()) < (9, 7):
        print("HyperTreeGrid reader checks require VTK 9.7 or later")
        raise SystemExit(77)
    def composite_blocks(output, count):
        if prefix == "pdc":
            assert output.IsA("vtkPartitionedDataSetCollection")
            assert output.GetNumberOfPartitionedDataSets() == count
            return [output.GetPartitionedDataSet(i) for i in range(count)]
        assert output.IsA("vtkMultiBlockDataSet") and output.GetNumberOfBlocks() == count
        return [output.GetBlock(i) for i in range(count)]

    def pieces(block):
        return ([block] if block.IsA("vtkHyperTreeGrid") else
                [block.GetPartitionAsDataObject(i) for i in range(block.GetNumberOfPartitions())])

    def check_interface_reader(part, values, indices):
        assert part.GetHasInterface()
        assert part.GetInterfaceNormalsName() == interface_names[0]
        assert part.GetInterfaceInterceptsName() == interface_names[1]
        for name, expected in zip(interface_names, interface_values(values)):
            array = part.GetCellData().GetArray(name)
            assert array is not None and array.GetNumberOfComponents() == 3
            np.testing.assert_array_equal(vtk_to_numpy(array)[np.asarray(indices, dtype=int)], expected)

    reader = vtkHDFReader()
    errors = []
    reader.AddObserver("ErrorEvent", lambda *args: errors.append(args))
    reader.SetFileName(str(path))
    reader.UpdateInformation()
    assert reader.GetNumberOfSteps() == 3
    for step in (0, 1, 2, 1, 0):
        reader.SetStep(step)
        reader.Update()
        assert not errors, errors
        output = reader.GetOutputDataObject(0)
        for mode, block in enumerate(composite_blocks(output, 3), 1):
            name = names[mode - 1]
            parts = pieces(block)
            assert len(parts) == nproc
            phase = int(mode == 3 and step > 0)
            shift = int(mode > 1 and step > 0)
            for part, ids in zip(parts, tree_ids(nproc, phase)):
                assert part.IsA("vtkHyperTreeGrid")
                np.testing.assert_array_equal(vtk_to_numpy(part.GetXCoordinates()), np.arange(3) + shift)
                data = part.GetCellData().GetArray("value")
                assert data is not None
                cursor = vtkHyperTreeGridNonOrientedGeometryCursor()
                visited = []
                indices = []
                for g in ids:
                    part.InitializeNonOrientedGeometryCursor(cursor, g, False)
                    assert cursor.HasTree()
                    refined = bool(phase and g == 0)
                    assert cursor.IsLeaf() != refined
                    assert not cursor.IsMasked()
                    bounds = [0.] * 6
                    cursor.GetBounds(bounds)
                    np.testing.assert_allclose(bounds, [g+shift, g+shift+1, 0, 0, 0, 0])
                    indices.append(cursor.GetGlobalNodeIndex())
                    visited.append(data.GetTuple1(cursor.GetGlobalNodeIndex()))
                    if refined:
                        assert cursor.GetNumberOfChildren() == 2
                        for child in range(2):
                            cursor.ToChild(child)
                            assert cursor.IsLeaf() and bool(cursor.IsMasked()) == (child == 1)
                            cursor.GetBounds(bounds)
                            np.testing.assert_allclose(bounds, [shift+.5*child, shift+.5*(child+1), 0, 0, 0, 0])
                            indices.append(cursor.GetGlobalNodeIndex())
                            visited.append(data.GetTuple1(cursor.GetGlobalNodeIndex()))
                            cursor.ToParent()
                assert part.GetNumberOfCells() == len(visited)
                np.testing.assert_array_equal(visited, cell_values(mode, step, ids), err_msg=f"{name}, step {step}")
                check_interface_reader(part, cell_values(mode, step, ids), indices)

    reader = vtkHDFReader()
    reader.AddObserver("ErrorEvent", lambda *args: errors.append(args))
    reader.SetFileName(str(static_path))
    reader.Update()
    assert not errors, errors
    parts = pieces(composite_blocks(reader.GetOutputDataObject(0), 1)[0])
    assert len(parts) == nproc
    for part, ids in zip(parts, tree_ids(nproc, 0)):
        assert part.GetNumberOfCells() == len(ids)
        cursor = vtkHyperTreeGridNonOrientedGeometryCursor()
        indices = []
        for g in ids:
            part.InitializeNonOrientedGeometryCursor(cursor, g, False)
            assert cursor.HasTree() and cursor.IsLeaf()
            indices.append(cursor.GetGlobalNodeIndex())
        check_interface_reader(part, cell_values(0, 0, ids), indices)
elif operation != "layout":
    raise ValueError(operation)

print(f"{prefix.upper()} HyperTreeGrid {operation} checks passed")
