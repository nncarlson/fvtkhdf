"""Check heterogeneous MB collection metadata and VTK block dispatch."""
import sys
from pathlib import Path

try:
    import h5py
    import numpy as np
except ImportError as exc:
    print(exc)
    sys.exit(77)

mode, directory = sys.argv[1:3]
prefix = sys.argv[3] if len(sys.argv) > 3 else "mb"
composite_type = b"PartitionedDataSetCollection" if prefix == "pdc" else b"MultiBlockDataSet"
types = ("ImageData", "RectilinearGrid", "StructuredGrid", "UnstructuredGrid", "HyperTreeGrid")
names = ("Repeated_Block", "Repeated_Block_1", "Repeated_Block_2", "Repeated_Block_3", "Repeated_Block_4")

if mode == "layout":
    for case in ("static", "temporal", "mixed"):
        with h5py.File(Path(directory) / f"{prefix}_collection_{case}.vtkhdf") as file:
            root = file["VTKHDF"]
            assert root.attrs["Type"] == composite_type
            np.testing.assert_array_equal(root.attrs["Version"], [2, 8])
            assert list(root["Assembly"]) == list(names)
            assert list(root) == ["Assembly"] + [f"Group_{i}" for i in range(5)]
            assert set(root["Assembly"]) == set(names)
            for i, (kind, name) in enumerate(zip(types, names)):
                block = root[f"Group_{i}"]
                assert block.attrs["Type"] == kind.encode()
                assert int(np.asarray(block.attrs["Index"]).flat[0]) == i
                assert block.attrs["Name"] == b"Repeated Block"
                link = root["Assembly"].get(name, getlink=True)
                assert isinstance(link, h5py.SoftLink)
                assert link.path == f"/VTKHDF/Group_{i}"
                np.testing.assert_array_equal(block["FieldData/id"][...].flat[:], [i + 1])
                temporal = case == "temporal" or (case == "mixed" and i != 0)
                assert ("Steps" in block) == temporal
                if temporal:
                    steps = block["Steps"]
                    assert int(np.asarray(steps.attrs["NSteps"]).flat[0]) == 3
                    np.testing.assert_array_equal(steps["Values"], [0, 1, 2])
                    np.testing.assert_array_equal(steps["FieldDataOffsets/clock"], [0, 1, 2])
                    np.testing.assert_array_equal(block["FieldData/clock"],
                                                  [[10 * (i + 1) + step] for step in range(3)])
elif mode == "reader":
    try:
        from vtkmodules.vtkIOHDF import vtkHDFReader
        from vtkmodules.vtkCommonCore import vtkVersion
        from vtkmodules.util.numpy_support import vtk_to_numpy
    except ImportError as exc:
        print(exc)
        sys.exit(77)
    if tuple(map(int, vtkVersion.GetVTKVersion().split(".")[:3])) < (9, 7, 0):
        print("RectilinearGrid and StructuredGrid require VTK 9.7 or later")
        sys.exit(77)
    from vtkhdf_reader_smoke import iter_leaf_datasets

    # Mixed static/temporal blocks are checked at the HDF5 level above because
    # VTK readers have trouble with that mixture (VTK issue 19963):
    # https://gitlab.kitware.com/vtk/vtk/-/work_items/19963
    # Both homogeneous timelines still contain all five dataset types.
    # The composite reader also omits leaf FieldData arrays; their values and
    # offsets are checked independently by the layout test.
    for case in ("static", "temporal"):
        reader = vtkHDFReader()
        errors = []
        reader.AddObserver("ErrorEvent", lambda *args: errors.append(args))
        reader.SetFileName(str(Path(directory) / f"{prefix}_collection_{case}.vtkhdf"))
        reader.UpdateInformation()
        assert not errors, errors
        if case == "temporal":
            assert reader.GetNumberOfSteps() == 3
        for step in range(3 if case == "temporal" else 1):
            if case == "temporal":
                reader.SetStep(step)
            reader.Update()
            assert not errors, errors
            output = reader.GetOutputDataObject(0)
            assert output.IsA("vtkMultiBlockDataSet")
            assert output.GetNumberOfBlocks() == len(types)
            for i, kind in enumerate(types):
                leaves = list(iter_leaf_datasets(output.GetBlock(i)))
                assert leaves and all(leaf.IsA("vtk" + kind) for leaf in leaves)
                if kind != "HyperTreeGrid":
                    assert sum(leaf.GetNumberOfPoints() for leaf in leaves) == 4
                else:
                    values = np.concatenate([vtk_to_numpy(leaf.GetCellData().GetArray("value")) for leaf in leaves])
                    np.testing.assert_array_equal(values, [5.])
                assert sum(leaf.GetNumberOfCells() for leaf in leaves) == 1
else:
    raise ValueError(mode)
