#!/usr/bin/env python3
"""Independent layout and VTK readback for RectilinearGrid/StructuredGrid."""
import sys
from pathlib import Path

try:
    import numpy as np
    import h5py
except ImportError as exc:
    print(exc)
    sys.exit(77)

mode, directory, *args = sys.argv[1:]
if mode == "reader":
    try:
        from vtkmodules.vtkIOHDF import vtkHDFReader
        from vtkmodules.vtkCommonCore import vtkVersion
        from vtkmodules.util.numpy_support import vtk_to_numpy
    except ImportError as exc:
        print(exc)
        sys.exit(77)
    version = tuple(map(int, vtkVersion.GetVTKVersion().split(".")[:3]))
    if version < (9, 7, 0):
        print("RectilinearGrid and StructuredGrid require VTK 9.7 or later")
        sys.exit(77)
    mask = int(args[0])
    if mask in (2, 4, 5, 6) and version < (9, 7, 2):
        print("VTK before 9.7.2 has a shared structured-data singleton-axis reader defect")
        sys.exit(77)
    masks = [mask]
else:
    masks = range(1, 8)

for mask in masks:
    dims = tuple(n if mask & (1 << axis) else 1 for axis, n in enumerate((5, 3, 2)))
    cdims = tuple(max(1, n - 1) for n in dims)
    k, j, i = np.indices(dims[::-1])
    value = i + 10*j + 100*k
    ck, cj, ci = np.indices(cdims[::-1])
    cell = ci + 10*cj + 100*ck
    xyz = np.stack((i + .125*j*j, 2*j + .25*k*k, 3*k + .5*i*j), axis=-1)
    axes = [np.arange(dims[0])**2, 2*np.arange(dims[1])**2+1, 3*np.arange(dims[2])**2-1]
    for grid in ("rectilinear", "structured"):
        for kind, width in (("real32", 4), ("real64", 8)):
            for temporal in (0, 1):
                path = Path(directory) / f"{grid}_{kind}_{mask}_{temporal}.vtkhdf"
                if mode == "layout":
                    with h5py.File(path) as f:
                        root = f["VTKHDF"]
                        assert root.attrs["Type"] == (b"RectilinearGrid" if grid == "rectilinear" else b"StructuredGrid")
                        assert root.attrs["Version"].shape == (2,)
                        np.testing.assert_array_equal(root.attrs["Dimensions"], dims)
                        assert "WholeExtent" not in root.attrs
                        if grid == "rectilinear":
                            for axis, expected in zip("XYZ", axes):
                                a = root[f"{axis}Coordinates"]
                                assert a.dtype.itemsize == width
                                np.testing.assert_array_equal(a, expected)
                        else:
                            a = root["Points"]
                            assert a.shape == (*dims[::-1], 3)
                            assert a.dtype.itemsize == width
                            assert np.prod(a.chunks) * width <= 64
                            np.testing.assert_array_equal(a, xyz)
                        # Data omits axes where the POINT dimension is singleton;
                        # active cell axes of length one must still be retained.
                        shape = tuple(n for n, p in zip(cdims[::-1], dims[::-1]) if p > 1)
                        expected_cell = cell.reshape(shape)
                        expected_value = np.squeeze(value)
                        if temporal:
                            assert root["Steps"].attrs["NSteps"] == 3
                            np.testing.assert_array_equal(root["Steps/Values"], [0, 1, 2])
                            np.testing.assert_array_equal(root["Steps/PointDataOffsets/value"], [0, 0, 1])
                            np.testing.assert_array_equal(root["Steps/CellDataOffsets/cell"], [0, 0, 0])
                            for offset in ("PointOffsets", "XCoordinatesOffsets", "YCoordinatesOffsets", "ZCoordinatesOffsets"):
                                assert offset not in root["Steps"]
                            expected_value = np.stack([expected_value, expected_value + 200])
                            expected_cell = expected_cell[None]
                        np.testing.assert_array_equal(root["PointData/value"], expected_value)
                        np.testing.assert_array_equal(root["CellData/cell"], expected_cell)
                else:
                    reader = vtkHDFReader()
                    errors = []
                    reader.AddObserver("ErrorEvent", lambda *event: errors.append(event))
                    reader.SetFileName(str(path))
                    reader.UpdateInformation()
                    assert reader.GetNumberOfSteps() == (3 if temporal else 1)
                    for step in range(3 if temporal else 1):
                        reader.SetStep(step)
                        reader.Update()
                        assert not errors, path
                        data = reader.GetOutputDataObject(0)
                        assert data.GetClassName() == ("vtkRectilinearGrid" if grid == "rectilinear" else "vtkStructuredGrid")
                        assert data.GetExtent() == tuple(v for n in dims for v in (0, n-1))
                        assert data.GetNumberOfPoints() == np.prod(dims)
                        assert data.GetNumberOfCells() == np.prod(cdims)
                        if grid == "rectilinear":
                            for a, expected in zip((data.GetXCoordinates(), data.GetYCoordinates(), data.GetZCoordinates()), axes):
                                np.testing.assert_array_equal(vtk_to_numpy(a), expected)
                        else:
                            np.testing.assert_array_equal(vtk_to_numpy(data.GetPoints().GetData()), xyz.reshape(-1, 3))
                        np.testing.assert_array_equal(vtk_to_numpy(data.GetPointData().GetArray("value")), value.ravel() + (200 if step == 2 else 0))
                        np.testing.assert_array_equal(vtk_to_numpy(data.GetCellData().GetArray("cell")), cell.ravel())
                        np.testing.assert_array_equal(vtk_to_numpy(data.GetFieldData().GetArray("energy")), [42+step])
                        assert data.GetPointData().GetScalars().GetName() == "value"
print(f"{mode} checks passed")
