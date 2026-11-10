#!/usr/bin/env python3
"""Independent checks for moving structured geometry and snapshot reuse."""
import sys
import warnings
from pathlib import Path

try:
    import h5py
    import numpy as np
except ImportError as exc:
    print(exc)
    sys.exit(77)

mode, directory, *args = sys.argv[1:]
directory = Path(directory)
if mode == "reader":
    try:
        from vtkmodules.vtkCommonCore import vtkVersion
        from vtkmodules.vtkIOHDF import vtkHDFReader
        from vtkmodules.util.numpy_support import vtk_to_numpy
    except ImportError as exc:
        print(exc)
        sys.exit(77)
    version = (vtkVersion.GetVTKMajorVersion(), vtkVersion.GetVTKMinorVersion(), vtkVersion.GetVTKBuildVersion())
    if version < (9, 7, 0):
        print("Moving structured grids require VTK 9.7 or later")
        sys.exit(77)
    mask = int(args[0])
    if mask in (2, 4, 5, 6) and version < (9, 7, 2):
        print("VTK before 9.7.2 has a structured-data singleton-axis reader defect")
        sys.exit(77)
    masks = [mask]
else:
    masks = range(1, 8)

GEOMETRY = np.array([0, 0, 1, 1, 2, 2, 2])
DATA = np.array([0, 0, 0, 1, 1, 2, 2])
PHASES = [0, 2, 4]
DATA_PHASES = [0, 3, 5]


def geometry(dims, phase):
    k, j, i = np.indices(dims[::-1])
    return np.stack((i + .125*j*j + .5*phase,
                     2*j + .25*k*k + .25*phase,
                     3*k + .5*i*j + .125*phase), axis=-1)


def axes(dims, phase):
    return [np.arange(dims[0])**2 + .5*phase,
            2*np.arange(dims[1])+1 + .25*phase,
            -np.arange(dims[2])-1 + .125*phase]


def reader_for(path, cache=True):
    reader = vtkHDFReader()
    if not cache and hasattr(reader, "SetUseCache"):
        # VTK 9.7.1 can leave the last data array attached on a cache hit
        # after a nonsequential request. This does not affect file offsets.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            reader.SetUseCache(False)
    errors = []
    reader.AddObserver("ErrorEvent", lambda *event: errors.append(event))
    reader.SetFileName(str(path))
    reader.UpdateInformation()
    assert not errors, (path, errors)
    return reader, errors


def check_coordinates(data, grid, dims, phase):
    if grid == "rectilinear":
        for actual, expected in zip((data.GetXCoordinates(), data.GetYCoordinates(), data.GetZCoordinates()), axes(dims, phase)):
            np.testing.assert_array_equal(vtk_to_numpy(actual), expected)
    else:
        np.testing.assert_array_equal(vtk_to_numpy(data.GetPoints().GetData()), geometry(dims, phase).reshape(-1, 3))


for mask in masks:
    dims = tuple(n if mask & (1 << a) else 1 for a, n in enumerate((5, 3, 2)))
    cdims = tuple(max(1, n-1) for n in dims)
    k, j, i = np.indices(dims[::-1])
    value = i + 10*j + 100*k
    ck, cj, ci = np.indices(cdims[::-1])
    cells = np.stack([ci+10*cj+100*ck+1000*c for c in (1, 2)], axis=-1)
    point_shape = tuple(n for n in dims[::-1] if n > 1)
    cell_shape = tuple(n for n, p in zip(cdims[::-1], dims[::-1]) if p > 1) + (2,)
    for grid in ("rectilinear", "structured"):
        for kind, width in (("real32", 4), ("real64", 8)):
            for before in (0, 1):
                path = directory / f"moving_{grid}_{kind}_{mask}_{before}.vtkhdf"
                if mode == "layout":
                    with h5py.File(path) as f:
                        root = f["VTKHDF"]
                        assert root.attrs["Type"] == (b"RectilinearGrid" if grid == "rectilinear" else b"StructuredGrid")
                        np.testing.assert_array_equal(root.attrs["Dimensions"], dims)
                        assert root["Steps"].attrs["NSteps"] == 7
                        np.testing.assert_array_equal(root["Steps/Values"], np.arange(7))
                        if grid == "structured":
                            a = root["Points"]
                            assert a.shape == (3, *dims[::-1], 3)
                            assert a.maxshape == (None, *dims[::-1], 3)
                            assert a.chunks[0] == 1
                            assert np.prod(a.chunks)*width <= 64
                            assert a.dtype.itemsize == width
                            np.testing.assert_array_equal(a, np.stack([geometry(dims, phase) for phase in PHASES]))
                            np.testing.assert_array_equal(root["Steps/PointOffsets"], GEOMETRY)
                        else:
                            for a, axis in enumerate("XYZ"):
                                payload = root[axis+"Coordinates"]
                                assert payload.shape == (3*dims[a],)
                                assert payload.maxshape == (None,)
                                assert payload.dtype.itemsize == width
                                np.testing.assert_array_equal(payload, np.concatenate([axes(dims, phase)[a] for phase in PHASES]))
                                offset = root[f"Steps/{axis}CoordinatesOffsets"]
                                assert offset.dtype.itemsize == 8
                                np.testing.assert_array_equal(offset, GEOMETRY*dims[a])
                        np.testing.assert_array_equal(root["PointData/value"], np.stack([(value+10*p).reshape(point_shape) for p in DATA_PHASES]))
                        np.testing.assert_array_equal(root["CellData/vector"], np.stack([(cells+10*p).reshape(cell_shape) for p in DATA_PHASES]))
                        np.testing.assert_array_equal(root["Steps/PointDataOffsets/value"], DATA)
                        np.testing.assert_array_equal(root["Steps/CellDataOffsets/vector"], DATA)
                        for group, name, expected in (
                                ("Point", "late_points", value+700),
                                ("Cell", "between_cells", ci+10*cj+100*ck+700),
                                ("Cell", "end_cells", ci+10*cj+100*ck+700)):
                            shape = point_shape if group == "Point" else cell_shape[:-1]
                            payload = root[f"{group}Data/{name}"]
                            assert payload.shape == (1, *shape)
                            assert payload.maxshape == payload.shape
                            np.testing.assert_array_equal(payload[0], expected.reshape(shape))
                            np.testing.assert_array_equal(root[f"Steps/{group}DataOffsets/{name}"], np.zeros(7))
                        np.testing.assert_array_equal(root["FieldData/energy"], np.array([42, 44, 46, 48])[:, None])
                        np.testing.assert_array_equal(root["Steps/FieldDataOffsets/energy"], [0, 0, 1, 1, 2, 2, 3])
                else:
                    for cache in (True, False):
                        reader, errors = reader_for(path, cache=cache)
                        assert reader.GetNumberOfSteps() == 7
                        for step in (range(7) if cache else (6, 0, 4, 1, 3, 2, 5, 6)):
                            reader.SetStep(step)
                            reader.Update()
                            assert not errors, (path, errors)
                            output = reader.GetOutputDataObject(0)
                            assert output.GetClassName() == ("vtkRectilinearGrid" if grid == "rectilinear" else "vtkStructuredGrid")
                            assert output.GetNumberOfPoints() == np.prod(dims)
                            assert output.GetNumberOfCells() == np.prod(cdims)
                            check_coordinates(output, grid, dims, PHASES[GEOMETRY[step]])
                            np.testing.assert_array_equal(vtk_to_numpy(output.GetPointData().GetArray("value")), value.ravel()+10*DATA_PHASES[DATA[step]])
                            np.testing.assert_array_equal(vtk_to_numpy(output.GetCellData().GetArray("vector")), cells.reshape(-1, 2)+10*DATA_PHASES[DATA[step]])
                            np.testing.assert_array_equal(vtk_to_numpy(output.GetFieldData().GetArray("energy")), [42+step-step % 2])
                            np.testing.assert_array_equal(vtk_to_numpy(output.GetPointData().GetArray("late_points")), value.ravel()+700)
                            for name in ("between_cells", "end_cells"):
                                np.testing.assert_array_equal(vtk_to_numpy(output.GetCellData().GetArray(name)),
                                                              (ci+10*cj+100*ck+700).ravel())
            if mask == 7:
                path = directory / f"moving_{grid}_{kind}_reopen.vtkhdf"
                if mode == "layout":
                    with h5py.File(path) as f:
                        g = f["VTKHDF"]
                        assert g["Steps"].attrs["NSteps"] == 1
                        names = ["PointOffsets"] if grid == "structured" else [a+"CoordinatesOffsets" for a in "XYZ"]
                        for name in names:
                            np.testing.assert_array_equal(g[f"Steps/{name}"], [0])
                        if grid == "structured":
                            np.testing.assert_array_equal(g["Points"], geometry(dims, 0)[None])
                        else:
                            for axis, expected in zip("XYZ", axes(dims, 0)):
                                np.testing.assert_array_equal(g[axis+"Coordinates"], expected)
                else:
                    reader, errors = reader_for(path)
                    reader.Update()
                    assert not errors, (path, errors)
                    check_coordinates(reader.GetOutputDataObject(0), grid, dims, 0)

if mode == "layout" or masks == [7]:
    for kind in ("real32", "real64"):
        path = directory / f"moving_structured_{kind}_empty.vtkhdf"
        if mode == "layout":
            with h5py.File(path) as f:
                g = f["VTKHDF"]
                assert g["Points"].shape == (2, 2, 2, 2, 3)
                np.testing.assert_array_equal(g["Steps/PointOffsets"], [0, 0, 1, 1])
                np.testing.assert_array_equal(g["Points"][0], np.ones((2, 2, 2, 3)))
                assert np.isnan(g["Points"][1]).all()
        else:
            reader, errors = reader_for(path)
            for step in (3, 0, 2, 1):
                reader.SetStep(step)
                reader.Update()
                assert not errors, (path, errors)
                points = vtk_to_numpy(reader.GetOutputDataObject(0).GetPoints().GetData())
                if step < 2:
                    np.testing.assert_array_equal(points, np.ones((8, 3)))
                else:
                    assert np.isnan(points).all()
print(f"moving geometry {mode} checks passed")
