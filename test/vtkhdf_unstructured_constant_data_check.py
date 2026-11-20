#!/usr/bin/env python3
"""Check constant unstructured arrays independently of the Fortran writer."""
from pathlib import Path
import sys

try:
    import h5py
    import numpy as np
except ImportError as exc:
    print(exc)
    sys.exit(77)

check, directory = sys.argv[1:]
directory = Path(directory)
if check == "reader":
    try:
        from vtkmodules.vtkIOHDF import vtkHDFReader
        from vtkmodules.util.numpy_support import vtk_to_numpy
    except ImportError as exc:
        print(exc)
        sys.exit(77)

kinds = {"int8": "uint8", "int32": "int32", "int64": "int64",
         "real32": "float32", "real64": "float64"}
prefixes = ("before", "during_first", "during_later", "after_data", "between", "after")
point_steps = [0, 1, 1, 3, 3]
cell_steps = [0, 0, 2, 3, 3]


def values(nproc, nw, association):
    parts = []
    for rank in range(nproc):
        nc = 0 if nproc > 1 and rank == 0 else rank + 1
        count = nc * (4 if association == "point" else 1)
        parts.append(50 * (nw - 1) + 10 * rank + np.arange(1, count + 1))
    return np.concatenate(parts)


def coordinates(nproc, nw, shift):
    tetra = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    parts = []
    for rank in range(nproc):
        nc = 0 if nproc > 1 and rank == 0 else rank + 1
        for cell in range(nc):
            parts.append(tetra + [100 * (nw - 1) + 10 * rank + 2 * cell, 0, shift])
    return np.concatenate(parts)


def constants(nproc, nw, association, mode):
    scalar = values(nproc, nw, association)
    result = {}
    for prefix in prefixes if mode else ("before", "after"):
        for kind, dtype in kinds.items():
            for suffix in ("scalar", "vector"):
                original = f"{prefix} {association} {kind} {suffix}"
                data = scalar if suffix == "scalar" else scalar[:, None] + [1, 2, 3]
                attribute = ("Scalars" if suffix == "scalar" else "Vectors") if prefix == "before" and kind == "real64" else None
                result[original.replace(" ", "_")] = (data.astype(dtype), original, attribute)
    if mode:
        first = "shared_name" if association == "point" else "shared_name_1"
        result[first] = (scalar.astype("float64"), "shared.name" if association == "point" else "shared name", None)
        result["shared_name_2"] = (scalar.astype("float64"), "shared name", None)
    return result


def leaves(data):
    if data.IsA("vtkDataSet"):
        return [data]
    return [data.GetPartition(i) for i in range(data.GetNumberOfPartitions())]


for mode, composite in [(0, False), (1, False), (2, False), (2, True)]:
    filename = "unstructured_constants_mb.vtkhdf" if composite else f"unstructured_constants_{mode}.vtkhdf"
    path = directory / filename
    roots = [f"VTKHDF/Group_{i}" for i in range(2)] if composite else ["VTKHDF"]
    with h5py.File(path) as f:
        nproc = int(f[f"{roots[0]}/FieldData/nproc"][...].ravel()[0])
        if composite:
            assert set(f["VTKHDF/Assembly"]) == {"solid", "solid_1"}
        if check == "layout":
            for nw, root in enumerate(roots, 1):
                g = f[root]
                assert ("Steps" in g) == bool(mode)
                npoint = len(values(nproc, nw, "point"))
                snapshots = 3 if mode == 2 else 1
                assert g["Points"].shape == (snapshots * npoint, 3)
                for i in range(snapshots):
                    np.testing.assert_array_equal(g["Points"][i * npoint:(i + 1) * npoint],
                                                  coordinates(nproc, nw, i))
                if mode:
                    assert g["Steps"].attrs["NSteps"] == 5
                    np.testing.assert_array_equal(g["Steps/Values"], np.arange(5))
                    offsets = np.array([0, 1, 1, 2, 2]) * npoint if mode == 2 else np.zeros(5)
                    np.testing.assert_array_equal(g["Steps/PointOffsets"], offsets)
                    np.testing.assert_array_equal(g["Steps/CellOffsets"], np.zeros(5))
                for association, group, steps in [("point", "PointData", point_steps),
                                                   ("cell", "CellData", cell_steps)]:
                    constant = constants(nproc, nw, association, mode)
                    changing = "shared_name_1" if association == "point" else "shared_name"
                    assert set(g[group]) == set(constant) | ({changing} if mode else set())
                    for name, (expected, original, attribute) in constant.items():
                        array = g[f"{group}/{name}"]
                        np.testing.assert_array_equal(array, expected)
                        assert array.dtype == expected.dtype
                        assert array.attrs["Name"] == original.encode()
                        if attribute:
                            assert array.attrs["Attribute"] == attribute.encode()
                        else:
                            assert "Attribute" not in array.attrs
                        offset_path = f"Steps/{group}Offsets/{name}"
                        if mode == 2:
                            offsets = g[offset_path]
                            assert offsets.dtype == np.dtype("int64")
                            assert offsets.maxshape == (None,)
                            np.testing.assert_array_equal(offsets, np.zeros(5))
                        else:
                            assert offset_path not in g
                    if mode:
                        base = values(nproc, nw, association)
                        written = [0, 1, 3] if association == "point" else [0, 2, 3]
                        np.testing.assert_array_equal(g[f"{group}/{changing}"],
                                                      np.concatenate([base + 1000 * s for s in written]))
                        np.testing.assert_array_equal(g[f"Steps/{group}Offsets/{changing}"],
                                                      [written.index(s) * len(base) for s in steps])
    if check == "reader":
        # A fresh reader at step 3 also checks offsets without a preceding cache hit.
        for sequence in ([0, 1, 2, 3, 4, 1], [3]) if mode else ([0],):
            reader = vtkHDFReader()
            errors = []
            reader.AddObserver("ErrorEvent", lambda *args: errors.append(args))
            reader.GetExecutive().AddObserver("ErrorEvent", lambda *args: errors.append(args))
            reader.SetFileName(str(path))
            reader.UpdateInformation()
            for step in sequence:
                reader.SetStep(step)
                reader.Update()
                assert not errors, (filename, step, errors)
                output = reader.GetOutputDataObject(0)
                for nw in range(1, len(roots) + 1):
                    data = output.GetBlock(nw - 1) if composite else output
                    parts = leaves(data)
                    actual = np.concatenate([vtk_to_numpy(p.GetPoints().GetData()) for p in parts])
                    shift = [0, 1, 1, 2, 2][step] if mode == 2 else 0
                    np.testing.assert_array_equal(actual, coordinates(nproc, nw, shift))
                    for association, method, steps in [("point", "GetPointData", point_steps),
                                                        ("cell", "GetCellData", cell_steps)]:
                        expected = {name: a for name, (a, _, _) in constants(nproc, nw, association, mode).items()}
                        if mode:
                            changing = "shared_name_1" if association == "point" else "shared_name"
                            expected[changing] = (values(nproc, nw, association) + 1000 * steps[step]).astype("float64")
                        for name, want in expected.items():
                            arrays = [getattr(p, method)().GetArray(name) for p in parts]
                            assert all(a is not None for a in arrays), (filename, step, name)
                            actual = np.concatenate([vtk_to_numpy(a) for a in arrays])
                            np.testing.assert_array_equal(actual, want, err_msg=f"{filename}, step {step}, {name}")
                            assert actual.dtype == want.dtype
                            assert arrays[0].GetNumberOfComponents() == (want.shape[1] if want.ndim == 2 else 1)
print(f"Unstructured constant-data {check} checks passed")
