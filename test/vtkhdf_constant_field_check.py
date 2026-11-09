#!/usr/bin/env python3
"""Check constant field values and complete temporal metadata independently."""
import sys
from pathlib import Path

try:
    import h5py
    import numpy as np
except ImportError as exc:
    print(exc)
    sys.exit(77)

mode, directory = sys.argv[1:]
directory = Path(directory)
if mode == "reader":
    try:
        from vtkmodules.vtkIOHDF import vtkHDFReader
        from vtkmodules.util.numpy_support import vtk_to_numpy
    except ImportError as exc:
        print(exc)
        sys.exit(77)

payloads = {
    "scalar": np.array([42]),
    "tuples": np.array([11, 12, 13, 14]),
    "vector": np.array([[11, 12, 13]]),
    "matrix": np.array([[31, 32], [33, 34], [35, 36]]),
    "empty": np.array([]),
}
kinds = {"int32": ("i", 4), "int64": ("i", 8), "real32": ("f", 4), "real64": ("f", 8)}
for grid, modes in ((1, range(4)), (2, range(2)), (3, range(3)), (4, range(3))):
    for file_mode in modes:
        path = directory / f"constant_fields_{grid}_{file_mode}.vtkhdf"
        temporal = file_mode != 0
        prefixes = ("before", "open", "between", "after") if temporal else ("before", "after")
        fields = {}
        for prefix in prefixes:
            for kind, dtype in kinds.items():
                for suffix, expected in payloads.items():
                    if suffix == "empty" and not temporal:
                        continue
                    fields[f"{prefix}_{kind}_{suffix}"] = (expected, dtype)
        if temporal:
            fields["before_real64_scalar_2"] = (np.array([222]), ("i", 4))
        if mode == "layout":
            with h5py.File(path) as f:
                g = f["VTKHDF"]
                assert ("Steps" in g) == temporal
                assert set(g["FieldData"]) == set(fields) | ({"before_real64_scalar_1"} if temporal else set())
                for name, (expected, dtype) in fields.items():
                    a = g[f"FieldData/{name}"]
                    nc = expected.shape[1] if expected.ndim == 2 else 1
                    nt = len(expected)
                    stored = expected.reshape(nt, nc) if temporal else expected
                    np.testing.assert_array_equal(a, stored)
                    assert (a.dtype.kind, a.dtype.itemsize) == dtype
                    original = "before_real64_scalar" if name.endswith("_2") else name
                    assert a.attrs["Name"] == original.encode()
                    if temporal:
                        np.testing.assert_array_equal(g[f"Steps/FieldDataOffsets/{name}"], np.zeros(4))
                        np.testing.assert_array_equal(g[f"Steps/FieldDataSizes/{name}"], np.tile([nc, nt], (4, 1)))
                if temporal:
                    name = "before_real64_scalar_1"
                    a = g[f"FieldData/{name}"][:]
                    np.testing.assert_array_equal(a[0, :1], [111])
                    np.testing.assert_array_equal(a[1], [301, 302])
                    np.testing.assert_array_equal(g[f"Steps/FieldDataOffsets/{name}"], [0, 0, 0, 1])
                    np.testing.assert_array_equal(g[f"Steps/FieldDataSizes/{name}"], [[1, 0], [1, 1], [1, 1], [2, 1]])
        else:
            reader = vtkHDFReader()
            errors = []
            reader.AddObserver("ErrorEvent", lambda *args: errors.append(args))
            reader.SetFileName(str(path))
            reader.UpdateInformation()
            for step in range(4) if temporal else [0]:
                if temporal:
                    reader.SetStep(step)
                reader.Update()
                assert not errors, (path, errors)
                fd = reader.GetOutputDataObject(0).GetFieldData()
                for name, (expected, dtype) in fields.items():
                    a = fd.GetArray(name)
                    assert a is not None, (path, step, name)
                    np.testing.assert_array_equal(vtk_to_numpy(a), expected)
                    assert a.GetNumberOfComponents() == (expected.shape[1] if expected.ndim == 2 else 1)
                    assert a.GetNumberOfTuples() == len(expected)
                    actual_dtype = vtk_to_numpy(a).dtype
                    assert (actual_dtype.kind, actual_dtype.itemsize) == dtype
                if temporal:
                    expected = np.array([]) if step == 0 else np.array([111]) if step < 3 else np.array([[301, 302]])
                    np.testing.assert_array_equal(vtk_to_numpy(fd.GetArray("before_real64_scalar_1")), expected)
print(f"Constant field {mode} checks passed for all writers and modes")
