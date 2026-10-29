#!/usr/bin/env python3
"""ImageData semantic readback and independent HDF5 layout checks."""

from __future__ import annotations

import argparse
from pathlib import Path

SKIP = 77
KINDS = {
    "int8": ("u", 1),
    "int32": ("i", 4),
    "int64": ("i", 8),
    "real32": ("f", 4),
    "real64": ("f", 8),
}
DIRECTION = [[1.0, 0.4, 0.1], [0.2, 2.0, 0.5], [0.0, 0.3, 1.5]]
FIELD_VALUES = [
    None,
    [42],
    [21, 22],
    [[16, 32]],
    [[1, 2, 3], [7, 8, 9]],
    [[1, 2, 3], [7, 8, 9]],
    [7],
    None,
    None,
    [[8, 9]],
    [[8, 9]],
]
FIELD_OFFSETS = [0, 0, 1, 3, 4, 4, 6, 7, 7, 7, 7]
FIELD_SIZES = [
    [1, 0],
    [1, 1],
    [1, 2],
    [2, 1],
    [3, 2],
    [3, 2],
    [1, 1],
    [1, 0],
    [1, 0],
    [2, 1],
    [2, 1],
]
REGISTRATION_CELLS = {
    "before": [91, 92],
    "before_1": [11, 11],
    "before_2": [12, 12],
    "before_3": [95, 96],
    "after": [13, 13],
    "after_1": [14, 14],
    "after_2": [93, 94],
}
REGISTRATION_POINTS = {
    "Unnamed_PointData": 21,
    "Unnamed_PointData_1": 22,
    "a_b": 23,
    "a_b_1": 24,
}


def naming_cases(group):
    # Explicit expected names also cover collisions with caller-supplied suffixes.
    return [
        ("Unnamed_" + group, "Unnamed_" + group),
        ("Unnamed_" + group + "_1", "Unnamed_" + group),
        ("a_b", "a/b"),
        ("a_b_1", "a.b"),
        ("a_b_1_1", "a_b_1"),
        ("a_b_2", "a b"),
        ("a_b_3", "a/b"),
        ("a_b_1_2", "a_b_1"),
        ("trimmed", "trimmed"),
        ("trimmed_1", "trimmed"),
        ("same", "same"),
        ("same_1", "same"),
    ]


def dependencies(reader=False):
    try:
        import numpy as np
        import h5py

        if reader:
            from vtkmodules.vtkIOHDF import vtkHDFReader
            from vtkmodules.vtkCommonCore import vtkVersion
            from vtkmodules.util.numpy_support import vtk_to_numpy

            return np, h5py, vtkHDFReader, vtkVersion, vtk_to_numpy
        return np, h5py
    except ImportError as exc:
        print(f"ImageData validation dependency unavailable: {exc}")
        raise SystemExit(SKIP)


def encoded(np, lower, dims, component=None, modulo=False, fraction=False):
    z, y, x = np.indices(tuple(dims[::-1]))
    values = (x + lower[0]) + 11 * (y + lower[1]) + 37 * (z + lower[2])
    if component is not None:
        values = values + 19 * component if modulo else values + 1000 * component
    if modulo:
        values %= 101
    if fraction:
        values = values + 0.25
    return values


def equal(np, actual, expected):
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-6)


def check_layout(directory):
    np, h5py = dependencies()
    with h5py.File(directory / "image_static.vtkhdf") as f:
        g = f["VTKHDF"]
        assert g.attrs["WholeExtent"].dtype.itemsize == 4
        equal(np, g.attrs["WholeExtent"], [2, 8, -3, 1, 4, 7])
        equal(np, g.attrs["Direction"], np.array(DIRECTION).ravel())
        equal(np, g.attrs["Origin"], [1, 2, 3])
        equal(np, g.attrs["Spacing"], [0.5, 0.25, 2])
        for association, dims in [("CellData", [6, 4, 3]), ("PointData", [7, 5, 4])]:
            for kind, (category, width) in KINDS.items():
                a = g[f"{association}/{kind}-scalar"]
                assert a.shape == tuple(dims[::-1])
                assert (a.dtype.kind, a.dtype.itemsize) == (category, width)
                equal(
                    np,
                    a[:],
                    encoded(
                        np,
                        [2, -3, 4],
                        dims,
                        modulo=True,
                        fraction=kind.startswith("real"),
                    ),
                )
                v = g[f"{association}/{kind}-vector"]
                assert v.shape == (*dims[::-1], 3)
                expected = np.stack(
                    [
                        encoded(np, [2, -3, 4], dims, c, True, kind.startswith("real"))
                        for c in range(1, 4)
                    ],
                    axis=-1,
                )
                equal(np, v[:], expected)
                assert np.prod(v.chunks) * width <= max(64, width)
        equal(np, g["CellData/a_b"][:], np.ones((3, 4, 6)))
        equal(np, g["CellData/a_b_1"][:], np.full((3, 4, 6), 2))
        assert g["CellData/a_b"].attrs["Name"] == b"a/b"
        assert g["CellData/a_b_1"].attrs["Name"] == b"a.b"
        check_static_fields(np, g["FieldData"])
    with h5py.File(directory / "image_temporal.vtkhdf") as f:
        g = f["VTKHDF"]
        assert g["Steps"].attrs["NSteps"] == 5
        equal(np, g["Steps/Values"][:], range(5))
        for name, location, offsets, phases, dims in [
            ("temperature", "CellData", [0, 0, 1, 1, 1], [0, 2], [6, 4, 3]),
            ("other", "CellData", [0, 0, 0, 0, 1], [0, 4], [6, 4, 3]),
            ("velocity", "PointData", [0, 0, 0, 1, 1], [0, 3], [7, 5, 4]),
        ]:
            a = g[f"{location}/{name}"]
            assert a.shape[0] == 2, (name, a.shape)
            assert a.chunks[0] == 1
            o = g[f"Steps/{location}Offsets/{name}"]
            assert o.dtype.itemsize == 8
            equal(np, o[:], offsets)
            for index, phase in enumerate(phases):
                expected = encoded(np, [0, 0, 0], dims) + phase * 1000
                if name == "velocity":
                    expected = np.stack(
                        [expected + c * 10000 for c in range(1, 4)], axis=-1
                    )
                equal(np, a[index], expected)
        for name in ["static_early", "static_late"]:
            a = g[f"CellData/{name}"]
            assert a.shape == (1, 3, 4, 6)
            equal(np, a[:], np.full(a.shape, 17))
            equal(np, g[f"Steps/CellDataOffsets/{name}"][:], [0] * 5)
        equal(np, g["FieldData/energy"][:], [[1], [3]])
        equal(np, g["Steps/FieldDataOffsets/energy"][:], [0, 0, 0, 1, 1])
        equal(
            np,
            g["Steps/FieldDataSizes/energy"][:],
            [[1, 0], [1, 1], [1, 1], [1, 1], [1, 1]],
        )
        equal(np, g["FieldData/history"][:], [[10, 20], [1, 2], [3, 4], [5, 6]])
        equal(np, g["Steps/FieldDataOffsets/history"][:], [0, 0, 0, 1, 1])
        equal(
            np,
            g["Steps/FieldDataSizes/history"][:],
            [[1, 0], [2, 1], [2, 1], [2, 3], [2, 3]],
        )
    for index, dims in enumerate(
        [[7, 5, 1], [1, 5, 4], [7, 1, 4], [7, 1, 1], [1, 5, 1], [1, 1, 4], [7, 5, 4]], 1
    ):
        lower = [-2, -3, -4] if index == 7 else [0, 0, 0]
        active = np.array(dims) > 1
        with h5py.File(directory / f"image_dim_{index}.vtkhdf") as f:
            for loc in ["PointData", "CellData"]:
                n = (
                    np.array(dims)
                    if loc == "PointData"
                    else np.maximum(1, np.array(dims) - 1)
                )
                shape = tuple(n[::-1][active[::-1]])
                for name in ["scalar", "vector"]:
                    a = f[f"VTKHDF/{loc}/{name}"]
                    assert a.shape == (shape if name == "scalar" else (*shape, 2))
                    expected = encoded(np, lower, n)
                    if name == "vector":
                        expected = np.stack([expected + 1000, expected + 2000], axis=-1)
                    equal(np, a[:], expected.reshape(a.shape))
                if loc == "CellData":
                    equal(
                        np,
                        f[f"VTKHDF/{loc}/noncontiguous"][:],
                        (encoded(np, lower, n) + 1000).reshape(shape),
                    )
    with h5py.File(directory / "image_empty.vtkhdf") as f:
        g = f["VTKHDF"]
        assert g["CellData/empty"].shape == (2, 3, 4, 6)
        equal(np, g["Steps/CellDataOffsets/empty"][:], [0, 1, 1])
    for mode in ["static", "temporal"]:
        image = directory / f"image_{mode}.vtkhdf"
        with (
            h5py.File(image) as f,
            h5py.File(directory / f"image_fields_{mode}_ug.vtkhdf") as reference,
        ):
            for location in ["FieldData"] + (
                ["Steps/FieldDataOffsets", "Steps/FieldDataSizes"]
                if mode == "temporal"
                else []
            ):
                actual = f[f"VTKHDF/{location}"]
                expected = reference[f"VTKHDF/{location}"]
                assert set(actual) == set(expected)
                for name in actual:
                    assert actual[name].dtype == expected[name].dtype
                    equal(np, actual[name][:], expected[name][:])
    with h5py.File(directory / "image_numeric_temporal.vtkhdf") as f:
        for loc, dims in [("CellData", [2, 1, 1]), ("PointData", [3, 2, 2])]:
            for kind in KINDS:
                for name in ["scalar", "vector"]:
                    path = f"{loc}/{kind}-{name}"
                    a = f[f"VTKHDF/{path}"]
                    expected = encoded(np, [0, 0, 0], dims, modulo=True)
                    if name == "vector":
                        expected = np.stack(
                            [
                                encoded(np, [0, 0, 0], dims, c, True)
                                for c in range(1, 3)
                            ],
                            axis=-1,
                        )
                    assert a.shape[0] == 1
                    equal(np, a[0], expected)
                    equal(np, f[f"VTKHDF/Steps/{loc}Offsets/{kind}-{name}"][:], [0, 0])
    for filename in ["image_fields_variable.vtkhdf", "image_fields_variable_ug.vtkhdf"]:
        with h5py.File(directory / filename) as f:
            g = f["VTKHDF"]
            for kind in list(KINDS)[1:]:
                a = g[f"FieldData/{kind}"]
                assert a.shape == (8, 3)
                assert a.maxshape == (None, None)
                assert a.chunks is not None
                assert (a.dtype.kind, a.dtype.itemsize) == KINDS[kind]
                offsets = g[f"Steps/FieldDataOffsets/{kind}"]
                assert offsets.dtype.itemsize == 8
                equal(np, offsets[:], FIELD_OFFSETS)
                equal(np, g[f"Steps/FieldDataSizes/{kind}"][:], FIELD_SIZES)
                for offset, (ncomp, ntuples), expected in zip(
                    FIELD_OFFSETS, FIELD_SIZES, FIELD_VALUES
                ):
                    if expected is not None:
                        equal(
                            np,
                            a[offset : offset + ntuples, :ncomp],
                            np.array(expected).reshape(ntuples, ncomp),
                        )
            assert g["FieldData/never_written"].shape == (0, 1)
            equal(np, g["Steps/FieldDataOffsets/never_written"][:], [0] * 11)
            equal(np, g["Steps/FieldDataSizes/never_written"][:], [[1, 0]] * 11)
    with h5py.File(directory / "image_nested.vtkhdf") as f:
        assert set(f) == {"container"}
        assert set(f["container"]) == {"leaf_static", "leaf_temporal"}
        static = f["container/leaf_static"]
        temporal = f["container/leaf_temporal"]
        for g in [static, temporal]:
            assert g.attrs["Version"].shape == (2,)
        assert static.attrs["Type"] == b"ImageData"
        assert temporal.attrs["Type"] == b"StructuredGrid"
        assert "Steps" not in static
        equal(np, static["CellData/value"][:], [[11]])
        equal(np, temporal["CellData/value"][:], [[[21]]])
        assert temporal["Steps"].attrs["NSteps"] == 2
        equal(np, temporal["Steps/Values"][:], [0, 1])
        equal(np, temporal["Steps/CellDataOffsets/value"][:], [0, 0])
    with h5py.File(directory / "image_registration_names.vtkhdf") as f:
        g = f["VTKHDF"]
        assert set(g["CellData"]) == set(REGISTRATION_CELLS)
        assert set(g["PointData"]) == set(REGISTRATION_POINTS)
        assert set(g["FieldData"]) == {"label", "label_1"}
        for name, expected in REGISTRATION_CELLS.items():
            a = g[f"CellData/{name}"]
            assert a.shape == (1, 1, 2)
            equal(np, a[:], np.array(expected).reshape(1, 1, 2))
            assert a.attrs["Name"] == (
                b"before" if name.startswith("before") else b"after"
            )
            equal(np, g[f"Steps/CellDataOffsets/{name}"][:], [0, 0])
        for name, expected in REGISTRATION_POINTS.items():
            a = g[f"PointData/{name}"]
            assert a.shape == (1, 2, 3)
            equal(np, a[:], np.full(a.shape, expected))
            equal(np, g[f"Steps/PointDataOffsets/{name}"][:], [0, 0])
        assert g["PointData/a_b"].attrs["Name"] == b"a/b"
        assert g["PointData/a_b_1"].attrs["Name"] == b"a.b"
        for name, expected in [("label", 31), ("label_1", 32)]:
            equal(np, g[f"FieldData/{name}"][:], [[expected]])
            assert g[f"FieldData/{name}"].attrs["Name"] == b"label"
            equal(np, g[f"Steps/FieldDataOffsets/{name}"][:], [0, 0])
            equal(np, g[f"Steps/FieldDataSizes/{name}"][:], [[1, 1], [1, 1]])
    with h5py.File(directory / "image_static_names.vtkhdf") as f:
        g = f["VTKHDF"]
        assert "Steps" not in g
        for group in ["CellData", "PointData", "FieldData"]:
            assert set(g[group]) == {"label", "label_1"}
            for name in ["label", "label_1"]:
                assert g[f"{group}/{name}"].attrs["Name"] == b"label"
        for path, shape, value in [
            ("CellData/label", (1, 2), 41),
            ("CellData/label_1", (1, 2), 42),
            ("PointData/label", (2, 3), 51),
            ("PointData/label_1", (2, 3, 2), 52),
            ("FieldData/label", (1,), 61),
            ("FieldData/label_1", (1,), 62),
        ]:
            assert g[path].shape == shape
            equal(np, g[path][:], np.full(shape, value))
        assert g["CellData/label"].dtype.kind == "i"
        assert g["CellData/label_1"].dtype.kind == "f"
    for mode in ["static", "temporal"]:
        with (
            h5py.File(directory / f"image_names_{mode}.vtkhdf") as image,
            h5py.File(directory / f"image_names_{mode}_ug.vtkhdf") as ug,
        ):
            for group in ["CellData", "PointData", "FieldData"]:
                expected_names = {name for name, _ in naming_cases(group)}
                for f in [image, ug]:
                    g = f["VTKHDF"]
                    assert set(g[group]) == expected_names
                    for j, (name, external) in enumerate(naming_cases(group), start=1):
                        a = g[f"{group}/{name}"]
                        assert a.attrs["Name"] == external.encode()
                        equal(np, a[:].ravel()[0], 100 + j)
                        if mode == "temporal":
                            path = f"Steps/{group}Offsets/{name}"
                            # UG omits offsets for static mesh data; ImageData
                            # stores an explicit zero offset for every step.
                            if f is ug and group != "FieldData" and j % 2 == 1:
                                assert path not in g
                            else:
                                equal(np, g[path][:], [0, 0])
                # FieldData payloads have the same representation in both writers.
                for name, _ in naming_cases(group):
                    if group == "FieldData":
                        equal(
                            np,
                            image[f"VTKHDF/{group}/{name}"][:],
                            ug[f"VTKHDF/{group}/{name}"][:],
                        )
    print("ImageData layout, values, UG FieldData parity, and snapshot reuse verified")


def check_static_fields(np, group):
    for kind in KINDS:
        if kind == "int8":
            continue
        equal(np, group[f"{kind}-value"][:], [7])
        equal(np, group[f"{kind}-tuples"][:], [1, 2, 3])
        equal(np, group[f"{kind}-vector"][:], [[4, 5, 6]])
        equal(np, group[f"{kind}-tensor"][:], [[1, 2], [3, 4], [5, 6]])


def read(path, step=None):
    _, _, Reader, _, to_numpy = dependencies(True)
    r = Reader()
    r.SetFileName(str(path))
    if step is not None:
        r.UpdateInformation()
        r.SetStep(step)
    errors = []
    r.AddObserver("ErrorEvent", lambda *_: errors.append(True))
    r.Update()
    assert not errors, f"VTK reader error for {path}"
    output = r.GetOutput()
    assert output is not None and output.IsA("vtkImageData")
    return r, output, to_numpy


def check_reader(directory, case):
    np, _, _, Version, _ = dependencies(True)
    version = (
        Version.GetVTKMajorVersion(),
        Version.GetVTKMinorVersion(),
        Version.GetVTKBuildVersion(),
    )
    # Both geometry reader fixes are queued for VTK 9.7.2. Require corrected
    # versions to pass these checks so a reader regression cannot cause a skip.
    if version < (9, 7, 2) and case == "static":
        print(
            f"VTK {Version.GetVTKVersion()}: positive nonzero extent checks require VTK 9.7.2; layout checks remain enabled"
        )
        raise SystemExit(SKIP)
    if (
        version < (9, 7, 2)
        and case.startswith("dim_")
        and int(case[4:]) in [2, 3, 5, 6]
    ):
        print(
            f"VTK {Version.GetVTKVersion()}: nontrailing singleton checks require VTK 9.7.2; layout checks remain enabled"
        )
        raise SystemExit(SKIP)
    print("VTK", Version.GetVTKVersion(), case)
    if case == "registration_names":
        for step in [0, 1]:
            _, im, to_numpy = read(directory / "image_registration_names.vtkhdf", step)
            assert im.GetCellData().GetNumberOfArrays() == len(REGISTRATION_CELLS)
            assert im.GetPointData().GetNumberOfArrays() == len(REGISTRATION_POINTS)
            for name, expected in REGISTRATION_CELLS.items():
                a = im.GetCellData().GetArray(name)
                assert a is not None and a.GetNumberOfComponents() == 1
                assert a.GetNumberOfTuples() == 2
                equal(np, to_numpy(a), expected)
            for name, expected in REGISTRATION_POINTS.items():
                a = im.GetPointData().GetArray(name)
                assert a is not None and a.GetNumberOfComponents() == 1
                assert a.GetNumberOfTuples() == 6
                equal(np, to_numpy(a), np.full(6, expected))
            for name, expected in [("label", 31), ("label_1", 32)]:
                a = im.GetFieldData().GetArray(name)
                assert a is not None and a.GetNumberOfComponents() == 1
                assert a.GetNumberOfTuples() == 1
                equal(np, to_numpy(a), [expected])
        _, im, to_numpy = read(directory / "image_static_names.vtkhdf")
        for data, tuples, fields in [
            (im.GetCellData(), 2, [("label", 1, 41), ("label_1", 1, 42)]),
            (im.GetPointData(), 6, [("label", 1, 51), ("label_1", 2, 52)]),
            (im.GetFieldData(), 1, [("label", 1, 61), ("label_1", 1, 62)]),
        ]:
            assert data.GetNumberOfArrays() == 2
            for name, components, value in fields:
                a = data.GetArray(name)
                assert a is not None and a.GetNumberOfComponents() == components
                assert a.GetNumberOfTuples() == tuples
                equal(np, to_numpy(a), np.full(to_numpy(a).shape, value))
        _, _, Reader, _, _ = dependencies(True)
        for mode, steps in [("static", [None]), ("temporal", [0, 1])]:
            for suffix in ["", "_ug"]:
                r = Reader()
                r.SetFileName(str(directory / f"image_names_{mode}{suffix}.vtkhdf"))
                errors = []
                r.AddObserver("ErrorEvent", lambda *_: errors.append(True))
                r.UpdateInformation()
                for step in steps:
                    if step is not None:
                        r.SetStep(step)
                    r.Update()
                    output = r.GetOutputDataObject(0)
                    leaf = (
                        output.GetPartition(0)
                        if output.IsA("vtkPartitionedDataSet")
                        else output
                    )
                    for group, data in [
                        ("CellData", leaf.GetCellData()),
                        ("PointData", leaf.GetPointData()),
                        ("FieldData", output.GetFieldData()),
                    ]:
                        names = {
                            data.GetArrayName(i)
                            for i in range(data.GetNumberOfArrays())
                        }
                        expected = {name for name, _ in naming_cases(group)}
                        if mode == "temporal" and group == "FieldData":
                            expected.add("Time")  # Added by the VTK reader.
                        assert names == expected
                        for j, (name, _) in enumerate(naming_cases(group), start=1):
                            a = data.GetArray(name)
                            assert a is not None and a.GetNumberOfComponents() == 1
                            equal(np, to_numpy(a)[0], 100 + j)
                assert not errors
        return
    if case == "temporal_fielddata":
        _, _, Reader, _, to_numpy = dependencies(True)

        def field_data(reader, step):
            reader.SetStep(step)
            reader.Update()
            output = reader.GetOutputDataObject(0)
            # FieldData belongs to the whole output, including a partitioned
            # UG's root rather than its individual partition datasets.
            return output.GetFieldData()

        for filename in ["image_temporal.vtkhdf", "image_fields_temporal_ug.vtkhdf"]:
            r = Reader()
            errors = []
            r.AddObserver("ErrorEvent", lambda *_: errors.append(True))
            r.SetFileName(str(directory / filename))
            r.UpdateInformation()
            for step in [1, 2, 3, 4]:
                a = field_data(r, step).GetArray("history")
                assert a is not None and a.GetNumberOfComponents() == 2
                assert a.GetNumberOfTuples() == (1 if step < 3 else 3)
                expected = [[10, 20]] if step < 3 else [[1, 2], [3, 4], [5, 6]]
                equal(np, to_numpy(a), expected)
            assert not errors, filename
        for filename in [
            "image_fields_variable.vtkhdf",
            "image_fields_variable_ug.vtkhdf",
        ]:
            r = Reader()
            errors = []
            r.AddObserver("ErrorEvent", lambda *_: errors.append(True))
            r.SetFileName(str(directory / filename))
            r.UpdateInformation()
            for step, ((ncomp, ntuples), expected) in enumerate(
                zip(FIELD_SIZES, FIELD_VALUES)
            ):
                fd = field_data(r, step)
                for kind in list(KINDS)[1:]:
                    a = fd.GetArray(kind)
                    if ntuples == 0:
                        assert a is None or a.GetNumberOfTuples() == 0
                    else:
                        assert a is not None and a.GetNumberOfComponents() == ncomp
                        assert a.GetNumberOfTuples() == ntuples
                        values = to_numpy(a)
                        assert (values.dtype.kind, values.dtype.itemsize) == KINDS[kind]
                        equal(np, values, expected)
                a = fd.GetArray("never_written")
                assert a is None or a.GetNumberOfTuples() == 0
            assert not errors, filename
        return
    if case == "static":
        _, im, to_numpy = read(directory / "image_static.vtkhdf")
        assert im.GetExtent() == (2, 8, -3, 1, 4, 7)
        for index in [[2, -3, 4], [8, 1, 7]]:
            point = im.GetPoint(im.ComputePointId(index))
            equal(
                np,
                point,
                np.array([1, 2, 3])
                + np.array(DIRECTION) @ (np.array(index) * [0.5, 0.25, 2]),
            )
        equal(
            np,
            [
                [im.GetDirectionMatrix().GetElement(i, j) for j in range(3)]
                for i in range(3)
            ],
            DIRECTION,
        )
        for loc, attrs, dims in [
            ("cell", im.GetCellData(), [6, 4, 3]),
            ("point", im.GetPointData(), [7, 5, 4]),
        ]:
            for kind in KINDS:
                for name in ["scalar", "vector"]:
                    a = attrs.GetArray(f"{kind}-{name}")
                    assert a is not None
                    expected = encoded(
                        np,
                        [2, -3, 4],
                        dims,
                        modulo=True,
                        fraction=kind.startswith("real"),
                    )
                    if name == "vector":
                        expected = np.stack(
                            [
                                encoded(
                                    np,
                                    [2, -3, 4],
                                    dims,
                                    c,
                                    True,
                                    kind.startswith("real"),
                                )
                                for c in range(1, 4)
                            ],
                            axis=-1,
                        )
                    equal(
                        np,
                        to_numpy(a),
                        (
                            expected.reshape((-1, 3))
                            if name == "vector"
                            else expected.ravel()
                        ),
                    )
        fd = im.GetFieldData()
        for kind in KINDS:
            if kind == "int8":
                continue
            for name, expected in [
                ("value", [7]),
                ("tuples", [1, 2, 3]),
                ("vector", [[4, 5, 6]]),
                ("tensor", [[1, 2], [3, 4], [5, 6]]),
            ]:
                equal(np, to_numpy(fd.GetArray(f"{kind}-{name}")), expected)
    elif case == "numeric":
        for step in [0, 1]:
            _, im, to_numpy = read(directory / "image_numeric_temporal.vtkhdf", step)
            for attrs, dims in [
                (im.GetCellData(), [2, 1, 1]),
                (im.GetPointData(), [3, 2, 2]),
            ]:
                for kind in KINDS:
                    equal(
                        np,
                        to_numpy(attrs.GetArray(f"{kind}-scalar")),
                        encoded(np, [0, 0, 0], dims, modulo=True).ravel(),
                    )
                    expected = np.stack(
                        [encoded(np, [0, 0, 0], dims, c, True) for c in range(1, 3)],
                        axis=-1,
                    )
                    equal(
                        np,
                        to_numpy(attrs.GetArray(f"{kind}-vector")),
                        expected.reshape((-1, 2)),
                    )
    elif case == "temporal":
        for step in [0, 1, 2, 3, 4]:
            _, im, to_numpy = read(directory / "image_temporal.vtkhdf", step)
            cp = 2 if step >= 2 else 0
            op = 4 if step == 4 else 0
            vp = 3 if step >= 3 else 0
            equal(
                np,
                to_numpy(im.GetCellData().GetArray("temperature")),
                (encoded(np, [0, 0, 0], [6, 4, 3]) + cp * 1000).ravel(),
            )
            equal(
                np,
                to_numpy(im.GetCellData().GetArray("other")),
                (encoded(np, [0, 0, 0], [6, 4, 3]) + op * 1000).ravel(),
            )
            expected = encoded(np, [0, 0, 0], [7, 5, 4]) + vp * 1000
            expected = np.stack([expected + c * 10000 for c in range(1, 4)], axis=-1)
            equal(
                np,
                to_numpy(im.GetPointData().GetArray("velocity")),
                expected.reshape((-1, 3)),
            )
            for name in ["static_early", "static_late"]:
                equal(np, to_numpy(im.GetCellData().GetArray(name)), np.full(72, 17))
            energy = im.GetFieldData().GetArray("energy")
            history = im.GetFieldData().GetArray("history")
            if step == 0:
                assert energy is None or energy.GetNumberOfTuples() == 0
                assert history is None or history.GetNumberOfTuples() == 0
            else:
                equal(np, to_numpy(energy), [1 if step < 3 else 3])
                expected = [[10, 20]] if step < 3 else [[1, 2], [3, 4], [5, 6]]
                equal(np, to_numpy(history), expected)
    else:
        index = int(case[4:])
        dims = [
            [7, 5, 1],
            [1, 5, 4],
            [7, 1, 4],
            [7, 1, 1],
            [1, 5, 1],
            [1, 1, 4],
            [7, 5, 4],
        ][index - 1]
        lo = [-2, -3, -4] if index == 7 else [0, 0, 0]
        _, im, to_numpy = read(directory / f"image_dim_{index}.vtkhdf")
        assert im.GetNumberOfPoints() == np.prod(dims)
        for attrs, n in [
            (im.GetPointData(), dims),
            (im.GetCellData(), np.maximum(1, np.array(dims) - 1)),
        ]:
            expected = encoded(np, lo, n)
            equal(np, to_numpy(attrs.GetArray("scalar")), expected.ravel())
            if attrs.GetArray("noncontiguous") is not None:
                equal(
                    np,
                    to_numpy(attrs.GetArray("noncontiguous")),
                    (expected + 1000).ravel(),
                )
            equal(
                np,
                to_numpy(attrs.GetArray("vector")),
                np.stack([expected + 1000, expected + 2000], axis=-1).reshape((-1, 2)),
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("check", choices=["layout", "reader"])
    parser.add_argument("directory", type=Path)
    parser.add_argument("case", nargs="?", default="static")
    args = parser.parse_args()
    if args.check == "layout":
        check_layout(args.directory)
    else:
        check_reader(args.directory, args.case)


if __name__ == "__main__":
    main()
