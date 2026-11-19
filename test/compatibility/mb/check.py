"""Independent expectations, HDF5 checks, and a pvpython reader worker."""
import argparse
import json
from pathlib import Path
import traceback

import numpy as np

TIMES = [0.0, 0.2, 0.7, 1.5, 3.0]
VTK_TYPES = dict(image="vtkImageData", rectilinear="vtkRectilinearGrid",
                 structured="vtkStructuredGrid", unstructured="vtkUnstructuredGrid")
ASSOCIATIONS = {"point": "PointData", "cell": "CellData", "field": "FieldData"}


def epoch(mode, step):
    return step // 2 if mode in (2, 3) else 0


def part_counts(mode, step, block, nproc):
    phase = int(mode == 3 and epoch(mode, step) == 1)
    return [0 if nproc > 1 and rank == (block - 1) % nproc else rank + 1 + phase
            for rank in range(nproc)]


def dimensions(block, ndim, nproc):
    dims = [3, 4 if ndim > 1 else 1, 2 if ndim > 2 else 1]
    if nproc > 1:
        dims[0] = sum(part_counts(0, 0, block, nproc)) + 1
    return dims


def counts(kind, mode, step, block=1, ndim=3, nproc=1):
    if kind == "unstructured":
        nc = sum(part_counts(mode, step, block, nproc))
        return 4 * nc, nc
    dims = dimensions(block, ndim, nproc)
    return int(np.prod(dims)), int(np.prod(np.maximum(1, np.array(dims) - 1)))


def coordinates(kind, mode, block, step, ndim=3, nproc=1):
    phase = epoch(mode, step)
    if kind == "unstructured":
        tetra = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
        return np.concatenate([tetra + [10 * block + 2 * c, 0, 0.25 * phase]
                               for c in range(counts(kind, mode, step, block, ndim, nproc)[1])])
    result = []
    nx, ny, nz = dimensions(block, ndim, nproc)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if kind == "image":
                    xyz = [10 * block + 0.5 * (i - 2), j, 2 + 2 * (k + 1)]
                elif kind == "rectilinear":
                    xyz = [10 * block + 0.5 * i + 0.1 * i * i, j + 0.1 * j * j,
                           2 * k + 0.1 * k * k + 0.25 * phase]
                else:
                    xyz = [10 * block + 0.5 * i + 0.05 * j * k,
                           j + 0.1 * i * k, 2 * k + 0.25 * phase]
                result.append(xyz)
    return np.array(result)


def arrays(kind, mode, block, step, association, ndim=3, nproc=1):
    if association == "field":
        result = {"block_id": np.array([block], dtype="int32"),
                  "direction": np.array([[1.0, 2.0, 3.0]])}
        if mode and mode != "fixed-constant":
            nt = 1 + step % 3
            result["clock"] = np.array([TIMES[step]])
            result["table"] = (1000 * block + 100 * step + np.arange(1, 2 * nt + 1)).reshape(nt, 2).astype(float)
        return result
    n = counts(kind, mode, step, block, ndim, nproc)[0 if association == "point" else 1]
    base = 1000 * block + np.arange(1, n + 1, dtype=float)
    result = {}
    if mode != 3:
        result[association + "_constant"] = base
        result[association + "_constant_vector"] = base[:, None] + [1, 2, 3]
    if mode and mode != "fixed-constant":
        phases = ([0, 1, 2, 3, 4] if mode == 3 else [0, 1, 1, 3, 3]) if association == "point" else [0, 0, 2, 2, 4]
        result[association + "_value"] = base + 100 * phases[step]
        result[association + "_vector"] = base[:, None] + 100 * phases[step] + [1, 2, 3]
    return result


def layout(case):
    import h5py

    ndim, nproc = case.get("dimensions", 3), case.get("writer_ranks", 1)
    with h5py.File(case["file"]) as file:
        root = file["VTKHDF"]
        if case["composite"]:
            assert root.attrs["Type"] == b"MultiBlockDataSet"
            assert "Steps" not in root  # Baseline files follow the composite specification.
            assert set(root["Assembly"]) == {f"block_{i}" for i in range(1, len(case["blocks"]) + 1)}
        for b, (kind, mode) in enumerate(case["blocks"], 1):
            dims = dimensions(b, ndim, nproc)
            g = root[f"Group_{b - 1}"] if case["composite"] else root
            assert g.attrs["Type"] == VTK_TYPES[kind][3:].encode()
            if case["composite"]:
                assert g.attrs["Index"] == b - 1
                assert g.attrs["Name"] == f"block_{b}".encode()
                link = root["Assembly"].get(f"block_{b}", getlink=True)
                assert isinstance(link, h5py.SoftLink) and link.path == f"/VTKHDF/Group_{b - 1}"
            assert ("Steps" in g) == bool(mode)
            if mode:
                assert g["Steps"].attrs["NSteps"] == 5
                np.testing.assert_array_equal(g["Steps/Values"], TIMES)
            for step in range(5 if mode else 1):
                want_xyz = coordinates(kind, mode, b, step, ndim, nproc)
                if kind == "image":
                    extent = np.array([[-2, 3, 1], np.array([-2, 3, 1]) + dims - 1]).T.ravel()
                    np.testing.assert_array_equal(g.attrs["WholeExtent"], extent)
                    np.testing.assert_array_equal(g.attrs["Origin"], [10 * b, -3, 2])
                    np.testing.assert_array_equal(g.attrs["Spacing"], [0.5, 1, 2])
                elif kind == "rectilinear":
                    axes = []
                    for name, length in zip(("XCoordinates", "YCoordinates", "ZCoordinates"), dims):
                        offset = int(g[f"Steps/{name}Offsets"][step]) if mode == 2 else 0
                        axes.append(np.asarray(g[name])[offset:offset + length])
                    xyz = np.array([[x, y, z] for z in axes[2] for y in axes[1] for x in axes[0]])
                    np.testing.assert_allclose(xyz, want_xyz, rtol=0, atol=1e-12)
                else:
                    data = g["Points"]
                    if kind == "unstructured":
                        offset = int(g["Steps/PointOffsets"][step]) if mode else 0
                        xyz = data[offset:offset + len(want_xyz)]
                    else:
                        xyz = data[epoch(mode, step)] if mode == 2 else data[...]
                    np.testing.assert_allclose(np.asarray(xyz).reshape(-1, 3), want_xyz, rtol=0, atol=1e-12)
                for association, group in ASSOCIATIONS.items():
                    expected = arrays(kind, mode, b, step, association, ndim, nproc)
                    assert set(g[group]) == set(expected)
                    for name, want in expected.items():
                        data = g[f"{group}/{name}"]
                        offset_path = f"Steps/{group}Offsets/{name}"
                        offset = int(g[offset_path][step]) if offset_path in g else 0
                        if association == "field":
                            actual = data[offset:offset + len(want)] if mode else data[...]
                            if mode:
                                shape = [want.shape[1] if want.ndim == 2 else 1, len(want)]
                                np.testing.assert_array_equal(g[f"Steps/FieldDataSizes/{name}"][step], shape)
                        elif kind == "unstructured":
                            actual = data[offset:offset + len(want)]
                        else:
                            actual = data[offset] if mode else data[...]
                        actual = np.asarray(actual).reshape(want.shape)
                        np.testing.assert_array_equal(actual, want)
                        assert actual.dtype == want.dtype
            if kind == "unstructured":
                parts = [n for step in ([0, 2, 4] if mode == 3 else [0])
                         for n in part_counts(mode, step, b, nproc) if n]
                np.testing.assert_array_equal(g["Connectivity"], np.concatenate([np.arange(4 * n) for n in parts]))
                np.testing.assert_array_equal(g["Offsets"], np.concatenate([4 * np.arange(n + 1) for n in parts]))
                np.testing.assert_array_equal(g["Types"], np.full(sum(parts), 10))
                np.testing.assert_array_equal(g["NumberOfCells"], parts)
                np.testing.assert_array_equal(g["NumberOfPoints"], np.array(parts) * 4)
    return {"status": "pass"}


def reader(case, cache="default", native_only=False):
    paraview_version = None
    if not native_only:
        import paraview
        from paraview import servermanager
        from paraview.simple import Delete, VTKHDFReader
        paraview_version = paraview.__version__
    from vtkmodules.vtkIOHDF import vtkHDFReader
    from vtkmodules.vtkCommonCore import vtkStringOutputWindow, vtkOutputWindow, vtkVersion
    from vtkmodules.vtkCommonExecutionModel import vtkStreamingDemandDrivenPipeline
    from vtkmodules.util.numpy_support import vtk_to_numpy

    result = {"paraview": paraview_version, "vtk": vtkVersion.GetVTKVersion(), "profiles": {}}
    temporal = any(mode for _, mode in case["blocks"])
    ndim, nproc = case.get("dimensions", 3), case.get("writer_ranks", 1)
    sequence = [0, 1, 2, 3, 4, 2, 0, 4] if temporal else [0]
    profiles = (("unchanged",) if native_only else
                ("unchanged", "geometry", "point", "cell", "field", "all", "late", "fetch"))
    for profile in profiles:
        failures = []
        checked = 0
        proxy = None
        capture = vtkStringOutputWindow()
        vtkOutputWindow.SetInstance(capture)

        def verify(label, actual, expected):
            nonlocal checked
            checked += 1
            try:
                assert np.shape(actual) == np.shape(expected), (np.shape(actual), np.shape(expected))
                np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)
            except (AssertionError, TypeError, ValueError) as exc:
                failures.append({"check": label, "error": str(exc)})

        try:
            available = {}
            if profile == "unchanged":
                native = vtkHDFReader()
                native.SetFileName(case["file"])
                if cache != "default":
                    native.SetUseCache(cache == "on")
                native.UpdateInformation()
                initial_mtime = native.GetMTime()
                verify("initial step count", native.GetNumberOfSteps(), len(TIMES) if temporal else 1)
                verify("initial temporal flag", bool(native.GetHasTemporalData()), temporal)
                info = native.GetOutputInformation(0)
                key = vtkStreamingDemandDrivenPipeline.TIME_STEPS()
                verify("time values", [info.Get(key, i) for i in range(info.Length(key))], TIMES if temporal else [])
                # Selection changes mark the reader modified and can hide the
                # composite temporal-state reset between information and data.
                # Keep the first read on a fresh, otherwise unchanged reader.
                steps = [0]
            else:
                proxy = VTKHDFReader(FileName=[case["file"]])
                if cache != "default":
                    proxy.UseCache = int(cache == "on")
                proxy.UpdatePipelineInformation()
                result["reader_options"] = {name: str(getattr(proxy, name)) for name in ("UseCache", "PieceDistribution")
                                            if name in proxy.ListProperties()}
                verify("time values", list(proxy.TimestepValues), TIMES if temporal else [])
                for association in ASSOCIATIONS:
                    prop = association.capitalize() + "ArrayStatus"
                    available[association] = list(getattr(proxy, prop).Available)
                    setattr(proxy, prop, available[association] if profile in (association, "all", "late", "fetch") else [])
                steps = ([3] if temporal else [0]) if profile == "late" else sequence
            for step in steps:
                if profile == "unchanged":
                    native.SetStep(step)
                    verify("reader unmodified before first data request", native.GetMTime(), initial_mtime)
                    native.Update()
                    verify("updated step count", native.GetNumberOfSteps(), len(TIMES) if temporal else 1)
                    verify("updated temporal flag", bool(native.GetHasTemporalData()), temporal)
                else:
                    proxy.UpdatePipeline(time=TIMES[step])
                # In this serial builtin session, inspect the server's actual
                # reader output. Fetch has its own transfer/append path and is
                # checked separately so its defects are not blamed on the file.
                if profile == "unchanged":
                    output = native.GetOutputDataObject(0)
                elif profile == "fetch":
                    output = servermanager.Fetch(proxy)
                else:
                    output = proxy.GetClientSideObject().GetOutputDataObject(0)
                if case["composite"]:
                    assert (output.IsA("vtkMultiBlockDataSet") or
                            output.IsA("vtkPartitionedDataSetCollection")), output.GetClassName()
                    size = (output.GetNumberOfBlocks() if output.IsA("vtkMultiBlockDataSet") else
                            output.GetNumberOfPartitionedDataSets())
                    assert size == len(case["blocks"])
                for b, (kind, mode) in enumerate(case["blocks"], 1):
                    data = output
                    if case["composite"]:
                        data = (output.GetBlock(b - 1) if output.IsA("vtkMultiBlockDataSet") else
                                output.GetPartitionedDataSet(b - 1))
                    label = f"step {step}, block {b}"
                    assert data is not None, label + ": missing block"
                    parts = ([data.GetPartition(i) for i in range(data.GetNumberOfPartitions())
                              if data.GetPartition(i) is not None] if data.IsA("vtkPartitionedDataSet") else [data])
                    assert parts and all(p.IsA(VTK_TYPES[kind]) for p in parts), label + ": wrong block type"
                    local_step = step if mode else 0
                    npnt, ncell = counts(kind, mode, local_step, b, ndim, nproc)
                    verify(label + " point/cell counts",
                           [sum(p.GetNumberOfPoints() for p in parts), sum(p.GetNumberOfCells() for p in parts)], [npnt, ncell])
                    xyz = np.array([p.GetPoint(i) for p in parts for i in range(p.GetNumberOfPoints())])
                    verify(label + " coordinates", xyz, coordinates(kind, mode, b, local_step, ndim, nproc))
                    if kind == "unstructured":
                        types = [p.GetCellType(i) for p in parts for i in range(p.GetNumberOfCells())]
                        verify(label + " cell types", types, [10] * ncell)
                        for index, part in enumerate(parts):
                            for c in range(part.GetNumberOfCells()):
                                cell = part.GetCell(c)
                                ids = [cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]
                                verify(label + f" part {index}, connectivity {c}", ids, np.arange(4 * c, 4 * c + 4))
                    for association, group in ASSOCIATIONS.items():
                        if profile not in (association, "unchanged", "all", "late", "fetch"):
                            continue
                        for name, want in arrays(kind, mode, b, local_step, association, ndim, nproc).items():
                            field = data.GetFieldData().GetArray(name) if association == "field" else None
                            actual_arrays = ([field] if field is not None else
                                             [getattr(p, "Get" + group)().GetArray(name) for p in parts])
                            if any(a is None for a in actual_arrays):
                                failures.append({"check": label + " " + group + "/" + name, "error": "missing array"})
                                continue
                            # Field data describes the logical dataset. It may
                            # reside on its container or be replicated on parts;
                            # distributed point/cell arrays instead concatenate.
                            values = ([vtk_to_numpy(a) for a in actual_arrays] if association == "field" else
                                      [np.concatenate([vtk_to_numpy(a) for a in actual_arrays])])
                            for actual in values:
                                verify(label + " " + name, actual, want)
                                if actual.dtype != want.dtype:
                                    failures.append({"check": label + " " + name, "error": "wrong numeric type"})
                            for array in actual_arrays:
                                verify(label + " " + name + " components", array.GetNumberOfComponents(),
                                       want.shape[1] if want.ndim == 2 else 1)
            log = capture.GetOutput()
            if "ERR|" in log or "ERROR:" in log:
                failures.append({"check": "reader diagnostics", "error": log})
            result["profiles"][profile] = {"status": "fail" if failures else "pass", "checks": checked,
                                            "failures": failures, "available": available}
        except Exception:
            result["profiles"][profile] = {"status": "fail", "checks": checked, "failures": failures,
                                            "error": traceback.format_exc(), "log": capture.GetOutput()}
        finally:
            if proxy is not None:
                Delete(proxy)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=("layout", "reader", "vtk"))
    parser.add_argument("case", type=Path)
    parser.add_argument("result", type=Path)
    parser.add_argument("--cache", choices=("default", "on", "off"), default="default")
    args = parser.parse_args()
    try:
        case = json.loads(args.case.read_text())
        value = (layout(case) if args.operation == "layout" else
                 reader(case, args.cache, native_only=args.operation == "vtk"))
    except Exception:
        value = {"status": "fail", "error": traceback.format_exc()}
    args.result.write_text(json.dumps(value, indent=2) + "\n")
    # Reader incompatibilities are results, not a failure of the harness itself.
    # RUN checks worker exit codes separately from the recorded observations.
