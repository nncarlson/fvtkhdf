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


def counts(kind, mode, step):
    if kind == "unstructured":
        nc = 2 if mode == 3 and epoch(mode, step) == 1 else 1
        return 4 * nc, nc
    return 24, 6


def coordinates(kind, mode, block, step):
    phase = epoch(mode, step)
    if kind == "unstructured":
        tetra = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
        return np.concatenate([tetra + [10 * block + 2 * c, 0, 0.25 * phase]
                               for c in range(counts(kind, mode, step)[1])])
    result = []
    for k in range(2):
        for j in range(4):
            for i in range(3):
                if kind == "image":
                    xyz = [10 * block + 0.5 * (i - 2), j, 2 + 2 * (k + 1)]
                elif kind == "rectilinear":
                    xyz = [10 * block + [0, 0.6, 1.4][i], [0, 1.1, 2.4, 3.9][j],
                           [0, 2.1][k] + 0.25 * phase]
                else:
                    xyz = [10 * block + 0.5 * i + 0.05 * j * k,
                           j + 0.1 * i * k, 2 * k + 0.25 * phase]
                result.append(xyz)
    return np.array(result)


def arrays(kind, mode, block, step, association):
    if association == "field":
        result = {"block_id": np.array([block], dtype="int32"),
                  "direction": np.array([[1.0, 2.0, 3.0]])}
        if mode and mode != "fixed-constant":
            nt = 1 + step % 3
            result["clock"] = np.array([TIMES[step]])
            result["table"] = (1000 * block + 100 * step + np.arange(1, 2 * nt + 1)).reshape(nt, 2).astype(float)
        return result
    n = counts(kind, mode, step)[0 if association == "point" else 1]
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

    with h5py.File(case["file"]) as file:
        root = file["VTKHDF"]
        if case["composite"]:
            assert root.attrs["Type"] == b"MultiBlockDataSet"
            assert "Steps" not in root  # Baseline files follow the composite specification.
            assert set(root["Assembly"]) == {f"block_{i}" for i in range(1, len(case["blocks"]) + 1)}
        for b, (kind, mode) in enumerate(case["blocks"], 1):
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
                want_xyz = coordinates(kind, mode, b, step)
                if kind == "image":
                    np.testing.assert_array_equal(g.attrs["WholeExtent"], [-2, 0, 3, 6, 1, 2])
                    np.testing.assert_array_equal(g.attrs["Origin"], [10 * b, -3, 2])
                    np.testing.assert_array_equal(g.attrs["Spacing"], [0.5, 1, 2])
                elif kind == "rectilinear":
                    axes = []
                    for name, length in zip(("XCoordinates", "YCoordinates", "ZCoordinates"), (3, 4, 2)):
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
                    expected = arrays(kind, mode, b, step, association)
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
                snapshots = [1, 2, 1] if mode == 3 else [1]
                np.testing.assert_array_equal(g["Connectivity"], np.concatenate([np.arange(4 * n) for n in snapshots]))
                np.testing.assert_array_equal(g["Offsets"], np.concatenate([4 * np.arange(n + 1) for n in snapshots]))
                np.testing.assert_array_equal(g["Types"], np.full(sum(snapshots), 10))
    return {"status": "pass"}


def reader(case):
    import paraview
    from paraview import servermanager
    from paraview.simple import Delete, VTKHDFReader
    from vtkmodules.vtkCommonCore import vtkStringOutputWindow, vtkOutputWindow, vtkVersion
    from vtkmodules.util.numpy_support import vtk_to_numpy

    result = {"paraview": paraview.__version__, "vtk": vtkVersion.GetVTKVersion(), "profiles": {}}
    temporal = any(mode for _, mode in case["blocks"])
    sequence = [0, 1, 2, 3, 4, 2, 0, 4] if temporal else [0]
    for profile in ("geometry", "point", "cell", "field", "all", "late", "fetch"):
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
            proxy = VTKHDFReader(FileName=[case["file"]])
            proxy.UpdatePipelineInformation()
            result["reader_options"] = {name: str(getattr(proxy, name)) for name in ("UseCache", "PieceDistribution")
                                        if name in proxy.ListProperties()}
            verify("time values", list(proxy.TimestepValues), TIMES if temporal else [])
            available = {}
            for association in ASSOCIATIONS:
                prop = association.capitalize() + "ArrayStatus"
                available[association] = list(getattr(proxy, prop).Available)
                setattr(proxy, prop, available[association] if profile in (association, "all", "late", "fetch") else [])
            for step in ([3] if temporal else [0]) if profile == "late" else sequence:
                proxy.UpdatePipeline(time=TIMES[step])
                # In this serial builtin session, inspect the server's actual
                # reader output. Fetch has its own transfer/append path and is
                # checked separately so its defects are not blamed on the file.
                if profile == "fetch":
                    output = servermanager.Fetch(proxy)
                else:
                    output = proxy.GetClientSideObject().GetOutputDataObject(0)
                if case["composite"]:
                    assert output.IsA("vtkMultiBlockDataSet"), output.GetClassName()
                    assert output.GetNumberOfBlocks() == len(case["blocks"])
                for b, (kind, mode) in enumerate(case["blocks"], 1):
                    data = output.GetBlock(b - 1) if case["composite"] else output
                    if data is not None and data.IsA("vtkPartitionedDataSet"):
                        assert data.GetNumberOfPartitions() == 1
                        data = data.GetPartition(0)
                    label = f"step {step}, block {b}"
                    assert data is not None, label + ": missing block"
                    assert data.IsA(VTK_TYPES[kind]), (label, data.GetClassName())
                    local_step = step if mode else 0
                    npnt, ncell = counts(kind, mode, local_step)
                    verify(label + " point/cell counts", [data.GetNumberOfPoints(), data.GetNumberOfCells()], [npnt, ncell])
                    xyz = np.array([data.GetPoint(i) for i in range(data.GetNumberOfPoints())])
                    verify(label + " coordinates", xyz, coordinates(kind, mode, b, local_step))
                    if kind == "unstructured":
                        types = [data.GetCellType(i) for i in range(data.GetNumberOfCells())]
                        verify(label + " cell types", types, [10] * ncell)
                        for c in range(min(ncell, data.GetNumberOfCells())):
                            cell = data.GetCell(c)
                            ids = [cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]
                            verify(label + f" connectivity {c}", ids, np.arange(4 * c, 4 * c + 4))
                    for association, group in ASSOCIATIONS.items():
                        if profile not in (association, "all", "late", "fetch"):
                            continue
                        attributes = getattr(data, "Get" + group)()
                        for name, want in arrays(kind, mode, b, local_step, association).items():
                            array = attributes.GetArray(name)
                            if array is None:
                                failures.append({"check": label + " " + group + "/" + name, "error": "missing array"})
                                continue
                            verify(label + " " + name, vtk_to_numpy(array), want)
                            verify(label + " " + name + " components", array.GetNumberOfComponents(),
                                   want.shape[1] if want.ndim == 2 else 1)
                            if vtk_to_numpy(array).dtype != want.dtype:
                                failures.append({"check": label + " " + name, "error": "wrong numeric type"})
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
    parser.add_argument("operation", choices=("layout", "reader"))
    parser.add_argument("case", type=Path)
    parser.add_argument("result", type=Path)
    args = parser.parse_args()
    try:
        value = globals()[args.operation](json.loads(args.case.read_text()))
    except Exception:
        value = {"status": "fail", "error": traceback.format_exc()}
    args.result.write_text(json.dumps(value, indent=2) + "\n")
    # Reader incompatibilities are results, not a failure of the harness itself.
    # RUN checks worker exit codes separately from the recorded observations.
