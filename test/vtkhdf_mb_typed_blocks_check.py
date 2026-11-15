#!/usr/bin/env python3
"""Independent layout and VTK readback checks for the typed UG block API."""

import argparse
from pathlib import Path

try:
    import h5py
    import numpy as np
except ImportError as exc:
    print(exc)
    raise SystemExit(77)


def ranks(nproc):
    return range(nproc - 1 if nproc > 1 else 1)


def counts(nproc, grow=False):
    return np.array([rank + 1 + int(grow) for rank in ranks(nproc)])


def coordinates(nproc, shift=0, grow=False):
    tetra = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    return np.concatenate(
        [tetra + [10 * rank + 2 * j, 0, shift]
         for rank, count in zip(ranks(nproc), counts(nproc, grow))
         for j in range(count)]
    )


def values(nproc, mode, step):
    phase = min(step, 1) if mode != 1 else 0
    return np.concatenate(
        [np.full(count, 100 * mode + 10 * phase + rank)
         for rank, count in zip(ranks(nproc), counts(nproc, mode == 3 and step > 0))]
    )


def nproc_from(root):
    return int(root["Group_0/FieldData/nproc"][...].flat[0])


def layout(directory):
    with h5py.File(directory / "mb_typed_static.vtkhdf") as file:
        root = file["VTKHDF"]
        nproc = nproc_from(root)
        assert len(root["Assembly"]) == 24
        for i in range(24):
            block = root[f"Group_{i}"]
            assert "Steps" not in block
            link = root["Assembly"].get(f"Block_{i + 1}", getlink=True)
            assert isinstance(link, h5py.SoftLink)
            assert link.path == f"/VTKHDF/Group_{i}"
            np.testing.assert_array_equal(block["Points"], coordinates(nproc))
            np.testing.assert_array_equal(block["CellData/value"], values(nproc, i + 1, 0))
        np.testing.assert_array_equal(root["Group_0/FieldData/alias"], [42])

    with h5py.File(directory / "mb_typed_modes.vtkhdf") as file:
        root = file["VTKHDF"]
        nproc = nproc_from(root)
        initial_cells = int(counts(nproc).sum())
        for mode in range(1, 4):
            block = root[f"Group_{mode - 1}"]
            steps = block["Steps"]
            assert int(np.asarray(steps.attrs["NSteps"]).flat[0]) == 3
            np.testing.assert_array_equal(steps["Values"], [0, 1, 2])
            np.testing.assert_array_equal(steps["NumberOfParts"], [len(ranks(nproc))] * 3)
            initial = coordinates(nproc)
            expected = initial if mode == 1 else np.concatenate(
                [initial, coordinates(nproc, shift=1, grow=mode == 3)]
            )
            np.testing.assert_array_equal(block["Points"], expected)
            offsets = [0, 0, 0] if mode == 1 else [0, len(initial), len(initial)]
            np.testing.assert_array_equal(steps["PointOffsets"], offsets)
            offsets = [0, initial_cells, initial_cells] if mode == 3 else [0, 0, 0]
            np.testing.assert_array_equal(steps["CellOffsets"], offsets)
            cells = counts(nproc)
            if mode == 3:
                cells = np.concatenate([cells, counts(nproc, grow=True)])
            np.testing.assert_array_equal(block["NumberOfCells"], cells)
            np.testing.assert_array_equal(block["NumberOfPoints"], 4 * cells)
            payload = values(nproc, mode, 0)
            if mode != 1:
                payload = np.concatenate([payload, values(nproc, mode, 1)])
            np.testing.assert_array_equal(block["CellData/value"], payload)
            offsets = [0, 0, 0] if mode == 1 else [0, initial_cells, initial_cells]
            np.testing.assert_array_equal(steps["CellDataOffsets/value"], offsets)
            np.testing.assert_array_equal(block["FieldData/clock"], [[mode], [20 + mode]])
            np.testing.assert_array_equal(steps["FieldDataOffsets/clock"], [0, 0, 1])
            for name, expected in (("constant", 70 + mode), ("late", 90 + mode), ("nproc", nproc)):
                np.testing.assert_array_equal(block[f"FieldData/{name}"], [[expected]])
                np.testing.assert_array_equal(steps[f"FieldDataOffsets/{name}"], [0, 0, 0])
            for name in ("clock", "constant", "late", "nproc"):
                np.testing.assert_array_equal(steps[f"FieldDataSizes/{name}"], np.ones((3, 2)))


def reader(directory):
    try:
        from vtkmodules.vtkIOHDF import vtkHDFReader
        from vtkmodules.util.numpy_support import vtk_to_numpy
    except ImportError as exc:
        print(exc)
        raise SystemExit(77)
    from vtkhdf_reader_smoke import iter_leaf_datasets

    for kind in ("static", "modes"):
        path = directory / f"mb_typed_{kind}.vtkhdf"
        with h5py.File(path) as file:
            nproc = nproc_from(file["VTKHDF"])
        vtk_reader = vtkHDFReader()
        errors = []
        vtk_reader.AddObserver("ErrorEvent", lambda *args: errors.append(args))
        vtk_reader.SetFileName(str(path))
        vtk_reader.UpdateInformation()
        assert not errors, errors
        if kind == "modes":
            assert vtk_reader.GetNumberOfSteps() == 3
        for step in range(3 if kind == "modes" else 1):
            if kind == "modes":
                vtk_reader.SetStep(step)
            vtk_reader.Update()
            assert not errors, errors
            output = vtk_reader.GetOutputDataObject(0)
            assert output.GetNumberOfBlocks() == (3 if kind == "modes" else 24)
            for i in range(output.GetNumberOfBlocks()):
                mode = i + 1
                leaves = list(iter_leaf_datasets(output.GetBlock(i)))
                actual_points = np.concatenate([vtk_to_numpy(leaf.GetPoints().GetData()) for leaf in leaves])
                actual_values = np.concatenate([vtk_to_numpy(leaf.GetCellData().GetArray("value")) for leaf in leaves])
                shift = min(step, 1) if kind == "modes" and mode != 1 else 0
                grow = kind == "modes" and mode == 3 and step > 0
                np.testing.assert_array_equal(actual_points, coordinates(nproc, shift, grow))
                np.testing.assert_array_equal(actual_values, values(nproc, mode, step if kind == "modes" else 0))
                assert sum(leaf.GetNumberOfCells() for leaf in leaves) == len(actual_values)
                assert all(leaf.GetCellType(j) == 10 for leaf in leaves for j in range(leaf.GetNumberOfCells()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("check", choices=("layout", "reader"))
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    (layout if args.check == "layout" else reader)(args.directory)
