#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path

SKIP_CODE = 77


def import_vtk_reader():
    errors = []

    try:
        from vtkmodules.vtkIOHDF import vtkHDFReader  # type: ignore

        return vtkHDFReader, None
    except Exception as exc:  # pragma: no cover - depends on local install
        errors.append(f"vtkmodules.vtkIOHDF import failed: {exc}")

    try:
        from vtk import vtkHDFReader  # type: ignore

        return vtkHDFReader, None
    except Exception as exc:  # pragma: no cover - depends on local install
        errors.append(f"vtk import failed: {exc}")

    return None, "; ".join(errors)


def fail(message: str) -> int:
    print(message, file=sys.stderr)
    return 1


def iter_leaf_datasets(data_object):
    if data_object is None:
        return

    if data_object.IsA("vtkCompositeDataSet"):
        iterator = data_object.NewIterator()
        iterator.VisitOnlyLeavesOn()
        iterator.SkipEmptyNodesOn()
        iterator.InitTraversal()
        while not iterator.IsDoneWithTraversal():
            leaf = iterator.GetCurrentDataObject()
            if leaf is not None:
                yield leaf
            iterator.GoToNextItem()
        return

    yield data_object


def summarize_output(data_object):
    class_name = data_object.GetClassName()
    leaf_classes = []
    non_empty_leaves = 0
    total_points = 0
    total_cells = 0

    for leaf in iter_leaf_datasets(data_object):
        leaf_classes.append(leaf.GetClassName())

        npoints = leaf.GetNumberOfPoints() if hasattr(leaf, "GetNumberOfPoints") else 0
        ncells = leaf.GetNumberOfCells() if hasattr(leaf, "GetNumberOfCells") else 0
        total_points += npoints
        total_cells += ncells

        if npoints > 0 or ncells > 0:
            non_empty_leaves += 1

    return {
        "class_name": class_name,
        "leaf_classes": leaf_classes,
        "non_empty_leaves": non_empty_leaves,
        "points": total_points,
        "cells": total_cells,
    }


def validate_temporal_read(reader):
    if not hasattr(reader, "GetNumberOfSteps") or not hasattr(reader, "SetStep"):
        return None

    steps = reader.GetNumberOfSteps()
    if steps < 2:
        raise RuntimeError(f"expected at least 2 time steps, got {steps}")

    return [0, 1]


def read_output(reader, step):
    if step is not None:
        reader.SetStep(step)
    reader.Update()
    output = reader.GetOutputDataObject(0)
    if output is None:
        output = reader.GetOutput()
    if output is None:
        raise RuntimeError("reader produced no output object")
    return output


def validate_ug_output(summary):
    if summary["points"] <= 0:
        raise RuntimeError("UG reader returned no points")
    if summary["cells"] <= 0:
        raise RuntimeError("UG reader returned no cells")
    if "vtkUnstructuredGrid" not in summary["leaf_classes"]:
        raise RuntimeError(
            "UG reader did not yield an unstructured-grid leaf "
            f"(root={summary['class_name']}, leaves={summary['leaf_classes']})"
        )


def validate_mb_output(output, summary):
    if not output.IsA("vtkCompositeDataSet"):
        raise RuntimeError(f"MB reader returned non-composite root {summary['class_name']}")
    if summary["non_empty_leaves"] < 2:
        raise RuntimeError(
            "MB reader returned fewer than 2 non-empty leaves "
            f"(root={summary['class_name']}, leaves={summary['leaf_classes']})"
        )
    if summary["points"] <= 0:
        raise RuntimeError("MB reader returned no points")
    if summary["cells"] <= 0:
        raise RuntimeError("MB reader returned no cells")


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Smoke-test VTK HDF reader availability and file readback.")
    parser.add_argument("kind", nargs="?", choices=("ug", "mb"))
    parser.add_argument("path", nargs="?")
    parser.add_argument("--check-import", action="store_true")
    args = parser.parse_args(argv)

    vtk_hdf_reader, error_message = import_vtk_reader()

    if args.check_import:
        if vtk_hdf_reader is None:
            return SKIP_CODE
        return 0

    if vtk_hdf_reader is None:
        print(error_message or "vtkHDFReader is not importable", file=sys.stderr)
        return SKIP_CODE
    if args.kind is None or args.path is None:
        return fail("expected dataset kind and file path")

    path = Path(args.path)
    if not path.is_file():
        return fail(f"missing input file: {path}")

    reader = vtk_hdf_reader()
    if hasattr(reader, "CanReadFile") and not reader.CanReadFile(str(path)):
        return fail(f"vtkHDFReader rejected input file: {path}")

    reader.SetFileName(str(path))
    if hasattr(reader, "UpdateInformation"):
        reader.UpdateInformation()

    try:
        steps_to_read = validate_temporal_read(reader)
        if steps_to_read is None:
            steps_to_read = [None]

        for step in steps_to_read:
            output = read_output(reader, step)
            summary = summarize_output(output)
            if args.kind == "ug":
                validate_ug_output(summary)
            else:
                validate_mb_output(output, summary)
    except Exception as exc:
        return fail(f"{path.name}: {exc}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
