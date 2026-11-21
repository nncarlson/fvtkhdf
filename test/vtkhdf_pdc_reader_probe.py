"""Record VTK reader compatibility separately from format conformance tests.

Only the VTK Python module is required. Each request sequence runs in a child
process because known reader defects can terminate Python. Reader failures are
findings, not failures of the writer's independently checked VTKHDF 2.8 output.
"""
import json
from pathlib import Path
import subprocess
import sys


def probe(directory, case, profile):
    import numpy as np
    from vtkmodules.vtkCommonCore import vtkVersion
    from vtkmodules.vtkIOHDF import vtkHDFReader
    from vtkmodules.util.numpy_support import vtk_to_numpy
    from vtkhdf_reader_smoke import iter_leaf_datasets

    reader = vtkHDFReader()
    errors = []
    reader.AddObserver("ErrorEvent", lambda *args: errors.append("reader ErrorEvent"))
    reader.SetFileName(str(Path(directory) / f"pdc_collection_{case}.vtkhdf"))
    reader.UpdateInformation()
    temporal = case != "static"
    expected_steps = 3 if temporal else 1
    assert reader.GetNumberOfSteps() == expected_steps, "initial step count"
    if profile != "unchanged":
        # Exercise the selection-modifying path without requiring ParaView.
        for selection in (reader.GetPointDataArraySelection(), reader.GetCellDataArraySelection(),
                          reader.GetFieldDataArraySelection()):
            selection.DisableAllArrays()
            selection.EnableAllArrays()
    kinds = ("ImageData", "RectilinearGrid", "StructuredGrid", "UnstructuredGrid", "HyperTreeGrid")
    for step in ([0, 2, 1, 0] if temporal else [0]):
        reader.SetStep(step)
        reader.Update()
        assert not errors, errors
        assert reader.GetNumberOfSteps() == expected_steps, "step count changed after read"
        output = reader.GetOutputDataObject(0)
        assert output.IsA("vtkPartitionedDataSetCollection"), output.GetClassName()
        assert output.GetNumberOfPartitionedDataSets() == 5, "block count"
        for i, kind in enumerate(kinds):
            block = output.GetPartitionedDataSet(i)
            leaves = list(iter_leaf_datasets(block))
            assert leaves and all(leaf.IsA("vtk" + kind) for leaf in leaves), kind
            assert sum(leaf.GetNumberOfCells() for leaf in leaves) == 1, kind + " cell count"
            if kind != "HyperTreeGrid":
                assert sum(leaf.GetNumberOfPoints() for leaf in leaves) == 4, kind + " point count"
            else:
                values = np.concatenate([vtk_to_numpy(p.GetCellData().GetArray("value")) for p in leaves])
                np.testing.assert_array_equal(values, [5.])
            # Record missing field arrays too, whether the reader places them on
            # the logical dataset or its parts. Layout tests verify their bytes.
            if profile != "fields":
                continue
            field = block.GetFieldData().GetArray("id")
            if field is None:
                field = leaves[0].GetFieldData().GetArray("id")
            assert field is not None, kind + " missing field data id"
            np.testing.assert_array_equal(vtk_to_numpy(field).ravel(), [i + 1])
    print(json.dumps({"vtk": vtkVersion.GetVTKVersion(), "status": "pass"}))


if __name__ == "__main__":
    if len(sys.argv) == 5 and sys.argv[1] == "--worker":
        probe(*sys.argv[2:])
    else:
        try:
            import numpy
            from vtkmodules.vtkIOHDF import vtkHDFReader
            from vtkmodules.vtkCommonCore import vtkVersion
        except ImportError as exc:
            print(exc)
            sys.exit(77)
        directory = Path(sys.argv[1])
        for case in ("static", "temporal", "mixed"):
            if not (directory / f"pdc_collection_{case}.vtkhdf").is_file():
                sys.exit(f"missing fixture: pdc_collection_{case}.vtkhdf")
        results = {"vtk": vtkVersion.GetVTKVersion(), "profiles": {}}
        for case in ("static", "temporal", "mixed"):
            for profile in ("unchanged", "selected", "fields"):
                command = [sys.executable, __file__, "--worker", str(directory), case, profile]
                try:
                    run = subprocess.run(command, capture_output=True, text=True, timeout=60)
                    result = {"status": "pass" if run.returncode == 0 else "reader failure",
                              "returncode": run.returncode, "stdout": run.stdout, "stderr": run.stderr}
                except subprocess.TimeoutExpired:
                    result = {"status": "reader timeout"}
                results["profiles"][f"{case}-{profile}"] = result
                print(f"{case}-{profile}: {result['status']}")
        report = directory / "pdc-reader-compatibility.json"
        report.write_text(json.dumps(results, indent=2) + "\n")
        print(f"Compatibility findings (not writer acceptance): {report}")
