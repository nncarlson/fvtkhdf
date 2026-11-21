"""Record VTK-only hierarchy readback separately from writer acceptance.

Fixtures have independent HDF5 conformance checks. Each reader runs in a child
process so a reader error, crash, or timeout becomes a compatibility finding.
"""
import json
from pathlib import Path
import subprocess
import sys

from vtkhdf_assembly_cases import trees


def probe(path, prefix, case, profile, nproc):
    import numpy as np
    from vtkmodules.vtkCommonDataModel import vtkCompositeDataSet
    from vtkmodules.vtkIOHDF import vtkHDFReader
    from vtkmodules.util.numpy_support import vtk_to_numpy
    from vtkhdf_reader_smoke import iter_leaf_datasets

    temporal = case in ("explicit", "implicit")
    reader = vtkHDFReader()
    errors = []
    reader.AddObserver("ErrorEvent", lambda *args: errors.append("reader ErrorEvent"))
    reader.SetFileName(path)
    reader.UpdateInformation()
    if profile == "selected":
        for selection in (reader.GetPointDataArraySelection(), reader.GetCellDataArraySelection()):
            selection.DisableAllArrays()
            selection.EnableAllArrays()

    def check_block(block, slot, step):
        leaves = list(iter_leaf_datasets(block))
        image = temporal and slot == 1
        kind = "vtkImageData" if image else "vtkUnstructuredGrid"
        assert leaves and all(leaf.IsA(kind) for leaf in leaves), (slot, kind)
        assert sum(leaf.GetNumberOfCells() for leaf in leaves) == nproc, slot
        if not image:
            assert sum(leaf.GetNumberOfPoints() for leaf in leaves) == 3 * nproc, slot
            actual = np.concatenate([vtk_to_numpy(leaf.GetCellData().GetArray("constant")) for leaf in leaves])
            np.testing.assert_array_equal(np.sort(actual.ravel()), slot + 1 + np.arange(nproc))
        if temporal and slot in (0, 1):
            name, offset = ("pressure", 100) if slot == 0 else ("temperature", 300)
            actual = np.concatenate([vtk_to_numpy(leaf.GetCellData().GetArray(name)) for leaf in leaves])
            np.testing.assert_array_equal(np.sort(actual.ravel()), offset + step + np.arange(nproc))

    def check_mb(output, expected, step):
        assert output.IsA("vtkMultiBlockDataSet"), output.GetClassName()
        assert output.GetNumberOfBlocks() == len(expected), list(expected)
        for i, (name, child) in enumerate(expected.items()):
            assert output.GetMetaData(i).Get(vtkCompositeDataSet.NAME()) == name, name
            data = output.GetBlock(i)
            assert data is not None, name
            if isinstance(child, int):
                check_block(data, child, step)
            else:
                check_mb(data, child, step)

    def check_assembly(assembly, node, expected):
        # VTK's PDC Assembly associates block indices directly with ordinary
        # nodes; soft-link leaves do not become separate vtkDataAssembly nodes.
        indices = {value for value in expected.values() if isinstance(value, int)}
        assert set(assembly.GetDataSetIndices(node, False)) == indices, node
        groups = {name: value for name, value in expected.items() if isinstance(value, dict)}
        children = assembly.GetChildNodes(node, False)
        assert [assembly.GetNodeName(child) for child in children] == list(groups), node
        for child, expected_child in zip(children, groups.values()):
            check_assembly(assembly, child, expected_child)

    for step in ([0, 1, 0] if temporal else [0]):
        reader.SetStep(step)
        reader.Update()
        assert not errors, errors
        assert reader.GetNumberOfSteps() == (2 if temporal else 1), "temporal state changed"
        output = reader.GetOutputDataObject(0)
        if prefix == "mb":
            check_mb(output, trees[case], step)
        else:
            assert output.IsA("vtkPartitionedDataSetCollection"), output.GetClassName()
            count = 4 if temporal else 0 if case == "empty" else 5 if case == "rootgrowth" else 2
            assert output.GetNumberOfPartitionedDataSets() == count
            for slot in range(count):
                check_block(output.GetPartitionedDataSet(slot), slot, step)
            assembly = output.GetDataAssembly()
            assert assembly is not None, "missing data assembly"
            check_assembly(assembly, 0, trees[case])
    print("hierarchy, geometry, and cell values verified")


if __name__ == "__main__":
    if sys.argv[1] == "--worker":
        probe(*sys.argv[2:6], int(sys.argv[6]))
    else:
        try:
            import numpy
            from vtkmodules.vtkIOHDF import vtkHDFReader
            from vtkmodules.vtkCommonCore import vtkVersion
        except ImportError as exc:
            print(exc)
            sys.exit(77)
        directory, prefix, nproc = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
        results = {"vtk": vtkVersion.GetVTKVersion(), "profiles": {}}
        for case in trees:
            path = directory / f"{prefix}_assembly_{case}.vtkhdf"
            if not path.is_file():
                sys.exit(f"missing fixture: {path}")
            for profile in (("unchanged", "selected") if case in ("explicit", "implicit") else ("unchanged",)):
                command = [sys.executable, __file__, "--worker", str(path), prefix, case, profile, nproc]
                try:
                    run = subprocess.run(command, capture_output=True, text=True, timeout=30)
                    result = {"status": "pass" if run.returncode == 0 else "reader failure",
                              "returncode": run.returncode, "stdout": run.stdout, "stderr": run.stderr}
                except subprocess.TimeoutExpired:
                    result = {"status": "reader timeout"}
                results["profiles"][f"{case}-{profile}"] = result
                print(f"{case}-{profile}: {result['status']}")
        report = directory / f"{prefix}-assembly-reader-compatibility.json"
        report.write_text(json.dumps(results, indent=2) + "\n")
        print(f"Compatibility findings (not writer acceptance): {report}")
