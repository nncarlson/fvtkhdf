"""Independently check shared Assembly semantics and unchanged block output."""
from pathlib import Path
import sys

try:
    import h5py
    import numpy as np
except ImportError as exc:
    print(exc)
    sys.exit(77)


def check_tree(group, expected):
    assert list(group) == list(expected), (group.name, list(group), list(expected))
    assert group.id.get_create_plist().get_link_creation_order() == 3, group.name
    assert not group.attrs, (group.name, dict(group.attrs))
    for name, child in expected.items():
        link = group.get(name, getlink=True)
        if isinstance(child, int):
            assert isinstance(link, h5py.SoftLink), (group.name, name)
            assert link.path == f"/VTKHDF/Group_{child}", (group.name, name, link.path)
        else:
            assert isinstance(link, h5py.HardLink), (group.name, name)
            check_tree(group[name], child)


directory, prefix = sys.argv[1:]
kind = b"MultiBlockDataSet" if prefix == "mb" else b"PartitionedDataSetCollection"
from vtkhdf_assembly_cases import trees

for case, tree in trees.items():
    path = Path(directory) / f"{prefix}_assembly_{case}.vtkhdf"
    with h5py.File(path) as file:
        root = file["VTKHDF"]
        assert root.attrs["Type"] == kind
        np.testing.assert_array_equal(root.attrs["Version"], [2, 8])
        assembly = root["Assembly"]
        check_tree(assembly, tree)
        nblocks = (4 if case in ("explicit", "implicit") else 0 if case == "empty"
                   else 5 if case == "rootgrowth" else 2)
        assert list(root) == ["Assembly"] + [f"Group_{i}" for i in range(nblocks)]
        for i in range(nblocks):
            block = root[f"Group_{i}"]
            assert int(np.asarray(block.attrs["Index"]).flat[0]) == i
            np.testing.assert_array_equal(block["FieldData/id"][...].ravel(), [i + 1])
        if case == "flat":
            # A rewrite would give these new creation-order values (2, 3).
            assert [assembly.id.links.get_info(name.encode()).corder
                    for name in ("A", "B")] == [0, 1]
            np.testing.assert_array_equal(root["Group_0/FieldData/after"][...].ravel(), [42])
        if case in ("explicit", "implicit"):
            fluid, solid = root["Group_0"], root["Group_1"]
            assert fluid.attrs["Type"] == b"UnstructuredGrid"
            assert solid.attrs["Type"] == b"ImageData"
            assert fluid.attrs["Name"] == b"Fluid.Name"
            assert solid.attrs["Name"] == b"Fluid Name"
            # Part sizes, mesh arrays, and values identify every MPI contribution.
            parts = fluid["NumberOfCells"][...].ravel()
            nproc = len(parts)  # fixed-mesh topology is stored once and reused
            np.testing.assert_array_equal(parts, np.ones(nproc))
            xyz = fluid["Points"][...].reshape(nproc, 3, 3)
            for rank in range(nproc):
                np.testing.assert_array_equal(xyz[rank], [[rank, 0, 0], [rank + 1, 0, 0], [rank, 1, 0]])
            np.testing.assert_array_equal(fluid["FieldData/after"][...].ravel(), [42])
            np.testing.assert_array_equal(fluid["CellData/constant"][...].ravel(), 1 + np.arange(nproc))
            np.testing.assert_array_equal(fluid["CellData/pressure"][...].reshape(2, nproc),
                                          [[100 + r for r in range(nproc)], [101 + r for r in range(nproc)]])
            np.testing.assert_array_equal(solid["CellData/temperature"][...].reshape(2, nproc),
                                          [[300 + r for r in range(nproc)], [301 + r for r in range(nproc)]])
            for i in range(4):
                steps = root[f"Group_{i}/Steps"]
                assert int(np.asarray(steps.attrs["NSteps"]).flat[0]) == 2
                np.testing.assert_array_equal(steps["Values"], [0, 1])
        print(f"{prefix} {case}: hierarchy and block output verified")
