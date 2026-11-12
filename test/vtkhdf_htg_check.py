#!/usr/bin/env python3
"""Independent layout and VTK reader checks for HyperTreeGrid fixtures.

The trees, their ownership, and the data follow the rules described in
vtkhdf_htg_test.F90. Expected values are derived here from the tree ids
stored in each part, so the checks apply to any number of MPI ranks.
"""
import sys
from pathlib import Path

try:
    import h5py
    import numpy as np
except ImportError as exc:
    print(exc)
    sys.exit(77)

mode, directory = sys.argv[1:3]
directory = Path(directory)
if mode == "reader":
    try:
        from vtkmodules.vtkCommonCore import vtkVersion
        from vtkmodules.vtkCommonDataModel import (
            vtkHyperTreeGrid,
            vtkHyperTreeGridNonOrientedGeometryCursor,
        )
        from vtkmodules.vtkIOHDF import vtkHDFReader
        from vtkmodules.util.numpy_support import vtk_to_numpy
    except ImportError as exc:
        print(exc)
        sys.exit(77)
    version = (vtkVersion.GetVTKMajorVersion(), vtkVersion.GetVTKMinorVersion(), vtkVersion.GetVTKBuildVersion())
    if version < (9, 7, 0):
        print("HyperTreeGrid reader checks require VTK 9.7 or later")
        sys.exit(77)

# name: dims, branch factor, topology phase per step, geometry phase per step,
# phases that supply a mask (None for no Mask dataset), static cell data,
# temporal cell data with the step whose values each step reads.
FIXTURES = {
    "htg_static": dict(dims=(4, 3, 2), bf=2, topology=None, geometry=None, masked=None,
                       static=["value", "vec", "ids", "big ids", "flag"], temporal={}, field=True),
    "htg_static_mask": dict(dims=(3, 1, 3), bf=3, topology=None, geometry=None, masked={0},
                            static=["value", "vec", "ids", "big ids", "flag"], temporal={}, field=True),
    "htg_fixed": dict(dims=(3, 4, 1), bf=2, topology=[0, 0, 0, 0], geometry=[0, 0, 0, 0], masked={0},
                      static=["ids"], temporal={"value": [0, 1, 1, 3], "vec": [0, 0, 0, 0]}, time=[0., 1., 2., 3.]),
    "htg_moving": dict(dims=(3, 3, 2), bf=2, topology=[0, 0, 0, 0], geometry=[0, 0, 1, 1], masked=None,
                       static=["ids"], temporal={"value": [0, 1, 2, 3]}, time=[0., 1., 2., 3.]),
    "htg_dynamic": dict(dims=(4, 3, 2), bf=2, topology=[0, 0, 1, 2], geometry=[0, 0, 0, 1], masked={1},
                        static=[], temporal={"value": [0, 1, 2, 3], "vec": [0, 0, 2, 3]}, time=[0., .25, .5, .75]),
}
DTYPES = {"value": "f8", "vec": "f4", "ids": "i4", "big ids": "i8", "flag": "u1"}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def axis(a, n, phase):
    i = np.arange(n, dtype=float)
    return (a + 1)*i + .25*i*i + .5*phase*(a + 1)


def pattern(g, phase, nchild):
    """Breadth-first refinement flags of tree g."""
    kind = (g + phase) % 3
    if kind == 0:
        return np.array([False])
    if kind == 1:
        return np.array([True] + [False]*nchild)
    flags = np.zeros(1 + 3*nchild, dtype=bool)
    flags[[0, 1, nchild]] = True
    return flags


def levels(flags, nchild):
    """Cell counts of each level implied by breadth-first flags."""
    counts, pos, n = [], 0, 1
    while n:
        counts.append(n)
        n, pos = int(flags[pos:pos+n].sum())*nchild, pos + n
    require(pos == len(flags), "inconsistent expected tree")
    return counts


def cell_values(name, g, n, step):
    lid = np.arange(n)
    if name == "value":
        return 1000.*g + lid + .5*step
    if name == "vec":
        return np.stack((np.full(n, g), lid, np.full(n, step)), axis=-1).astype("f4")
    if name == "flag":
        return lid % 2
    return 1000*g + lid


def expected_part(fx, tree_ids, tphase):
    """Expected topology datasets of a part owning tree_ids, in order."""
    nchild = fx["bf"]**sum(d > 1 for d in fx["dims"])
    flags = [pattern(g, tphase, nchild) for g in tree_ids]
    depths = [levels(f, nchild) for f in flags]
    bits = [f[:len(f)-d[-1]] for f, d in zip(flags, depths)]
    bits = np.concatenate(bits) if bits else np.zeros(0, bool)
    cells = np.concatenate(flags) if flags else np.zeros(0, bool)
    lid = np.concatenate([np.arange(len(f)) for f in flags]) if flags else np.zeros(0, int)
    return dict(
        flags=flags,
        depth_per_tree=[len(d) for d in depths],
        cells_per_depth=[n for d in depths for n in d],
        descriptor_bits=len(bits),
        descriptors=np.packbits(bits),
        ncell=len(cells),
        mask=np.packbits(~cells & (lid % 5 == 4)),
    )


def nbytes(bits):
    return (bits + 7)//8


def check_layout(name, fx):
    with h5py.File(directory/f"{name}.vtkhdf") as f:
        root = f["VTKHDF"]
        dims = fx["dims"]
        require(root.attrs["Type"] == b"HyperTreeGrid", f"{name}: Type")
        require(tuple(root.attrs["Version"]) == (2, 8), f"{name}: Version")
        require(tuple(root.attrs["Dimensions"]) == dims, f"{name}: Dimensions")
        require(root.attrs["BranchFactor"] == fx["bf"], f"{name}: BranchFactor")
        require(root.attrs["TransposedRootIndexing"] == 0, f"{name}: TransposedRootIndexing")
        require("PointData" not in root, f"{name}: PointData must be absent")
        for dataset in ("TreeIds", "DepthPerTree", "NumberOfCellsPerTreeDepth", "NumberOfTrees",
                        "NumberOfDepths", "NumberOfCells", "DescriptorsSize"):
            require(root[dataset].dtype == np.int64 and root[dataset].ndim == 1, f"{name}: {dataset} type")
        require(root["Descriptors"].dtype == np.uint8, f"{name}: Descriptors type")
        for a in "XYZ":
            require(root[f"{a}Coordinates"].dtype == np.float64, f"{name}: {a}Coordinates type")
        has_mask = fx["masked"] is not None
        require(("Mask" in root) == has_mask, f"{name}: Mask presence")

        temporal = fx["topology"] is not None
        require(("Steps" in root) == temporal, f"{name}: Steps presence")
        if temporal:
            steps = root["Steps"]
            nsteps = len(fx["topology"])
            require(steps.attrs["NSteps"] == nsteps, f"{name}: NSteps")
            require(np.array_equal(steps["Values"][:], fx["time"]), f"{name}: Values")
            offsets = {key: steps[key][:] for key in steps if key.endswith("Offsets") and key not in
                       ("CellDataOffsets", "FieldDataOffsets")}
            for key, values in offsets.items():
                require(values.shape == (nsteps,) and values.dtype == np.int64, f"{name}: {key} shape")
            nparts = steps["NumberOfParts"][:]
            topology, geometry = fx["topology"], fx["geometry"]
        else:
            nsteps = 1
            offsets = {}
            nparts = [len(root["NumberOfTrees"])]
            topology, geometry = [0], [0]

        def offset(key, step):
            return int(offsets[key][step]) if temporal else 0

        ntrees_total = int(np.prod([max(d - 1, 1) for d in dims]))
        sizes = dict(TreeIds=0, NumberOfCellsPerTreeDepth=0, Descriptors=0, Mask=0, Parts=0, X=0, Y=0, Z=0)
        data_sizes = {key: 0 for key in fx["temporal"]}
        for step in range(nsteps):
            npart = int(nparts[step])
            require(npart >= 1, f"{name}: step {step} has no parts")
            part0 = offset("PartOffsets", step)
            ntree = root["NumberOfTrees"][part0:part0+npart]
            ncell = root["NumberOfCells"][part0:part0+npart]
            ndepth = root["NumberOfDepths"][part0:part0+npart]
            nbits = root["DescriptorsSize"][part0:part0+npart]
            require(len(ntree) == npart, f"{name}: step {step} part metadata")
            require(npart == int(nparts[0]), f"{name}: step {step} part count changed")
            tree0 = offset("TreeIdsOffsets", step)
            require(offset("DepthPerTreeOffsets", step) == tree0, f"{name}: DepthPerTreeOffsets")
            all_ids = root["TreeIds"][tree0:tree0+ntree.sum()]
            require(sorted(all_ids) == list(range(ntrees_total)), f"{name}: step {step} tree ids")
            depth0 = offset("NumberOfCellsPerTreeDepthOffsets", step)
            desc0 = offset("DescriptorsOffsets", step)
            mask0 = offset("MaskOffsets", step)
            cell0 = 0
            for part in range(npart):
                ids = root["TreeIds"][tree0:tree0+ntree[part]]
                want = expected_part(fx, ids, topology[step])
                where = f"{name}: step {step} part {part}"
                require(ncell[part] == want["ncell"], f"{where} NumberOfCells")
                require(list(root["DepthPerTree"][tree0:tree0+ntree[part]]) == want["depth_per_tree"],
                        f"{where} DepthPerTree")
                require(ndepth[part] == len(want["cells_per_depth"]), f"{where} NumberOfDepths")
                require(list(root["NumberOfCellsPerTreeDepth"][depth0:depth0+ndepth[part]]) ==
                        want["cells_per_depth"], f"{where} NumberOfCellsPerTreeDepth")
                require(nbits[part] == want["descriptor_bits"], f"{where} DescriptorsSize")
                require(np.array_equal(root["Descriptors"][desc0:desc0+nbytes(nbits[part])], want["descriptors"]),
                        f"{where} Descriptors")
                if has_mask:
                    mask = want["mask"] if topology[step] in fx["masked"] else np.zeros(nbytes(ncell[part]), "u1")
                    require(np.array_equal(root["Mask"][mask0:mask0+nbytes(ncell[part])], mask), f"{where} Mask")
                for a, axis_name in enumerate("XYZ"):
                    start = offset(f"{axis_name}CoordinatesOffsets", step) + part*dims[a]
                    require(np.array_equal(root[f"{axis_name}Coordinates"][start:start+dims[a]],
                                           axis(a, dims[a], geometry[step])), f"{where} {axis_name}Coordinates")
                for key, source in fx["temporal"].items():
                    start = int(root[f"Steps/CellDataOffsets/{key}"][step]) + cell0
                    got = root[f"CellData/{key}"][start:start+ncell[part]]
                    lid = 0
                    for g, flags in zip(ids, want["flags"]):
                        require(np.array_equal(got[lid:lid+len(flags)], cell_values(key, g, len(flags), source[step])),
                                f"{where} CellData/{key}")
                        lid += len(flags)
                for key in fx["static"]:
                    got = root[f"CellData/{key.replace(' ', '_')}"][cell0:cell0+ncell[part]]
                    lid = 0
                    for g, flags in zip(ids, want["flags"]):
                        require(np.array_equal(got[lid:lid+len(flags)], cell_values(key, g, len(flags), 0)),
                                f"{where} CellData/{key}")
                        lid += len(flags)
                tree0 += ntree[part]
                depth0 += ndepth[part]
                desc0 += nbytes(nbits[part])
                mask0 += nbytes(ncell[part])
                cell0 += ncell[part]

            # Storage is appended only for a new topology, geometry, or data.
            new_topology = step == 0 or topology[step] != topology[step-1]
            new_parts = step == 0 or npart != int(nparts[step-1])
            new_geometry = step == 0 or geometry[step] != geometry[step-1] or new_parts
            if new_topology:
                require(part0 == sizes["Parts"], f"{name}: step {step} PartOffsets")
                require(offset("TreeIdsOffsets", step) == sizes["TreeIds"], f"{name}: step {step} TreeIdsOffsets")
                sizes["Parts"] += npart
                sizes["TreeIds"] += ntree.sum()
                sizes["NumberOfCellsPerTreeDepth"] += ndepth.sum()
                sizes["Descriptors"] += sum(nbytes(n) for n in nbits)
                sizes["Mask"] += sum(nbytes(n) for n in ncell)
            elif temporal:
                for key in ("PartOffsets", "TreeIdsOffsets", "DescriptorsOffsets", "MaskOffsets",
                            "NumberOfCellsPerTreeDepthOffsets"):
                    require(offsets[key][step] == offsets[key][step-1], f"{name}: step {step} {key} reuse")
            for a, axis_name in enumerate("XYZ"):
                if new_geometry:
                    require(offset(f"{axis_name}CoordinatesOffsets", step) == sizes[axis_name],
                            f"{name}: step {step} {axis_name}CoordinatesOffsets")
                    sizes[axis_name] += npart*dims[a]
                elif temporal:
                    key = f"{axis_name}CoordinatesOffsets"
                    require(offsets[key][step] == offsets[key][step-1], f"{name}: step {step} {key} reuse")
            for key, source in fx["temporal"].items():
                if source[step] == step:
                    data_sizes[key] += ncell.sum()

        require(len(root["NumberOfTrees"]) == sizes["Parts"], f"{name}: NumberOfTrees size")
        require(len(root["NumberOfCells"]) == sizes["Parts"], f"{name}: NumberOfCells size")
        require(len(root["NumberOfDepths"]) == sizes["Parts"], f"{name}: NumberOfDepths size")
        require(len(root["DescriptorsSize"]) == sizes["Parts"], f"{name}: DescriptorsSize size")
        require(len(root["TreeIds"]) == sizes["TreeIds"], f"{name}: TreeIds size")
        require(len(root["DepthPerTree"]) == sizes["TreeIds"], f"{name}: DepthPerTree size")
        require(len(root["NumberOfCellsPerTreeDepth"]) == sizes["NumberOfCellsPerTreeDepth"],
                f"{name}: NumberOfCellsPerTreeDepth size")
        require(len(root["Descriptors"]) == sizes["Descriptors"], f"{name}: Descriptors size")
        if has_mask:
            require(len(root["Mask"]) == sizes["Mask"], f"{name}: Mask size")
        for a in "XYZ":
            require(len(root[f"{a}Coordinates"]) == sizes[a], f"{name}: {a}Coordinates size")
        for key, size in data_sizes.items():
            dataset = root[f"CellData/{key}"]
            require(dataset.shape[0] == size and dataset.dtype == np.dtype(DTYPES[key]), f"{name}: CellData/{key}")
        for key in fx["static"]:
            dataset = root[f"CellData/{key.replace(' ', '_')}"]
            require(dataset.dtype == np.dtype(DTYPES[key]), f"{name}: CellData/{key} type")
            require(dataset.attrs["Name"] == key.encode(), f"{name}: CellData/{key} Name")
        value = root["CellData/value"]
        require(value.attrs.get("Attribute") == (None if name in ("htg_moving", "htg_dynamic") else b"Scalars"),
                f"{name}: value Attribute")


def check_coarse_layout():
    with h5py.File(directory/"htg_coarse.vtkhdf") as f:
        root = f["VTKHDF"]
        require(tuple(root.attrs["Dimensions"]) == (5, 1, 1), "htg_coarse: Dimensions")
        require(root["Descriptors"].shape == (1,), "htg_coarse: Descriptors holds one unused byte")
        require((root["DescriptorsSize"][:] == 0).all(), "htg_coarse: DescriptorsSize")
        require(sorted(root["TreeIds"][:]) == [0, 1, 2, 3], "htg_coarse: TreeIds")
        require((root["DepthPerTree"][:] == 1).all(), "htg_coarse: DepthPerTree")
        require((root["NumberOfCellsPerTreeDepth"][:] == 1).all(), "htg_coarse: NumberOfCellsPerTreeDepth")
        require(np.array_equal(root["CellData/value"][:], 1000.*root["TreeIds"][:]), "htg_coarse: value")


def pieces(output):
    if isinstance(output, vtkHyperTreeGrid):
        return [output]
    return [output.GetPartitionAsDataObject(i) for i in range(output.GetNumberOfPartitions())]


def cell_bounds(root_bounds, dims, bf, flags):
    """Bounds of every cell of a tree, in breadth-first order."""
    active = [a for a in range(3) if dims[a] > 1]
    nchild = bf**len(active)
    bounds = [np.array(root_bounds, dtype=float)]
    for cell, refined in enumerate(flags):
        parent = bounds[cell]
        if not refined:
            continue
        for child in range(nchild):
            b = parent.copy()
            c = child
            for a in active:
                lo, hi = parent[2*a], parent[2*a+1]
                index, c = c % bf, c // bf
                b[2*a] = lo + (hi - lo)*index/bf
                b[2*a+1] = lo + (hi - lo)*(index + 1)/bf
            bounds.append(b)
    return bounds


def check_tree(htg, where, g, flags, dims, bf, coords, masked):
    """Compare a tree's structure, geometry, and mask; return its global cell indices."""
    nx, ny = max(dims[0] - 1, 1), max(dims[1] - 1, 1)
    ijk = (g % nx, (g // nx) % ny, g // (nx*ny))
    root_bounds = []
    for a in range(3):
        lo = coords[a][ijk[a]]
        hi = coords[a][ijk[a]+1] if dims[a] > 1 else lo
        root_bounds += [lo, hi]
    want_bounds = cell_bounds(root_bounds, dims, bf, flags)
    nchild = bf**sum(d > 1 for d in dims)

    cursor = vtkHyperTreeGridNonOrientedGeometryCursor()
    htg.InitializeNonOrientedGeometryCursor(cursor, int(g), False)
    require(cursor.HasTree(), f"{where} tree {g} missing")
    start = cursor.GetGlobalNodeIndex()
    seen = {}

    def visit():
        local = cursor.GetGlobalNodeIndex() - start
        bounds = [0.]*6
        cursor.GetBounds(bounds)
        seen[local] = (cursor.IsLeaf(), bounds, cursor.IsMasked())
        if not cursor.IsLeaf():
            require(cursor.GetNumberOfChildren() == nchild, f"{where} tree {g} child count")
            for child in range(nchild):
                cursor.ToChild(child)
                visit()
                cursor.ToParent()

    visit()
    require(sorted(seen) == list(range(len(flags))), f"{where} tree {g} cell indices")
    for local, refined in enumerate(flags):
        leaf, bounds, is_masked = seen[local]
        require(leaf == (not refined), f"{where} tree {g} cell {local} refinement")
        require(np.allclose(bounds, want_bounds[local], rtol=0, atol=1e-12), f"{where} tree {g} cell {local} bounds")
        want_masked = masked and not refined and local % 5 == 4
        require(bool(is_masked) == bool(want_masked), f"{where} tree {g} cell {local} mask")
    return start + np.arange(len(flags))


def check_reader(name, fx):
    path = directory/f"{name}.vtkhdf"
    with h5py.File(path) as f:
        root = f["VTKHDF"]
        temporal = "Steps" in root
        nsteps = int(root["Steps"].attrs["NSteps"]) if temporal else 1
        part_ids = []
        for step in range(nsteps):
            npart = int(root["Steps/NumberOfParts"][step]) if temporal else len(root["NumberOfTrees"])
            part0 = int(root["Steps/PartOffsets"][step]) if temporal else 0
            tree0 = int(root["Steps/TreeIdsOffsets"][step]) if temporal else 0
            ids = []
            for n in root["NumberOfTrees"][part0:part0+npart]:
                ids.append(root["TreeIds"][tree0:tree0+n])
                tree0 += n
            part_ids.append(ids)

    reader = vtkHDFReader()
    reader.SetFileName(str(path))
    reader.UpdateInformation()
    require(reader.GetNumberOfSteps() == nsteps, f"{name}: reader step count")
    dims, bf = fx["dims"], fx["bf"]
    nchild = bf**sum(d > 1 for d in dims)
    topology = fx["topology"] or [0]
    geometry = fx["geometry"] or [0]
    for step in list(range(nsteps)) + list(range(nsteps - 1, -1, -1)):
        reader.SetStep(step)
        reader.Update()
        # Field data belongs to the output as a whole, a partitioned dataset if there are several parts.
        field_data = reader.GetOutputDataObject(0).GetFieldData()
        output = pieces(reader.GetOutputDataObject(0))
        require(len(output) == len(part_ids[step]), f"{name}: step {step} piece count")
        for part, (htg, ids) in enumerate(zip(output, part_ids[step])):
            where = f"{name}: step {step} part {part}"
            require(isinstance(htg, vtkHyperTreeGrid), f"{where} is not a HyperTreeGrid")
            require(tuple(htg.GetDimensions()) == dims, f"{where} dimensions")
            require(htg.GetBranchFactor() == bf, f"{where} branch factor")
            coords = [vtk_to_numpy(c) for c in (htg.GetXCoordinates(), htg.GetYCoordinates(), htg.GetZCoordinates())]
            for a in range(3):
                require(np.array_equal(coords[a], axis(a, dims[a], geometry[step])), f"{where} coordinates {a}")
            flags = [pattern(g, topology[step], nchild) for g in ids]
            require(htg.GetNumberOfCells() == sum(len(f) for f in flags), f"{where} cell count")
            masked = fx["masked"] is not None and topology[step] in fx["masked"]
            require((htg.GetMask() is not None) == (fx["masked"] is not None), f"{where} mask presence")
            cell_data = htg.GetCellData()
            names = sorted(cell_data.GetArrayName(i) for i in range(cell_data.GetNumberOfArrays()))
            # The reader labels HyperTreeGrid arrays by their HDF5 dataset names.
            keys = {key.replace(" ", "_"): key for key in fx["static"] + list(fx["temporal"])}
            require(names == sorted(keys), f"{where} cell arrays {names}")
            for g, tree_flags in zip(ids, flags):
                index = check_tree(htg, where, g, tree_flags, dims, bf, coords, masked)
                for array_name, key in keys.items():
                    source = fx["temporal"][key][step] if key in fx["temporal"] else 0
                    got = vtk_to_numpy(cell_data.GetArray(array_name))[index]
                    require(np.array_equal(got, cell_values(key, g, len(tree_flags), source)),
                            f"{where} tree {g} {key}")
            if name not in ("htg_moving", "htg_dynamic"):
                scalars = cell_data.GetScalars()
                require(scalars is not None and scalars.GetName() == "value", f"{where} active scalars")
        if fx.get("field"):
            field = field_data.GetArray("field-value")
            require(field is not None and field.GetValue(0) == 42.0, f"{name}: field data")
        if "time" in fx and name == "htg_fixed":
            field = field_data.GetArray("time")
            require(field is not None and field.GetValue(0) == fx["time"][step], f"{name}: step {step} time field")


def check_coarse_reader():
    reader = vtkHDFReader()
    reader.SetFileName(str(directory/"htg_coarse.vtkhdf"))
    reader.Update()
    total = 0
    for htg in pieces(reader.GetOutputDataObject(0)):
        require(htg.GetMask() is None, "htg_coarse: unexpected mask")
        values = vtk_to_numpy(htg.GetCellData().GetArray("value"))
        cursor = vtkHyperTreeGridNonOrientedGeometryCursor()
        for g in range(4):
            htg.InitializeNonOrientedGeometryCursor(cursor, g, False)
            if not cursor.HasTree():
                continue
            require(cursor.IsLeaf(), f"htg_coarse: tree {g} is refined")
            require(values[cursor.GetGlobalNodeIndex()] == 1000.*g, f"htg_coarse: tree {g} value")
            total += 1
    require(total == 4, "htg_coarse: tree count")


for fixture_name, fixture in FIXTURES.items():
    if mode == "layout":
        check_layout(fixture_name, fixture)
    else:
        check_reader(fixture_name, fixture)
if mode == "layout":
    check_coarse_layout()
else:
    check_coarse_reader()
