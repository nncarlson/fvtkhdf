# ParaView compatibility including HyperTreeGrid

Checked 2026-11-19 using `283ef55` plus the uncommitted MultiBlock HyperTreeGrid
integration. Serial files used GCC 16.1.1; four-rank files used NAG 7.2 and
MPICH 4.3.2. Both installations used HDF5 1.14.6.

The ParaView 6.2.1 (VTK 9.7.2) run covers 120 distinct files: 49 baseline,
48 reduced-dimensional, and 23 MPI-written fixtures. All passed independent
HDF5 checks before readback. Eighteen baseline files were reread with caching
enabled, plus six temporal MPI-written HTG controls with caching enabled.
All reading used a serial builtin session.

- HTG geometry and cell data pass in all four modes, as standalone output and
  MB blocks, in 1D/2D/3D, and with unequal MPI contributions and empty ranks.
  Cursor checks verify every root, refined parent and child, cell bounds,
  masks, and values. Dynamic ownership changes retain one stored part per rank.
- Five-type collections with homogeneous static or temporal timelines pass
  geometry and point/cell checks, including different temporal modes and both
  block orders. All five MB block types lose field data in the reader.
- Mixed static/temporal HTG collections terminate the reader with a floating-
  point exception in both block orders. Diagnostics first report missing
  `Steps/TreeIdsOffsets` in the static block. Completed `unchanged` observations
  are preserved; termination occurs in the `geometry` profile, so later
  profiles are unobserved. Both baseline driver runs exit 1 for these workers;
  other fixtures still run and retain their observations.
- Promoting the static HTG block to fixed-mesh mode restores geometry and cell
  readback without registering temporal data or rewriting unchanged values.
  This does not restore MB field data.
- With `UseCache=1`, serial standalone temporal HTG files lose constant field
  arrays after the initial read. A fresh reader starting at a later step works;
  four-rank standalone HTG controls retain their fields with caching enabled.
  This matches [VTK #20251](https://gitlab.kitware.com/vtk/vtk/-/work_items/20251)
  and [fVTKHDF #67](https://github.com/nncarlson/fvtkhdf/issues/67).
- The existing unchanged-reader temporal-state reset also affects HTG blocks.
  The visible structured-array failure remains distinct from ParaView's
  selection-modifying profiles. See
  [VTK #20246](https://gitlab.kitware.com/vtk/vtk/-/work_items/20246) and
  [fVTKHDF #68](https://github.com/nncarlson/fvtkhdf/issues/68).
- The standalone RectilinearGrid Fetch-only failure remains observable; no
  additional HTG Fetch failure was found beyond field loss or a mixed timeline.

Seventeen existing baseline files containing HTG were also read with ParaView
6.2.0 (VTK 9.7.1). HTG-only homogeneous-timeline cases have the same geometry,
cell, and field results, and mixed HTG timelines terminate in the same way.
The older reader additionally fails structured arrays in heterogeneous
collections. Two diagnostic controls confirm this is not introduced by HTG:
a four-type static file without HTG still fails, whereas standalone
RectilinearGrid arrays read correctly (its Fetch-only failure remains).
These older-reader observations do not change the ParaView 6.2.1 findings.

`not-applicable` means the HTG-only case has no point data. The All, Fresh later
step, and Fetch profiles include field checks: their failures in homogeneous
MB files do not imply geometry or cell-data failures. `worker failed` means
that profile was not completed after native process termination. Detailed
JSON and logs are generated outside the repository by the harness.

## Baseline, ParaView 6.2.1

Writer: 283ef55-plus-hypertree-mb
pvpython: /opt/paraview/6.2.1/bin/pvpython
ParaView: 6.2.1; VTK: 9.7.2
Cache: default

| Case | Unchanged reader | Geometry | Point | Cell | Field | All | Fresh later step | Fetch |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| standalone-image-static | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-image-static | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-image-fixed | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-image-fixed | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-rectilinear-static | pass | pass | pass | pass | pass | pass | pass | fail |
| mb-rectilinear-static | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-rectilinear-fixed | pass | pass | pass | pass | pass | pass | pass | fail |
| mb-rectilinear-fixed | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-rectilinear-moving | pass | pass | pass | pass | pass | pass | pass | fail |
| mb-rectilinear-moving | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-structured-static | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-structured-static | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-structured-fixed | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-structured-fixed | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-structured-moving | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-structured-moving | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-unstructured-static | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-unstructured-static | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-unstructured-fixed | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-unstructured-fixed | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-unstructured-moving | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-unstructured-moving | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-unstructured-dynamic | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-unstructured-dynamic | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-hypertree-static | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-static | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-fixed | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-fixed | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-moving | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-moving | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-dynamic | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-dynamic | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| all-static | fail | pass | pass | pass | fail | fail | fail | fail |
| all-fixed | fail | pass | pass | pass | fail | fail | fail | fail |
| different-modes | fail | pass | pass | pass | fail | fail | fail | fail |
| different-modes-reversed | fail | pass | pass | pass | fail | fail | fail | fail |
| mixed-modes | fail | pass | fail | fail | fail | fail | fail | fail |
| mixed-modes-reversed | fail | fail | fail | fail | fail | fail | fail | fail |
| unstructured-mixed | fail | fail | fail | fail | fail | fail | fail | fail |
| unstructured-mixed-reversed | fail | fail | fail | fail | fail | fail | fail | fail |
| rectilinear-mixed | fail | pass | fail | fail | fail | fail | fail | fail |
| rectilinear-mixed-reversed | fail | fail | fail | fail | fail | fail | fail | fail |
| structured-mixed | fail | pass | fail | fail | fail | fail | fail | fail |
| structured-mixed-reversed | fail | fail | fail | fail | fail | fail | fail | fail |
| hypertree-mixed | fail | worker failed | worker failed | worker failed | worker failed | worker failed | worker failed | worker failed |
| hypertree-mixed-reversed | fail | worker failed | worker failed | worker failed | worker failed | worker failed | worker failed | worker failed |
| unstructured-workaround | fail | pass | pass | pass | fail | fail | fail | fail |
| rectilinear-workaround | fail | pass | pass | pass | fail | fail | fail | fail |
| hypertree-workaround | fail | pass | not-applicable | pass | fail | fail | fail | fail |

## Reduced dimensions, ParaView 6.2.1

Writer: 283ef55-plus-hypertree-mb
pvpython: /opt/paraview/6.2.1/bin/pvpython
ParaView: 6.2.1; VTK: 9.7.2
Cache: default

| Case | Unchanged reader | Geometry | Point | Cell | Field | All | Fresh later step | Fetch |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| standalone-image-static-1d | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-image-static-1d | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-image-fixed-1d | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-image-fixed-1d | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-rectilinear-static-1d | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-rectilinear-static-1d | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-rectilinear-fixed-1d | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-rectilinear-fixed-1d | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-rectilinear-moving-1d | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-rectilinear-moving-1d | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-structured-static-1d | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-structured-static-1d | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-structured-fixed-1d | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-structured-fixed-1d | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-structured-moving-1d | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-structured-moving-1d | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-hypertree-static-1d | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-static-1d | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-fixed-1d | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-fixed-1d | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-moving-1d | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-moving-1d | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-dynamic-1d | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-dynamic-1d | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-image-static-2d | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-image-static-2d | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-image-fixed-2d | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-image-fixed-2d | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-rectilinear-static-2d | pass | pass | pass | pass | pass | pass | pass | fail |
| mb-rectilinear-static-2d | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-rectilinear-fixed-2d | pass | pass | pass | pass | pass | pass | pass | fail |
| mb-rectilinear-fixed-2d | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-rectilinear-moving-2d | pass | pass | pass | pass | pass | pass | pass | fail |
| mb-rectilinear-moving-2d | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-structured-static-2d | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-structured-static-2d | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-structured-fixed-2d | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-structured-fixed-2d | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-structured-moving-2d | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-structured-moving-2d | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-hypertree-static-2d | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-static-2d | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-fixed-2d | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-fixed-2d | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-moving-2d | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-moving-2d | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-dynamic-2d | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-dynamic-2d | fail | pass | not-applicable | pass | fail | fail | fail | fail |

## Four-rank output, ParaView 6.2.1

Writer: 283ef55-plus-hypertree-mb
pvpython: /opt/paraview/6.2.1/bin/pvpython
ParaView: 6.2.1; VTK: 9.7.2
Cache: default

| Case | Unchanged reader | Geometry | Point | Cell | Field | All | Fresh later step | Fetch |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| standalone-image-fixed-mpi | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-image-fixed-mpi | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-rectilinear-moving-mpi | pass | pass | pass | pass | pass | pass | pass | fail |
| mb-rectilinear-moving-mpi | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-structured-moving-mpi | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-structured-moving-mpi | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-unstructured-fixed-mpi | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-unstructured-fixed-mpi | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-unstructured-moving-mpi | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-unstructured-moving-mpi | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-unstructured-dynamic-mpi | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-unstructured-dynamic-mpi | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-hypertree-static-mpi | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-static-mpi | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-fixed-mpi | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-fixed-mpi | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-moving-mpi | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-moving-mpi | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-dynamic-mpi | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-dynamic-mpi | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| all-fixed-mpi | fail | pass | pass | pass | fail | fail | fail | fail |
| different-modes-mpi | fail | pass | pass | pass | fail | fail | fail | fail |
| different-modes-reversed-mpi | fail | pass | pass | pass | fail | fail | fail | fail |

## Serial reader cache enabled, ParaView 6.2.1

Writer: 283ef55-plus-hypertree-mb
pvpython: /opt/paraview/6.2.1/bin/pvpython
ParaView: 6.2.1; VTK: 9.7.2
Cache: on

| Case | Unchanged reader | Geometry | Point | Cell | Field | All | Fresh later step | Fetch |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| standalone-image-fixed | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-image-fixed | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-rectilinear-moving | pass | pass | pass | pass | pass | pass | pass | fail |
| mb-rectilinear-moving | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-structured-moving | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-structured-moving | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-unstructured-moving | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-unstructured-moving | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-unstructured-dynamic | pass | pass | pass | pass | pass | pass | pass | pass |
| mb-unstructured-dynamic | fail | pass | pass | pass | fail | fail | fail | fail |
| standalone-hypertree-fixed | pass | pass | not-applicable | pass | fail | fail | pass | fail |
| mb-hypertree-fixed | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-moving | pass | pass | not-applicable | pass | fail | fail | pass | fail |
| mb-hypertree-moving | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-dynamic | pass | pass | not-applicable | pass | fail | fail | pass | fail |
| mb-hypertree-dynamic | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| different-modes | fail | pass | pass | pass | fail | fail | fail | fail |
| different-modes-reversed | fail | pass | pass | pass | fail | fail | fail | fail |

## MPI-written HTG cache controls, ParaView 6.2.1

Writer: 283ef55-plus-hypertree-mb
pvpython: /opt/paraview/6.2.1/bin/pvpython
ParaView: 6.2.1; VTK: 9.7.2
Cache: on

| Case | Unchanged reader | Geometry | Point | Cell | Field | All | Fresh later step | Fetch |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| standalone-hypertree-fixed-mpi | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-fixed-mpi | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-moving-mpi | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-moving-mpi | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-dynamic-mpi | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-dynamic-mpi | fail | pass | not-applicable | pass | fail | fail | fail | fail |

## HTG-containing baseline subset, ParaView 6.2.0

Writer: 283ef55-plus-hypertree-mb
pvpython: /opt/paraview/6.2.0/bin/pvpython
ParaView: 6.2.0; VTK: 9.7.1
Cache: default

| Case | Unchanged reader | Geometry | Point | Cell | Field | All | Fresh later step | Fetch |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| standalone-hypertree-static | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-static | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-fixed | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-fixed | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-moving | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-moving | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| standalone-hypertree-dynamic | pass | pass | not-applicable | pass | pass | pass | pass | pass |
| mb-hypertree-dynamic | fail | pass | not-applicable | pass | fail | fail | fail | fail |
| all-static | fail | pass | fail | fail | fail | fail | fail | fail |
| all-fixed | fail | pass | fail | fail | fail | fail | fail | fail |
| different-modes | fail | pass | fail | fail | fail | fail | fail | fail |
| different-modes-reversed | fail | pass | fail | fail | fail | fail | fail | fail |
| mixed-modes | fail | pass | fail | fail | fail | fail | fail | fail |
| mixed-modes-reversed | fail | fail | fail | fail | fail | fail | fail | fail |
| hypertree-mixed | fail | worker failed | worker failed | worker failed | worker failed | worker failed | worker failed | worker failed |
| hypertree-mixed-reversed | fail | worker failed | worker failed | worker failed | worker failed | worker failed | worker failed | worker failed |
| hypertree-workaround | fail | pass | not-applicable | pass | fail | fail | fail | fail |

## Controls without HTG, ParaView 6.2.0

Writer: 283ef55-plus-hypertree-mb; controls without HTG
pvpython: /opt/paraview/6.2.0/bin/pvpython
ParaView: 6.2.0; VTK: 9.7.1
Cache: default

| Case | Unchanged reader | Geometry | Point | Cell | Field | All | Fresh later step | Fetch |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| four-type-static | fail | pass | fail | fail | fail | fail | fail | fail |
| standalone-rectilinear-static | pass | pass | pass | pass | pass | pass | pass | fail |
