# ParaView 6.2.1 baseline

Historical four-type results. See the [extended five-type results](results-hypertree-6.2.1.md)
for the current matrix including HyperTreeGrid.

Checked 2026-11-17 with serial fVTKHDF at commit `7967b14`, ParaView 6.2.1,
and bundled VTK 9.7.2. All 38 fixtures passed independent layout checks.
The initial scope is serial 3D output; parallel and reduced-dimensional cases
were not investigated. Default reader options were used: `UseCache=0` and
`PieceDistribution=Block`. This run does not establish behavior with caching
enabled.

- All standalone controls passed the server-side reader checks, including field
  data. Fetching standalone RectilinearGrid output corrupted point/cell values
  even though the server-side output was correct. This is recorded separately
  in the Fetch column; it was not observed for the MB counterparts.
- All-static and all-temporal MB collections passed geometry and point/cell
  checks. Static and temporal field arrays were absent for every MB block type,
  and the reader advertised no selectable field arrays.
- Mixed static/temporal collections failed. Block order changed the symptoms:
  missing time information, missing blocks, and array read errors occurred.
- Promoting the otherwise static block to fixed-mesh mode restored geometry and
  point/cell checks, with constant data and no temporal-array registration or
  repeated writes on that block. This did not restore field data.
- The All, Fresh later step, and Fetch profiles include field-data checks, so
  those profiles fail even for MB cases with correct geometry and point/cell
  output. Consult the individual columns before interpreting a failure.
- The separate direct `vtkHDFReader` issue documented in the manual was
  reproduced with the earlier minimal report file, but that same file read
  correctly through ParaView's reader pipeline. It is not a general failure of
  temporal structured blocks in ParaView 6.2.1.

The table is a snapshot of `summary.md` produced by `run.py`. Detailed results
include per-profile failures, selected arrays, reader options, and versions.
Regenerate the fixtures and rerun the checker to investigate a changed reader.

Writer: 7967b14
pvpython: /opt/paraview/6.2.1/bin/pvpython
ParaView: 6.2.1; VTK: 9.7.2

| Case | Geometry | Point | Cell | Field | All | Fresh later step | Fetch |
| --- | --- | --- | --- | --- | --- | --- | --- |
| standalone-image-static | pass | pass | pass | pass | pass | pass | pass |
| mb-image-static | pass | pass | pass | fail | fail | fail | fail |
| standalone-image-fixed | pass | pass | pass | pass | pass | pass | pass |
| mb-image-fixed | pass | pass | pass | fail | fail | fail | fail |
| standalone-rectilinear-static | pass | pass | pass | pass | pass | pass | fail |
| mb-rectilinear-static | pass | pass | pass | fail | fail | fail | fail |
| standalone-rectilinear-fixed | pass | pass | pass | pass | pass | pass | fail |
| mb-rectilinear-fixed | pass | pass | pass | fail | fail | fail | fail |
| standalone-rectilinear-moving | pass | pass | pass | pass | pass | pass | fail |
| mb-rectilinear-moving | pass | pass | pass | fail | fail | fail | fail |
| standalone-structured-static | pass | pass | pass | pass | pass | pass | pass |
| mb-structured-static | pass | pass | pass | fail | fail | fail | fail |
| standalone-structured-fixed | pass | pass | pass | pass | pass | pass | pass |
| mb-structured-fixed | pass | pass | pass | fail | fail | fail | fail |
| standalone-structured-moving | pass | pass | pass | pass | pass | pass | pass |
| mb-structured-moving | pass | pass | pass | fail | fail | fail | fail |
| standalone-unstructured-static | pass | pass | pass | pass | pass | pass | pass |
| mb-unstructured-static | pass | pass | pass | fail | fail | fail | fail |
| standalone-unstructured-fixed | pass | pass | pass | pass | pass | pass | pass |
| mb-unstructured-fixed | pass | pass | pass | fail | fail | fail | fail |
| standalone-unstructured-moving | pass | pass | pass | pass | pass | pass | pass |
| mb-unstructured-moving | pass | pass | pass | fail | fail | fail | fail |
| standalone-unstructured-dynamic | pass | pass | pass | pass | pass | pass | pass |
| mb-unstructured-dynamic | pass | pass | pass | fail | fail | fail | fail |
| all-static | pass | pass | pass | fail | fail | fail | fail |
| all-fixed | pass | pass | pass | fail | fail | fail | fail |
| different-modes | pass | pass | pass | fail | fail | fail | fail |
| different-modes-reversed | pass | pass | pass | fail | fail | fail | fail |
| mixed-modes | pass | fail | fail | fail | fail | fail | fail |
| mixed-modes-reversed | fail | fail | fail | fail | fail | fail | fail |
| unstructured-mixed | fail | fail | fail | fail | fail | fail | fail |
| unstructured-mixed-reversed | fail | fail | fail | fail | fail | fail | fail |
| rectilinear-mixed | pass | fail | fail | fail | fail | fail | fail |
| rectilinear-mixed-reversed | fail | fail | fail | fail | fail | fail | fail |
| structured-mixed | pass | fail | fail | fail | fail | fail | fail |
| structured-mixed-reversed | fail | fail | fail | fail | fail | fail | fail |
| unstructured-workaround | pass | pass | pass | fail | fail | fail | fail |
| rectilinear-workaround | pass | pass | pass | fail | fail | fail | fail |
