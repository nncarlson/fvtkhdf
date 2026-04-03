# Test Notes

The test suite is organized into a few categories with different goals.

## Writer Contract Tests

These are the small Fortran executables that exercise one writer mode or one
specific API contract.

- `vtkhdf_ug_static_file_type`: UG static mesh and static data writes.
- `vtkhdf_ug_fixed_file_type`: UG fixed-mesh temporal writes.
- `vtkhdf_moving_ug_file_type`: UG moving-mesh geometry updates.
- `vtkhdf_dynamic_ug_file_type`: UG dynamic-mesh topology and geometry updates.
- `vtkhdf_mb_basic_blocks_file_type`: MB ordinary two-block static and temporal writes.
- `vtkhdf_chunk_settings`: chunk-layout configuration checks.
- `test_leak`: repeated HDF5 init sanity check.

## Interop Fixtures

The broad UG and MB programs are no longer the primary writer coverage. They now
mainly produce representative files for outer checks.

- `vtkhdf_ug_fixture_file`: broad UG fixture producer.
- `vtkhdf_mb_fixture_file`: broad MB fixture producer.
- `vtkhdf_ug_reader_smoke`: VTK readback of the UG fixture.
- `vtkhdf_mb_reader_smoke`: VTK readback of the MB fixture.

The broad fixtures are intentionally narrower than before. They should only
carry the datasets needed by the readback and attribute checks layered on top of
them.

## Structural Checks

These are `h5dump`-based checks layered on top of fixture or regression files.
They protect on-disk metadata that reader smoke alone would miss.

- UG fixture attribute checks for static and temporal point/cell data.
- Dynamic UG temporal mesh checks for `NumberOfPoints`, `NumberOfCells`,
  `PointOffsets`, `CellOffsets`, and `ConnectivityIdOffsets`.
- MB fixture checks for assembly links and block-name-to-`Group_i` mapping.
- MB zero-part checks for sparse-block `NumberOfParts`,
  `NumberOfPoints`, and temporal point-data offsets.
- Chunk-layout checks for default and custom chunk sizing.

## Targeted MPI Regressions

These are MPI-specific tests where the distributed layout is the behavior under
test rather than just the execution mode.

- `vtkhdf_mb_zero_part_file_type`: one MB block includes a 0-sized part on one rank.
- `vtkhdf_mb_zero_part_reader_smoke`: VTK readback of that targeted zero-part file.

## Negative Contract Tests

These are expected-failure tests registered with `WILL_FAIL`. They protect API
ordering and handle-usage assertions that intentionally abort.

- `vtkhdf_ug_fail_late_register`
- `vtkhdf_ug_fail_wrong_cell_count`
- `vtkhdf_ug_fail_data_before_mesh`
- `vtkhdf_mb_fail_wrong_block_handle`
- `vtkhdf_mb_fail_late_block`

## Optional VTK Python Reader Checks

The reader smoke tests use VTK's `vtkHDFReader`. They are enabled only when the
configured Python interpreter can import the `vtk` package.

A repo-local virtual environment is a convenient way to keep this optional
dependency isolated:

```bash
python3 -m venv .venv
.venv/bin/pip install vtk
```

Configure CMake to use that interpreter:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DENABLE_MPI=OFF \
  -DPython3_EXECUTABLE=$PWD/.venv/bin/python
```

If `vtk` is not available, the build still succeeds and the reader checks are
skipped.

## fpm Notes

The pure Fortran executable tests are the best candidates for eventual `fpm`
coverage. The Python reader smoke tests, `h5dump`-based checks, and `CTest`
`WILL_FAIL` negative tests are CMake-first by design.
