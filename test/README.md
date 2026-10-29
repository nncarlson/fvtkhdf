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

## ImageData Checks

`vtkhdf_image_file_type` produces static, temporal, dimensionality, numeric-kind,
and all-empty snapshot fixtures in serial and MPI configurations. It includes
halo selections, differing Fortran actual bounds, repeated field names,
sanitization collisions, changed decomposition, unequal block counts padded
with zero-sized contributions, late static fields, and every generated numeric
kind for both scalar and multicomponent temporal data.
Registration-name checks exercise duplicate and blank labels, sanitization
collisions, and static/temporal label collisions in both creation orders.
Readback verifies that each returned handle writes its own dataset and that
name-based static writes create distinct fields. Static cell/point handles
assemble additional blocks, including empty first contributions, for every
numeric kind and both scalar and multicomponent data. A UG/ImageData parity fixture
checks the same labels in every association, with blank names, trimming,
sanitization collisions, caller-supplied suffixes, and mixed static/temporal
fields; HDF5 metadata and VTK readback must agree. Repeated static names
may have different numeric kinds or component counts; static FieldData writes
remain complete, independent values.

`vtkhdf_image_layout` requires NumPy and h5py. It independently checks values,
dataset rank/order, chunks, 32-bit WholeExtent, 64-bit offsets, and snapshot
counts. It compares static and temporal FieldData datasets and metadata with
UG reference files, including variable shape, initial omission, and reuse.
An internal structured-core fixture creates independently named static and
temporal leaves under an existing HDF5 group, checks their data and temporal
offsets, and verifies that closing the leaves leaves the parent handles open.
Its temporal leaf exercises moving-mode attribute storage independently of
the ImageData wrapper's more restrictive mode policy.
Independent temporal FieldData checks cover all four numeric kinds, growth
and reduction of component counts, varying tuple counts, noncontiguous input,
empty writes, fields never written, and non-root ranks supplying different
values. Both dataset axes must be extendible and offsets must count tuples.

`vtkhdf_image_reader_*` checks geometry, physical coordinates, every static
numeric kind, and temporal point/cell values with `vtkHDFReader`. Positive-extent
and nontrailing-singleton checks are skipped for VTK versions earlier than
9.7.2; they are mandatory for 9.7.2 and later so reader regressions fail the
tests. Layout checks remain enabled for all versions. Complete mesh-associated
validation was performed with a VTK source
build reporting 9.7.20261022. The temporal FieldData reader test strictly checks
component counts, tuple counts, types, and every value for both ImageData and
UG files; it does not skip a malformed FieldData result.

The expected-failure ImageData tests cover local selection validity, type and
component consistency, ordering, missing initial snapshots, foreign handles,
FieldData restrictions and repeated writes to one FieldData handle in a step,
and MPI metadata disagreement.
Creation rejects moving, dynamic, and invalid ImageData modes, including
disagreement between MPI ranks about otherwise supported modes.
Tests intentionally do not diagnose global coverage or overlap.

UG, MultiBlock, and ImageData file wrappers own public field handles; their cores
use registry indices. Foreign-file handle checks cover cell, point, and field
data in UG and MultiBlock, including matching field and block indices, alongside
the existing ImageData ownership check.

See [the writer implementation guide](../doc/implementation.md) for the common
organization of the UG and structured cores and the responsibilities of their
shared data modules.
