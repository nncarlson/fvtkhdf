# Test Notes

The test suite is organized into a few categories with different goals.

## Writer Contract Tests

These are the small Fortran executables that exercise one writer mode or one
specific API contract.

- `vtkhdf_unstructured_static_file_type`: UnstructuredGrid static mesh and static data writes.
- `vtkhdf_unstructured_fixed_file_type`: UnstructuredGrid fixed-mesh temporal writes.
- `vtkhdf_moving_unstructured_file_type`: UnstructuredGrid moving-mesh geometry updates.
- `vtkhdf_dynamic_unstructured_file_type`: UnstructuredGrid dynamic-mesh topology and geometry updates.
- `vtkhdf_mb_basic_blocks_file_type`: MB ordinary two-block static and temporal writes.
- `vtkhdf_mb_collections`: all five simple block types share ownership, naming,
  and static, temporal, or mixed lifecycles; layout and reader checks verify
  assembly links, block types, metadata, and field forwarding.
- `vtkhdf_mb_hypertree`: fixed, moving, and dynamic HyperTreeGrid blocks share
  an MB timeline. Layout and VTK cursor checks verify custom chunk settings,
  per-block geometry and topology, masks, cell values and omitted-write reuse,
  explicit and implicit step finalization, and empty parts with changing tree
  ownership. MPI builds run both two- and four-rank fixtures.
- `vtkhdf_chunk_settings`: chunk-layout configuration checks.
- `test_leak`: repeated HDF5 init sanity check.
- `vtkhdf_simple_dataset_lifecycle`: all five concrete block types borrow one
  parent file and communicator, write static and temporal field data, and
  close/reinitialize through the abstract lifecycle interface. It checks close
  with an open step, repeated close, reset of geometry state, and preservation
  of the parent resources. Separate standalone siblings exercise failed
  creation, reuse, flushing, and owned-resource cleanup.

## Interop Fixtures

The broad UnstructuredGrid and MB programs are no longer the primary writer coverage. They now
mainly produce representative files for outer checks.

- `vtkhdf_unstructured_fixture_file`: broad UnstructuredGrid fixture producer.
- `vtkhdf_mb_fixture_file`: broad MB fixture producer.
- `vtkhdf_unstructured_reader_smoke`: VTK readback of the UnstructuredGrid fixture.
- `vtkhdf_mb_reader_smoke`: VTK readback of the MB fixture.

The broad fixtures are intentionally narrower than before. They should only
carry the datasets needed by the readback and attribute checks layered on top of
them.

## Structural Checks

These are `h5dump`-based checks layered on top of fixture or regression files.
They protect on-disk metadata that reader smoke alone would miss.
The checks remain registered when `h5dump` is unavailable at configuration
time; CTest reports them as skipped. Run them with
`ctest --test-dir build -L h5dump`.

- UnstructuredGrid fixture attribute checks for static and temporal point/cell data.
- Dynamic UnstructuredGrid temporal mesh checks for `NumberOfPoints`, `NumberOfCells`,
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

These tests protect API ordering and handle-usage assertions that intentionally
abort. Most use `WILL_FAIL`. `vtkhdf_image_fail_unstructured_duplicate` instead uses
`vtkhdf_expect_assertion.py` to require both a nonzero exit and the specific
duplicate-write assertion diagnostic. Its UnstructuredGrid setup runs independently of the
ImageData setup, and the expected source location is resolved from the template
so unrelated setup errors cannot satisfy the test.

Repeatable `--assertion` options allow a test to accept any of several precise
assertion locations. The interface `operation_agreement` test accepts either
the write-side or registration-side agreement assertion: any rank may abort
first. `vtkhdf_expect_assertion_checks` verifies both alternatives with isolated
diagnostics, including a surviving nonzero-rank log, while still rejecting
unrelated errors, successful commands, missing diagnostics, and stale logs.

`vtkhdf_assert_returning_mpi` compiles the actual fatal handler with isolated,
test-only MPI and HDF5 modules. Its `MPI_Abort` deliberately returns; the checker
requires a failed `INSIST` to call `MPI_Abort` and exit unsuccessfully, while a
successful `INSIST` must exit normally without calling `MPI_Abort`. The driver
ends normally if the fatal handler returns, so the exit status detects a missing
local termination. The test runs as one process,
without real MPI/HDF5 libraries or an MPI launcher, in both serial and MPI
configurations, including Release. It tests the application's fallback, not
the native MPI runtime's abort behavior.

The MPI `geometry_mismatch`, `name_mismatch`, and `mode_mismatch` checks are
registered only for Debug configurations: they exercise debug-only `ASSERT`
diagnostics. Local argument and lifecycle checks use unconditional `INSIST`
and remain tested in Release. See the [runtime checking policy](../doc/developer.md).

- `vtkhdf_unstructured_fail_late_register`
- `vtkhdf_unstructured_fail_wrong_cell_count`
- `vtkhdf_unstructured_fail_data_before_mesh`
- `vtkhdf_mb_fail_wrong_block_handle`
- `vtkhdf_mb_fail_late_block`

## Optional VTK Python Reader Checks

The reader smoke tests use VTK's `vtkHDFReader`. They are enabled only when the
configured Python interpreter can import the `vtk` package.

`vtkhdf_mb_reader_smoke` always performs its ordinary composite readback checks
when VTK is available. The separate `vtkhdf_mb_reader_temporal_state` also
requires the step count and temporal flag to survive each data read. It skips
released VTK versions through 9.7.2 because of VTK issue 20246 (MR !14030).
Development builds and newer versions run this check strictly; their failures
are not converted to skips. The independent compatibility harness continues
to exercise the defect on all available reader versions.

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

The independent [MB compatibility project](compatibility/mb/README.md) checks
an installed library through a selected ParaView version's `pvpython`. It is
run explicitly, outside the top-level build and CTest suite, and records reader
limitations separately from writer/layout failures. Fixtures can be reused
across reader versions.

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
numeric kind and both scalar and multicomponent data. A UnstructuredGrid/ImageData parity fixture
checks the same labels in every association, with blank names, trimming,
sanitization collisions, caller-supplied suffixes, and mixed static/temporal
fields; HDF5 metadata and VTK readback must agree. Repeated static names
may have different numeric kinds or component counts; static FieldData writes
remain complete, independent values.

`vtkhdf_image_layout` requires NumPy and h5py. It independently checks values,
dataset rank/order, chunks, 32-bit WholeExtent, 64-bit offsets, and snapshot
counts. It compares static and temporal FieldData datasets and metadata with
UnstructuredGrid reference files, including variable shape, initial omission, and reuse.
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
UnstructuredGrid files; it does not skip a malformed FieldData result.

The expected-failure ImageData tests cover local selection validity, type and
component consistency, ordering, missing initial snapshots, foreign handles,
field-data ordering and repeated writes to one field-data handle in a step,
and MPI metadata disagreement.
Creation rejects moving, dynamic, and invalid ImageData modes, including
disagreement between MPI ranks about otherwise supported modes.
Tests intentionally do not diagnose global coverage or overlap.

`vtkhdf_constant_field_test.F90.fypp` covers by-name constant field data for all
simple writers and supported modes, including changing UnstructuredGrid topology. It writes
all four numeric kinds as scalars, scalar tuples, vectors, and noncontiguous
multicomponent arrays. Temporal files also include empty constants. Fields are
created before mesh definition, during an open step, between steps, and after
the last finalized step. Mixed constant/registered names exercise disambiguation
in both creation orders. Independent layout and VTK checks verify complete
offset/size histories, values at every step, rank-zero payload selection, and
unchanged registered-field omission, reuse, and shape changes. Non-temporal
files retain their original layouts without a Steps group. VTK versions before
9.7 read the UnstructuredGrid and ImageData fixtures; the RectilinearGrid and StructuredGrid
reader cases run only with VTK 9.7 or later. Layout checks cover all types.

The simple writers own their registries and construct public data handles.
Each handle records its owning simple dataset, including a MultiBlock leaf.
Writes go directly to typed block pointers; the composite does not forward data
handles or add a block index to them.
Foreign-file handle checks cover cell, point, and field data in UnstructuredGrid and
MultiBlock, including matching registry indices, alongside the existing
ImageData ownership check.

See the [developer guide](../doc/developer.md) for writer organization,
resource ownership, shared data helpers, and the runtime checking policy.

## HyperTreeGrid

The shared lifecycle fixture (`vtkhdf_simple_dataset_test.F90`) includes HTG
file and block siblings. It checks repeated close and reuse, changed root grids
and branch factors, automatic finalization of an open step, inherited field
data methods, failed standalone creation and reopening, and preservation of
the borrowed parent file, group, and communicator. This exercises the block
directly without MultiBlock integration.

`vtkhdf_hypertree_test.F90` writes static, fixed, moving, and dynamic fixtures in
one, two, and three dimensions with branch factors 2 and 3. Trees follow three
refinement patterns chosen by tree id and topology phase, and are dealt to the
MPI ranks by a phase-dependent rule that leaves the last rank without trees in
one phase, producing an empty part while retaining one stored part per rank.
In the two-rank dynamic fixture, cell counts change from `[35, 35]` to
`[70, 0]` with `NumberOfParts` remaining 2. Both layout and reader checks cover
this transition and returning to the earlier populated state. The fixtures
cover both index kinds, both coordinate kinds, either order of
geometry and topology, every cell data kind, masks supplied for all, some, or
none of the topologies, a grid without any refined tree, static cell data
written during a later step, and temporal omission and reuse.

`vtkhdf_hypertree_check.py` derives the expected trees from the tree ids of each
stored part, so the same checks apply to any number of ranks. The layout check
uses h5py to verify attributes, types, per-part metadata, packed descriptor
and mask bytes, coordinates, cell data, every temporal offset, and that
storage grows only for a new topology, geometry, or data write. The reader
check requires VTK 9.7 or later. It traverses every tree with a cursor,
comparing refinement, cell bounds, masking, and cell data against the
expected values at every step, read in ascending and then descending order.

`vtkhdf_hypertree_fail_test.F90` checks grid definition, coordinate lengths, tree
ids, incomplete and oversized trees, masked refined cells, ordering of mesh
and data writes, mode restrictions, late registration, and reuse of temporal
data after a change of cell count.

The HTG fixtures additionally cover all seven active-axis combinations, both
branch factors, and both root-indexing orders. Independent VTK cursor checks
verify the resulting tree positions, including transposed root indexing.

`vtkhdf_hypertree_api_test.F90.fypp` exercises every cell-data numeric kind as a
scalar and vector in static and temporal files, using noncontiguous arrays,
empty MPI contributions, changing ownership, and temporal omission/reuse.
Its layout checker verifies values, types, offsets, and custom versus default
chunk sizes. Cell-name fixtures cover blank names, trimming, sanitization and
suffix collisions, and static/temporal collisions in both creation orders;
reader checks verify that each temporal handle selects its own dataset.

The field-data fixture tests only HTG's integration with the shared component:
correct parent groups, rank-zero values, completed-step counts for late
constants, write-flag resets, and finalization of the last open step. Field
readback uses a fresh reader per step because VTK 9.7.2 can drop reused constant
field arrays during repeated reads. The shared component's full kind/shape
matrix is not repeated here.

Additional negative checks cover foreign cell and field handles, duplicate cell
writes, writes outside an open step, late field registration, missing initial
cell values, and invalid modes. Foreign-handle and duplicate-write tests require
the specific assertion diagnostic, so unrelated setup failures cannot satisfy
them. MPI rank disagreements in mode, dimensions, and branch factor run only in
Debug builds, where the corresponding diagnostics are enabled.

## Constant Point and Cell Data on Moving Unstructured Meshes

`vtkhdf_unstructured_constant_data_test.F90.fypp` writes standalone static,
fixed-mesh, and moving-mesh controls, plus two moving UnstructuredGrid blocks
in one MB file. All five numeric kinds have scalar and noncontiguous vector
payloads. Constants are created before stepping, in early and later open
steps, after temporal writes, between steps, and after finalization. The MB
owner closes an open final step. MPI partitions have unequal sizes and rank
zero is empty.

Independent layout and VTK checks verify one stored constant payload, complete
zero-offset histories, geometry movement and reuse, unchanged temporal-array
reuse, names and attributes, and identical per-block behavior. Static and
fixed-mesh controls retain their existing offset layouts. Expected-failure
cases retain dynamic-mesh rejection, collective entity-count checks, and
first-write requirements for registered temporal data.

## Moving RectilinearGrid and StructuredGrid

`vtkhdf_moving_grid_test.F90.fypp` writes both real kinds, all seven active-axis
combinations, and initial geometry before or during the first step. Seven time
steps reference three geometry snapshots, while point/cell data uses a different
replacement schedule. Geometry ownership changes across snapshots; tests include
halo sections, empty contributions, all-empty geometry replacements, writer
reuse, and both two- and four-rank execution. The shared lifecycle test also
uses moving geometry in borrowed leaves.
Static point/cell data is introduced during an open step, between steps, and
after the last step (open or finalized). Incremental static writes span different
steps; layout and reader checks verify their final values at every time, including
earlier steps, with a single snapshot and complete zero-offset histories.

`vtkhdf_moving_grid_check.py` independently verifies payloads and per-step offsets
with h5py, then checks geometry and data with VTK. Sequential reads use the
default cache. Nonsequential reads disable the legacy cache where that API is
available: VTK 9.7.1 can retain a previously attached point/cell array on a cache
hit. This defect is fixed in VTK 9.7.2, and ParaView 6.2.1 is expected to include
the fix. The file's offsets and uncached values are correct. Nontrailing singleton
reader cases use the same existing version-dependent skips as fixed grids;
all layout checks run.

`vtkhdf_moving_grid_fail_test.F90` checks mode, geometry ordering, immutable
bounds, coordinate kind/size consistency, and the required initial temporal-data
write before reuse.

The `vtkhdf_hypertree_interface_*` tests check paired material-interface
writes and registration with both real kinds. Independent HDF5 checks and
VTK readback verify actual names, every component, storage kinds, paired
offsets, replacement and reuse, and singleton reset on close. MPI fixtures
run with two and four ranks, including empty parts, uneven ownership,
noncontiguous inputs, and topology changes. Contract failures check the
intended assertion diagnostic, including first-step data on empty ranks,
kind mismatches, lifecycle rules, and Debug collective agreement.

The `vtkhdf_mb_hypertree_*` and `vtkhdf_pdc_hypertree_*` tests also check static
and temporal interface pairs on composite HTG blocks. Layout and VTK readback
checks cover block-level attributes and component values, replacement after
topology changes, omission reuse, and step finalization by the owning file,
including empty MPI contributions.
