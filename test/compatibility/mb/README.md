# MultiBlock reader compatibility

This is an independent CMake project, like the consumers in `test/packaging`.
It is not included by the library's top-level CMake project or its CTest suite.
It generates small files using an installed fVTKHDF library and checks them
through ParaView's reader pipeline. Reader incompatibilities are recorded as
observations so that one failure does not stop the investigation.

## Build and run

Use an installed fVTKHDF library from the branch under investigation and the
same Fortran compiler used to build it. CMake detects whether the installation
supports MPI and builds the generator accordingly. The default matrix uses
serial writing and reading, full 3D grids, and double-precision point/cell
arrays. The additional suites below extend writing, dimensions, and caching;
all reader checks still use a serial builtin ParaView session.

From the repository root, with the compiler and HDF5 environment loaded:

```sh
cmake -S test/compatibility/mb -B build-mb-compat \
  -DCMAKE_PREFIX_PATH=/path/to/fvtkhdf-install -DCMAKE_BUILD_TYPE=Debug
cmake --build build-mb-compat --parallel

python3 test/compatibility/mb/run.py \
  --generator build-mb-compat/generate \
  --pvpython /opt/paraview/6.2.1/bin/pvpython \
  --work-dir build-mb-fixtures \
  --results-dir build-mb-results-6.2.1 \
  --writer-label COMMIT_OR_VERSION
```

The controlling Python needs NumPy and h5py. ParaView's `pvpython` needs its
bundled NumPy; it does not need h5py or access to the controlling Python's
packages. No rendering or GUI is required. `--disable-registry` avoids reading
or modifying ParaView's user settings.

To generate and validate fixtures without a reader, use `--generate-only`
and omit `--pvpython`. To check another ParaView version against exactly the
same files, use:

```sh
python3 test/compatibility/mb/run.py --read-only \
  --pvpython /path/to/another/paraview/bin/pvpython \
  --work-dir build-mb-fixtures \
  --results-dir build-mb-results-another-version
```

Choose a different results directory for each version. Each directory contains
`summary.md`, aggregate `results.json`, per-file JSON observations, and logs.
Results identify the actual ParaView and VTK versions. Generated VTKHDF files,
case descriptions, writer/layout logs, and the matrix manifest remain in the
work directory. Keep these generated artifacts outside the source tree or in
directories beginning with `build`, which the repository ignores.

The driver exits nonzero for generation/layout failures, worker crashes, or
timeouts (`--timeout`, default 90 seconds per worker). A completed reader
investigation exits zero even when compatibility checks fail: consult the
summary and detailed observations. A profile-level Python exception is recorded
as a failed observation with its traceback, never as a pass or a skipped case.
When a worker terminates during native reader execution, its last active
profile and completed observations are retained. The mixed static/temporal
HTG cases currently cause such worker failures; this does not prevent the
remaining files from being checked.

## Selective expansion

`--suite dimensions` generates 48 files: all eight supported structured and
four HyperTreeGrid type/mode configurations in 1D and 2D, each as MB and standalone output.
Only trailing singleton axes are used: `(nx,1,1)` and `(nx,ny,1)`.
Nontrailing singleton arrangements are deliberately excluded from this
black-box comparison with standalone controls.

```sh
python3 test/compatibility/mb/run.py --suite dimensions \
  --generator build-mb-compat/generate \
  --pvpython /opt/paraview/6.2.1/bin/pvpython \
  --work-dir build-mb-dimensions --results-dir build-mb-dimensions-6.2.1
```

For MPI writing, configure another build of this project against an **MPI
fVTKHDF installation** and load its matching compiler/MPI/HDF5 environment.
`--suite parallel --ranks 4` generates 23 files: six representative structured/
unstructured configurations and all four HyperTreeGrid modes with standalone
controls, plus three five-type collections.
Every block has unequal rank contributions and a block-dependent empty rank,
including rank zero. Structured grids use disjoint cell slabs and explicitly
owned shared points; unstructured parts grow and shrink in dynamic mode.
HyperTreeGrid keeps one partition per rank, including empty ranks, and rotates
tree ownership between ranks in dynamic mode while retaining the global root grid.
Coordinate vectors and field data still follow their rank-zero API contracts.

```sh
python3 test/compatibility/mb/run.py --suite parallel --ranks 4 \
  --generator build-mb-compat-mpi/generate --generate-only \
  --work-dir build-mb-parallel

python3 test/compatibility/mb/run.py --read-only \
  --pvpython /opt/paraview/6.2.1/bin/pvpython \
  --work-dir build-mb-parallel --results-dir build-mb-parallel-6.2.1
```

The driver uses `mpiexec -n N` for an MPI generator (including a one-rank run).
Use `--mpiexec` to select the launcher matching the installed writer. A serial
generator is rejected when multiple ranks are requested. MPI runs require
elevated execution permissions in the development sandbox. Reader execution
remains serial: this suite checks MPI-written files, not distributed reading.

`--suite cache` selects 18 representative temporal files from the existing
baseline manifest and enables `UseCache=1`. These include both standalone and
MB output, geometry/data reuse, dynamic meshes, forward/backward requests,
and heterogeneous collections in both orders. No files need regenerating:

```sh
python3 test/compatibility/mb/run.py --suite cache --read-only \
  --pvpython /opt/paraview/6.2.1/bin/pvpython \
  --work-dir build-mb-fixtures --results-dir build-mb-cache-6.2.1
```

`--cache on|off|default` can also be used with any suite. The actual reader
options are recorded in the results. When reading an existing manifest,
`--suite` selects its corresponding cases and rejects missing fixtures;
omitting it reads every manifest case. The 256 ordered-pair matrix remains
outside this selective expansion.

## Matrix

`run.py` defines the cases as configuration data:

- All 16 supported type/mode combinations as single-block MB files.
- The same 16 configurations as standalone controls.
- Five-type collections with all-static, all-fixed, and different temporal
  modes; the latter is also written in reverse block order.
- Mixed static/temporal collections in both block orders: a five-type
  collection, two unstructured blocks, static rectilinear with temporal
  structured, static structured with temporal image, and two HyperTreeGrid blocks.
- Three workaround cases that replace the static block by fixed-mesh mode,
  leaving that block's arrays constant with no registered temporal data.

This gives 49 files. The larger ordered-pair matrix is not generated by default.
There are no per-case Fortran programs or checked-in binary fixtures.

Each temporal file has five nonuniform time values. Moving geometry changes at
steps 2 and 4. Point data changes at steps 0, 1, and 3; cell data changes at steps
0, 2, and 4. Dynamic unstructured meshes grow and shrink, with point data also
rewritten when its entity count changes. Constant arrays retain their values.

Every block has scalar and three-component cell arrays, static field
data (`block_id`, `direction`), and, in temporal modes, a scalar `clock` and a
two-component `table` whose tuple count changes. Dynamic unstructured blocks
and HyperTreeGrid blocks have only temporal cell data in dynamic mode, as required
by the API. The other types also have point arrays; HyperTreeGrid has no point
data, so point profiles for HTG-only files are marked `not-applicable`. Values identify
the block, entity, component, and data step. All blocks use the same array names
and compatible types/component counts.
The fixture-only `fixed-constant` argument selects ordinary
`VTKHDF_FIXED_MESH` mode and suppresses temporal-data registration/writes for
that block; it is not a new library mode.

The layout worker independently checks values, offsets, field sizes, geometry,
unstructured connectivity, HTG root ownership and packed topology/mask bits,
and assembly links before invoking ParaView. It
does not add a composite-level `Steps` group or change any file to accommodate
a reader. Such workarounds belong in separately identified diagnostic copies.

The reader worker uses `paraview.simple.VTKHDFReader`, enables array selections
explicitly, updates the pipeline, and inspects the actual server-side output
in the serial builtin session. Separate
profiles read geometry alone, point data, cell data, field data, and all data.
Each profile starts with a fresh reader and checks forward and nonsequential
time requests. A final fresh-reader profile starts directly at a later time.
A separate `unchanged` profile uses `vtkHDFReader` directly and reads step zero
immediately after `UpdateInformation()`. It leaves array selections unchanged,
asserts that the reader's modification time has not changed before the data
request, and checks that the step count and temporal flag survive the read.
It also verifies geometry and array values at step zero. Changing an array
selection after discovering metadata marks the reader modified and can mask
the composite temporal-state reset reported in VTK issue 20246; the other
profiles retain their selection changes to test those paths independently.

The `unchanged` profile also runs with ordinary Python and VTK, without ParaView:

```sh
python3 test/compatibility/mb/check.py vtk \
  build-mb-fixtures/mb-image-fixed.json build-mb-results-vtk.json
```

As with the ParaView worker, inspect the recorded `profiles.unchanged.status`;
reader incompatibilities are recorded in JSON rather than the process exit code.
A separate `fetch` profile checks transfer through `servermanager.Fetch`.
Its transfer/append path can introduce errors even when the reader output is
correct; these are distinguished from reader failures.
Checks cover time discovery, block types, geometry, counts, unstructured
connectivity, HTG cursor bounds/refinement/masks, values, component counts, and numeric types. Reader errors and
missing arrays are failures; they are not replaced by successful metadata-only
checks. Extra reader-generated arrays such as `Time` are permitted.

The independent expected-value functions are in `check.py`. A passing
standalone control and failing MB case help isolate composite-specific
behavior. An array appearing in ParaView's selection list is not sufficient:
its values must be present on the expected output leaf.

The current findings and practical workarounds are documented in the
[composite dataset reference](../../../doc/sphinx/composite.rst).
The [ParaView 6.2.1 baseline](results-6.2.1.md) retains the initial case summary;
detailed JSON observations and logs are generated by the harness.
The [selective expansion results](results-expanded-6.2.1.md) record the
additional dimensionality, MPI-writing, and cache-enabled checks.
The [extended five-type results](results-hypertree-6.2.1.md) include HyperTreeGrid
in all four suites and supersede those historical four-type summaries.

## Shared Assembly probes

The main CTest suite also generates nested MB and PDC assembly fixtures,
including repeated placements, empty groups, default root references, and
commitment during an open time step. Independent h5py layout checks determine
VTKHDF 2.8 conformance. Reader probes use only VTK's Python module and record
results in `test/mb-assembly-reader-compatibility.json` and
`test/pdc-assembly-reader-compatibility.json` under the build directory.
Reader crashes and timeouts are captured in isolated workers and do not fail
writer acceptance. Both unchanged-reader and selection-changing sequences
exercise the nested temporal fixtures.
