# Developer Guide

This guide describes the source organization, resource ownership, runtime
checking policy, and validation practices used by fVTKHDF.

## Source Organization

The public file writers implement the API documented in the
[reference manual](sphinx/index.rst). The UnstructuredGrid writer contains its
mesh and data implementation directly. The three structured writers share
`vtkhdf_structured_base`; each concrete writer supplies its geometry methods,
VTK dataset type name, and supported modes.

`vtkhdf_simple_dataset` defines the common lifecycle contract used by composite
writers: cleanup, time stepping, and dataset-state queries. It does not define
geometry or data-writing interfaces. MultiBlock stores polymorphic simple
writers and dispatches format-specific operations to its supported block types.
Its public block interface currently supports UnstructuredGrid only.

| Module or type | Responsibility |
| --- | --- |
| `vtkhdf_ug_file_type` | Unstructured mesh, data registration and writes, and lifecycle policy |
| `vtkhdf_structured_base_class` | Logical grids, rectangular region writes, shared registration, and structured lifecycle policy |
| `vtkhdf_image_file_type`, `vtkhdf_rectilinear_file_type`, `vtkhdf_structured_file_type` | Geometry definitions and format-specific geometry storage |
| `vtkhdf_mb_file_type` | Composite assembly, block handles, and a shared timeline |
| `vtkhdf_field_data_type` | Shared field data registration, writes, temporal metadata, and cleanup |
| `vtkhdf_data_support` | Name resolution, registry records, data groups, and chunk sizing |
| `vtkhdf_file_context_type`, `vtkhdf_context_type` | File operations and communicator management |
| `vtkhdf_h5`, `vtkhdf_h5_c_binding`, `vtkhdf_h5_c_shim.c` | HDF5 operations and C bindings |

Edit `.F90.fypp` templates rather than generated `.F90` build products. When
adding a template, update both `src/CMakeLists.txt` and
`fpm/generate_fpm_src.py`. Keep generated overload families straightforward,
with the same organization for point and cell data.

### Resource and Handle Ownership

For simple dataset writers, `create` owns a standalone file and duplicated MPI
communicator through a `vtkhdf_file_context`. Internal `init` methods create
simple VTK datasets under a parent group and borrow the composite writer's
parallel context. Closing a leaf must not release its parent's file or
communicator. Implicit Fortran finalization cannot perform collective cleanup;
callers must use `close`.

Each writer owns its data registry and constructs the public handles it returns.
A data handle records the root HDF5 group ID of its owning VTK dataset; handle
access validates that owner. MultiBlock additionally records and checks the
block index before forwarding a handle to its leaf writer. Components remain
private so callers cannot construct indices independently of the owning writer.

`named_record` represents an entry in a namespace, including composite block
names and point, cell, or field dataset names. `dataset_record` adds temporal
bookkeeping. Association-specific records add the information needed for their
storage and reuse rules. Name resolution is shared so normalization and
collision handling remain consistent across writer types.

### Registration and Time Stepping

Data creation validates the local contract, reserves a unique name and registry
entry, creates the payload storage, and records the original name and optional
VTK attribute role. Each write by name creates a new dataset. Structured static
writes may return a handle for subsequent region writes to the same dataset.

Starting a step finalizes an open step, records the time, and resets temporal
write flags. Finalization checks required mesh data and appends the per-step
metadata. Omitted temporal writes reuse the preceding value where permitted;
structured region writes within one step assemble a single replacement.

The UnstructuredGrid and structured implementations use different mesh and
point/cell layouts. UnstructuredGrid appends entity tuples and records tuple
offsets. Structured writers use whole-grid snapshots, with several region writes
able to contribute to one snapshot. Repeated offsets provide reuse without
allocating another snapshot. Geometry snapshots and data snapshots are
independent of one another and of the number of time steps.

Field data uses the same implementation in all simple writers. In temporal
files, complete values are stored as tuple rows with a component dimension;
per-step offsets and sizes describe the current value. A static field in a
temporal file stores one value and repeats its offset and sizes for every step,
including earlier finalized steps when the field is created late. Only rank
zero supplies field data values.

## Runtime Checks and Collective-Call Contracts

This is the developer policy for deciding whether a condition needs a check,
and whether that check belongs in production builds. It applies to serial and
MPI builds. Existing code is not entirely consistent with this policy; use the
reason for a check, rather than a nearby assertion, to choose its form.

### Error-handling philosophy

fVTKHDF makes a best effort to complete a request, resolving incidental problems
when it can preserve the intended meaning. For example, conflicting dataset
names are disambiguated automatically. Even if a clash originates in erroneous
application input, returning an error would give the application no useful
recovery decision when the library can already accommodate the names.

User-facing procedures return errors only when the application can meaningfully
act on them during execution. Currently, file creation provides such a status
return for failure to create the requested file: the application might choose
another path or continue without that output. Programming errors, such as an
invalid handle or incorrect call order, require correcting the program and
are not treated as recoverable API outcomes.

Detected programming errors and unrecoverable library failures follow a
**fail-fast** policy: report the failure and terminate through the assertion
machinery. Continuing could corrupt output or leave MPI ranks waiting in
incompatible operations. A new status return should therefore have a concrete,
supported recovery use; an error code alone does not establish recoverability.
The checks below provide diagnostics and enforce contracts within this model.

### Principles

Collective participation and identical values for non-distributed arguments are
API preconditions. Callers must satisfy them. Checking every aspect of collective
usage is expensive and cannot reliably diagnose arbitrary violations: a rank
that never enters the operation cannot participate in its validation either.

Use inexpensive, unconditional checks for argument validity and lifecycle
constraints. Use debug-only communication to diagnose disagreements between
ranks. Keep communication required to compute results or establish shared state,
and retain collective checks on library-operation success where needed.

Communication cost includes latency and synchronization, not just data volume.
Several reductions or broadcasts around each small region write can be costly
even when they exchange only a few integers.

### Choosing a check

| Situation | Usual choice | Examples |
| --- | --- | --- |
| Already established by construction, an earlier check, or the language interface | No additional check | Rechecking a validated private helper's arguments; comparing an association tag solely to recognize which method other ranks entered |
| Useful diagnostic that shared arguments agree across ranks | `ASSERT(ctx%identical(...))` | Global bounds, names, numeric kinds, component counts, modes, optional-argument presence |
| Local argument or lifecycle validity | Local `INSIST(...)` | Supported mode, valid handle ownership/index, array shape, selection within the array and global grid, open time step |
| Success of an HDF5/MPI operation | Unconditional status handling; local or collective `INSIST` as appropriate | A valid HDF5 ID, successful transfer, successful MPI operation |
| Result required by the algorithm or shared control flow | Unconditional computation and communication | Global entity counts, offsets, whether any rank changed its mesh, whether an existing mesh remains globally valid |
| Actionable failure with a supported recovery path | Return the documented status and message | Failure to create the requested output file |

These are defaults, not a requirement to add a check for every conceivable
misuse. Choose checks that protect a concrete operation or provide a useful
diagnostic. Avoid repeating an established invariant throughout internal helpers.

### `ASSERT`: optional diagnostics

[`vtkhdf_assert.inc`](../src/vtkhdf_assert.inc) removes `ASSERT` when `NDEBUG`
is defined, including evaluation of its argument. The CMake build defines
`NDEBUG` for configurations other than Debug. Other build systems must select
their assertion configuration consistently across all ranks.

For example:

```fortran
ASSERT(ctx%identical([mode]))
INSIST(this%supports_mode(mode))
```

The first check diagnoses rank disagreement in a debug build. The second
unconditionally rejects a locally unsupported mode.

Expressions in `ASSERT` must have no effects needed for correct execution.
Never put an HDF5 operation, communicator initialization, required reduction,
or state update inside it. Perform required work first, then check its result.

Agreement checks are diagnostic aids under an otherwise compatible collective
sequence. They do not make arbitrary mismatched API calls safe. Do not build a
protocol of matching constant tags or padding values across different methods
in an attempt to recognize every such mismatch.

### `INSIST`: local contracts

`INSIST` remains active with `NDEBUG`. Use it for inexpensive conditions that
must hold before proceeding:

```fortran
INSIST(this%step_open)
INSIST(allocated(vars))
INSIST(index >= 1 .and. index <= size(vars))
```

Check prerequisites separately before evaluating expressions that depend on
them. Fortran does not guarantee short-circuit evaluation; combining allocation
and array access in one logical expression does not protect that access.

For a nonempty rectangular selection, check locally that its logical bounds
fit both the global grid and the array's mapped footprint. These are local
argument checks; they do not need a `global_all` merely to make every rank
report the same assertion outcome. Empty contributions skip those checks and
do not need matching `global_all(.true.)` calls.

Similarly, if only rank zero's coordinate vectors are referenced, only rank
zero needs to check their lengths. Other ranks' ignored arguments need not be
validated. A locally computed invalid condition is still a fatal contract error.

### Library failures and required reductions

Keep library operations outside assertions and handle their status
unconditionally. For example, on a collective dataset-creation path:

```fortran
dset_id = H5Dcreate(loc_id, name, type_id, space_id)
INSIST(ctx%global_all(dset_id >= 0))
```

Retain collective success checks that prevent participants from continuing
after an operation failed on any rank. This does not require adding a reduction
after every local HDF5 or MPI call; local status checks can be sufficient for
local operations. Follow the operation's semantics and the surrounding error
handling. Treat an ordinary MPI error as fatal unless an explicit recovery path
has been designed; further collective communication may itself be impossible.

Do not classify every `global_all` as validation overhead. If its result changes
shared writer state or controls subsequent collective operations, it is part of
the algorithm. For example:

```fortran
geometry_valid = ctx%global_all(local_geometry_valid)
```

If later code uses `geometry_valid` to decide whether geometry must be supplied
again, the reduction is required in every build. It is distinct from reducing a
local argument condition solely to assert that all arguments are valid.

### Fatal failure and verification

A failed local contract must terminate the MPI job, so other ranks are not left
waiting in later collectives. The current
[`vtkhdf_assert`](../src/vtkhdf_assert.F90) reports the rank, source location,
and HDF5 error stack, then calls `MPI_Abort(MPI_COMM_WORLD,...)` in MPI builds;
serial builds use `error stop`. This assumes an active MPI runtime. Preserve
job-wide termination when changing the assertion machinery; cleanup requiring
collective participation does not belong in that failure path.

Test the behavior protected by a check, rather than duplicating its expression
in a test. Tests that expect an `ASSERT` failure belong only in builds where
assertions are enabled. MPI failure tests should have timeouts and exercise
specific diagnostics, without assuming arbitrary collective-call mismatches
will be diagnosed cleanly. Include empty-rank contributions when changing
validation around distributed writes, and check that production builds do not
depend on any computation removed with `ASSERT`.

## Validation

Run serial and MPI builds when changing shared writer or HDF5 code. Check the
installed-package consumers when changing public modules or exported build
requirements, and exercise the fpm source-generation path when changing source
lists or preprocessing.

For new formats or storage changes, use both direct HDF5 layout checks and VTK
readback. The layout checks verify shapes, offsets, and stored snapshot counts;
reader checks verify geometry and the interpreted point, cell, and field values.
Use a VTK version supporting the feature under test, and retain older-reader
coverage for supported formats. Skip only known unsupported cases, rather than
skipping a complete check that also covers supported types.

Include empty contributions, multiple regions, noncontiguous actual arguments,
nonzero logical bounds, temporal reuse, and changed rank ownership where they
apply. Region overlap and incomplete coverage are caller-contract violations;
the library does not attempt global geometric validation.

For example, a serial Debug build including examples can be checked with:

```sh
cmake -B build-debug -DCMAKE_BUILD_TYPE=Debug -DENABLE_MPI=OFF -DBUILD_EXAMPLES=ON
cmake --build build-debug --parallel
ctest --test-dir build-debug --output-on-failure
```

Repeat with `-DENABLE_MPI=ON` and the appropriate MPI compiler environment.
Compiler-sensitive changes should be validated with the supported compiler
families. Runtime checking changes also require a production build to establish
that correctness does not depend on debug-only expressions.

For documentation changes, use a strict Sphinx build:

```sh
sphinx-build -b html -n -W --keep-going doc/sphinx build-doc/html
```
