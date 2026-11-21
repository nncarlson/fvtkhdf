# HyperTreeGrid material-interface API specification

Status: proposed API, agreed design from the PR #57 review discussion.
Date: 2026-11-16.

This specifies two type-bound methods on `vtkhdf_hypertree_file` for writing the
normals and intercepts arrays that describe a HyperTreeGrid material
interface. The writer manages their dataset names, root attributes, and
temporal handles. Implementation is separate from this specification.

## Public interface

```fortran
call file%write_interface_data(normals, intercepts)
call file%register_temporal_interface_data(mold)
```

Both are subroutines. Registration returns no public handle. The writer
contains the two underlying cell-data handles privately.

`write_interface_data` performs a static write if no interface has been
defined, or a temporal write if a temporal interface has been registered.
The method has no names, handles, attribute-role, or static/temporal selector
arguments; registration determines temporal behavior.

This API replaces the proposed `interface_normals_name` and
`interface_intercepts_name` arguments to `write_mesh`. Root-grid definition
returns to:

```fortran
call file%write_mesh(dims [, branch_factor] [, transposed_root_indexing])
```

There is no separate normals-only or intercepts-only write method.

## Arguments and numeric kinds

### `write_interface_data(normals, intercepts)`

- Both arguments are input rank-2 arrays of shape `(3,ncell)`, where `ncell`
  is the calling rank's current number of HTG cells.
- Both are real arrays of the same kind: either `real32` or `real64` from
  `iso_fortran_env`.
- The first dimension indexes components. The second indexes every local
  cell, refined or leaf, in the same tree-by-tree, breadth-first order as
  ordinary HTG cell data.
- The supplied kind is preserved in storage. These arrays are not promoted
  to `real64`.
- In temporal use, their kind must match the registered mold's kind. A
  `real64` write to an interface registered as `real32`, or the reverse, is
  an error.
- Their kind need not match the kind of coordinate arguments or coordinate
  storage.
- Noncontiguous actual arrays follow the same supported behavior as
  ordinary cell-data writes.

The component count is always three. Integer arrays, mixed-kind pairs,
flattened arrays, and different component counts are not supported by this
API. The writer stores the supplied values; it does not compute the
interface or normalize the normals. Interpretation of the values follows
VTK's HTG interface-array conventions and is the caller's responsibility.

### `register_temporal_interface_data(mold)`

- `mold` is an input scalar of kind `real32` or `real64`.
- Its type and kind select the storage type of both interface arrays.
- Its value is never referenced.
- A scalar is sufficient because the three-component layout is prescribed.

Examples:

```fortran
call file%register_temporal_interface_data(0.0_real32)
call file%register_temporal_interface_data(0.0_real64)
```

These illustrate alternatives: only one registration is permitted for a
given open writer.

## Names and VTKHDF representation

The conventional dataset base names are:

| Array | Base name in `CellData` |
| --- | --- |
| Normals | `vtkInterfaceNormals` |
| Intercepts | `vtkInterfaceIntercepts` |

The caller does not select these names. The writer reserves them using the
existing cell-data naming machinery. Collisions with ordinary cell arrays
are resolved by its normal suffix rules, independently for each array.
The names are not globally reserved: ordinary cell data can use these
labels before or after interface definition without replacing the interface
datasets.

When the pair is defined, the writer writes these root string attributes:

| Root attribute | Value |
| --- | --- |
| `InterfaceNormalsName` | Actual stored name of the normals dataset |
| `InterfaceInterceptsName` | Actual stored name of the intercepts dataset |

For example, if an ordinary cell array already owns `vtkInterfaceNormals`, the
interface normals might be stored as `vtkInterfaceNormals_1`, and
`InterfaceNormalsName` must contain `vtkInterfaceNormals_1`.

Both datasets receive the same naming metadata as ordinary cell datasets.
The special interface attributes do not imply an ordinary active-array role
such as `Normals` or `Vectors`.

Static data uses the existing static cell-data layout. Temporal registration
creates two three-component temporal cell datasets and their normal
`Steps/CellDataOffsets` datasets. Chunk sizing and other storage policy
follow the existing cell-data implementation.

The root attributes are written at static definition or temporal
registration. Their values remain fixed for the open writer's lifetime;
temporal writes update the referenced datasets, not these names. A writer
that never defines an interface does not acquire either attribute or either
interface dataset.

## Lifecycle and singleton rule

An open writer has one of three interface states:

| State | Meaning |
| --- | --- |
| Undefined | No interface pair has been defined |
| Static | A static pair has been written |
| Temporal | A temporal pair has been registered |

The only definition transitions are:

| Call | Required state | Result |
| --- | --- | --- |
| First static `write_interface_data` | Undefined | Write both arrays and attributes; enter Static |
| `register_temporal_interface_data` | Undefined | Create both temporal datasets and attributes; enter Temporal |

**Definition is permitted at most once across both operations combined.**
In particular:

- A second static write is an error, even if its values or kind differ.
- A second registration is an error, even if its mold is identical.
- Registration after a static write is an error.
- A temporal interface cannot be changed to static or renamed.

This is the additional lifecycle rule relative to ordinary by-name cell-data
writes, which create another disambiguated dataset on each call. Repeated
interface definitions must not silently create a new pair or repoint the
attributes.

Once in Temporal state, `write_interface_data` supplies values to the same
registered pair, subject to the temporal rules below. It never creates a
new pair.

Explicit `close` follows the ordinary writer cleanup and open-step
finalization rules. It clears the interface state and private handles.
Creating a new file with the same Fortran object begins with an Undefined
interface state.

## Writer modes

| Writer mode | Static interface definition | Temporal interface registration |
| --- | --- | --- |
| `VTKHDF_STATIC` | Allowed | Not allowed |
| `VTKHDF_FIXED_MESH` | Allowed | Allowed |
| `VTKHDF_MOVING_MESH` | Allowed | Allowed |
| `VTKHDF_DYNAMIC_MESH` | Not allowed | Allowed |

These are the existing HTG cell-data mode restrictions. The interface API
adds no different mesh policy.

### Static writes

A static write requires an open writer and a complete mesh: root-grid
definition, geometry, and topology must already have been supplied.

In a temporal writer that permits static cell data, the static pair may be
introduced at the same points in the lifecycle as ordinary static cell
data, including during or between steps. It applies to all steps, including
earlier steps, under the ordinary static cell-data behavior.

Writing the pair establishes the same mesh-change barrier as any other
cell-data write. In particular, after writing it during an open step, the
mesh cannot be changed again until the next step.

### Temporal registration and writes

Registration requires a temporal writer and must precede the first
`start_time_step` call. It does not require a mesh to have been written and
does not establish a mesh-change barrier.

After registration:

1. A write requires an open time step and a complete current mesh.
2. Both arrays must be written in the first time step.
3. The pair may be written at most once in each step.
4. Writing the pair establishes the same mesh-change barrier as an ordinary
   temporal cell-data write.
5. Omitting the write in a later step reuses both arrays' most recently
   written values. Reuse requires that every rank's local cell count remain
   compatible with that write, under the existing cell-data reuse checks.
6. If cell counts change incompatibly, both arrays must be rewritten before
   the step is finalized.

Equal cell counts do not establish that values still describe the same
cells. After migration, refinement changes, or cell reordering, the caller
is responsible for rewriting the pair whenever the old values no longer
correspond to the current cell order.

Both arrays are written or reused together. Independent updates, mixed
static/temporal treatment of the two arrays, and partial writes are not
provided by this convenience API. This is an API choice, not a restriction
of the VTKHDF format. Supplying unchanged normals along with updated
intercepts is valid.

## MPI behavior and validation

The methods are collective over the writer's communicator and follow the
same collective-call contract as ordinary cell data:

- All ranks make the same registration/write calls in the same order and
  choose the same numeric kind.
- Registration mold values may differ because they are ignored.
- Values and local cell counts may differ between ranks.
- A rank with no cells supplies arrays of shape `(3,0)` and still
  participates in the calls.
- The singleton state and static/temporal choice must agree across ranks.

The implementation must check singleton state, mode, registration timing,
open-step requirements, kind, and array shapes as appropriate. Both arrays'
structural arguments should be validated before either payload is written.

Local contract violations use the existing unconditional assertion/fatal
error policy. Cross-rank agreement diagnostics follow the existing Debug
checking policy. These subroutines add no recoverable status return and no
transactional guarantee for failures during HDF5 operations.

The private temporal handles must participate in the ordinary registry,
step-flag reset, reuse validation, offset emission, and cleanup paths. The
absence of a public handle must not bypass those checks.

## Usage examples

The fragments below use serial `create` signatures. MPI callers additionally
supply their communicator and participate collectively. Array construction
and ordinary mesh calls are abbreviated.

### Static interface

```fortran
use, intrinsic :: iso_fortran_env, only: real32
use vtkhdf_hypertree_file_type

! file, stat, errmsg, dims, normals, and intercepts are declared by the caller.
call file%create('static.vtkhdf', stat, errmsg, mode=VTKHDF_STATIC)
if (stat /= 0) error stop errmsg
call file%write_mesh(dims)
! Write geometry and topology; populate real32 arrays of shape (3,ncell).
call file%write_interface_data(normals, intercepts)
call file%close()
```

### Temporal interface

```fortran
use, intrinsic :: iso_fortran_env, only: real64
use vtkhdf_hypertree_file_type

call file%create('dynamic.vtkhdf', stat, errmsg, mode=VTKHDF_DYNAMIC_MESH)
if (stat /= 0) error stop errmsg
call file%register_temporal_interface_data(0.0_real64)
call file%write_mesh(dims)

do step = 1, nsteps
  call file%start_time_step(times(step))
  ! Supply initial mesh, or replacements when needed.
  ! Populate real64 normals and intercepts of shape (3,current_ncell).
  ! The first step requires a write. Later omission reuses the whole pair.
  if (step == 1 .or. interface_changed) then
    call file%write_interface_data(normals, intercepts)
  end if
  call file%finalize_time_step()
end do
call file%close()
```

The caller's `interface_changed` decision must include any mesh changes
that make reuse inappropriate.

## Acceptance checks for implementation

- Static and temporal pairs work with both real kinds. Layout checks verify
  that storage preserves the chosen kind, and VTK readback verifies names,
  enabled interfaces, and every component value.
- Collisions with ordinary cell-data names before and after definition are
  disambiguated without changing which arrays the attributes identify.
- No interface datasets or attributes appear when this API is unused.
- Repeated static definition, repeated registration, and registration after
  static definition are rejected. After registration, writes always take
  the temporal path.
- Wrong mode, late registration, writes outside a step, missing initial
  data, repeated writes within a step, invalid shapes, and kind mismatches
  are rejected through the intended contract checks.
- Paired omission reuses both offset histories; paired replacement appends
  both arrays. Changed cell counts require replacement. Root attribute
  values remain unchanged.
- Two- and four-rank runs cover empty contributions, uneven ownership,
  noncontiguous arrays, and topology changes using the normal HTG fixtures.
- Closing and creating a new file with the same object resets the singleton
  state. No new public handle type is introduced.

Reading, XML output, compression, and inclusion in composite writers are
outside this API change.
