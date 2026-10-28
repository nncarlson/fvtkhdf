# Writer implementation

The UnstructuredGrid and structured cores use the same organization:
initialization and cleanup, geometry, static data, temporal registration,
writes by index, time stepping, and registry/storage helpers. Their Fypp
families generate separate cell and point procedures with the same naming
conventions. There is no runtime cell/point/field selector in these interfaces.

## Responsibilities

| Layer | Responsibilities |
| --- | --- |
| File wrappers | Public overloads, file/communicator lifetime, HDF5 error context, public handles and ownership checks |
| Mesh cores | Geometry, mode/lifecycle policy, mesh-specific cell/point storage and validation |
| `vtkhdf_data_support` | Common records, name/designation attributes, write-flag reset, group cleanup, chunk sizing |
| `vtkhdf_field_data_type` | FieldData records, static and temporal writes, registration, metadata, cleanup |
| `vtkhdf_named_record_type` | Common named records, external-name defaults, sanitization, and collision resolution for data fields and MultiBlock Assembly links |
| `vtkhdf_h5` | HDF5 dataset creation, append, and region operations |

Both cores contain `cell_data`, `point_data`, and `field_data`. Each association
has `vars`, `group_id`, and `offsets_id`; FieldData also has `sizes_id`.
`data_var_item` holds the common identity and temporal bookkeeping.
Mesh records extend it with UG's last local entity count or structured
dataset IDs and scalar/component interpretation. FieldData uses the same
record and list implementation in both cores.

## Creation and registration

Both mesh cores follow the same sequence:

1. Check the call's programming contracts and lifecycle.
2. Reserve the name and registry entry using `reserve_cell_data_name` or
   `reserve_point_data_name`. Each returns external and internal names and an
   integer index. Static writes reserve a fresh entry on every named call.
3. Create or write the mesh-specific payload dataset.
4. Call the common `write_data_name_attr` to preserve the external label as
   `Name` and any designation as `Attribute`.
5. Create applicable temporal metadata.

The core returns an index; the file wrapper constructs the public typed
handle. The handle's owner is the root group of its VTK dataset, including an
individual MultiBlock leaf. Wrappers validate ownership before passing an
index to the core. MultiBlock additionally checks its block index.

FieldData wrappers check lifecycle policy and delegate creation, registration,
and writes to `vtkhdf_field_data_type`. Static FieldData writes supply the complete
value. Temporal FieldData uses tuple rows, component/tuple sizes, and omission
reuse identically in both cores. Only rank 0 supplies the stored value.

## Time stepping and storage

Starting a step finalizes an open step, records its time, resets temporal
write flags, and opens the new step. Finalization checks required mesh data,
appends mesh-specific offsets and common FieldData metadata, and closes the
step. Flags are reset at the start of a step only.

Cell/point storage differs because the file layouts differ. UG appends entity
tuples and records tuple offsets. Structured cores write rectangular regions
into whole-grid snapshots and record snapshot indices. Multiple structured
block writes can therefore contribute to one field in a step. Named static
structured writes optionally return handles for additional blocks.
Structured static fields also have explicit per-step zero offsets; UG static
fields have no temporal offset dataset. These storage differences are local
to the mesh cores.

Name collisions are handled on a best-effort basis by choosing unique internal
names. They are not programming errors. Assertions enforce programming
contracts, such as handle ownership, lifecycle, and compatible block types.
