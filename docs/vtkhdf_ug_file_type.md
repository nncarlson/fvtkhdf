# VTKHDF_UG_FILE_TYPE

This module defines the `vtkhdf_ug_file` derived type for exporting mesh-based solution data to a VTKHDF file readable by ParaView. It produces an UnstructuredGrid (UG) dataset and supports both static and time-dependent cell and point data on a single static mesh.

## Type Definition

```fortran
type :: vtkhdf_ug_file
```

## Methods

### create
Initializes a new VTKHDF unstructured grid file.

**MPI Build:**
```fortran
call file%create(filename, comm, stat, errmsg, is_temporal)
```
- `filename` (character(*), intent(in)): The name of the file to create.
- `comm` (integer, intent(in)): MPI communicator.
- `stat` (integer, intent(out)): Status code (0 for success).
- `errmsg` (character(:), allocatable, intent(out)): Error message if `stat` is non-zero.
- `is_temporal` (logical, intent(in), optional): Set to `.true.` if the file will contain time-dependent data.

**Serial Build:**
```fortran
call file%create(filename, stat, errmsg, is_temporal)
```
- `filename` (character(*), intent(in)): The name of the file to create.
- `stat` (integer, intent(out)): Status code (0 for success).
- `errmsg` (character(:), allocatable, intent(out)): Error message if `stat` is non-zero.
- `is_temporal` (logical, intent(in), optional): Set to `.true.` if the file will contain time-dependent data.

---

### close
Cleanly closes HDF5 identifiers and the file, and default-initializes the object.
```fortran
call file%close()
```

---

### flush
Flushes the file to disk.
```fortran
call file%flush()
```

---

### write_mesh
Writes the unstructured mesh geometry and topology to the file. This procedure **must** be called before any data writing or registration procedures.

```fortran
call file%write_mesh(x, cnode, xcnode, types)
```
- `x` (real(real32/real64), intent(in)): Node coordinates, rank-2 array of shape `(3, npoints)`.
- `cnode` (integer(int32/int64), intent(in)): Connectivity array (1-based indices).
- `xcnode` (integer(int32/int64), intent(in)): Connectivity offset array (1-based indices).
- `types` (integer(int8), intent(in)): VTK cell types (e.g., `VTK_TETRA`, `VTK_HEXAHEDRON`).

---

### write_cell_data / write_point_data
Writes static cell-based or point-based data to a new named dataset.

```fortran
call file%write_cell_data(name, array)
call file%write_point_data(name, array)
```
- `name` (character(*), intent(in)): The name of the dataset.
- `array` (numeric, intent(in)): The data array. Can be scalar (rank-1) or vector (rank-2). For temporal files, this dataset is static and not associated with any time step.

---

### write_time_step
Marks the start of a new time step.

```fortran
call file%write_time_step(time)
```
- `time` (real(real64), intent(in)): The time value for the new step. Subsequent temporal data writes will be associated with this step.

---

### register_temporal_cell_data / register_temporal_point_data
Registers a name as a time-dependent dataset. This must be called before `write_time_step` and before writing any temporal data.

```fortran
call file%register_temporal_cell_data(name, mold)
call file%register_temporal_point_data(name, mold)
```
- `name` (character(*), intent(in)): The name of the dataset.
- `mold` (numeric, intent(in)): Defines the data type and shape (scalar or vector). Values are ignored.

---

### write_temporal_cell_data / write_temporal_point_data
Writes cell-based or point-based data for the current time step.

```fortran
call file%write_temporal_cell_data(name, array)
call file%write_temporal_point_data(name, array)
```
- `name` (character(*), intent(in)): The name of the registered temporal dataset.
- `array` (numeric, intent(in)): The data array for the current time step.

## Example Usage

```fortran
program example
  use, intrinsic :: iso_fortran_env, only: r8 => real64
  use vtkhdf_ug_file_type
  use vtkhdf_vtk_cell_types
  implicit none
  type(vtkhdf_ug_file) :: viz
  integer :: stat
  character(:), allocatable :: errmsg

  ! Mesh data...
  real(r8) :: x(3,8)
  integer  :: cnode(8), xcnode(2)
  integer(1) :: types(1)

  ! Initialize mesh data here...

  call viz%create('output.vtkhdf', stat, errmsg, is_temporal=.true.)
  if (stat /= 0) error stop errmsg

  call viz%write_mesh(x, cnode, xcnode, types)

  ! Register temporal data
  call viz%register_temporal_point_data('velocity', [0.0_r8, 0.0_r8, 0.0_r8])

  ! Time step 1
  call viz%write_time_step(0.0_r8)
  call viz%write_temporal_point_data('velocity', initial_velocity)

  ! Time step 2
  call viz%write_time_step(1.0_r8)
  call viz%write_temporal_point_data('velocity', next_velocity)

  call viz%close()
end program example
```
