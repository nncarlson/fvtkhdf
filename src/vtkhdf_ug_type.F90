!!
!! VTKHDF_UG_TYPE
!!
!! This module defines an auxiliary derived type for managing an HDF5 tree that
!! stores an UnstructuredGrid VTKHDF dataset. It serves as a core component of
!! the VTKHDF_FILE type.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! January 2026
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! The up-to-date specification for VTKHDF is at
!! https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html
!!
!! This module was written for version 2.5 of the format.
!!
!! This implementation supports static and time dependent data (point and
!! cell), but both assume a single static mesh.
!!

#include "f90_assert.fpp"

module vtkhdf_ug_type

  use,intrinsic :: iso_fortran_env
  use vtkhdf_h5_c_binding
  use vtkhdf_h5
  use mpi
  implicit none
  private

  type, public :: vtkhdf_ug
    integer :: comm = MPI_COMM_NULL
    integer(hid_t) :: root_id=0, cgrp_id=0, pgrp_id=0
    integer :: nproc, npart
    integer :: nnode, ncell, nnode_tot, ncell_tot
    logical :: temporal = .false.
    integer :: nsteps = -1
    integer(hid_t) :: steps_id=0, cogrp_id=0, pogrp_id=0
    type(temporal_data), pointer :: temporal_point_data => null()
    type(temporal_data), pointer :: temporal_cell_data => null()
  contains
    procedure :: init
    procedure :: write_mesh
    procedure :: write_time_step
    procedure :: write_cell_data_real32,  write_cell_data_real64
    procedure :: write_cell_data_int32,   write_cell_data_int64
    procedure :: write_point_data_real32, write_point_data_real64
    procedure :: write_point_data_int32,  write_point_data_int64
    procedure :: register_temporal_cell_data_real32,  register_temporal_cell_data_real64
    procedure :: register_temporal_cell_data_int32,   register_temporal_cell_data_int64
    procedure :: register_temporal_point_data_real32, register_temporal_point_data_real64
    procedure :: register_temporal_point_data_int32,  register_temporal_point_data_int64
    procedure :: write_temporal_cell_data_real32,  write_temporal_cell_data_real64
    procedure :: write_temporal_cell_data_int32,   write_temporal_cell_data_int64
    procedure :: write_temporal_point_data_real32, write_temporal_point_data_real64
    procedure :: write_temporal_point_data_int32,  write_temporal_point_data_int64
    procedure :: close
  end type

  type :: temporal_data
    character(:), allocatable :: name
    integer :: next_offset = 0
    logical :: flag = .false.
    type(temporal_data), pointer :: next => null()
  contains
    final :: temporal_data_delete
  end type

  integer, parameter :: vtkhdf_version(*) = [2,5]

contains

  !! Create a new VTKHDF UnstructuredGrid type group. To create a group
  !! supporting time-dependent data, specify the optional argument
  !! TEMPORAL to true; otherwise a static group is created by default.

  subroutine init(this, loc_id, name, comm, stat, errmsg, temporal)

    class(vtkhdf_ug), intent(inout) :: this
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer, intent(in) :: comm
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    logical, intent(in), optional :: temporal

    integer :: ierr

    this%comm = comm
    call MPI_Comm_size(this%comm, this%nproc, ierr)
    INSIST(ierr == MPI_SUCCESS)

    this%root_id = H5Gcreate(loc_id, name)
    if (global_any(this%root_id < 0, this%comm)) then
      stat = 1
      errmsg = 'unable to create group'
      return
    end if

    call h5_write_attr(this%root_id, 'Version', vtkhdf_version, this%comm, stat, errmsg)
    if (stat /= 0) return
    call h5_write_attr(this%root_id, 'Type', 'UnstructuredGrid', this%comm, stat, errmsg)
    if (stat /= 0) return

    this%cgrp_id = H5Gcreate(this%root_id, 'CellData')
    if (global_any(this%cgrp_id < 0, this%comm)) then
      stat = 1
      errmsg = 'unable to create "CellData" group'
      return
    end if

    this%pgrp_id = H5Gcreate(this%root_id, 'PointData')
    if (global_any(this%pgrp_id < 0, this%comm)) then
      stat = 1
      errmsg = 'unable to create "PointData" group'
      return
    end if

    if (present(temporal)) this%temporal = temporal
    if (this%temporal) then

      this%steps_id = H5Gcreate(this%root_id, 'Steps')
      if (global_any(this%steps_id < 0, this%comm)) then
        stat = 1
        errmsg = 'unable to create "Steps" group'
        return
      end if

      this%cogrp_id = H5Gcreate(this%steps_id, 'CellDataOffsets')
      if (global_any(this%cogrp_id < 0, this%comm)) then
        stat = 1
        errmsg = 'unable to create "CellDataOffsets" group'
        return
      end if

      this%pogrp_id = H5Gcreate(this%steps_id, 'PointDataOffsets')
      if (global_any(this%pogrp_id < 0, this%comm)) then
        stat = 1
        errmsg = 'unable to create "PointDataOffsets" group'
        return
      end if

      this%nsteps = 0
      call h5_write_attr(this%steps_id, 'NSteps', this%nsteps, this%comm, stat, errmsg)
      if (stat /= 0) return

      associate (imold => [1], rmold => [1.0_real64], chunk_size => 100)
        call h5_create_unlimited_dataset(this%steps_id, 'Values', rmold, chunk_size, this%comm, stat, errmsg)
        if (stat /= 0) return
        call h5_create_unlimited_dataset(this%steps_id, 'PointOffsets', imold, chunk_size, this%comm, stat, errmsg)
        if (stat /= 0) return
        call h5_create_unlimited_dataset(this%steps_id, 'CellOffsets', imold, chunk_size, this%comm, stat, errmsg)
        if (stat /= 0) return
        call h5_create_unlimited_dataset(this%steps_id, 'ConnectivityIdOffsets', imold, chunk_size, this%comm, stat, errmsg)
        if (stat /= 0) return
        call h5_create_unlimited_dataset(this%steps_id, 'NumberOfParts', imold, chunk_size, this%comm, stat, errmsg)
        if (stat /= 0) return
        call h5_create_unlimited_dataset(this%steps_id, 'PartOffsets', imold, chunk_size, this%comm, stat, errmsg)
        if (stat /= 0) return
      end associate

    end if

    stat = 0

  end subroutine init

  subroutine close(this)
    class(vtkhdf_ug), intent(inout) :: this
    integer :: ierr
    if (this%cogrp_id > 0) ierr = H5Gclose(this%cogrp_id)
    if (this%pogrp_id > 0) ierr = H5Gclose(this%pogrp_id)
    if (this%steps_id > 0) ierr = H5Gclose(this%steps_id)
    if (this%cgrp_id > 0) ierr = H5Gclose(this%cgrp_id)
    if (this%pgrp_id > 0) ierr = H5Gclose(this%pgrp_id)
    if (this%root_id > 0) ierr = H5Gclose(this%root_id)
    if (associated(this%temporal_cell_data)) deallocate(this%temporal_cell_data)
    if (associated(this%temporal_point_data)) deallocate(this%temporal_point_data)
    call default_initialize(this)
  contains
    subroutine default_initialize(this)
      type(vtkhdf_ug), intent(out) :: this
    end subroutine
  end subroutine

  recursive subroutine temporal_data_delete(this)
    type(temporal_data), intent(inout) :: this
    if (associated(this%next)) deallocate(this%next)
  end subroutine

  !! Write the unstructured mesh to the file. The mesh is described in the
  !! conventional manner by the X, CNODE, and XCNODE arrays. The additional
  !! array TYPES provides an unambiguous specification of the cell types,
  !! that we would otherwise infer from the OFFSETS array. This procedure
  !! must be called before any of the following procedures.

  !! NB: Paraview and/or the VTKHDF reader currently has a problem handling
  !! 0-sized parts (see https://gitlab.kitware.com/vtk/vtk/-/issues/19923).
  !! As a workaround, we omit any empty part by not contributing anything
  !! at all to the NumberOf* and Offsets datasets instead of the 0 we
  !! normally would. Other datasets like Connectivity already contribute
  !! nothing for empty parts.

  subroutine write_mesh(this, x, cnode, xcnode, types, stat, errmsg)

    class(vtkhdf_ug), intent(inout) :: this
    real(real64), intent(in) :: x(:,:)
    integer, intent(in) :: cnode(:), xcnode(:)
    integer(int8), intent(in) :: types(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: ierr, idum(0)

    INSIST(this%root_id > 0)

    this%nnode = size(x, dim=2)
    this%ncell = size(types)

    ASSERT(size(xcnode) == this%ncell+1)
    ASSERT(size(cnode) == xcnode(size(xcnode))-1)
    ASSERT(minval(cnode) >= 1)
    ASSERT(maxval(cnode) <= this%nnode)

    call MPI_Allreduce(this%nnode, this%nnode_tot, 1, MPI_INTEGER, MPI_SUM, this%comm, ierr)
    INSIST(ierr == MPI_SUCCESS)
    call MPI_Allreduce(this%ncell, this%ncell_tot, 1, MPI_INTEGER, MPI_SUM, this%comm, ierr)
    INSIST(ierr == MPI_SUCCESS)
    call MPI_Allreduce(merge(1, 0, this%ncell>0), this%npart, 1, MPI_INTEGER, MPI_SUM, this%comm, ierr) ! see NB above
    INSIST(ierr == MPI_SUCCESS)

    if (this%ncell > 0) then
      call h5_write_dataset(this%root_id, 'NumberOfPoints', this%nnode, this%comm, stat, errmsg)
    else ! see NB above
      call h5_write_dataset(this%root_id, 'NumberOfPoints', idum, this%comm, stat, errmsg)
    end if
    if (stat /= 0) return

    if (this%ncell > 0) then
      call h5_write_dataset(this%root_id, 'NumberOfCells',  this%ncell, this%comm, stat, errmsg)
    else ! see NB above
      call h5_write_dataset(this%root_id, 'NumberOfCells',  idum, this%comm, stat, errmsg)
    end if
    if (stat /= 0) return

    if (this%ncell > 0) then
      call h5_write_dataset(this%root_id, 'NumberOfConnectivityIds', size(cnode), this%comm, stat, errmsg)
    else ! see NB above
      call h5_write_dataset(this%root_id, 'NumberOfConnectivityIds', idum, this%comm, stat, errmsg)
    end if
    if (stat /= 0) return

    if (this%ncell > 0) then
      call h5_write_dataset(this%root_id, 'Offsets', xcnode-1, this%comm, stat, errmsg) ! offsets instead of starting indices
    else ! see NB above
      call h5_write_dataset(this%root_id, 'Offsets', idum, this%comm, stat, errmsg)
    end if
    if (stat /= 0) return

    if (this%ncell > 0) then
      call h5_write_dataset(this%root_id, 'Connectivity', cnode-1, this%comm, stat, errmsg)  ! 0-based indexing
    else ! workaround for gfortran bug https://gcc.gnu.org/bugzilla/show_bug.cgi?id=123899
      call h5_write_dataset(this%root_id, 'Connectivity', idum, this%comm, stat, errmsg)
    endif
    if (stat /= 0) return

    call h5_write_dataset(this%root_id, 'Types', types, this%comm, stat, errmsg)
    if (stat /= 0) return

    call h5_write_dataset(this%root_id, 'Points', x, this%comm, stat, errmsg)
    if (stat /= 0) return

  end subroutine

  !! Writes the cell DATA to a new CellData group dataset NAME. Scalar,
  !! vector, and tensor cell-based data are supported. In the case of a
  !! temporal file supporting time-dependent datasets, this dataset is
  !! static and not associated with any time step.

  subroutine write_cell_data_real32(this, name, data, stat, errmsg)
    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, allocatable :: dims(:)
    dims = shape(data)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%ncell)
    call h5_write_dataset(this%cgrp_id, name, data, this%comm, stat, errmsg)
  end subroutine

  subroutine write_cell_data_real64(this, name, data, stat, errmsg)
    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, allocatable :: dims(:)
    dims = shape(data)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%ncell)
    call h5_write_dataset(this%cgrp_id, name, data, this%comm, stat, errmsg)
  end subroutine

  subroutine write_cell_data_int32(this, name, data, stat, errmsg)
    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, allocatable :: dims(:)
    dims = shape(data)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%ncell)
    call h5_write_dataset(this%cgrp_id, name, data, this%comm, stat, errmsg)
  end subroutine

  subroutine write_cell_data_int64(this, name, data, stat, errmsg)
    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, allocatable :: dims(:)
    dims = shape(data)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%ncell)
    call h5_write_dataset(this%cgrp_id, name, data, this%comm, stat, errmsg)
  end subroutine

  !! Writes the point DATA to a new PointData group dataset NAME. Scalar,
  !! vector, and tensor point-based data are supported. In the case of a
  !! temporal file supporting time-dependent datasets, this dataset is
  !! static and not associated with any time step.

  subroutine write_point_data_real32(this, name, data, stat, errmsg)
    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, allocatable :: dims(:)
    dims = shape(data)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%nnode)
    call h5_write_dataset(this%pgrp_id, name, data, this%comm, stat, errmsg)
  end subroutine

  subroutine write_point_data_real64(this, name, data, stat, errmsg)
    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, allocatable :: dims(:)
    dims = shape(data)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%nnode)
    call h5_write_dataset(this%pgrp_id, name, data, this%comm, stat, errmsg)
  end subroutine

  subroutine write_point_data_int32(this, name, data, stat, errmsg)
    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, allocatable :: dims(:)
    dims = shape(data)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%nnode)
    call h5_write_dataset(this%pgrp_id, name, data, this%comm, stat, errmsg)
  end subroutine

  subroutine write_point_data_int64(this, name, data, stat, errmsg)
    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, allocatable :: dims(:)
    dims = shape(data)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%nnode)
    call h5_write_dataset(this%pgrp_id, name, data, this%comm, stat, errmsg)
  end subroutine

  !! Register the specified NAME as a time-dependent CellData group dataset.
  !! This writes no data, but only configures some necessary metadata.
  !! The MOLD array argument shall have the same type, kind, and rank as the
  !! actual cell data, and the same extent in all but the last dimension,
  !! but the array values themselves are not accessed. Scalar, vector, and
  !! tensor-valued cell data are supported (rank-1, 2, and 3 MOLD).

  subroutine register_temporal_cell_data_real32(this, name, mold, stat, errmsg)

    class(vtkhdf_ug), intent(inout) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: mold(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(temporal_data), pointer :: new

    INSIST(this%nsteps == 0)

    !! Dataset for the cell data
    associate (chunk_size => this%ncell_tot)
      call h5_create_unlimited_dataset(this%cgrp_id, name, mold, chunk_size, this%comm, stat, errmsg)
    end associate
    if (stat /= 0) return

    !! Dataset for its time step offsets
    associate (mold => [1], chunk_size => 100)
      call h5_create_unlimited_dataset(this%cogrp_id, name, mold, chunk_size, this%comm, stat, errmsg)
    end associate
    if (stat /= 0) return

    !! Add to list CellData datasets
    allocate(new)
    new%name = name
    new%next => this%temporal_cell_data
    this%temporal_cell_data => new

  end subroutine register_temporal_cell_data_real32

  subroutine register_temporal_cell_data_real64(this, name, mold, stat, errmsg)

    class(vtkhdf_ug), intent(inout) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: mold(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(temporal_data), pointer :: new

    INSIST(this%nsteps == 0)

    !! Dataset for the cell data
    associate (chunk_size => this%ncell_tot)
      call h5_create_unlimited_dataset(this%cgrp_id, name, mold, chunk_size, this%comm, stat, errmsg)
    end associate
    if (stat /= 0) return

    !! Dataset for its time step offsets
    associate (mold => [1], chunk_size => 100)
      call h5_create_unlimited_dataset(this%cogrp_id, name, mold, chunk_size, this%comm, stat, errmsg)
    end associate
    if (stat /= 0) return

    !! Add to list CellData datasets
    allocate(new)
    new%name = name
    new%next => this%temporal_cell_data
    this%temporal_cell_data => new

  end subroutine register_temporal_cell_data_real64

  subroutine register_temporal_cell_data_int32(this, name, mold, stat, errmsg)

    class(vtkhdf_ug), intent(inout) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: mold(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(temporal_data), pointer :: new

    INSIST(this%nsteps == 0)

    !! Dataset for the cell data
    associate (chunk_size => this%ncell_tot)
      call h5_create_unlimited_dataset(this%cgrp_id, name, mold, chunk_size, this%comm, stat, errmsg)
    end associate
    if (stat /= 0) return

    !! Dataset for its time step offsets
    associate (mold => [1], chunk_size => 100)
      call h5_create_unlimited_dataset(this%cogrp_id, name, mold, chunk_size, this%comm, stat, errmsg)
    end associate
    if (stat /= 0) return

    !! Add to list CellData datasets
    allocate(new)
    new%name = name
    new%next => this%temporal_cell_data
    this%temporal_cell_data => new

  end subroutine register_temporal_cell_data_int32

  subroutine register_temporal_cell_data_int64(this, name, mold, stat, errmsg)

    class(vtkhdf_ug), intent(inout) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: mold(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(temporal_data), pointer :: new

    INSIST(this%nsteps == 0)

    !! Dataset for the cell data
    associate (chunk_size => this%ncell_tot)
      call h5_create_unlimited_dataset(this%cgrp_id, name, mold, chunk_size, this%comm, stat, errmsg)
    end associate
    if (stat /= 0) return

    !! Dataset for its time step offsets
    associate (mold => [1], chunk_size => 100)
      call h5_create_unlimited_dataset(this%cogrp_id, name, mold, chunk_size, this%comm, stat, errmsg)
    end associate
    if (stat /= 0) return

    !! Add to list CellData datasets
    allocate(new)
    new%name = name
    new%next => this%temporal_cell_data
    this%temporal_cell_data => new

  end subroutine register_temporal_cell_data_int64

  !! Register the specified NAME as a time-dependent PointData group dataset.
  !! This writes no data, but only configures some necessary internal metadata.
  !! The MOLD array argument shall have the same type, kind, and rank as the
  !! actual point data, and the same extent in all but the last dimension,
  !! but the array values themselves are not accessed. Scalar, vector, and
  !! tensor-valued point data are supported (rank-1, 2, and 3 MOLD).

  subroutine register_temporal_point_data_real32(this, name, mold, stat, errmsg)

    class(vtkhdf_ug), intent(inout) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: mold(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(temporal_data), pointer :: new

    INSIST(this%nsteps == 0)

    !! Dataset for the point data
    associate (chunk_size => this%nnode_tot)
      call h5_create_unlimited_dataset(this%pgrp_id, name, mold, chunk_size, this%comm, stat, errmsg)
    end associate
    if (stat /= 0) return

    !! Dataset for its time step offsets
    associate (mold => [1], chunk_size => 100)
      call h5_create_unlimited_dataset(this%pogrp_id, name, mold, chunk_size, this%comm, stat, errmsg)
    end associate
    if (stat /= 0) return

    !! Add to list PointData datasets
    allocate(new)
    new%name = name
    new%next => this%temporal_point_data
    this%temporal_point_data => new

  end subroutine register_temporal_point_data_real32

  subroutine register_temporal_point_data_real64(this, name, mold, stat, errmsg)

    class(vtkhdf_ug), intent(inout) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: mold(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(temporal_data), pointer :: new

    INSIST(this%nsteps == 0)

    !! Dataset for the point data
    associate (chunk_size => this%nnode_tot)
      call h5_create_unlimited_dataset(this%pgrp_id, name, mold, chunk_size, this%comm, stat, errmsg)
    end associate
    if (stat /= 0) return

    !! Dataset for its time step offsets
    associate (mold => [1], chunk_size => 100)
      call h5_create_unlimited_dataset(this%pogrp_id, name, mold, chunk_size, this%comm, stat, errmsg)
    end associate
    if (stat /= 0) return

    !! Add to list PointData datasets
    allocate(new)
    new%name = name
    new%next => this%temporal_point_data
    this%temporal_point_data => new

  end subroutine register_temporal_point_data_real64

  subroutine register_temporal_point_data_int32(this, name, mold, stat, errmsg)

    class(vtkhdf_ug), intent(inout) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: mold(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(temporal_data), pointer :: new

    INSIST(this%nsteps == 0)

    !! Dataset for the point data
    associate (chunk_size => this%nnode_tot)
      call h5_create_unlimited_dataset(this%pgrp_id, name, mold, chunk_size, this%comm, stat, errmsg)
    end associate
    if (stat /= 0) return

    !! Dataset for its time step offsets
    associate (mold => [1], chunk_size => 100)
      call h5_create_unlimited_dataset(this%pogrp_id, name, mold, chunk_size, this%comm, stat, errmsg)
    end associate
    if (stat /= 0) return

    !! Add to list PointData datasets
    allocate(new)
    new%name = name
    new%next => this%temporal_point_data
    this%temporal_point_data => new

  end subroutine register_temporal_point_data_int32

  subroutine register_temporal_point_data_int64(this, name, mold, stat, errmsg)

    class(vtkhdf_ug), intent(inout) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: mold(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(temporal_data), pointer :: new

    INSIST(this%nsteps == 0)

    !! Dataset for the point data
    associate (chunk_size => this%nnode_tot)
      call h5_create_unlimited_dataset(this%pgrp_id, name, mold, chunk_size, this%comm, stat, errmsg)
    end associate
    if (stat /= 0) return

    !! Dataset for its time step offsets
    associate (mold => [1], chunk_size => 100)
      call h5_create_unlimited_dataset(this%pogrp_id, name, mold, chunk_size, this%comm, stat, errmsg)
    end associate
    if (stat /= 0) return

    !! Add to list PointData datasets
    allocate(new)
    new%name = name
    new%next => this%temporal_point_data
    this%temporal_point_data => new

  end subroutine register_temporal_point_data_int64

  !! Mark the start of a new time step with time value TIME. Subsequent output
  !! of time-dependent datasets will be associated with this time step.

  subroutine write_time_step(this, time)

    class(vtkhdf_ug), intent(inout) :: this
    real(real64), intent(in) :: time

    type(temporal_data), pointer :: tmp
    integer :: stat
    character(:), allocatable :: errmsg

    INSIST(this%nsteps >= 0)

    this%nsteps = this%nsteps + 1
    call h5_write_attr(this%steps_id, 'NSteps', this%nsteps, this%comm, stat, errmsg)
    INSIST(stat == 0)

    call h5_append_to_dataset(this%steps_id, 'Values', time, this%comm, stat, errmsg, root=0)
    INSIST(stat == 0)

    !! A single mesh is used for all time steps so there are no offsets
    call h5_append_to_dataset(this%steps_id, 'PointOffsets', 0, this%comm, stat, errmsg, root=0)
    INSIST(stat == 0)
    call h5_append_to_dataset(this%steps_id, 'CellOffsets', 0, this%comm, stat, errmsg, root=0)
    INSIST(stat == 0)
    call h5_append_to_dataset(this%steps_id, 'ConnectivityIdOffsets', 0, this%comm, stat, errmsg, root=0)
    INSIST(stat == 0)
    call h5_append_to_dataset(this%steps_id, 'NumberOfParts', this%npart, this%comm, stat, errmsg, root=0)
    INSIST(stat == 0)
    call h5_append_to_dataset(this%steps_id, 'PartOffsets', 0, this%comm, stat, errmsg, root=0)
    INSIST(stat == 0)

    !! Cell and point data offsets for each of the temporal datasets for this timestep

    tmp => this%temporal_point_data
    do while (associated(tmp))
      tmp%flag = .false.  ! dataset not yet written for this time step
      call h5_append_to_dataset(this%pogrp_id, tmp%name, tmp%next_offset, this%comm, stat, errmsg, root=0)
      INSIST(stat == 0)
      tmp => tmp%next
    end do

    tmp => this%temporal_cell_data
    do while (associated(tmp))
      tmp%flag = .false.  ! dataset not yet written for this time step
      call h5_append_to_dataset(this%cogrp_id, tmp%name, tmp%next_offset, this%comm, stat, errmsg, root=0)
      INSIST(stat == 0)
      tmp => tmp%next
    end do

  end subroutine write_time_step

  !! Write the cell DATA to the time-dependent CellData group dataset NAME.
  !! The data is associated with the current time step. The type and shape
  !! of DATA must be consistent with the MOLD argument passed to the call
  !! of REGISTER_TEMPORAL_CELL_DATA that created this dataset.

  subroutine write_temporal_cell_data_real32(this, name, data, stat, errmsg)

    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg

    integer, allocatable :: dims(:)
    type(temporal_data), pointer :: dset

    INSIST(this%nsteps >= 0)

    dims = shape(data)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%ncell)

    dset => this%temporal_cell_data
    do while (associated(dset))
      if (dset%name == name) exit
      dset => dset%next
    end do

    if (associated(dset)) then
      if (dset%flag) then
        stat = 1
        errmsg = 'dataset "' // name // '" already written for this time step'
        return
      else
        call h5_append_to_dataset(this%cgrp_id, name, data, this%comm, stat, errmsg)
        if (stat /= 0) return
        dset%next_offset = dset%next_offset + this%ncell_tot
        dset%flag = .true. ! dataset has been written for this time step
      end if
    else
      stat = 1
      errmsg = '"' // name // '" is not a valid temporal cell dataset'
      return
    end if

  end subroutine write_temporal_cell_data_real32

  subroutine write_temporal_cell_data_real64(this, name, data, stat, errmsg)

    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg

    integer, allocatable :: dims(:)
    type(temporal_data), pointer :: dset

    INSIST(this%nsteps >= 0)

    dims = shape(data)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%ncell)

    dset => this%temporal_cell_data
    do while (associated(dset))
      if (dset%name == name) exit
      dset => dset%next
    end do

    if (associated(dset)) then
      if (dset%flag) then
        stat = 1
        errmsg = 'dataset "' // name // '" already written for this time step'
        return
      else
        call h5_append_to_dataset(this%cgrp_id, name, data, this%comm, stat, errmsg)
        if (stat /= 0) return
        dset%next_offset = dset%next_offset + this%ncell_tot
        dset%flag = .true. ! dataset has been written for this time step
      end if
    else
      stat = 1
      errmsg = '"' // name // '" is not a valid temporal cell dataset'
      return
    end if

  end subroutine write_temporal_cell_data_real64

  subroutine write_temporal_cell_data_int32(this, name, data, stat, errmsg)

    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg

    integer, allocatable :: dims(:)
    type(temporal_data), pointer :: dset

    INSIST(this%nsteps >= 0)

    dims = shape(data)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%ncell)

    dset => this%temporal_cell_data
    do while (associated(dset))
      if (dset%name == name) exit
      dset => dset%next
    end do

    if (associated(dset)) then
      if (dset%flag) then
        stat = 1
        errmsg = 'dataset "' // name // '" already written for this time step'
        return
      else
        call h5_append_to_dataset(this%cgrp_id, name, data, this%comm, stat, errmsg)
        if (stat /= 0) return
        dset%next_offset = dset%next_offset + this%ncell_tot
        dset%flag = .true. ! dataset has been written for this time step
      end if
    else
      stat = 1
      errmsg = '"' // name // '" is not a valid temporal cell dataset'
      return
    end if

  end subroutine write_temporal_cell_data_int32

  subroutine write_temporal_cell_data_int64(this, name, data, stat, errmsg)

    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg

    integer, allocatable :: dims(:)
    type(temporal_data), pointer :: dset

    INSIST(this%nsteps >= 0)

    dims = shape(data)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%ncell)

    dset => this%temporal_cell_data
    do while (associated(dset))
      if (dset%name == name) exit
      dset => dset%next
    end do

    if (associated(dset)) then
      if (dset%flag) then
        stat = 1
        errmsg = 'dataset "' // name // '" already written for this time step'
        return
      else
        call h5_append_to_dataset(this%cgrp_id, name, data, this%comm, stat, errmsg)
        if (stat /= 0) return
        dset%next_offset = dset%next_offset + this%ncell_tot
        dset%flag = .true. ! dataset has been written for this time step
      end if
    else
      stat = 1
      errmsg = '"' // name // '" is not a valid temporal cell dataset'
      return
    end if

  end subroutine write_temporal_cell_data_int64

  !! Write the point DATA to the time-dependent PointData group dataset NAME.
  !! The data is associated with the current time step. The type and shape
  !! of DATA must be consistent with the MOLD argument passed to the call
  !! of REGISTER_TEMPORAL_POINT_DATA that created this dataset.

  subroutine write_temporal_point_data_real32(this, name, data, stat, errmsg)

    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg

    integer, allocatable :: dims(:)
    type(temporal_data), pointer :: dset

    INSIST(this%nsteps >= 0)

    dims = shape(data)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%nnode)

    dset => this%temporal_point_data
    do while (associated(dset))
      if (dset%name == name) exit
      dset => dset%next
    end do

    if (associated(dset)) then
      if (dset%flag) then
        stat = 1
        errmsg = 'dataset "' // name // '" already written for this time step'
        return
      else
        call h5_append_to_dataset(this%pgrp_id, name, data, this%comm, stat, errmsg)
        if (stat /= 0) return
        dset%next_offset = dset%next_offset + this%nnode_tot
        dset%flag = .true. ! dataset has been written for this time step
      end if
    else
      stat = 1
      errmsg = '"' // name // '" is not a valid temporal point dataset'
      return
    end if

  end subroutine write_temporal_point_data_real32

  subroutine write_temporal_point_data_real64(this, name, data, stat, errmsg)

    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg

    integer, allocatable :: dims(:)
    type(temporal_data), pointer :: dset

    INSIST(this%nsteps >= 0)

    dims = shape(data)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%nnode)

    dset => this%temporal_point_data
    do while (associated(dset))
      if (dset%name == name) exit
      dset => dset%next
    end do

    if (associated(dset)) then
      if (dset%flag) then
        stat = 1
        errmsg = 'dataset "' // name // '" already written for this time step'
        return
      else
        call h5_append_to_dataset(this%pgrp_id, name, data, this%comm, stat, errmsg)
        if (stat /= 0) return
        dset%next_offset = dset%next_offset + this%nnode_tot
        dset%flag = .true. ! dataset has been written for this time step
      end if
    else
      stat = 1
      errmsg = '"' // name // '" is not a valid temporal point dataset'
      return
    end if

  end subroutine write_temporal_point_data_real64

  subroutine write_temporal_point_data_int32(this, name, data, stat, errmsg)

    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg

    integer, allocatable :: dims(:)
    type(temporal_data), pointer :: dset

    INSIST(this%nsteps >= 0)

    dims = shape(data)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%nnode)

    dset => this%temporal_point_data
    do while (associated(dset))
      if (dset%name == name) exit
      dset => dset%next
    end do

    if (associated(dset)) then
      if (dset%flag) then
        stat = 1
        errmsg = 'dataset "' // name // '" already written for this time step'
        return
      else
        call h5_append_to_dataset(this%pgrp_id, name, data, this%comm, stat, errmsg)
        if (stat /= 0) return
        dset%next_offset = dset%next_offset + this%nnode_tot
        dset%flag = .true. ! dataset has been written for this time step
      end if
    else
      stat = 1
      errmsg = '"' // name // '" is not a valid temporal point dataset'
      return
    end if

  end subroutine write_temporal_point_data_int32

  subroutine write_temporal_point_data_int64(this, name, data, stat, errmsg)

    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg

    integer, allocatable :: dims(:)
    type(temporal_data), pointer :: dset

    INSIST(this%nsteps >= 0)

    dims = shape(data)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%nnode)

    dset => this%temporal_point_data
    do while (associated(dset))
      if (dset%name == name) exit
      dset => dset%next
    end do

    if (associated(dset)) then
      if (dset%flag) then
        stat = 1
        errmsg = 'dataset "' // name // '" already written for this time step'
        return
      else
        call h5_append_to_dataset(this%pgrp_id, name, data, this%comm, stat, errmsg)
        if (stat /= 0) return
        dset%next_offset = dset%next_offset + this%nnode_tot
        dset%flag = .true. ! dataset has been written for this time step
      end if
    else
      stat = 1
      errmsg = '"' // name // '" is not a valid temporal point dataset'
      return
    end if

  end subroutine write_temporal_point_data_int64

end module vtkhdf_ug_type
