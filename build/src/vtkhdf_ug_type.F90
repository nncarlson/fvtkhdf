# 1 "/app/src/vtkhdf_ug_type.F90.fypp"
!!
!! VTKHDF_UG_TYPE
!!
!! This module defines an internal derived type used to manage an HDF5 tree
!! storing an UnstructuredGrid VTKHDF dataset. It serves as a core component of
!! the VTKHDF_MB_FILE and VTKHDF_UG_file types. This implementation supports
!! static and time-dependent point and cell data, but assumes a single static
!! mesh.
!!
!! Copyright (c) 2026 Neil Carlson <neil.n.carlson@gmail.com>
!! SPDX-License-Identifier: BSD-2-Clause
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! The latest specification for the VTKHDF file format is available at:
!! https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html#
!!

#include "vtkhdf_assert.inc"

# 25 "/app/src/vtkhdf_ug_type.F90.fypp"
# 27 "/app/src/vtkhdf_ug_type.F90.fypp"
# 28 "/app/src/vtkhdf_ug_type.F90.fypp"
# 29 "/app/src/vtkhdf_ug_type.F90.fypp"

module vtkhdf_ug_type

  use,intrinsic :: iso_fortran_env
  use vtkhdf_h5_c_binding
  use vtkhdf_h5
  use vtkhdf_ctx_type
  use vtkhdf_h5e_context_type
  implicit none
  private

  type :: temporal_cache_wrapper
    type(temporal_data), pointer :: p => null()
  end type

  type, public :: vtkhdf_ug
    type(vtkhdf_ctx) :: ctx ! unowned copy
    integer(hid_t) :: root_id=H5I_INVALID_HID
    integer(hid_t) :: cgrp_id=H5I_INVALID_HID, pgrp_id=H5I_INVALID_HID
    integer :: nnode, ncell, nnode_tot, ncell_tot, npart
    logical :: is_temporal = .false.
    integer :: nsteps = -1
    integer(hid_t) :: steps_id=H5I_INVALID_HID
    integer(hid_t) :: cogrp_id=H5I_INVALID_HID, pogrp_id=H5I_INVALID_HID
    type(temporal_data), pointer :: temporal_point_data => null()
    type(temporal_data), pointer :: temporal_cell_data => null()
    type(temporal_cache_wrapper), pointer :: cell_cache => null()
    type(temporal_cache_wrapper), pointer :: point_cache => null()
  contains
    procedure :: init
    procedure :: close
# 61 "/app/src/vtkhdf_ug_type.F90.fypp"
# 62 "/app/src/vtkhdf_ug_type.F90.fypp"
    generic :: write_mesh => write_mesh_real32_int32
    procedure, private :: write_mesh_real32_int32
# 62 "/app/src/vtkhdf_ug_type.F90.fypp"
    generic :: write_mesh => write_mesh_real32_int64
    procedure, private :: write_mesh_real32_int64
# 65 "/app/src/vtkhdf_ug_type.F90.fypp"
# 61 "/app/src/vtkhdf_ug_type.F90.fypp"
# 62 "/app/src/vtkhdf_ug_type.F90.fypp"
    generic :: write_mesh => write_mesh_real64_int32
    procedure, private :: write_mesh_real64_int32
# 62 "/app/src/vtkhdf_ug_type.F90.fypp"
    generic :: write_mesh => write_mesh_real64_int64
    procedure, private :: write_mesh_real64_int64
# 65 "/app/src/vtkhdf_ug_type.F90.fypp"
# 66 "/app/src/vtkhdf_ug_type.F90.fypp"
    procedure :: write_time_step
# 68 "/app/src/vtkhdf_ug_type.F90.fypp"
# 69 "/app/src/vtkhdf_ug_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_int32
    procedure, private :: write_cell_data_int32
    generic :: register_temporal_cell_data => register_temporal_cell_data_int32
    procedure, private :: register_temporal_cell_data_int32
    generic :: write_temporal_cell_data => write_temporal_cell_data_int32
    procedure, private :: write_temporal_cell_data_int32
# 69 "/app/src/vtkhdf_ug_type.F90.fypp"
    generic :: write_point_data => write_point_data_int32
    procedure, private :: write_point_data_int32
    generic :: register_temporal_point_data => register_temporal_point_data_int32
    procedure, private :: register_temporal_point_data_int32
    generic :: write_temporal_point_data => write_temporal_point_data_int32
    procedure, private :: write_temporal_point_data_int32
# 76 "/app/src/vtkhdf_ug_type.F90.fypp"
# 68 "/app/src/vtkhdf_ug_type.F90.fypp"
# 69 "/app/src/vtkhdf_ug_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_int64
    procedure, private :: write_cell_data_int64
    generic :: register_temporal_cell_data => register_temporal_cell_data_int64
    procedure, private :: register_temporal_cell_data_int64
    generic :: write_temporal_cell_data => write_temporal_cell_data_int64
    procedure, private :: write_temporal_cell_data_int64
# 69 "/app/src/vtkhdf_ug_type.F90.fypp"
    generic :: write_point_data => write_point_data_int64
    procedure, private :: write_point_data_int64
    generic :: register_temporal_point_data => register_temporal_point_data_int64
    procedure, private :: register_temporal_point_data_int64
    generic :: write_temporal_point_data => write_temporal_point_data_int64
    procedure, private :: write_temporal_point_data_int64
# 76 "/app/src/vtkhdf_ug_type.F90.fypp"
# 68 "/app/src/vtkhdf_ug_type.F90.fypp"
# 69 "/app/src/vtkhdf_ug_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_real32
    procedure, private :: write_cell_data_real32
    generic :: register_temporal_cell_data => register_temporal_cell_data_real32
    procedure, private :: register_temporal_cell_data_real32
    generic :: write_temporal_cell_data => write_temporal_cell_data_real32
    procedure, private :: write_temporal_cell_data_real32
# 69 "/app/src/vtkhdf_ug_type.F90.fypp"
    generic :: write_point_data => write_point_data_real32
    procedure, private :: write_point_data_real32
    generic :: register_temporal_point_data => register_temporal_point_data_real32
    procedure, private :: register_temporal_point_data_real32
    generic :: write_temporal_point_data => write_temporal_point_data_real32
    procedure, private :: write_temporal_point_data_real32
# 76 "/app/src/vtkhdf_ug_type.F90.fypp"
# 68 "/app/src/vtkhdf_ug_type.F90.fypp"
# 69 "/app/src/vtkhdf_ug_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_real64
    procedure, private :: write_cell_data_real64
    generic :: register_temporal_cell_data => register_temporal_cell_data_real64
    procedure, private :: register_temporal_cell_data_real64
    generic :: write_temporal_cell_data => write_temporal_cell_data_real64
    procedure, private :: write_temporal_cell_data_real64
# 69 "/app/src/vtkhdf_ug_type.F90.fypp"
    generic :: write_point_data => write_point_data_real64
    procedure, private :: write_point_data_real64
    generic :: register_temporal_point_data => register_temporal_point_data_real64
    procedure, private :: register_temporal_point_data_real64
    generic :: write_temporal_point_data => write_temporal_point_data_real64
    procedure, private :: write_temporal_point_data_real64
# 76 "/app/src/vtkhdf_ug_type.F90.fypp"
# 77 "/app/src/vtkhdf_ug_type.F90.fypp"
    final :: vtkhdf_ug_delete
  end type

  type :: temporal_data
    character(:), allocatable :: name
    integer :: next_offset = 0
    logical :: flag = .false.
    type(temporal_data), pointer :: next => null()
  contains
    final :: temporal_data_delete
  end type

contains

  !! Finalizer for VTKHDF_UG_FILE objects. We free heap memory we own but avoid
  !! doing things that may require syncronization with other ranks (MPI/PHDF5)
  !! because implicit finalization occurs it is not guaranteed to be collective
  !! or ordered with respect to other ranks when it occurs. This can leak HDF5
  !! IDs, but that is unavoidable. Users should always use CLOSE to do a proper
  !! collective cleanup.

  subroutine vtkhdf_ug_delete(this)
    type(vtkhdf_ug), intent(inout) :: this
    if (associated(this%temporal_cell_data)) deallocate(this%temporal_cell_data)
    if (associated(this%temporal_point_data)) deallocate(this%temporal_point_data)
    if (associated(this%cell_cache)) deallocate(this%cell_cache)
    if (associated(this%point_cache)) deallocate(this%point_cache)
    this%cogrp_id = H5I_INVALID_HID
    this%pogrp_id = H5I_INVALID_HID
    this%steps_id = H5I_INVALID_HID
    this%cgrp_id  = H5I_INVALID_HID
    this%pgrp_id  = H5I_INVALID_HID
    this%root_id  = H5I_INVALID_HID
  end subroutine

  subroutine close(this)
    class(vtkhdf_ug), intent(inout) :: this
    integer :: ierr
    if (H5Iis_valid(this%cogrp_id) > 0) ierr = H5Gclose(this%cogrp_id)
    if (H5Iis_valid(this%pogrp_id) > 0) ierr = H5Gclose(this%pogrp_id)
    if (H5Iis_valid(this%steps_id) > 0) ierr = H5Gclose(this%steps_id)
    if (H5Iis_valid(this%cgrp_id) > 0) ierr = H5Gclose(this%cgrp_id)
    if (H5Iis_valid(this%pgrp_id) > 0) ierr = H5Gclose(this%pgrp_id)
    if (H5Iis_valid(this%root_id) > 0) ierr = H5Gclose(this%root_id)
    call finalization(this) ! free local memory
  contains
    subroutine finalization(this)
      class(vtkhdf_ug), intent(out) :: this
    end subroutine
  end subroutine

  recursive subroutine temporal_data_delete(this)
    type(temporal_data), intent(inout) :: this
    if (associated(this%next)) deallocate(this%next)
  end subroutine

  !! Create a new VTKHDF UnstructuredGrid type group. To create a group
  !! supporting time-dependent data, specify the optional argument
  !! TEMPORAL to true; otherwise a static group is created by default.

  subroutine init(this, loc_id, name, ctx, temporal)

    use vtkhdf_version_param

    class(vtkhdf_ug), intent(inout) :: this
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    type(vtkhdf_ctx), intent(in) :: ctx
    logical, intent(in), optional :: temporal

    integer :: ierr

    this%ctx = ctx

    INSIST(this%ctx%global_all(loc_id >= 0))

    this%root_id = H5Gcreate(loc_id, name)
    INSIST(this%ctx%global_all(this%root_id >= 0))

    call h5_write_attr(this%ctx, this%root_id, 'Version', vtkhdf_version)
    call h5_write_attr(this%ctx, this%root_id, 'Type', 'UnstructuredGrid')

    this%cgrp_id = H5Gcreate(this%root_id, 'CellData')
    INSIST(this%ctx%global_all(this%cgrp_id >= 0))

    this%pgrp_id = H5Gcreate(this%root_id, 'PointData')
    INSIST(this%ctx%global_all(this%pgrp_id >= 0))

    if (present(temporal)) this%is_temporal = temporal
    if (this%is_temporal) then
      this%steps_id = H5Gcreate(this%root_id, 'Steps')
      INSIST(this%ctx%global_all(this%steps_id >= 0))

      this%cogrp_id = H5Gcreate(this%steps_id, 'CellDataOffsets')
      INSIST(this%ctx%global_all(this%cogrp_id >= 0))

      this%pogrp_id = H5Gcreate(this%steps_id, 'PointDataOffsets')
      INSIST(this%ctx%global_all(this%pogrp_id >= 0))

      this%nsteps = 0
      call h5_write_attr(this%ctx, this%steps_id, 'NSteps', this%nsteps)

      associate (mold => 1.0_real64, chunk_size => 100)
        call h5_create_unlimited_dataset(this%ctx, this%steps_id, 'Values', mold, chunk_size)
      end associate

      ! NB: Only static meshes are supported. ParaView 6.1 uses this attribute
      ! setting to signal it should use mesh caching -- read once and reuse.
      call h5_write_attr(this%ctx, this%steps_id, 'StaticMesh', 1)

      ! NB: These are no longer needed in ParaView 6.1 when StaticMesh=1, but
      ! we retain them for compatibility with earlier versions. They will be
      ! filled with 0 offsets and the constant number of parts.
      associate (mold => 1, chunk_size => 100)
        call h5_create_unlimited_dataset(this%ctx, this%steps_id, 'PointOffsets', mold, chunk_size)
        call h5_create_unlimited_dataset(this%ctx, this%steps_id, 'CellOffsets', mold, chunk_size)
        call h5_create_unlimited_dataset(this%ctx, this%steps_id, 'ConnectivityIdOffsets', mold, chunk_size)
        call h5_create_unlimited_dataset(this%ctx, this%steps_id, 'NumberOfParts', mold, chunk_size)
        call h5_create_unlimited_dataset(this%ctx, this%steps_id, 'PartOffsets', mold, chunk_size)
      end associate

      allocate(this%cell_cache)
      allocate(this%point_cache)
    end if

  end subroutine init

  !! Write the unstructured mesh to the file. The mesh is described in the
  !! conventional manner by the X, CNODE, and XCNODE arrays. The additional
  !! array TYPES specifies the VTK cell types. This procedure must be called
  !! before any of the following procedures.

  !! NB: Paraview and/or the VTKHDF reader currently has a problem handling
  !! 0-sized parts (see https://gitlab.kitware.com/vtk/vtk/-/issues/19923).
  !! As a workaround, we omit any empty part by not contributing anything
  !! at all to the NumberOf* and Offsets datasets instead of the 0 we
  !! normally would. Other datasets like Connectivity already contribute
  !! nothing for empty parts.
# 216 "/app/src/vtkhdf_ug_type.F90.fypp"
# 217 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_mesh_real32_int32(this, x, cnode, xcnode, types)

    class(vtkhdf_ug), intent(inout) :: this
    real(real32), intent(in) :: x(:,:)
    integer(int32), intent(in) :: cnode(:), xcnode(:)
    integer(int8), intent(in) :: types(:)

    integer :: ierr, idum(0)

    INSIST(this%root_id >= 0)

    this%nnode = size(x, dim=2)
    this%ncell = size(types)

    INSIST(size(xcnode) == this%ncell+1)
    INSIST(size(cnode) == xcnode(size(xcnode))-1)
    INSIST(minval(cnode) >= 1)
    INSIST(maxval(cnode) <= this%nnode)

    call this%ctx%global_sum(this%nnode, this%nnode_tot)
    call this%ctx%global_sum(this%ncell, this%ncell_tot)
    call this%ctx%global_sum(merge(1, 0, this%ncell>0), this%npart) ! see NB above

    if (this%ncell > 0) then ! see NB above
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfPoints', this%nnode)
    else
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfPoints', idum)
    end if

    if (this%ncell > 0) then ! see NB above
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfCells',  this%ncell)
    else
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfCells',  idum)
    end if

    if (this%ncell > 0) then ! see NB above
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfConnectivityIds', size(cnode))
    else
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfConnectivityIds', idum)
    end if

    if (this%ncell > 0) then ! see NB above
      call h5_write_dataset(this%ctx, this%root_id, 'Offsets', xcnode-1) ! offsets instead of starting indices
    else
#ifdef GNU_PR123899
      block
        integer(int32) :: idum(0)
        call h5_write_dataset(this%ctx, this%root_id, 'Offsets', idum)
      end block
#else
      call h5_write_dataset(this%ctx, this%root_id, 'Offsets', [integer(int32)::])
#endif
    end if

#ifdef GNU_PR123899
    if (this%ncell > 0) then
      call h5_write_dataset(this%ctx, this%root_id, 'Connectivity', cnode-1)  ! 0-based indexing
    else
      block
        integer(int32) :: idum(0)
        call h5_write_dataset(this%ctx, this%root_id, 'Connectivity', idum)
      end block
    endif
#else
    call h5_write_dataset(this%ctx, this%root_id, 'Connectivity', cnode-1)  ! 0-based indexing
#endif

    call h5_write_dataset(this%ctx, this%root_id, 'Types', types)
    call h5_write_dataset(this%ctx, this%root_id, 'Points', x)

  end subroutine write_mesh_real32_int32

# 217 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_mesh_real32_int64(this, x, cnode, xcnode, types)

    class(vtkhdf_ug), intent(inout) :: this
    real(real32), intent(in) :: x(:,:)
    integer(int64), intent(in) :: cnode(:), xcnode(:)
    integer(int8), intent(in) :: types(:)

    integer :: ierr, idum(0)

    INSIST(this%root_id >= 0)

    this%nnode = size(x, dim=2)
    this%ncell = size(types)

    INSIST(size(xcnode) == this%ncell+1)
    INSIST(size(cnode) == xcnode(size(xcnode))-1)
    INSIST(minval(cnode) >= 1)
    INSIST(maxval(cnode) <= this%nnode)

    call this%ctx%global_sum(this%nnode, this%nnode_tot)
    call this%ctx%global_sum(this%ncell, this%ncell_tot)
    call this%ctx%global_sum(merge(1, 0, this%ncell>0), this%npart) ! see NB above

    if (this%ncell > 0) then ! see NB above
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfPoints', this%nnode)
    else
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfPoints', idum)
    end if

    if (this%ncell > 0) then ! see NB above
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfCells',  this%ncell)
    else
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfCells',  idum)
    end if

    if (this%ncell > 0) then ! see NB above
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfConnectivityIds', size(cnode))
    else
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfConnectivityIds', idum)
    end if

    if (this%ncell > 0) then ! see NB above
      call h5_write_dataset(this%ctx, this%root_id, 'Offsets', xcnode-1) ! offsets instead of starting indices
    else
#ifdef GNU_PR123899
      block
        integer(int64) :: idum(0)
        call h5_write_dataset(this%ctx, this%root_id, 'Offsets', idum)
      end block
#else
      call h5_write_dataset(this%ctx, this%root_id, 'Offsets', [integer(int64)::])
#endif
    end if

#ifdef GNU_PR123899
    if (this%ncell > 0) then
      call h5_write_dataset(this%ctx, this%root_id, 'Connectivity', cnode-1)  ! 0-based indexing
    else
      block
        integer(int64) :: idum(0)
        call h5_write_dataset(this%ctx, this%root_id, 'Connectivity', idum)
      end block
    endif
#else
    call h5_write_dataset(this%ctx, this%root_id, 'Connectivity', cnode-1)  ! 0-based indexing
#endif

    call h5_write_dataset(this%ctx, this%root_id, 'Types', types)
    call h5_write_dataset(this%ctx, this%root_id, 'Points', x)

  end subroutine write_mesh_real32_int64

# 291 "/app/src/vtkhdf_ug_type.F90.fypp"
# 216 "/app/src/vtkhdf_ug_type.F90.fypp"
# 217 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_mesh_real64_int32(this, x, cnode, xcnode, types)

    class(vtkhdf_ug), intent(inout) :: this
    real(real64), intent(in) :: x(:,:)
    integer(int32), intent(in) :: cnode(:), xcnode(:)
    integer(int8), intent(in) :: types(:)

    integer :: ierr, idum(0)

    INSIST(this%root_id >= 0)

    this%nnode = size(x, dim=2)
    this%ncell = size(types)

    INSIST(size(xcnode) == this%ncell+1)
    INSIST(size(cnode) == xcnode(size(xcnode))-1)
    INSIST(minval(cnode) >= 1)
    INSIST(maxval(cnode) <= this%nnode)

    call this%ctx%global_sum(this%nnode, this%nnode_tot)
    call this%ctx%global_sum(this%ncell, this%ncell_tot)
    call this%ctx%global_sum(merge(1, 0, this%ncell>0), this%npart) ! see NB above

    if (this%ncell > 0) then ! see NB above
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfPoints', this%nnode)
    else
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfPoints', idum)
    end if

    if (this%ncell > 0) then ! see NB above
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfCells',  this%ncell)
    else
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfCells',  idum)
    end if

    if (this%ncell > 0) then ! see NB above
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfConnectivityIds', size(cnode))
    else
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfConnectivityIds', idum)
    end if

    if (this%ncell > 0) then ! see NB above
      call h5_write_dataset(this%ctx, this%root_id, 'Offsets', xcnode-1) ! offsets instead of starting indices
    else
#ifdef GNU_PR123899
      block
        integer(int32) :: idum(0)
        call h5_write_dataset(this%ctx, this%root_id, 'Offsets', idum)
      end block
#else
      call h5_write_dataset(this%ctx, this%root_id, 'Offsets', [integer(int32)::])
#endif
    end if

#ifdef GNU_PR123899
    if (this%ncell > 0) then
      call h5_write_dataset(this%ctx, this%root_id, 'Connectivity', cnode-1)  ! 0-based indexing
    else
      block
        integer(int32) :: idum(0)
        call h5_write_dataset(this%ctx, this%root_id, 'Connectivity', idum)
      end block
    endif
#else
    call h5_write_dataset(this%ctx, this%root_id, 'Connectivity', cnode-1)  ! 0-based indexing
#endif

    call h5_write_dataset(this%ctx, this%root_id, 'Types', types)
    call h5_write_dataset(this%ctx, this%root_id, 'Points', x)

  end subroutine write_mesh_real64_int32

# 217 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_mesh_real64_int64(this, x, cnode, xcnode, types)

    class(vtkhdf_ug), intent(inout) :: this
    real(real64), intent(in) :: x(:,:)
    integer(int64), intent(in) :: cnode(:), xcnode(:)
    integer(int8), intent(in) :: types(:)

    integer :: ierr, idum(0)

    INSIST(this%root_id >= 0)

    this%nnode = size(x, dim=2)
    this%ncell = size(types)

    INSIST(size(xcnode) == this%ncell+1)
    INSIST(size(cnode) == xcnode(size(xcnode))-1)
    INSIST(minval(cnode) >= 1)
    INSIST(maxval(cnode) <= this%nnode)

    call this%ctx%global_sum(this%nnode, this%nnode_tot)
    call this%ctx%global_sum(this%ncell, this%ncell_tot)
    call this%ctx%global_sum(merge(1, 0, this%ncell>0), this%npart) ! see NB above

    if (this%ncell > 0) then ! see NB above
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfPoints', this%nnode)
    else
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfPoints', idum)
    end if

    if (this%ncell > 0) then ! see NB above
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfCells',  this%ncell)
    else
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfCells',  idum)
    end if

    if (this%ncell > 0) then ! see NB above
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfConnectivityIds', size(cnode))
    else
      call h5_write_dataset(this%ctx, this%root_id, 'NumberOfConnectivityIds', idum)
    end if

    if (this%ncell > 0) then ! see NB above
      call h5_write_dataset(this%ctx, this%root_id, 'Offsets', xcnode-1) ! offsets instead of starting indices
    else
#ifdef GNU_PR123899
      block
        integer(int64) :: idum(0)
        call h5_write_dataset(this%ctx, this%root_id, 'Offsets', idum)
      end block
#else
      call h5_write_dataset(this%ctx, this%root_id, 'Offsets', [integer(int64)::])
#endif
    end if

#ifdef GNU_PR123899
    if (this%ncell > 0) then
      call h5_write_dataset(this%ctx, this%root_id, 'Connectivity', cnode-1)  ! 0-based indexing
    else
      block
        integer(int64) :: idum(0)
        call h5_write_dataset(this%ctx, this%root_id, 'Connectivity', idum)
      end block
    endif
#else
    call h5_write_dataset(this%ctx, this%root_id, 'Connectivity', cnode-1)  ! 0-based indexing
#endif

    call h5_write_dataset(this%ctx, this%root_id, 'Types', types)
    call h5_write_dataset(this%ctx, this%root_id, 'Points', x)

  end subroutine write_mesh_real64_int64

# 291 "/app/src/vtkhdf_ug_type.F90.fypp"
# 292 "/app/src/vtkhdf_ug_type.F90.fypp"
  !! Writes the cell data ARRAY to a new CellData group dataset NAME. In the
  !! case of a temporal file supporting time-dependent datasets, this dataset
  !! is static and not associated with any time step. Note that VTK only
  !! supports scalar and vector-valued mesh data (rank-1 and 2 ARRAY); this
  !! should be a compile time constraint implemented by the public API.
# 298 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_cell_data_int32(this, name, array)
    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: array(..)
    integer, allocatable :: dims(:)
    dims = shape(array)
    INSIST(this%ctx%global_all(dims(size(dims)) == this%ncell))
    call h5_write_dataset(this%ctx, this%cgrp_id, name, array)
  end subroutine
# 298 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_cell_data_int64(this, name, array)
    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: array(..)
    integer, allocatable :: dims(:)
    dims = shape(array)
    INSIST(this%ctx%global_all(dims(size(dims)) == this%ncell))
    call h5_write_dataset(this%ctx, this%cgrp_id, name, array)
  end subroutine
# 298 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_cell_data_real32(this, name, array)
    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: array(..)
    integer, allocatable :: dims(:)
    dims = shape(array)
    INSIST(this%ctx%global_all(dims(size(dims)) == this%ncell))
    call h5_write_dataset(this%ctx, this%cgrp_id, name, array)
  end subroutine
# 298 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_cell_data_real64(this, name, array)
    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(..)
    integer, allocatable :: dims(:)
    dims = shape(array)
    INSIST(this%ctx%global_all(dims(size(dims)) == this%ncell))
    call h5_write_dataset(this%ctx, this%cgrp_id, name, array)
  end subroutine
# 309 "/app/src/vtkhdf_ug_type.F90.fypp"

  !! Writes the point data ARRAY to a new PointData group dataset NAME. In the
  !! case of a temporal file supporting time-dependent datasets, this dataset
  !! is static and not associated with any time step. Note that VTK only
  !! supports scalar and vector-valued mesh data (rank-1 and 2 ARRAY); this
  !! should be a compile time constraint implemented by the public API.
# 316 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_point_data_int32(this, name, array)
    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: array(..)
    integer, allocatable :: dims(:)
    dims = shape(array)
    INSIST(this%ctx%global_all(dims(size(dims)) == this%nnode))
    call h5_write_dataset(this%ctx, this%pgrp_id, name, array)
  end subroutine
# 316 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_point_data_int64(this, name, array)
    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: array(..)
    integer, allocatable :: dims(:)
    dims = shape(array)
    INSIST(this%ctx%global_all(dims(size(dims)) == this%nnode))
    call h5_write_dataset(this%ctx, this%pgrp_id, name, array)
  end subroutine
# 316 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_point_data_real32(this, name, array)
    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: array(..)
    integer, allocatable :: dims(:)
    dims = shape(array)
    INSIST(this%ctx%global_all(dims(size(dims)) == this%nnode))
    call h5_write_dataset(this%ctx, this%pgrp_id, name, array)
  end subroutine
# 316 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_point_data_real64(this, name, array)
    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(..)
    integer, allocatable :: dims(:)
    dims = shape(array)
    INSIST(this%ctx%global_all(dims(size(dims)) == this%nnode))
    call h5_write_dataset(this%ctx, this%pgrp_id, name, array)
  end subroutine
# 327 "/app/src/vtkhdf_ug_type.F90.fypp"

  !! Register NAME as a time-dependent CellData group dataset. This writes
  !! no data, but configures a new extendable HDF5 dataset for future writing.
  !!
  !! The rank-1 array MOLD determines the dataset's fixed leading dimensions:
  !! - For scalar-valued point data: MOLD has size 0;
  !! - For vector-valued point data: MOLD size equals the number of components.
  !!
  !! The type and kind of MOLD defines the HDF5 datatype, and its value is
  !! ignored. The last dimension of written data will correspond to the cell
  !! index, making it the extendable (slowest-varying) dimension of the dataset.
  !! All MPI ranks must provide identical MOLD characteristics.

  !! Register NAME as a time-dependent CellData group dataset. This writes
  !! no data, but configures a new extendable HDF5 dataset for future writing.
  !!
  !! MOLD determines the dataset's fixed leading dimensions:
  !! - For scalar-valued cell data: MOLD is a scalar;
  !! - For vector-valued cell data: MOLD is a 1D array whose size equals the
  !!   number of components.
  !!
  !! The type and kind of MOLD defines the HDF5 datatype, and its value is
  !! ignored. The last dimension of written data will correspond to the cell
  !! index, making it the extendable (slowest-varying) dimension of the dataset.
  !! All MPI ranks must provide identical MOLD characteristics.
# 353 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine register_temporal_cell_data_int32(this, name, mold)

    class(vtkhdf_ug), intent(inout) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: mold(..)

    type(temporal_data), pointer :: new, tail

    INSIST(this%nsteps >= 0) ! must be configured for temporal
    INSIST(this%nsteps == 0) ! cannot register new dataset after time stepping has begun

    !! Dataset for the cell data
    associate (chunk_size => this%ncell_tot)
      call h5_create_unlimited_dataset(this%ctx, this%cgrp_id, name, mold, chunk_size)
    end associate

    !! Dataset for its time step offsets
    associate (mold => 1, chunk_size => 100)
      call h5_create_unlimited_dataset(this%ctx, this%cogrp_id, name, mold, chunk_size)
    end associate

    !! Add to list CellData datasets (append to end)
    allocate(new)
    new%name = name
    new%next => null()
    if (.not. associated(this%temporal_cell_data)) then
       this%temporal_cell_data => new
    else
       tail => this%temporal_cell_data
       do while (associated(tail%next))
          tail => tail%next
       end do
       tail%next => new
    end if

  end subroutine register_temporal_cell_data_int32
# 353 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine register_temporal_cell_data_int64(this, name, mold)

    class(vtkhdf_ug), intent(inout) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: mold(..)

    type(temporal_data), pointer :: new, tail

    INSIST(this%nsteps >= 0) ! must be configured for temporal
    INSIST(this%nsteps == 0) ! cannot register new dataset after time stepping has begun

    !! Dataset for the cell data
    associate (chunk_size => this%ncell_tot)
      call h5_create_unlimited_dataset(this%ctx, this%cgrp_id, name, mold, chunk_size)
    end associate

    !! Dataset for its time step offsets
    associate (mold => 1, chunk_size => 100)
      call h5_create_unlimited_dataset(this%ctx, this%cogrp_id, name, mold, chunk_size)
    end associate

    !! Add to list CellData datasets (append to end)
    allocate(new)
    new%name = name
    new%next => null()
    if (.not. associated(this%temporal_cell_data)) then
       this%temporal_cell_data => new
    else
       tail => this%temporal_cell_data
       do while (associated(tail%next))
          tail => tail%next
       end do
       tail%next => new
    end if

  end subroutine register_temporal_cell_data_int64
# 353 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine register_temporal_cell_data_real32(this, name, mold)

    class(vtkhdf_ug), intent(inout) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: mold(..)

    type(temporal_data), pointer :: new, tail

    INSIST(this%nsteps >= 0) ! must be configured for temporal
    INSIST(this%nsteps == 0) ! cannot register new dataset after time stepping has begun

    !! Dataset for the cell data
    associate (chunk_size => this%ncell_tot)
      call h5_create_unlimited_dataset(this%ctx, this%cgrp_id, name, mold, chunk_size)
    end associate

    !! Dataset for its time step offsets
    associate (mold => 1, chunk_size => 100)
      call h5_create_unlimited_dataset(this%ctx, this%cogrp_id, name, mold, chunk_size)
    end associate

    !! Add to list CellData datasets (append to end)
    allocate(new)
    new%name = name
    new%next => null()
    if (.not. associated(this%temporal_cell_data)) then
       this%temporal_cell_data => new
    else
       tail => this%temporal_cell_data
       do while (associated(tail%next))
          tail => tail%next
       end do
       tail%next => new
    end if

  end subroutine register_temporal_cell_data_real32
# 353 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine register_temporal_cell_data_real64(this, name, mold)

    class(vtkhdf_ug), intent(inout) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: mold(..)

    type(temporal_data), pointer :: new, tail

    INSIST(this%nsteps >= 0) ! must be configured for temporal
    INSIST(this%nsteps == 0) ! cannot register new dataset after time stepping has begun

    !! Dataset for the cell data
    associate (chunk_size => this%ncell_tot)
      call h5_create_unlimited_dataset(this%ctx, this%cgrp_id, name, mold, chunk_size)
    end associate

    !! Dataset for its time step offsets
    associate (mold => 1, chunk_size => 100)
      call h5_create_unlimited_dataset(this%ctx, this%cogrp_id, name, mold, chunk_size)
    end associate

    !! Add to list CellData datasets (append to end)
    allocate(new)
    new%name = name
    new%next => null()
    if (.not. associated(this%temporal_cell_data)) then
       this%temporal_cell_data => new
    else
       tail => this%temporal_cell_data
       do while (associated(tail%next))
          tail => tail%next
       end do
       tail%next => new
    end if

  end subroutine register_temporal_cell_data_real64
# 391 "/app/src/vtkhdf_ug_type.F90.fypp"

  !! Register NAME as a time-dependent PointData group dataset. This writes
  !! no data, but configures a new extendable HDF5 dataset for future writing.
  !!
  !! MOLD determines the dataset's fixed leading dimensions:
  !! - For scalar-valued point data: MOLD is a scalar;
  !! - For vector-valued point data: MOLD is a 1D array whose size equals the
  !!   number of components.
  !!
  !! The type and kind of MOLD defines the HDF5 datatype, and its value is
  !! ignored. The last dimension of written data will correspond to the point
  !! index, making it the extendable (slowest-varying) dimension of the dataset.
  !! All MPI ranks must provide identical MOLD characteristics.
# 405 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine register_temporal_point_data_int32(this, name, mold)

    class(vtkhdf_ug), intent(inout) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: mold(..)

    type(temporal_data), pointer :: new, tail

    INSIST(this%nsteps >= 0) ! must be configured for temporal
    INSIST(this%nsteps == 0) ! cannot register new dataset after time stepping has begun

    !! Dataset for the point data
    associate (chunk_size => this%nnode_tot)
      call h5_create_unlimited_dataset(this%ctx, this%pgrp_id, name, mold, chunk_size)
    end associate

    !! Dataset for its time step offsets
    associate (mold => 1, chunk_size => 100)
      call h5_create_unlimited_dataset(this%ctx, this%pogrp_id, name, mold, chunk_size)
    end associate

    !! Add to list PointData datasets (append to end)
    allocate(new)
    new%name = name
    new%next => null()
    if (.not. associated(this%temporal_point_data)) then
       this%temporal_point_data => new
    else
       tail => this%temporal_point_data
       do while (associated(tail%next))
          tail => tail%next
       end do
       tail%next => new
    end if

  end subroutine register_temporal_point_data_int32
# 405 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine register_temporal_point_data_int64(this, name, mold)

    class(vtkhdf_ug), intent(inout) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: mold(..)

    type(temporal_data), pointer :: new, tail

    INSIST(this%nsteps >= 0) ! must be configured for temporal
    INSIST(this%nsteps == 0) ! cannot register new dataset after time stepping has begun

    !! Dataset for the point data
    associate (chunk_size => this%nnode_tot)
      call h5_create_unlimited_dataset(this%ctx, this%pgrp_id, name, mold, chunk_size)
    end associate

    !! Dataset for its time step offsets
    associate (mold => 1, chunk_size => 100)
      call h5_create_unlimited_dataset(this%ctx, this%pogrp_id, name, mold, chunk_size)
    end associate

    !! Add to list PointData datasets (append to end)
    allocate(new)
    new%name = name
    new%next => null()
    if (.not. associated(this%temporal_point_data)) then
       this%temporal_point_data => new
    else
       tail => this%temporal_point_data
       do while (associated(tail%next))
          tail => tail%next
       end do
       tail%next => new
    end if

  end subroutine register_temporal_point_data_int64
# 405 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine register_temporal_point_data_real32(this, name, mold)

    class(vtkhdf_ug), intent(inout) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: mold(..)

    type(temporal_data), pointer :: new, tail

    INSIST(this%nsteps >= 0) ! must be configured for temporal
    INSIST(this%nsteps == 0) ! cannot register new dataset after time stepping has begun

    !! Dataset for the point data
    associate (chunk_size => this%nnode_tot)
      call h5_create_unlimited_dataset(this%ctx, this%pgrp_id, name, mold, chunk_size)
    end associate

    !! Dataset for its time step offsets
    associate (mold => 1, chunk_size => 100)
      call h5_create_unlimited_dataset(this%ctx, this%pogrp_id, name, mold, chunk_size)
    end associate

    !! Add to list PointData datasets (append to end)
    allocate(new)
    new%name = name
    new%next => null()
    if (.not. associated(this%temporal_point_data)) then
       this%temporal_point_data => new
    else
       tail => this%temporal_point_data
       do while (associated(tail%next))
          tail => tail%next
       end do
       tail%next => new
    end if

  end subroutine register_temporal_point_data_real32
# 405 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine register_temporal_point_data_real64(this, name, mold)

    class(vtkhdf_ug), intent(inout) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: mold(..)

    type(temporal_data), pointer :: new, tail

    INSIST(this%nsteps >= 0) ! must be configured for temporal
    INSIST(this%nsteps == 0) ! cannot register new dataset after time stepping has begun

    !! Dataset for the point data
    associate (chunk_size => this%nnode_tot)
      call h5_create_unlimited_dataset(this%ctx, this%pgrp_id, name, mold, chunk_size)
    end associate

    !! Dataset for its time step offsets
    associate (mold => 1, chunk_size => 100)
      call h5_create_unlimited_dataset(this%ctx, this%pogrp_id, name, mold, chunk_size)
    end associate

    !! Add to list PointData datasets (append to end)
    allocate(new)
    new%name = name
    new%next => null()
    if (.not. associated(this%temporal_point_data)) then
       this%temporal_point_data => new
    else
       tail => this%temporal_point_data
       do while (associated(tail%next))
          tail => tail%next
       end do
       tail%next => new
    end if

  end subroutine register_temporal_point_data_real64
# 443 "/app/src/vtkhdf_ug_type.F90.fypp"

  !! Mark the start of a new time step with time value TIME. Subsequent output
  !! of time-dependent datasets will be associated with this time step.

  subroutine write_time_step(this, time)

    class(vtkhdf_ug), intent(inout) :: this
    real(real64), intent(in) :: time

    type(temporal_data), pointer :: tmp

    INSIST(this%nsteps >= 0)

    this%nsteps = this%nsteps + 1
    call h5_write_attr(this%ctx, this%steps_id, 'NSteps', this%nsteps)

    call h5_append_to_dataset(this%ctx, this%steps_id, 'Values', time, root=0)

    !! A single mesh is used for all time steps so there are no offsets
    call h5_append_to_dataset(this%ctx, this%steps_id, 'PointOffsets', 0, root=0)
    call h5_append_to_dataset(this%ctx, this%steps_id, 'CellOffsets', 0, root=0)
    call h5_append_to_dataset(this%ctx, this%steps_id, 'ConnectivityIdOffsets', 0, root=0)
    call h5_append_to_dataset(this%ctx, this%steps_id, 'NumberOfParts', this%npart, root=0)
    call h5_append_to_dataset(this%ctx, this%steps_id, 'PartOffsets', 0, root=0)

    !! Cell and point data offsets for each of the temporal datasets for this timestep

    tmp => this%temporal_point_data
    do while (associated(tmp))
      tmp%flag = .false.  ! dataset not yet written for this time step
      call h5_append_to_dataset(this%ctx, this%pogrp_id, tmp%name, tmp%next_offset, root=0)
      tmp => tmp%next
    end do

    tmp => this%temporal_cell_data
    do while (associated(tmp))
      tmp%flag = .false.  ! dataset not yet written for this time step
      call h5_append_to_dataset(this%ctx, this%cogrp_id, tmp%name, tmp%next_offset, root=0)
      tmp => tmp%next
    end do

  end subroutine write_time_step

  !! Write the cell-based ARRAY to the time-dependent CellData group dataset
  !! NAME. The data is associated with the current time step. The type and
  !! leading dimensions of ARRAY must match the MOLD provided during the
  !! registration of NAME. The last (or only) dimension of ARRAY must
  !  correspond to the local number of cells on this MPI rank.
# 492 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_temporal_cell_data_int32(this, name, array)

    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: array(..)

    integer, allocatable :: dims(:)
    type(temporal_data), pointer :: dset

    INSIST(this%nsteps > 0)

    dims = shape(array)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%ncell)

    dset => null()
    if (associated(this%cell_cache)) then
      ! Check current
      if (associated(this%cell_cache%p)) then
        if (this%cell_cache%p%name == name) then
          dset => this%cell_cache%p
        ! Check next (optimization for sequential access)
        elseif (associated(this%cell_cache%p%next)) then
          if (this%cell_cache%p%next%name == name) then
             dset => this%cell_cache%p%next
          end if
        end if
      end if
    end if

    if (.not. associated(dset)) then
      dset => this%temporal_cell_data
      do while (associated(dset))
        if (dset%name == name) exit
        dset => dset%next
      end do
      INSIST(associated(dset))
    end if

    ! Update cache
    if (associated(this%cell_cache)) then
       this%cell_cache%p => dset
    end if

    INSIST(.not.dset%flag) ! no double writing in a single step
    call h5_append_to_dataset(this%ctx, this%cgrp_id, name, array)
    dset%next_offset = dset%next_offset + this%ncell_tot
    dset%flag = .true. ! dataset has been written for this step

  end subroutine write_temporal_cell_data_int32
# 492 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_temporal_cell_data_int64(this, name, array)

    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: array(..)

    integer, allocatable :: dims(:)
    type(temporal_data), pointer :: dset

    INSIST(this%nsteps > 0)

    dims = shape(array)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%ncell)

    dset => null()
    if (associated(this%cell_cache)) then
      ! Check current
      if (associated(this%cell_cache%p)) then
        if (this%cell_cache%p%name == name) then
          dset => this%cell_cache%p
        ! Check next (optimization for sequential access)
        elseif (associated(this%cell_cache%p%next)) then
          if (this%cell_cache%p%next%name == name) then
             dset => this%cell_cache%p%next
          end if
        end if
      end if
    end if

    if (.not. associated(dset)) then
      dset => this%temporal_cell_data
      do while (associated(dset))
        if (dset%name == name) exit
        dset => dset%next
      end do
      INSIST(associated(dset))
    end if

    ! Update cache
    if (associated(this%cell_cache)) then
       this%cell_cache%p => dset
    end if

    INSIST(.not.dset%flag) ! no double writing in a single step
    call h5_append_to_dataset(this%ctx, this%cgrp_id, name, array)
    dset%next_offset = dset%next_offset + this%ncell_tot
    dset%flag = .true. ! dataset has been written for this step

  end subroutine write_temporal_cell_data_int64
# 492 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_temporal_cell_data_real32(this, name, array)

    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: array(..)

    integer, allocatable :: dims(:)
    type(temporal_data), pointer :: dset

    INSIST(this%nsteps > 0)

    dims = shape(array)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%ncell)

    dset => null()
    if (associated(this%cell_cache)) then
      ! Check current
      if (associated(this%cell_cache%p)) then
        if (this%cell_cache%p%name == name) then
          dset => this%cell_cache%p
        ! Check next (optimization for sequential access)
        elseif (associated(this%cell_cache%p%next)) then
          if (this%cell_cache%p%next%name == name) then
             dset => this%cell_cache%p%next
          end if
        end if
      end if
    end if

    if (.not. associated(dset)) then
      dset => this%temporal_cell_data
      do while (associated(dset))
        if (dset%name == name) exit
        dset => dset%next
      end do
      INSIST(associated(dset))
    end if

    ! Update cache
    if (associated(this%cell_cache)) then
       this%cell_cache%p => dset
    end if

    INSIST(.not.dset%flag) ! no double writing in a single step
    call h5_append_to_dataset(this%ctx, this%cgrp_id, name, array)
    dset%next_offset = dset%next_offset + this%ncell_tot
    dset%flag = .true. ! dataset has been written for this step

  end subroutine write_temporal_cell_data_real32
# 492 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_temporal_cell_data_real64(this, name, array)

    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(..)

    integer, allocatable :: dims(:)
    type(temporal_data), pointer :: dset

    INSIST(this%nsteps > 0)

    dims = shape(array)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%ncell)

    dset => null()
    if (associated(this%cell_cache)) then
      ! Check current
      if (associated(this%cell_cache%p)) then
        if (this%cell_cache%p%name == name) then
          dset => this%cell_cache%p
        ! Check next (optimization for sequential access)
        elseif (associated(this%cell_cache%p%next)) then
          if (this%cell_cache%p%next%name == name) then
             dset => this%cell_cache%p%next
          end if
        end if
      end if
    end if

    if (.not. associated(dset)) then
      dset => this%temporal_cell_data
      do while (associated(dset))
        if (dset%name == name) exit
        dset => dset%next
      end do
      INSIST(associated(dset))
    end if

    ! Update cache
    if (associated(this%cell_cache)) then
       this%cell_cache%p => dset
    end if

    INSIST(.not.dset%flag) ! no double writing in a single step
    call h5_append_to_dataset(this%ctx, this%cgrp_id, name, array)
    dset%next_offset = dset%next_offset + this%ncell_tot
    dset%flag = .true. ! dataset has been written for this step

  end subroutine write_temporal_cell_data_real64
# 544 "/app/src/vtkhdf_ug_type.F90.fypp"

  !! Write the point-based ARRAY to the time-dependent PointData group dataset
  !! NAME. The data is associated with the current time step. The type and
  !! leading dimensions of ARRAY must match the MOLD provided during the
  !! registration of NAME. The last (or only) dimension of ARRAY must
  !  correspond to the local number of points on this MPI rank.
# 551 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_temporal_point_data_int32(this, name, array)

    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: array(..)

    integer, allocatable :: dims(:)
    type(temporal_data), pointer :: dset

    INSIST(this%nsteps > 0)

    dims = shape(array)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%nnode)

    dset => null()
    if (associated(this%point_cache)) then
      ! Check current
      if (associated(this%point_cache%p)) then
        if (this%point_cache%p%name == name) then
          dset => this%point_cache%p
        ! Check next (optimization for sequential access)
        elseif (associated(this%point_cache%p%next)) then
          if (this%point_cache%p%next%name == name) then
             dset => this%point_cache%p%next
          end if
        end if
      end if
    end if

    if (.not. associated(dset)) then
      dset => this%temporal_point_data
      do while (associated(dset))
        if (dset%name == name) exit
        dset => dset%next
      end do
      INSIST(associated(dset))
    end if

    ! Update cache
    if (associated(this%point_cache)) then
       this%point_cache%p => dset
    end if

    INSIST(.not.dset%flag) ! no double writing in a single step
    call h5_append_to_dataset(this%ctx, this%pgrp_id, name, array)
    dset%next_offset = dset%next_offset + this%nnode_tot
    dset%flag = .true. ! dataset has been written for this step

  end subroutine write_temporal_point_data_int32
# 551 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_temporal_point_data_int64(this, name, array)

    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: array(..)

    integer, allocatable :: dims(:)
    type(temporal_data), pointer :: dset

    INSIST(this%nsteps > 0)

    dims = shape(array)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%nnode)

    dset => null()
    if (associated(this%point_cache)) then
      ! Check current
      if (associated(this%point_cache%p)) then
        if (this%point_cache%p%name == name) then
          dset => this%point_cache%p
        ! Check next (optimization for sequential access)
        elseif (associated(this%point_cache%p%next)) then
          if (this%point_cache%p%next%name == name) then
             dset => this%point_cache%p%next
          end if
        end if
      end if
    end if

    if (.not. associated(dset)) then
      dset => this%temporal_point_data
      do while (associated(dset))
        if (dset%name == name) exit
        dset => dset%next
      end do
      INSIST(associated(dset))
    end if

    ! Update cache
    if (associated(this%point_cache)) then
       this%point_cache%p => dset
    end if

    INSIST(.not.dset%flag) ! no double writing in a single step
    call h5_append_to_dataset(this%ctx, this%pgrp_id, name, array)
    dset%next_offset = dset%next_offset + this%nnode_tot
    dset%flag = .true. ! dataset has been written for this step

  end subroutine write_temporal_point_data_int64
# 551 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_temporal_point_data_real32(this, name, array)

    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: array(..)

    integer, allocatable :: dims(:)
    type(temporal_data), pointer :: dset

    INSIST(this%nsteps > 0)

    dims = shape(array)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%nnode)

    dset => null()
    if (associated(this%point_cache)) then
      ! Check current
      if (associated(this%point_cache%p)) then
        if (this%point_cache%p%name == name) then
          dset => this%point_cache%p
        ! Check next (optimization for sequential access)
        elseif (associated(this%point_cache%p%next)) then
          if (this%point_cache%p%next%name == name) then
             dset => this%point_cache%p%next
          end if
        end if
      end if
    end if

    if (.not. associated(dset)) then
      dset => this%temporal_point_data
      do while (associated(dset))
        if (dset%name == name) exit
        dset => dset%next
      end do
      INSIST(associated(dset))
    end if

    ! Update cache
    if (associated(this%point_cache)) then
       this%point_cache%p => dset
    end if

    INSIST(.not.dset%flag) ! no double writing in a single step
    call h5_append_to_dataset(this%ctx, this%pgrp_id, name, array)
    dset%next_offset = dset%next_offset + this%nnode_tot
    dset%flag = .true. ! dataset has been written for this step

  end subroutine write_temporal_point_data_real32
# 551 "/app/src/vtkhdf_ug_type.F90.fypp"

  subroutine write_temporal_point_data_real64(this, name, array)

    class(vtkhdf_ug), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(..)

    integer, allocatable :: dims(:)
    type(temporal_data), pointer :: dset

    INSIST(this%nsteps > 0)

    dims = shape(array)
    INSIST(size(dims) >= 1 .and. size(dims) <= 3)
    INSIST(dims(size(dims)) == this%nnode)

    dset => null()
    if (associated(this%point_cache)) then
      ! Check current
      if (associated(this%point_cache%p)) then
        if (this%point_cache%p%name == name) then
          dset => this%point_cache%p
        ! Check next (optimization for sequential access)
        elseif (associated(this%point_cache%p%next)) then
          if (this%point_cache%p%next%name == name) then
             dset => this%point_cache%p%next
          end if
        end if
      end if
    end if

    if (.not. associated(dset)) then
      dset => this%temporal_point_data
      do while (associated(dset))
        if (dset%name == name) exit
        dset => dset%next
      end do
      INSIST(associated(dset))
    end if

    ! Update cache
    if (associated(this%point_cache)) then
       this%point_cache%p => dset
    end if

    INSIST(.not.dset%flag) ! no double writing in a single step
    call h5_append_to_dataset(this%ctx, this%pgrp_id, name, array)
    dset%next_offset = dset%next_offset + this%nnode_tot
    dset%flag = .true. ! dataset has been written for this step

  end subroutine write_temporal_point_data_real64
# 603 "/app/src/vtkhdf_ug_type.F90.fypp"

end module vtkhdf_ug_type
