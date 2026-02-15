# 1 "/app/src/vtkhdf_ug_file_type.F90.fypp"
!!
!! VTKHDF_UG_FILE_TYPE
!!
!! This module defines a derived type for exporting mesh-based solution data
!! to a VTKHDF file readable by the ParaView visualization tool. The type
!! produces a UnstructuredGrid-type file (UG) dataset, and supports static
!! and time-dependent cell and point data, assuming a single static mesh.
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

# 24 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 26 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 27 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 28 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 29 "/app/src/vtkhdf_ug_file_type.F90.fypp"

module vtkhdf_ug_file_type

  use,intrinsic :: iso_fortran_env
  use vtkhdf_h5_c_binding
  use vtkhdf_h5
  use vtkhdf_ug_type
  use vtkhdf_ctx_type
  use vtkhdf_h5e_context_type
  use vtkhdf_version_param
  implicit none
  private

  public :: vtkhdf_version

  type, public :: vtkhdf_ug_file
    private
    type(vtkhdf_ctx) :: ctx
    integer(hid_t)  :: file_id=H5I_INVALID_HID
    type(vtkhdf_ug) :: ug
  contains
    procedure :: create
    procedure :: close
# 53 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 54 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_mesh => write_mesh_real32_int32
    procedure, private :: write_mesh_real32_int32
# 54 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_mesh => write_mesh_real32_int64
    procedure, private :: write_mesh_real32_int64
# 57 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 53 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 54 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_mesh => write_mesh_real64_int32
    procedure, private :: write_mesh_real64_int32
# 54 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_mesh => write_mesh_real64_int64
    procedure, private :: write_mesh_real64_int64
# 57 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 58 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    procedure :: write_time_step
# 60 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 61 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 62 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_int32_scalar
    procedure, private :: write_cell_data_int32_scalar
    generic :: register_temporal_cell_data => register_temporal_cell_data_int32_scalar
    procedure, private :: register_temporal_cell_data_int32_scalar
    generic :: write_temporal_cell_data => write_temporal_cell_data_int32_scalar
    procedure, private :: write_temporal_cell_data_int32_scalar
# 62 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_int32_vector
    procedure, private :: write_cell_data_int32_vector
    generic :: register_temporal_cell_data => register_temporal_cell_data_int32_vector
    procedure, private :: register_temporal_cell_data_int32_vector
    generic :: write_temporal_cell_data => write_temporal_cell_data_int32_vector
    procedure, private :: write_temporal_cell_data_int32_vector
# 69 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 61 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 62 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_point_data => write_point_data_int32_scalar
    procedure, private :: write_point_data_int32_scalar
    generic :: register_temporal_point_data => register_temporal_point_data_int32_scalar
    procedure, private :: register_temporal_point_data_int32_scalar
    generic :: write_temporal_point_data => write_temporal_point_data_int32_scalar
    procedure, private :: write_temporal_point_data_int32_scalar
# 62 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_point_data => write_point_data_int32_vector
    procedure, private :: write_point_data_int32_vector
    generic :: register_temporal_point_data => register_temporal_point_data_int32_vector
    procedure, private :: register_temporal_point_data_int32_vector
    generic :: write_temporal_point_data => write_temporal_point_data_int32_vector
    procedure, private :: write_temporal_point_data_int32_vector
# 69 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 70 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 60 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 61 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 62 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_int64_scalar
    procedure, private :: write_cell_data_int64_scalar
    generic :: register_temporal_cell_data => register_temporal_cell_data_int64_scalar
    procedure, private :: register_temporal_cell_data_int64_scalar
    generic :: write_temporal_cell_data => write_temporal_cell_data_int64_scalar
    procedure, private :: write_temporal_cell_data_int64_scalar
# 62 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_int64_vector
    procedure, private :: write_cell_data_int64_vector
    generic :: register_temporal_cell_data => register_temporal_cell_data_int64_vector
    procedure, private :: register_temporal_cell_data_int64_vector
    generic :: write_temporal_cell_data => write_temporal_cell_data_int64_vector
    procedure, private :: write_temporal_cell_data_int64_vector
# 69 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 61 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 62 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_point_data => write_point_data_int64_scalar
    procedure, private :: write_point_data_int64_scalar
    generic :: register_temporal_point_data => register_temporal_point_data_int64_scalar
    procedure, private :: register_temporal_point_data_int64_scalar
    generic :: write_temporal_point_data => write_temporal_point_data_int64_scalar
    procedure, private :: write_temporal_point_data_int64_scalar
# 62 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_point_data => write_point_data_int64_vector
    procedure, private :: write_point_data_int64_vector
    generic :: register_temporal_point_data => register_temporal_point_data_int64_vector
    procedure, private :: register_temporal_point_data_int64_vector
    generic :: write_temporal_point_data => write_temporal_point_data_int64_vector
    procedure, private :: write_temporal_point_data_int64_vector
# 69 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 70 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 60 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 61 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 62 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_real32_scalar
    procedure, private :: write_cell_data_real32_scalar
    generic :: register_temporal_cell_data => register_temporal_cell_data_real32_scalar
    procedure, private :: register_temporal_cell_data_real32_scalar
    generic :: write_temporal_cell_data => write_temporal_cell_data_real32_scalar
    procedure, private :: write_temporal_cell_data_real32_scalar
# 62 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_real32_vector
    procedure, private :: write_cell_data_real32_vector
    generic :: register_temporal_cell_data => register_temporal_cell_data_real32_vector
    procedure, private :: register_temporal_cell_data_real32_vector
    generic :: write_temporal_cell_data => write_temporal_cell_data_real32_vector
    procedure, private :: write_temporal_cell_data_real32_vector
# 69 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 61 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 62 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_point_data => write_point_data_real32_scalar
    procedure, private :: write_point_data_real32_scalar
    generic :: register_temporal_point_data => register_temporal_point_data_real32_scalar
    procedure, private :: register_temporal_point_data_real32_scalar
    generic :: write_temporal_point_data => write_temporal_point_data_real32_scalar
    procedure, private :: write_temporal_point_data_real32_scalar
# 62 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_point_data => write_point_data_real32_vector
    procedure, private :: write_point_data_real32_vector
    generic :: register_temporal_point_data => register_temporal_point_data_real32_vector
    procedure, private :: register_temporal_point_data_real32_vector
    generic :: write_temporal_point_data => write_temporal_point_data_real32_vector
    procedure, private :: write_temporal_point_data_real32_vector
# 69 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 70 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 60 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 61 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 62 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_real64_scalar
    procedure, private :: write_cell_data_real64_scalar
    generic :: register_temporal_cell_data => register_temporal_cell_data_real64_scalar
    procedure, private :: register_temporal_cell_data_real64_scalar
    generic :: write_temporal_cell_data => write_temporal_cell_data_real64_scalar
    procedure, private :: write_temporal_cell_data_real64_scalar
# 62 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_real64_vector
    procedure, private :: write_cell_data_real64_vector
    generic :: register_temporal_cell_data => register_temporal_cell_data_real64_vector
    procedure, private :: register_temporal_cell_data_real64_vector
    generic :: write_temporal_cell_data => write_temporal_cell_data_real64_vector
    procedure, private :: write_temporal_cell_data_real64_vector
# 69 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 61 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 62 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_point_data => write_point_data_real64_scalar
    procedure, private :: write_point_data_real64_scalar
    generic :: register_temporal_point_data => register_temporal_point_data_real64_scalar
    procedure, private :: register_temporal_point_data_real64_scalar
    generic :: write_temporal_point_data => write_temporal_point_data_real64_scalar
    procedure, private :: write_temporal_point_data_real64_scalar
# 62 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    generic :: write_point_data => write_point_data_real64_vector
    procedure, private :: write_point_data_real64_vector
    generic :: register_temporal_point_data => register_temporal_point_data_real64_vector
    procedure, private :: register_temporal_point_data_real64_vector
    generic :: write_temporal_point_data => write_temporal_point_data_real64_vector
    procedure, private :: write_temporal_point_data_real64_vector
# 69 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 70 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 71 "/app/src/vtkhdf_ug_file_type.F90.fypp"
    final :: vtkhdf_ug_file_delete
  end type

contains

  !! Finalizer for VTKHDF_UG_FILE objects. We free heap memory we own but avoid
  !! doing things that may require syncronization with other ranks (MPI/PHDF5)
  !! because where implicit finalization occurs it is not guaranteed to be
  !! collective or ordered with respect to other ranks. This can leak HDF5
  !! IDs and the MPI communicator, but that is unavoidable. Users should
  !! always use CLOSE to do a proper collective cleanup and close of the file.

  subroutine vtkhdf_ug_file_delete(this)
    type(vtkhdf_ug_file), intent(inout) :: this
    this%file_id = H5I_INVALID_HID
  end subroutine

  !! Cleanly closes H5 identifiers and the file, and default initializes.
  subroutine close(this)
    class(vtkhdf_ug_file), intent(inout) :: this
    integer :: ierr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%close
    if (H5Iis_valid(this%file_id) > 0) ierr = H5Fclose(this%file_id)
    call err_ctx%restore(0)
    call this%ctx%free  ! frees duplicated comm in MPI build
    call finalize(this) ! free local memory
  contains
    subroutine finalize(this)
      class(vtkhdf_ug_file), intent(out) :: this
    end subroutine
  end subroutine

#ifdef USE_MPI
  subroutine create(this, filename, comm, stat, errmsg, is_temporal)
#else
  subroutine create(this, filename, stat, errmsg, is_temporal)
#endif

    class(vtkhdf_ug_file), intent(out) :: this
    character(*), intent(in) :: filename
#ifdef USE_MPI
    integer, intent(in) :: comm
#endif
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    logical, intent(in), optional :: is_temporal

    integer(hid_t) :: fapl
    integer :: ierr
    type(vtkhdf_h5e_context) :: err_ctx

    ! MPI build: duplicates comm, caches rank/size
    ! Serial build: no comm, sets rank=0, size=1
#ifdef USE_MPI
    call this%ctx%init(comm)
#else
    call this%ctx%init()
#endif

    call init_hdf5
    call err_ctx%capture

    !! Open the file (read/write), overwrite any existing file
    call this%ctx%get_fapl(fapl) ! configures for parallel I/O in MPI build
    this%file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl)
    ierr = H5Pclose(fapl)
    if (this%ctx%global_any(this%file_id < 0)) then
      stat = 1
      errmsg = 'failed to open file'
      if (this%file_id > 0) ierr = H5Fclose(this%file_id)
      this%file_id = H5I_INVALID_HID
      call this%ctx%free
      call err_ctx%restore(stat)
      return
    end if

    call this%ug%init(this%file_id, 'VTKHDF', this%ctx, is_temporal)

    stat = 0
    call err_ctx%restore(stat)

  end subroutine create

  !! Write the unstructured mesh to the file. The mesh is described in the
  !! conventional manner by the X, CNODE, and XCNODE arrays. The additional
  !! array TYPES specifies the VTK cell types. This procedure must be called
  !! before any of the following procedures.
# 161 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 162 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_mesh_real32_int32(this, x, cnode, xcnode, types)
    class(vtkhdf_ug_file), intent(inout) :: this
    real(real32), intent(in) :: x(:,:)
    integer(int32), intent(in) :: cnode(:), xcnode(:)
    integer(int8), intent(in) :: types(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_mesh(x, cnode, xcnode, types)
    call err_ctx%restore(0)
  end subroutine
# 162 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_mesh_real32_int64(this, x, cnode, xcnode, types)
    class(vtkhdf_ug_file), intent(inout) :: this
    real(real32), intent(in) :: x(:,:)
    integer(int64), intent(in) :: cnode(:), xcnode(:)
    integer(int8), intent(in) :: types(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_mesh(x, cnode, xcnode, types)
    call err_ctx%restore(0)
  end subroutine
# 174 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 161 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 162 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_mesh_real64_int32(this, x, cnode, xcnode, types)
    class(vtkhdf_ug_file), intent(inout) :: this
    real(real64), intent(in) :: x(:,:)
    integer(int32), intent(in) :: cnode(:), xcnode(:)
    integer(int8), intent(in) :: types(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_mesh(x, cnode, xcnode, types)
    call err_ctx%restore(0)
  end subroutine
# 162 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_mesh_real64_int64(this, x, cnode, xcnode, types)
    class(vtkhdf_ug_file), intent(inout) :: this
    real(real64), intent(in) :: x(:,:)
    integer(int64), intent(in) :: cnode(:), xcnode(:)
    integer(int8), intent(in) :: types(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_mesh(x, cnode, xcnode, types)
    call err_ctx%restore(0)
  end subroutine
# 174 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 175 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  !! Writes the cell-based data ARRAY to a new named cell dataset. Scalar,
  !! vector, and tensor data are supported. For temporal files, this dataset
  !! is static and not associated with any time step.
# 180 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 181 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 182 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_cell_data_int32_scalar(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: array(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 182 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_cell_data_int32_vector(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: array(:,:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 193 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 181 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 182 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_cell_data_int64_scalar(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: array(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 182 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_cell_data_int64_vector(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: array(:,:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 193 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 181 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 182 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_cell_data_real32_scalar(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: array(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 182 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_cell_data_real32_vector(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: array(:,:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 193 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 181 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 182 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_cell_data_real64_scalar(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 182 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_cell_data_real64_vector(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(:,:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 193 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 194 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 180 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 181 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 182 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_point_data_int32_scalar(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: array(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 182 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_point_data_int32_vector(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: array(:,:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 193 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 181 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 182 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_point_data_int64_scalar(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: array(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 182 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_point_data_int64_vector(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: array(:,:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 193 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 181 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 182 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_point_data_real32_scalar(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: array(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 182 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_point_data_real32_vector(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: array(:,:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 193 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 181 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 182 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_point_data_real64_scalar(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 182 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_point_data_real64_vector(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(:,:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 193 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 194 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 195 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  !! Register NAME as a time-dependent PointData dataset for the VTKHDF file.
  !! This prepares the file to receive point-based data at each time step
  !! but performs no actual data write.
  !!
  !! The MOLD argument defines the nature of the point-based data:
  !! - For scalar-valued data: MOLD is a scalar (e.g., a single real value);
  !! - For vector-valued data: MOLD is a 1D array where size equals the number
  !!   of components (e.g., size 3 for velocity).
  !!
  !! The type and kind of MOLD determine the field's data type in the file.
  !! While the values in MOLD are ignored, its characteristics must be
  !! identical across all MPI ranks.
  !!
  !! Subsequent writes to this field assume the last (or only) dimension
  !! of the data array corresponds to the point index.
# 212 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 213 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 214 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine register_temporal_cell_data_int32_scalar(this, name, mold)
    class(vtkhdf_ug_file), intent(inout) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: mold
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture !TODO? push this down?
    call this%ug%register_temporal_cell_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 214 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine register_temporal_cell_data_int32_vector(this, name, mold)
    class(vtkhdf_ug_file), intent(inout) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: mold(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture !TODO? push this down?
    call this%ug%register_temporal_cell_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 225 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 213 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 214 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine register_temporal_cell_data_int64_scalar(this, name, mold)
    class(vtkhdf_ug_file), intent(inout) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: mold
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture !TODO? push this down?
    call this%ug%register_temporal_cell_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 214 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine register_temporal_cell_data_int64_vector(this, name, mold)
    class(vtkhdf_ug_file), intent(inout) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: mold(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture !TODO? push this down?
    call this%ug%register_temporal_cell_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 225 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 213 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 214 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine register_temporal_cell_data_real32_scalar(this, name, mold)
    class(vtkhdf_ug_file), intent(inout) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: mold
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture !TODO? push this down?
    call this%ug%register_temporal_cell_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 214 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine register_temporal_cell_data_real32_vector(this, name, mold)
    class(vtkhdf_ug_file), intent(inout) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: mold(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture !TODO? push this down?
    call this%ug%register_temporal_cell_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 225 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 213 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 214 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine register_temporal_cell_data_real64_scalar(this, name, mold)
    class(vtkhdf_ug_file), intent(inout) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: mold
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture !TODO? push this down?
    call this%ug%register_temporal_cell_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 214 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine register_temporal_cell_data_real64_vector(this, name, mold)
    class(vtkhdf_ug_file), intent(inout) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: mold(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture !TODO? push this down?
    call this%ug%register_temporal_cell_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 225 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 226 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 212 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 213 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 214 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine register_temporal_point_data_int32_scalar(this, name, mold)
    class(vtkhdf_ug_file), intent(inout) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: mold
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture !TODO? push this down?
    call this%ug%register_temporal_point_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 214 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine register_temporal_point_data_int32_vector(this, name, mold)
    class(vtkhdf_ug_file), intent(inout) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: mold(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture !TODO? push this down?
    call this%ug%register_temporal_point_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 225 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 213 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 214 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine register_temporal_point_data_int64_scalar(this, name, mold)
    class(vtkhdf_ug_file), intent(inout) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: mold
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture !TODO? push this down?
    call this%ug%register_temporal_point_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 214 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine register_temporal_point_data_int64_vector(this, name, mold)
    class(vtkhdf_ug_file), intent(inout) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: mold(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture !TODO? push this down?
    call this%ug%register_temporal_point_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 225 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 213 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 214 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine register_temporal_point_data_real32_scalar(this, name, mold)
    class(vtkhdf_ug_file), intent(inout) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: mold
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture !TODO? push this down?
    call this%ug%register_temporal_point_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 214 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine register_temporal_point_data_real32_vector(this, name, mold)
    class(vtkhdf_ug_file), intent(inout) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: mold(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture !TODO? push this down?
    call this%ug%register_temporal_point_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 225 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 213 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 214 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine register_temporal_point_data_real64_scalar(this, name, mold)
    class(vtkhdf_ug_file), intent(inout) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: mold
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture !TODO? push this down?
    call this%ug%register_temporal_point_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 214 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine register_temporal_point_data_real64_vector(this, name, mold)
    class(vtkhdf_ug_file), intent(inout) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: mold(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture !TODO? push this down?
    call this%ug%register_temporal_point_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 225 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 226 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 227 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  !! Mark the start of a new time step with time value TIME. Subsequent output
  !! of time-dependent datasets will be associated with this time step.

  subroutine write_time_step(this, time)
    class(vtkhdf_ug_file), intent(inout) :: this
    real(real64), intent(in) :: time
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_time_step(time)
    call err_ctx%restore(0)
  end subroutine

  !! Write the cell/point-based data ARRAY to the named time-dependent dataset.
  !! The data is associated with the current time step.
# 243 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 244 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 245 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_temporal_cell_data_int32_scalar(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: array(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_temporal_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 245 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_temporal_cell_data_int32_vector(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: array(:,:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_temporal_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 256 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 244 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 245 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_temporal_cell_data_int64_scalar(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: array(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_temporal_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 245 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_temporal_cell_data_int64_vector(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: array(:,:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_temporal_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 256 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 244 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 245 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_temporal_cell_data_real32_scalar(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: array(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_temporal_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 245 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_temporal_cell_data_real32_vector(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: array(:,:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_temporal_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 256 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 244 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 245 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_temporal_cell_data_real64_scalar(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_temporal_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 245 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_temporal_cell_data_real64_vector(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(:,:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_temporal_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 256 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 257 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 243 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 244 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 245 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_temporal_point_data_int32_scalar(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: array(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_temporal_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 245 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_temporal_point_data_int32_vector(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int32), intent(in) :: array(:,:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_temporal_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 256 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 244 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 245 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_temporal_point_data_int64_scalar(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: array(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_temporal_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 245 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_temporal_point_data_int64_vector(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    integer(int64), intent(in) :: array(:,:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_temporal_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 256 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 244 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 245 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_temporal_point_data_real32_scalar(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: array(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_temporal_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 245 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_temporal_point_data_real32_vector(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real32), intent(in) :: array(:,:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_temporal_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 256 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 244 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 245 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_temporal_point_data_real64_scalar(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_temporal_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 245 "/app/src/vtkhdf_ug_file_type.F90.fypp"

  subroutine write_temporal_point_data_real64_vector(this, name, array)
    class(vtkhdf_ug_file), intent(in) :: this
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(:,:)
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%ug%write_temporal_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 256 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 257 "/app/src/vtkhdf_ug_file_type.F90.fypp"
# 258 "/app/src/vtkhdf_ug_file_type.F90.fypp"

end module vtkhdf_ug_file_type
