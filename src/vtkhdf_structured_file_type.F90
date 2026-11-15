!!
!! VTKHDF_STRUCTURED_FILE_TYPE
!!
!! Standalone file ownership for the shared structured dataset writer.
!! Dataset writing and temporal state are inherited from the writer; only
!! CREATE, FLUSH, and file cleanup belong here. Explicit CLOSE is collective.
!! Implicit finalization frees local allocations without performing I/O.
!!
!! Copyright (c) 2026 Neil Carlson <neil.n.carlson@gmail.com>
!! SPDX-License-Identifier: BSD-2-Clause
!!

module vtkhdf_structured_file_type

  use,intrinsic :: iso_fortran_env, only: int64
  use vtkhdf_structured_writer_type, only: vtkhdf_structured_writer
  use vtkhdf_file_context_type
  use vtkhdf_h5_c_binding, only: init_hdf5
  use vtkhdf_h5e_context_type
  use vtkhdf_data_handle_types
  use vtkhdf_version_param
  use vtkhdf_file_modes, only: VTKHDF_STATIC, VTKHDF_FIXED_MESH, VTKHDF_MOVING_MESH
  implicit none
  private

  public :: vtkhdf_version, VTKHDF_STATIC, VTKHDF_FIXED_MESH, VTKHDF_MOVING_MESH
  public :: vtkhdf_cell_data_handle, vtkhdf_point_data_handle, vtkhdf_field_data_handle

  type, extends(vtkhdf_structured_writer), public :: vtkhdf_structured_file
    private
    type(vtkhdf_file_context) :: file
  contains
    procedure :: close, flush
#ifdef USE_MPI_F08
    generic :: create => create_f08, create_int
    procedure, private :: create_f08, create_int
#elif defined(USE_MPI)
    generic :: create => create_int
    procedure, private :: create_int
#else
    procedure :: create
#endif
  end type

contains

#ifdef USE_MPI_F08
  subroutine create_f08(this, filename, comm, stat, errmsg, mode, data_chunk_bytes)
    use mpi_f08, only: MPI_Comm
    class(vtkhdf_structured_file), intent(out) :: this
    character(*), intent(in) :: filename
    type(MPI_Comm), intent(in) :: comm
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: mode
    integer(int64), intent(in), optional :: data_chunk_bytes
    type(vtkhdf_h5e_context) :: err_ctx
    call init_hdf5
    call err_ctx%capture
    call this%create(filename, comm%mpi_val, stat, errmsg, mode, data_chunk_bytes)
    call err_ctx%restore(stat)
  end subroutine
#endif
#ifdef USE_MPI
  subroutine create_int(this, filename, comm, stat, errmsg, mode, data_chunk_bytes)
    class(vtkhdf_structured_file), intent(out) :: this
    character(*), intent(in) :: filename
    integer, intent(in) :: comm
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: mode
    integer(int64), intent(in), optional :: data_chunk_bytes
    type(vtkhdf_h5e_context) :: err_ctx
    call init_hdf5
    call err_ctx%capture
    call this%file%open(filename, comm, stat, errmsg)
    if (stat == 0) call this%init(this%file%id, 'VTKHDF', this%file%ctx, mode, data_chunk_bytes)
    call err_ctx%restore(stat)
  end subroutine
#else
  subroutine create(this, filename, stat, errmsg, mode, data_chunk_bytes)
    class(vtkhdf_structured_file), intent(out) :: this
    character(*), intent(in) :: filename
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: mode
    integer(int64), intent(in), optional :: data_chunk_bytes
    type(vtkhdf_h5e_context) :: err_ctx
    call init_hdf5
    call err_ctx%capture
    call this%file%open(filename, stat, errmsg)
    if (stat == 0) call this%init(this%file%id, 'VTKHDF', this%file%ctx, mode, data_chunk_bytes)
    call err_ctx%restore(stat)
  end subroutine
#endif

  subroutine flush(this)
    class(vtkhdf_structured_file), intent(in) :: this
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%file%flush
    call err_ctx%restore(0)
  end subroutine

  subroutine close(this)
    class(vtkhdf_structured_file), intent(inout) :: this
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%close_dataset
    call this%file%close
    call reset(this)
    call err_ctx%restore(0)
  contains
    subroutine reset(this)
      class(vtkhdf_structured_file), intent(out) :: this
    end subroutine
  end subroutine close

end module vtkhdf_structured_file_type
