!!
!! VTKHDF_CTX_TYPE
!!
!! An internal abstraction layer that encapsulates the collective operations
!! used by the VTKHDF modules. The public interface is identical in serial and
!! MPI-parallel builds; this file provides the MPI-parallel implementation.
!!
!! Copyright (c) 2026 Neil Carlson <neil.n.carlson@gmail.com>
!! SPDX-License-Identifier: BSD-2-Clause
!!

#include "vtkhdf_assert.inc"

module vtkhdf_ctx_type

  use,intrinsic :: iso_fortran_env
  use mpi
  implicit none
  private

  type, public :: vtkhdf_ctx
    integer :: comm = MPI_COMM_NULL
    integer :: rank = 0
    integer :: size = 1
  contains
    procedure :: init
    procedure :: free
    procedure :: global_any, global_all
    generic :: global_sum => global_sum_int32_0, global_sum_int64_0
    generic :: global_min => global_min_int32_1, global_min_int64_1
    procedure, private :: global_sum_int32_0, global_sum_int64_0
    procedure, private :: global_min_int32_1, global_min_int64_1
    generic :: scan_sum => scan_sum_int64_0
    procedure, private :: scan_sum_int64_0
    procedure :: get_dxpl
    procedure :: get_fapl
  end type

contains

  subroutine init(this, comm)
    class(vtkhdf_ctx), intent(out) :: this
    integer, intent(in) :: comm
    integer :: ierr
    call MPI_Comm_dup(comm, this%comm, ierr)
    INSIST(ierr == 0)
    call MPI_Comm_rank(this%comm, this%rank, ierr)
    call MPI_Comm_size(this%comm, this%size, ierr)
  end subroutine

  subroutine free(this)
    class(vtkhdf_ctx), intent(inout) :: this
    integer :: ierr
    if (this%comm /= MPI_COMM_NULL) call MPI_Comm_free(this%comm, ierr)
    this%comm = MPI_COMM_NULL
    this%rank = 0
    this%size = 1
  end subroutine

  !! Returns true if MASK is true on at least one rank; otherwise returns false.
  logical function global_any(this, mask)
    class(vtkhdf_ctx), intent(in) :: this
    logical, intent(in) :: mask
    integer :: local, global, ierr
    local = merge(1, 0, mask) ! 1 is T, 0 is F
    call MPI_Allreduce(local, global, 1, MPI_INTEGER, MPI_MAX, this%comm, ierr)
    INSIST(ierr == MPI_SUCCESS)
    global_any = (global == 1)
  end function

  !! Returns true if MASK is true on all ranks; otherwise returns false.
  logical function global_all(this, mask)
    class(vtkhdf_ctx), intent(in) :: this
    logical, intent(in) :: mask
    integer :: local, global, ierr
    local = merge(1, 0, mask) ! 1 is T, 0 is F
    call MPI_Allreduce(local, global, 1, MPI_INTEGER, MPI_MIN, this%comm, ierr)
    INSIST(ierr == MPI_SUCCESS)
    global_all = (global == 1)
  end function

  subroutine global_sum_int32_0(this, send, recv)
    class(vtkhdf_ctx), intent(in) :: this
    integer(int32), intent(in)  :: send
    integer(int32), intent(out) :: recv
    integer :: ierr
    call MPI_Allreduce(send, recv, 1, MPI_INT32_T, MPI_SUM, this%comm, ierr)
  end subroutine

  subroutine global_sum_int64_0(this, send, recv)
    class(vtkhdf_ctx), intent(in) :: this
    integer(int64), intent(in)  :: send
    integer(int64), intent(out) :: recv
    integer :: ierr
    call MPI_Allreduce(send, recv, 1, MPI_INT64_T, MPI_SUM, this%comm, ierr)
  end subroutine

  subroutine global_min_int32_1(this, send, recv)
    class(vtkhdf_ctx), intent(in) :: this
    integer(int32), intent(in)  :: send(:)
    integer(int32), intent(out) :: recv(:)
    integer :: ierr
    INSIST(size(send) == size(recv))
    call MPI_Allreduce(send, recv, size(send), MPI_INT32_T, MPI_MIN, this%comm, ierr)
    INSIST(ierr == MPI_SUCCESS)
  end subroutine

  subroutine global_min_int64_1(this, send, recv)
    class(vtkhdf_ctx), intent(in) :: this
    integer(int64), intent(in)  :: send(:)
    integer(int64), intent(out) :: recv(:)
    integer :: ierr
    INSIST(size(send) == size(recv))
    call MPI_Allreduce(send, recv, size(send), MPI_INT64_T, MPI_MIN, this%comm, ierr)
    INSIST(ierr == MPI_SUCCESS)
  end subroutine

  subroutine scan_sum_int64_0(this, send, recv)
    class(vtkhdf_ctx), intent(in) :: this
    integer(int64), intent(in)  :: send
    integer(int64), intent(out) :: recv
    integer :: ierr
    call MPI_Scan(send, recv, 1, MPI_INT64_T, MPI_SUM, this%comm, ierr)
    INSIST(ierr == MPI_SUCCESS)
  end subroutine

  !! Returns a data transfer property list configured for collective I/O
  subroutine get_dxpl(this, dxpl)
    use vtkhdf_h5_c_binding
    class(vtkhdf_ctx), intent(in) :: this
    integer(hid_t), intent(out) :: dxpl
    integer :: ierr
    dxpl = H5Pcreate(H5P_DATASET_XFER)
    INSIST(dxpl > 0)
    ierr = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE)
    INSIST(ierr == 0)
  end subroutine

  !! Returns a file access property list configured for parallel I/0
  subroutine get_fapl(this, fapl)
    use vtkhdf_h5_c_binding
    class(vtkhdf_ctx), intent(in) :: this
    integer(hid_t), intent(out) :: fapl
    integer :: stat
    fapl = H5Pcreate(H5P_FILE_ACCESS)
    stat = H5Pset_fapl_mpio(fapl, this%comm)
    stat = H5Pset_all_coll_metadata_ops(fapl, is_collective=.true._c_bool)
    stat = H5Pset_coll_metadata_write(fapl, is_collective=.true._c_bool)
  end subroutine

end module vtkhdf_ctx_type
