!!
!! VTKHDF_CTX_TYPE
!!
!! An internal abstraction layer that encapsulates the collective operations
!! used by the VTKHDF modules. The public interface is identical in serial and
!! MPI-parallel builds; this file provides the serial implementation with
!! semantically equivalent, non-MPI behavior.
!!
!! Copyright (c) 2026 Neil Carlson <neil.n.carlson@gmail.com>
!! SPDX-License-Identifier: BSD-2-Clause
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! In serial builds, collective operations degenerate to local operations that
!! preserve the same semantics as the MPI implementation.
!!

#include "vtkhdf_assert.inc"

module vtkhdf_ctx_type

  use,intrinsic :: iso_fortran_env
  implicit none
  private

  type, public :: vtkhdf_ctx
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

  subroutine init(this)
    class(vtkhdf_ctx), intent(out) :: this
  end subroutine

  subroutine free(this)
    class(vtkhdf_ctx), intent(inout) :: this
    this%rank = 0
    this%size = 1
  end subroutine

  !! Returns true if MASK is true on at least one rank; otherwise returns false.
  logical function global_any(this, mask)
    class(vtkhdf_ctx), intent(in) :: this
    logical, intent(in) :: mask
    global_any = mask
  end function

  !! Returns true if MASK is true on all ranks; otherwise returns false.
  logical function global_all(this, mask)
    class(vtkhdf_ctx), intent(in) :: this
    logical, intent(in) :: mask
    global_all = mask
  end function

  subroutine global_sum_int32_0(this, send, recv)
    class(vtkhdf_ctx), intent(in) :: this
    integer(int32), intent(in)  :: send
    integer(int32), intent(out) :: recv
    recv = send
  end subroutine

  subroutine global_sum_int64_0(this, send, recv)
    class(vtkhdf_ctx), intent(in) :: this
    integer(int64), intent(in)  :: send
    integer(int64), intent(out) :: recv
    recv = send
  end subroutine

  subroutine global_min_int32_1(this, send, recv)
    class(vtkhdf_ctx), intent(in) :: this
    integer(int32), intent(in)  :: send(:)
    integer(int32), intent(out) :: recv(:)
    recv = send
  end subroutine

  subroutine global_min_int64_1(this, send, recv)
    class(vtkhdf_ctx), intent(in) :: this
    integer(int64), intent(in)  :: send(:)
    integer(int64), intent(out) :: recv(:)
    recv = send
  end subroutine

  subroutine scan_sum_int64_0(this, send, recv)
    class(vtkhdf_ctx), intent(in) :: this
    integer(int64), intent(in)  :: send
    integer(int64), intent(out) :: recv
    recv = send
  end subroutine

  !! Returns a default data transfer property list (INDEPENDENT)
  subroutine get_dxpl(this, dxpl)
    use vtkhdf_h5_c_binding
    class(vtkhdf_ctx), intent(in) :: this
    integer(hid_t), intent(out) :: dxpl
    dxpl = H5Pcreate(H5P_DATASET_XFER)
  end subroutine

  !! Returns a default file access property lis
  subroutine get_fapl(this, fapl)
    use vtkhdf_h5_c_binding
    class(vtkhdf_ctx), intent(in) :: this
    integer(hid_t), intent(out) :: fapl
    fapl = H5Pcreate(H5P_FILE_ACCESS)
  end subroutine

end module vtkhdf_ctx_type
