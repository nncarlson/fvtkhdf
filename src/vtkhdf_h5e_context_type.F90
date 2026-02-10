!!
!! VTKHDF_H5E_CONTEXT_TYPE
!!
!! This module defines a small internal derived type that acts as a context
!! guard for the HDF5 error state.
!!
!! * CAPTURE saves the current HDF5 error handler, disables automatic error
!!   printing, and clears the current error stack.
!! * RESTORE reinstates the saved handler and, when STAT == 0, clears the error
!!   stack so successful returns do not leave stale errors.
!! * CLEAR explicitly discards expected errors while the context is active.
!!
!! Intended usage is to manage the HDF5 error context at public API boundaries:
!! call CAPTURE on entry and call RESTORE before each return. Internal
!! procedures inherit the active context and need not manage it. On unsuccessful
!! returns (STAT /= 0), the error stack is preserved for the caller to inspect.
!!
!! Copyright (c) 2026 Neil Carlson <neil.n.carlson@gmail.com>
!! SPDX-License-Identifier: BSD-2-Clause
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! Typical usage:
!!
!!   type(vtkhdf_h5e_context) :: err_ctx
!!   integer :: stat
!!
!!   call err_ctx%capture()
!!   ...
!!   call try_open_attribute(stat)
!!   if (stat /= 0) then
!!     call err_ctx%clear()      ! expected failure
!!     call create_attribute(stat)
!!   end if
!!   ...
!!   if (stat /= 0) then
!!     call err_ctx%restore(stat)
!!     return
!!   end if
!!   ...
!!   call err_ctx%restore(0)
!!   return
!!

module vtkhdf_h5e_context_type

  use vtkhdf_h5_c_binding
  implicit none
  private

  type, public :: vtkhdf_h5e_context
    private
    type(c_funptr) :: old_func = c_null_funptr
    type(c_ptr) :: old_data = c_null_ptr
  contains
    procedure :: capture
    procedure :: restore
    procedure :: clear
  end type

contains

  subroutine capture(this)
    class(vtkhdf_h5e_context), intent(inout) :: this
    integer :: ierr
    ierr = H5Eget_auto2(H5E_DEFAULT, this%old_func, this%old_data)
    ierr = H5Eset_auto2(H5E_DEFAULT, C_NULL_FUNPTR, C_NULL_PTR)
    ierr = H5Eclear2(H5E_DEFAULT)
  end subroutine

  subroutine restore(this, stat)
    class(vtkhdf_h5e_context), intent(inout) :: this
    integer, intent(in) :: stat
    integer :: ierr
    if (stat == 0) ierr = H5Eclear2(H5E_DEFAULT)
    ierr = H5Eset_auto2(H5E_DEFAULT, this%old_func, this%old_data)
  end subroutine

  subroutine clear(this)
    class(vtkhdf_h5e_context), intent(inout) :: this
    integer :: ierr
    ierr = H5Eclear2(H5E_DEFAULT)
  end subroutine

end module vtkhdf_h5e_context_type
