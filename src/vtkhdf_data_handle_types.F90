!!
!! VTKHDF_DATA_HANDLE_TYPES
!!
!! This module defines the cell, point, and field data handle types shared by
!! the file interfaces. Handles identify fields; block placement is supplied
!! separately at each write.
!!
!! Copyright (c) 2026 Neil Carlson <neil.n.carlson@gmail.com>
!! SPDX-License-Identifier: BSD-2-Clause
!!

module vtkhdf_data_handle_types

  use,intrinsic :: iso_fortran_env, only: int64
  implicit none
  private

  type, abstract :: vtkhdf_data_handle
    private
    integer :: index = -1
    ! BLOCK_INDEX is interpreted only by the MultiBlock front end.
    integer :: block_index = -1
    ! Root HDF5 group ID of the owning VTK dataset (an individual leaf block
    ! within a MultiBlock file). File wrappers use this only to check
    ! handle ownership during the open writer's lifetime.
    integer(int64) :: owner = -1
  contains
    procedure :: set, get_index, set_block_index, get_block_index
  end type

  type, extends(vtkhdf_data_handle), public :: vtkhdf_cell_data_handle
  end type

  type, extends(vtkhdf_data_handle), public :: vtkhdf_point_data_handle
  end type

  type, extends(vtkhdf_data_handle), public :: vtkhdf_field_data_handle
  end type

contains

  subroutine set_block_index(this, block_index)
    class(vtkhdf_data_handle), intent(inout) :: this
    integer, intent(in) :: block_index
    this%block_index = block_index
  end subroutine

  pure function get_block_index(this) result(block_index)
    class(vtkhdf_data_handle), intent(in) :: this
    integer :: block_index
    block_index = this%block_index
  end function

  subroutine set(this, index, owner)
    class(vtkhdf_data_handle), intent(out) :: this
    integer, intent(in) :: index
    integer(int64), intent(in), optional :: owner
    this%index = index
    if (present(owner)) this%owner = owner
  end subroutine

  pure function get_index(this, owner) result(index)
    class(vtkhdf_data_handle), intent(in) :: this
    integer(int64), intent(in), optional :: owner
    integer :: index
    index = this%index
    if (present(owner)) then
      if (owner /= this%owner) index = -1
    end if
  end function

end module vtkhdf_data_handle_types
