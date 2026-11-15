!!
!! VTKHDF_DATA_HANDLE_TYPES
!!
!! This module defines the cell, point, and field data handle types shared by
!! simple writers and composite blocks. Handles identify data arrays within
!! their owning VTK dataset.
!!
!! SPDX-License-Identifier: BSD-2-Clause
!!

module vtkhdf_data_handle_types

  use,intrinsic :: iso_fortran_env, only: int64
  implicit none
  private

  type, abstract :: vtkhdf_data_handle
    private
    integer :: index_ = -1
    ! Root HDF5 group ID of the owning VTK dataset (an individual leaf block
    ! within a MultiBlock file). Writers use this only to check
    ! handle ownership during the open writer's lifetime.
    integer(int64) :: owner = -1
  contains
    procedure :: set, index
  end type

  type, extends(vtkhdf_data_handle), public :: vtkhdf_cell_data_handle
  end type

  type, extends(vtkhdf_data_handle), public :: vtkhdf_point_data_handle
  end type

  type, extends(vtkhdf_data_handle), public :: vtkhdf_field_data_handle
  end type

contains

  subroutine set(this, index, owner)
    class(vtkhdf_data_handle), intent(out) :: this
    integer, intent(in) :: index
    integer(int64), intent(in) :: owner
    this%index_ = index
    this%owner = owner
  end subroutine

  pure function index(this, owner) result(value)
    class(vtkhdf_data_handle), intent(in) :: this
    integer(int64), intent(in) :: owner
    integer :: value
    value = this%index_
    if (owner /= this%owner) value = -1
  end function

end module vtkhdf_data_handle_types
