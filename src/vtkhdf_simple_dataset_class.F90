!! Common lifecycle contract for simple VTKHDF writers and composite leaves.
!! Geometry and data interfaces belong to the UG and structured branches.
!!
!! Copyright (c) 2026 Neil Carlson <neil.n.carlson@gmail.com>
!! SPDX-License-Identifier: BSD-2-Clause
module vtkhdf_simple_dataset_class
  use,intrinsic :: iso_fortran_env, only: real32, real64
  use vtkhdf_h5_c_binding, only: hid_t
  implicit none
  private

  type, abstract, public :: vtkhdf_simple_dataset
  contains
    procedure(close_interface), deferred :: close
    procedure(root_interface), deferred :: root
    procedure(query_interface), deferred :: temporal_enabled, step_is_open
    generic :: start_time_step => start_time_step_real32, start_time_step_real64
    procedure, private :: start_time_step_real32
    procedure(step_interface), deferred :: start_time_step_real64
    procedure(close_interface), deferred :: finalize_time_step
  end type

  abstract interface
    subroutine close_interface(this)
      import vtkhdf_simple_dataset
      class(vtkhdf_simple_dataset), intent(inout) :: this
    end subroutine
    function root_interface(this) result(id)
      import vtkhdf_simple_dataset, hid_t
      class(vtkhdf_simple_dataset), intent(in) :: this
      integer(hid_t) :: id
    end function
    pure function query_interface(this) result(value)
      import vtkhdf_simple_dataset
      class(vtkhdf_simple_dataset), intent(in) :: this
      logical :: value
    end function
    subroutine step_interface(this, time)
      import vtkhdf_simple_dataset, real64
      class(vtkhdf_simple_dataset), intent(inout) :: this
      real(real64), intent(in) :: time
    end subroutine
  end interface

contains

  subroutine start_time_step_real32(this, time)
    class(vtkhdf_simple_dataset), intent(inout) :: this
    real(real32), intent(in) :: time
    call this%start_time_step_real64(real(time,real64))
  end subroutine

end module
