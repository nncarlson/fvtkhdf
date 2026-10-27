!!
!! VTKHDF_FILE_MODES
!!
!! This module defines the shared writer mode constants. Each dataset format
!! supports an appropriate subset of the static, fixed-mesh, moving-mesh, and
!! dynamic-mesh modes.
!!
!! Copyright (c) 2026 Neil Carlson <neil.n.carlson@gmail.com>
!! SPDX-License-Identifier: BSD-2-Clause
!!

module vtkhdf_file_modes
  implicit none
  private
  integer, parameter, public :: VTKHDF_STATIC       = 0
  integer, parameter, public :: VTKHDF_FIXED_MESH   = 1
  integer, parameter, public :: VTKHDF_MOVING_MESH  = 2
  integer, parameter, public :: VTKHDF_DYNAMIC_MESH = 3
end module vtkhdf_file_modes
