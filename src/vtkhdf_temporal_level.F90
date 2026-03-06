!!
!! VTKHDF_TEMPORAL_LEVEL
!!
!! This module defines an enumeration for VTKHDF temporal mesh levels.
!!
!! Copyright (c) 2026 Neil Carlson <neil.n.carlson@gmail.com>
!! SPDX-License-Identifier: BSD-2-Clause
!!

module vtkhdf_temporal_level
  use,intrinsic :: iso_fortran_env, only: int8
  implicit none
  private

  public :: VTKHDF_STATIC_MESH
  public :: VTKHDF_DEFORMED_MESH
  public :: VTKHDF_TEMPORAL_MESH

  integer(int8), parameter :: VTKHDF_STATIC_MESH = 1_int8
  integer(int8), parameter :: VTKHDF_DEFORMED_MESH = 2_int8
  integer(int8), parameter :: VTKHDF_TEMPORAL_MESH = 3_int8

end module vtkhdf_temporal_level
