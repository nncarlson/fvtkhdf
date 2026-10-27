!!
!! VTKHDF_DATA_NAMES
!!
!! This module provides the common naming policy for VTK data fields:
!! external-name defaults, sanitization, and collision resolution. A common
!! parent type allows name searches in the mesh-specific field registries.
!!
!! Copyright (c) 2026 Neil Carlson <neil.n.carlson@gmail.com>
!! SPDX-License-Identifier: BSD-2-Clause
!!

module vtkhdf_data_names

  implicit none
  private

  public :: vtkhdf_named_item, resolve_data_name, sanitized_name

  !! Common parent of field records, allowing name searches in either registry.
  !! Each writer retains ownership of its records and their storage metadata.
  type :: vtkhdf_named_item
    character(:), allocatable :: dataset_name
  end type

contains

  !! Normalize the external label and choose an unused internal dataset name.
  !! An unallocated actual array for ITEMS is treated as an absent argument.
  subroutine resolve_data_name(name, default_name, external_name, internal_name, items)

    character(*), intent(in) :: name, default_name
    character(:), allocatable, intent(out) :: external_name, internal_name
    class(vtkhdf_named_item), intent(in), optional :: items(:)

    character(:), allocatable :: base
    character(32) :: string
    integer :: j, suffix

    external_name = trim(name)
    if (external_name == '') external_name = default_name
    base = sanitized_name(external_name)
    internal_name = base
    if (.not.present(items)) return
    j = 1
    suffix = 1
    do while (j <= size(items))
      if (internal_name == items(j)%dataset_name) then
        write(string,'(i0)') suffix
        internal_name = base//'_'//trim(string)
        suffix = suffix + 1
        j = 1
      else
        j = j + 1
      end if
    end do

  end subroutine resolve_data_name

  !! Replace blacklisted characters with underscore, following vtkHDFWriter.
  pure function sanitized_name(name) result(clean_name)
    character(*), intent(in) :: name
    character(len_trim(name)) :: clean_name
    integer :: j
    clean_name = trim(name)
    do j = 1, len(clean_name)
      if (scan(clean_name(j:j),'/. ') > 0) clean_name(j:j) = '_'
    end do
  end function

end module vtkhdf_data_names
