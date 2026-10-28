!!
!! VTKHDF_NAMED_RECORD_TYPE
!!
!! This module provides the common naming policy for VTK data fields and
!! MultiBlock Assembly links: external-name defaults, sanitization, and
!! collision resolution within an HDF5 group.
!!
!! Copyright (c) 2026 Neil Carlson <neil.n.carlson@gmail.com>
!! SPDX-License-Identifier: BSD-2-Clause
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! VTKHDF_NAMED_RECORD represents an object participating in an HDF5 naming
!! namespace. Its INTERNAL_NAME component identifies a dataset or link
!! within a group. Extensions add the metadata needed by each kind of object.
!! UnstructuredGrid and structured writers keep arrays of records for their
!! cell, point, and field datasets; MultiBlock keeps an array of block records.
!! The common parent lets the naming routines search any of these arrays.
!! Each array also serves as the catalog of names already used in its group:
!! CellData, PointData, or FieldData for data fields, or Assembly for blocks.
!! Data-field catalogs belong to an individual VTK dataset (a leaf in a
!! MultiBlock file). Block internal names identify Assembly soft links; the
!! corresponding block groups have separate generated names such as Group_0.
!!
!! Naming has two parts. First, normalize the user-supplied label by trimming
!! it and supplying a default if empty. This external label is preserved by
!! the caller in the dataset's Name attribute, or the block's Name attribute
!! for MultiBlock. The private SANITIZED_NAME helper replaces '/', '.', and
!! spaces with underscores to form a candidate internal HDF5 name.
!! Second, disambiguate that candidate against the group's existing names
!! by appending '_1', '_2', etc. until an unused name is found. This leaves
!! the external label unchanged.
!!
!! RESOLVE_RECORD_NAME combines normalization and disambiguation. It examines
!! the supplied records without modifying them; the caller reserves the name
!! by appending a record to its own array after resolution.
!!

module vtkhdf_named_record_type

  implicit none
  private

  !! Common parent of records whose internal names share an HDF5 namespace.
  !! Each writer retains ownership of its records and their metadata.
  type, public :: vtkhdf_named_record
    character(:), allocatable :: internal_name
  end type

  public :: resolve_record_name

contains

  !! Normalize the external label and choose an unused internal HDF5 name.
  !! An unallocated actual array for ITEMS is treated as an absent argument.
  subroutine resolve_record_name(name, default_name, external_name, internal_name, items)

    character(*), intent(in) :: name, default_name
    character(:), allocatable, intent(out) :: external_name, internal_name
    class(vtkhdf_named_record), intent(in), optional :: items(:)

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
      if (internal_name == items(j)%internal_name) then
        write(string,'(i0)') suffix
        internal_name = base//'_'//trim(string)
        suffix = suffix + 1
        j = 1
      else
        j = j + 1
      end if
    end do

  end subroutine resolve_record_name

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

end module vtkhdf_named_record_type
