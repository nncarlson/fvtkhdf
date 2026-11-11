!!
!! VTKHDF_DATA_SUPPORT
!!
!! This module defines common record and HDF5 group types shared by the
!! UnstructuredGrid, structured, MultiBlock, and FieldData implementations.
!! It also provides name resolution and dataset chunk sizing.
!!
!! Copyright (c) 2026 Neil Carlson <neil.n.carlson@gmail.com>
!! SPDX-License-Identifier: BSD-2-Clause
!!

module vtkhdf_data_support

  use,intrinsic :: iso_fortran_env
  use vtkhdf_h5_c_binding
  implicit none
  private

  public :: STEPS_CHUNK_SIZE
  public :: resolve_record_name, data_chunk_size

  integer(int64), parameter :: STEPS_CHUNK_SIZE = 512

  !! Base record for members of one naming namespace: the blocks of a composite
  !! VTK dataset, or one kind of dataset (cell, point, or field) of a simple
  !! VTK dataset. Users keep an array of records, possibly of an extended type,
  !! which RESOLVE_RECORD_NAME searches to choose an unused internal name.

  type, public :: named_record
    character(:), allocatable :: internal_name
  end type

  !! Record for a cell, point, or field HDF5 dataset, adding temporal write
  !! tracking, offsets, and shape information.

  type, extends(named_record), public :: dataset_record
    logical :: is_temporal = .false., wrote_this_step = .false.
    integer(int64) :: offset = 0, next_offset = 0
    integer(int32) :: ncomp = 1, tuples = 0
  contains
    procedure :: reset_step_flags
  end type

  !! Common base for the pair of HDF5 group handles of one data association.
  !! GROUP_ID identifies CellData, PointData, or FieldData; OFFSETS_ID identifies
  !! its temporal-offset group under Steps, and remains invalid for static files.
  !! CLOSE releases only these handles; extensions manage additional resources.

  type, public :: data_groups
    integer(hid_t) :: group_id = H5I_INVALID_HID, offsets_id = H5I_INVALID_HID
  contains
    procedure :: close
  end type

contains

  !! Normalize the external label and choose an unused internal HDF5 name.
  !! Trim NAME and use DEFAULT_NAME if blank. Replace '/', '.', and spaces with
  !! underscores, then append '_1', '_2', etc. as needed to avoid existing names.
  !! ITEMS is not modified; the caller reserves the name by appending a record
  !! and preserves the external label in the dataset or block's Name attribute.
  !! An unallocated actual array for ITEMS is treated as an absent argument.

  subroutine resolve_record_name(name, default_name, external_name, internal_name, items)

    character(*), intent(in) :: name, default_name
    character(:), allocatable, intent(out) :: external_name, internal_name
    class(named_record), intent(in), optional :: items(:)

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

  subroutine close(this)
    class(data_groups), intent(inout) :: this
    integer :: ierr
    if (H5Iis_valid(this%offsets_id) > 0) ierr = H5Gclose(this%offsets_id)
    if (H5Iis_valid(this%group_id) > 0) ierr = H5Gclose(this%group_id)
    this%offsets_id = H5I_INVALID_HID
    this%group_id = H5I_INVALID_HID
  end subroutine

  !! Clear temporal write tracking; retain offsets and shapes for omission reuse.
  elemental subroutine reset_step_flags(this)
    class(dataset_record), intent(inout) :: this
    if (this%is_temporal) this%wrote_this_step = .false.
  end subroutine

  pure function data_chunk_size(data_chunk_bytes, value_bytes) result(chunk_size)
    integer(int64), intent(in) :: data_chunk_bytes, value_bytes
    integer(int64) :: chunk_size
    chunk_size = max(1_int64, data_chunk_bytes/value_bytes)
  end function

end module vtkhdf_data_support
