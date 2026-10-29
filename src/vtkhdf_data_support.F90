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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! NAMED_RECORD represents an object participating in an HDF5 naming
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
!! Cell and point payload storage belongs to the mesh-specific cores.
!! FieldData records, storage, and temporal metadata are implemented in
!! vtkhdf_field_data_type.
!!
!! DATA_GROUPS is the common base for the pair of HDF5 group handles of one
!! data association. GROUP_ID identifies CellData, PointData, or FieldData
!! beneath the owning VTK dataset's root. OFFSETS_ID identifies the corresponding
!! group beneath Steps (CellDataOffsets, PointDataOffsets, or FieldDataOffsets);
!! it remains invalid for static files. Each group contains datasets named for
!! the association's individual variables.
!!
!! UG and structured each define a MESH_DATA_LIST extension containing their
!! own variable records. FIELD_DATA_LIST in vtkhdf_field_data_type adds its
!! records and the FieldDataSizes group handle. CLOSE releases and invalidates
!! only the two base handles; extensions or their owners close any additional
!! resources separately. FIELD_DATA_LIST overrides CLOSE for that purpose.
!! Its override invokes the parent binding after closing the sizes group.
!!

module vtkhdf_data_support

  use,intrinsic :: iso_fortran_env
  use vtkhdf_h5_c_binding
  implicit none
  private

  public :: STEPS_CHUNK_SIZE
  public :: resolve_record_name, data_chunk_size

  integer(int64), parameter :: STEPS_CHUNK_SIZE = 512

  !! Common parent of records whose internal names share an HDF5 namespace.
  !! Each writer retains ownership of its records and their metadata.
  type, public :: named_record
    character(:), allocatable :: internal_name
  end type

  type, extends(named_record), public :: data_var_record
    logical :: is_temporal = .false., wrote_this_step = .false.
    integer(int64) :: offset = 0, next_offset = 0
    integer(int32) :: ncomp = 1, tuples = 0
  contains
    procedure :: reset_step_flags
  end type

  !! Common base for the pair of HDF5 group handles of one data association.
  type, public :: data_groups
    integer(hid_t) :: group_id = H5I_INVALID_HID, offsets_id = H5I_INVALID_HID
  contains
    procedure :: close
  end type

contains

  !! Normalize the external label and choose an unused internal HDF5 name.
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
    class(data_var_record), intent(inout) :: this
    if (this%is_temporal) this%wrote_this_step = .false.
  end subroutine

  pure function data_chunk_size(data_chunk_bytes, value_bytes) result(chunk_size)
    integer(int64), intent(in) :: data_chunk_bytes, value_bytes
    integer(int64) :: chunk_size
    chunk_size = max(1_int64, data_chunk_bytes/value_bytes)
  end function

end module vtkhdf_data_support
