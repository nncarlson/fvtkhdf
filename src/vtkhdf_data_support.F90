!!
!! VTKHDF_DATA_SUPPORT
!!
!! This module provides common data records and helpers shared by the
!! UnstructuredGrid, structured, and FieldData implementations: dataset name
!! attributes, temporal write-flag reset, HDF5 group cleanup, and chunk sizing.
!!
!! Copyright (c) 2026 Neil Carlson <neil.n.carlson@gmail.com>
!! SPDX-License-Identifier: BSD-2-Clause
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! Cell and point payload storage belongs to the mesh-specific cores.
!! FieldData records, storage, and temporal metadata are implemented in
!! vtkhdf_field_data_type.
!!

#include "vtkhdf_assert.inc"

module vtkhdf_data_support

  use,intrinsic :: iso_fortran_env
  use vtkhdf_h5_c_binding
  use vtkhdf_h5
  use vtkhdf_ctx_type
  use vtkhdf_named_record_type, only: vtkhdf_named_record
  implicit none
  private

  public :: data_var_item, data_groups, STEPS_CHUNK_SIZE
  public :: write_data_name_attr, reset_data_step_flags, close_data_groups
  public :: data_chunk_size

  integer(int64), parameter :: STEPS_CHUNK_SIZE = 512

  type, extends(vtkhdf_named_record) :: data_var_item
    logical :: is_temporal = .false., wrote_this_step = .false.
    integer(int64) :: offset = 0, next_offset = 0
    integer(int32) :: ncomp = 1, tuples = 0
  end type

  !! HDF5 groups for one data association; each child supplies its records.
  type :: data_groups
    integer(hid_t) :: group_id = H5I_INVALID_HID, offsets_id = H5I_INVALID_HID
  end type

contains

  subroutine close_data_groups(data)
    class(data_groups), intent(inout) :: data
    integer :: ierr
    if (H5Iis_valid(data%offsets_id) > 0) ierr = H5Gclose(data%offsets_id)
    if (H5Iis_valid(data%group_id) > 0) ierr = H5Gclose(data%group_id)
    data%offsets_id = H5I_INVALID_HID
    data%group_id = H5I_INVALID_HID
  end subroutine

  subroutine reset_data_step_flags(vars)
    class(data_var_item), intent(inout), optional :: vars(:)
    integer :: j
    if (.not.present(vars)) return
    do j = 1, size(vars)
      if (vars(j)%is_temporal) vars(j)%wrote_this_step = .false.
    end do
  end subroutine

  subroutine write_data_name_attr(ctx, grp_id, dataset_name, external_name, attribute)
    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: grp_id
    character(*), intent(in) :: dataset_name, external_name
    character(*), intent(in), optional :: attribute
    integer(hid_t) :: dset_id
    integer :: ierr
    dset_id = H5Dopen(grp_id, dataset_name)
    INSIST(ctx%global_all(dset_id >= 0))
    call h5_write_attr(ctx, dset_id, 'Name', external_name)
    if (present(attribute)) then
      INSIST(ctx%identical(trim(attribute)))
      if (attribute /= '') call h5_write_attr(ctx, dset_id, 'Attribute', trim(attribute))
    end if
    ierr = H5Dclose(dset_id)
    INSIST(ctx%global_all(ierr >= 0))
  end subroutine

  pure function data_chunk_size(data_chunk_bytes, value_bytes) result(chunk_size)
    integer(int64), intent(in) :: data_chunk_bytes, value_bytes
    integer(int64) :: chunk_size
    chunk_size = max(1_int64, data_chunk_bytes/value_bytes)
  end function
end module vtkhdf_data_support
