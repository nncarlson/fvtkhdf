# 1 "/app/src/vtkhdf_h5.F90.fypp"
!!
!! VTKHDF_H5
!!
!! An application-specific abstraction over parallel HDF5 that provides
!! single-call procedures for writing attributes and datasets, creating and
!! appending extendable datasets, and other utilities for private use by
!! VTKHDF modules.
!!
!! Copyright (c) 2026 Neil Carlson <neil.n.carlson@gmail.com>
!! SPDX-License-Identifier: BSD-2-Clause
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! All public procedures are collective operations. This module follows a
!! "fail-fast" philosophy: internal logic errors, schema mismatches, or
!! unrecoverable HDF5/filesystem failures trigger an immediate global
!! MPI_Abort (via f90_assert) to prevent deadlocks and ensure file integrity.
!!

#include "vtkhdf_assert.inc"

# 31 "/app/src/vtkhdf_h5.F90.fypp"
# 32 "/app/src/vtkhdf_h5.F90.fypp"

module vtkhdf_h5

  use,intrinsic :: iso_fortran_env
  use vtkhdf_ctx_type
  use vtkhdf_h5_c_binding
  implicit none
  private

  public :: h5_write_attr, h5_create_unlimited_dataset, h5_write_dataset, h5_append_to_dataset

  interface h5_write_attr
# 45 "/app/src/vtkhdf_h5.F90.fypp"
    procedure write_attr_int8
# 45 "/app/src/vtkhdf_h5.F90.fypp"
    procedure write_attr_int32
# 45 "/app/src/vtkhdf_h5.F90.fypp"
    procedure write_attr_int64
# 45 "/app/src/vtkhdf_h5.F90.fypp"
    procedure write_attr_real32
# 45 "/app/src/vtkhdf_h5.F90.fypp"
    procedure write_attr_real64
# 47 "/app/src/vtkhdf_h5.F90.fypp"
    procedure write_attr_string
  end interface

  interface h5_write_dataset
# 52 "/app/src/vtkhdf_h5.F90.fypp"
    procedure write_dataset_int8
# 52 "/app/src/vtkhdf_h5.F90.fypp"
    procedure write_dataset_int32
# 52 "/app/src/vtkhdf_h5.F90.fypp"
    procedure write_dataset_int64
# 52 "/app/src/vtkhdf_h5.F90.fypp"
    procedure write_dataset_real32
# 52 "/app/src/vtkhdf_h5.F90.fypp"
    procedure write_dataset_real64
# 54 "/app/src/vtkhdf_h5.F90.fypp"
  end interface

  interface h5_create_unlimited_dataset
# 58 "/app/src/vtkhdf_h5.F90.fypp"
    procedure create_unlimited_dataset_int32
# 58 "/app/src/vtkhdf_h5.F90.fypp"
    procedure create_unlimited_dataset_int64
# 58 "/app/src/vtkhdf_h5.F90.fypp"
    procedure create_unlimited_dataset_real32
# 58 "/app/src/vtkhdf_h5.F90.fypp"
    procedure create_unlimited_dataset_real64
# 60 "/app/src/vtkhdf_h5.F90.fypp"
  end interface

  interface h5_append_to_dataset
# 64 "/app/src/vtkhdf_h5.F90.fypp"
    procedure append_to_dataset_int32
# 64 "/app/src/vtkhdf_h5.F90.fypp"
    procedure append_to_dataset_int64
# 64 "/app/src/vtkhdf_h5.F90.fypp"
    procedure append_to_dataset_real32
# 64 "/app/src/vtkhdf_h5.F90.fypp"
    procedure append_to_dataset_real64
# 66 "/app/src/vtkhdf_h5.F90.fypp"
  end interface

contains
# 70 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine write_attr_int8(ctx, obj_id, name, value)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: name
    integer(int8), intent(in) :: value(..)

    integer(hid_t) :: type_id, attr_id
    integer :: ierr

    type_id = H5T_NATIVE_UINT8

    call get_attr_id(ctx, obj_id, name, type_id, shape(value, hsize_t), attr_id)

    select rank (value)
    rank (0)
      ierr = H5Awrite(attr_id, type_id, [value])
    rank (1)
      ierr = H5Awrite(attr_id, type_id, value)
    rank (2)
      ierr = H5Awrite(attr_id, type_id, value)
    rank default
      INSIST(.false.)
    end select
    INSIST(ctx%global_all(ierr >= 0))

    ierr = H5Aclose(attr_id)

  end subroutine write_attr_int8

# 70 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine write_attr_int32(ctx, obj_id, name, value)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: name
    integer(int32), intent(in) :: value(..)

    integer(hid_t) :: type_id, attr_id
    integer :: ierr

    type_id = H5T_NATIVE_INT32

    call get_attr_id(ctx, obj_id, name, type_id, shape(value, hsize_t), attr_id)

    select rank (value)
    rank (0)
      ierr = H5Awrite(attr_id, type_id, [value])
    rank (1)
      ierr = H5Awrite(attr_id, type_id, value)
    rank (2)
      ierr = H5Awrite(attr_id, type_id, value)
    rank default
      INSIST(.false.)
    end select
    INSIST(ctx%global_all(ierr >= 0))

    ierr = H5Aclose(attr_id)

  end subroutine write_attr_int32

# 70 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine write_attr_int64(ctx, obj_id, name, value)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: name
    integer(int64), intent(in) :: value(..)

    integer(hid_t) :: type_id, attr_id
    integer :: ierr

    type_id = H5T_NATIVE_INT64

    call get_attr_id(ctx, obj_id, name, type_id, shape(value, hsize_t), attr_id)

    select rank (value)
    rank (0)
      ierr = H5Awrite(attr_id, type_id, [value])
    rank (1)
      ierr = H5Awrite(attr_id, type_id, value)
    rank (2)
      ierr = H5Awrite(attr_id, type_id, value)
    rank default
      INSIST(.false.)
    end select
    INSIST(ctx%global_all(ierr >= 0))

    ierr = H5Aclose(attr_id)

  end subroutine write_attr_int64

# 70 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine write_attr_real32(ctx, obj_id, name, value)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: name
    real(real32), intent(in) :: value(..)

    integer(hid_t) :: type_id, attr_id
    integer :: ierr

    type_id = H5T_NATIVE_FLOAT

    call get_attr_id(ctx, obj_id, name, type_id, shape(value, hsize_t), attr_id)

    select rank (value)
    rank (0)
      ierr = H5Awrite(attr_id, type_id, [value])
    rank (1)
      ierr = H5Awrite(attr_id, type_id, value)
    rank (2)
      ierr = H5Awrite(attr_id, type_id, value)
    rank default
      INSIST(.false.)
    end select
    INSIST(ctx%global_all(ierr >= 0))

    ierr = H5Aclose(attr_id)

  end subroutine write_attr_real32

# 70 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine write_attr_real64(ctx, obj_id, name, value)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: name
    real(real64), intent(in) :: value(..)

    integer(hid_t) :: type_id, attr_id
    integer :: ierr

    type_id = H5T_NATIVE_DOUBLE

    call get_attr_id(ctx, obj_id, name, type_id, shape(value, hsize_t), attr_id)

    select rank (value)
    rank (0)
      ierr = H5Awrite(attr_id, type_id, [value])
    rank (1)
      ierr = H5Awrite(attr_id, type_id, value)
    rank (2)
      ierr = H5Awrite(attr_id, type_id, value)
    rank default
      INSIST(.false.)
    end select
    INSIST(ctx%global_all(ierr >= 0))

    ierr = H5Aclose(attr_id)

  end subroutine write_attr_real64

# 102 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine write_attr_string(ctx, obj_id, name, value)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: name
#ifdef INTEL_BUG20240327
    character(*), intent(in) :: value
#else
    character(*), intent(in) :: value(..)
#endif

    integer(hid_t) :: type_id, attr_id
    integer :: ierr

    type_id = H5Tcopy(H5T_NATIVE_CHARACTER)
    INSIST(type_id > 0)
    ierr = H5Tset_size(type_id, len(value, kind=c_size_t))
    INSIST(ierr >= 0)

    call get_attr_id(ctx, obj_id, name, type_id, shape(value, hsize_t), attr_id)

#ifdef INTEL_BUG20240327
    ierr = H5Awrite(attr_id, type_id, [value])
#else
    select rank (value)
    rank (0)
      ierr = H5Awrite(attr_id, type_id, [value])
    rank (1)
      ierr = H5Awrite(attr_id, type_id, value)
    rank (2)
      ierr = H5Awrite(attr_id, type_id, value)
    rank default
      INSIST(.false.)
    end select
#endif
    INSIST(ctx%global_all(ierr >= 0))

    ierr = H5Aclose(attr_id)
    ierr = H5Tclose(type_id)

  end subroutine write_attr_string


  subroutine get_attr_id(ctx, obj_id, name, type_id, dims, attr_id)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: name
    integer(hid_t), intent(in) :: type_id
    integer(hsize_t), intent(in) :: dims(:) ! the Fortran shape
    integer(hid_t), intent(out) :: attr_id

    integer(hid_t) :: space_id
    integer :: htri, ierr

    htri = H5Aexists(obj_id, name)
    if (ctx%global_all(htri > 0)) then
      attr_id = H5Aopen(obj_id, name)
      INSIST(attr_id > 0)
      block ! Check that the input type matches the file type
        integer(hid_t) :: file_type_id
        file_type_id = H5Aget_type(attr_id)
        INSIST(ctx%global_all(H5Tequal(file_type_id, type_id) > 0))
        ierr = H5Tclose(file_type_id)
      end block
      block ! Check that the input shape matches the file type
        integer(hid_t) :: file_space_id
        integer :: file_ndims
        integer(hsize_t), allocatable :: file_dims(:)
        file_space_id = H5Aget_space(attr_id)
        file_ndims = H5Sget_simple_extent_ndims(file_space_id)
        ! Check rank compatibility (scalar = 0)
        INSIST(ctx%global_all(file_ndims == size(dims)))
        if (file_ndims > 0) then ! Check dimensions if not scalar
          allocate(file_dims(file_ndims))
          ierr = H5Sget_simple_extent_dims(file_space_id, file_dims)
          ! Note: file_dims is C-ordered, dims is Fortran-ordered
          INSIST(ctx%global_all(all(file_dims == dims(size(dims):1:-1))))
        end if
        ierr = H5Sclose(file_space_id)
      end block
      return
    else if (ctx%global_any(htri /= 0)) then
      INSIST(.false.) ! ranks disagree on existence of the attribute
    end if

    !! It doesn't exist, so create it
    select case (size(dims))
    case (0)
      space_id = H5Screate(H5S_SCALAR)
    case (1)
      space_id = H5Screate(dims)
    case (2:)
      block
        integer(hsize_t) :: c_dims(size(dims))
        c_dims = dims(size(dims):1:-1) ! reverse order for the C shape
        space_id = H5Screate(c_dims)
      end block
    end select
    INSIST(space_id >= 0)

    attr_id = H5Acreate(obj_id, name, type_id, space_id)
    INSIST(ctx%global_all(attr_id >= 0))

    ierr = H5Sclose(space_id)

  end subroutine get_attr_id

  !! Write a distributed ARRAY to a new dataset NAME at location LOC_ID in the
  !! HDF5 file tree structure. The written dataset is the concatenation of the
  !! arrays from all the MPI ranks. The rank of ARRAY must be the same on each
  !! MPI rank, and in the multi-rank case, the extent must be the same in all
  !! but the distributed last dimension. If the optional argument ROOT is
  !! present, the behavior is different. In this case ROOT specifies the MPI
  !! rank of the only rank that writes ARRAY, and the written dataset consists
  !! of that data alone. The content of ARRAY on all other processes is ignored,
  !! however the previous constraints on the rank and extents continues to
  !! apply. Moreover this must still be called collectively as creation of the
  !! dataset itself is a collective operation.
# 223 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine write_dataset_int8(ctx, loc_id, name, array, root)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(int8), intent(in) :: array(..)
    integer, intent(in), optional :: root

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id, dxpl_id
    integer :: ierr

    type_id = H5T_NATIVE_UINT8

    call write_dataset_common(ctx, loc_id, name, type_id, shape(array,hsize_t), &
        dset_id, mem_space_id, data_space_id, root)

    call ctx%get_dxpl(dxpl_id) ! collective data transfer in MPI build

    select rank (array)
    rank (0)
      ierr = H5Dwrite(dset_id, type_id, [array], mem_space_id, data_space_id, dxpl_id)
    rank (1)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl_id)
    rank (2)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl_id)
    rank default
      INSIST(.false.)
    end select
    INSIST(ctx%global_any(ierr >= 0))

    !! Cleanup
    ierr = H5Pclose(dxpl_id)
    ierr = H5Sclose(mem_space_id)
    ierr = H5Sclose(data_space_id)
    ierr = H5Dclose(dset_id)

  end subroutine write_dataset_int8

# 223 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine write_dataset_int32(ctx, loc_id, name, array, root)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(int32), intent(in) :: array(..)
    integer, intent(in), optional :: root

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id, dxpl_id
    integer :: ierr

    type_id = H5T_NATIVE_INT32

    call write_dataset_common(ctx, loc_id, name, type_id, shape(array,hsize_t), &
        dset_id, mem_space_id, data_space_id, root)

    call ctx%get_dxpl(dxpl_id) ! collective data transfer in MPI build

    select rank (array)
    rank (0)
      ierr = H5Dwrite(dset_id, type_id, [array], mem_space_id, data_space_id, dxpl_id)
    rank (1)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl_id)
    rank (2)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl_id)
    rank default
      INSIST(.false.)
    end select
    INSIST(ctx%global_any(ierr >= 0))

    !! Cleanup
    ierr = H5Pclose(dxpl_id)
    ierr = H5Sclose(mem_space_id)
    ierr = H5Sclose(data_space_id)
    ierr = H5Dclose(dset_id)

  end subroutine write_dataset_int32

# 223 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine write_dataset_int64(ctx, loc_id, name, array, root)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(int64), intent(in) :: array(..)
    integer, intent(in), optional :: root

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id, dxpl_id
    integer :: ierr

    type_id = H5T_NATIVE_INT64

    call write_dataset_common(ctx, loc_id, name, type_id, shape(array,hsize_t), &
        dset_id, mem_space_id, data_space_id, root)

    call ctx%get_dxpl(dxpl_id) ! collective data transfer in MPI build

    select rank (array)
    rank (0)
      ierr = H5Dwrite(dset_id, type_id, [array], mem_space_id, data_space_id, dxpl_id)
    rank (1)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl_id)
    rank (2)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl_id)
    rank default
      INSIST(.false.)
    end select
    INSIST(ctx%global_any(ierr >= 0))

    !! Cleanup
    ierr = H5Pclose(dxpl_id)
    ierr = H5Sclose(mem_space_id)
    ierr = H5Sclose(data_space_id)
    ierr = H5Dclose(dset_id)

  end subroutine write_dataset_int64

# 223 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine write_dataset_real32(ctx, loc_id, name, array, root)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(real32), intent(in) :: array(..)
    integer, intent(in), optional :: root

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id, dxpl_id
    integer :: ierr

    type_id = H5T_NATIVE_FLOAT

    call write_dataset_common(ctx, loc_id, name, type_id, shape(array,hsize_t), &
        dset_id, mem_space_id, data_space_id, root)

    call ctx%get_dxpl(dxpl_id) ! collective data transfer in MPI build

    select rank (array)
    rank (0)
      ierr = H5Dwrite(dset_id, type_id, [array], mem_space_id, data_space_id, dxpl_id)
    rank (1)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl_id)
    rank (2)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl_id)
    rank default
      INSIST(.false.)
    end select
    INSIST(ctx%global_any(ierr >= 0))

    !! Cleanup
    ierr = H5Pclose(dxpl_id)
    ierr = H5Sclose(mem_space_id)
    ierr = H5Sclose(data_space_id)
    ierr = H5Dclose(dset_id)

  end subroutine write_dataset_real32

# 223 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine write_dataset_real64(ctx, loc_id, name, array, root)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(..)
    integer, intent(in), optional :: root

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id, dxpl_id
    integer :: ierr

    type_id = H5T_NATIVE_DOUBLE

    call write_dataset_common(ctx, loc_id, name, type_id, shape(array,hsize_t), &
        dset_id, mem_space_id, data_space_id, root)

    call ctx%get_dxpl(dxpl_id) ! collective data transfer in MPI build

    select rank (array)
    rank (0)
      ierr = H5Dwrite(dset_id, type_id, [array], mem_space_id, data_space_id, dxpl_id)
    rank (1)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl_id)
    rank (2)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl_id)
    rank default
      INSIST(.false.)
    end select
    INSIST(ctx%global_any(ierr >= 0))

    !! Cleanup
    ierr = H5Pclose(dxpl_id)
    ierr = H5Sclose(mem_space_id)
    ierr = H5Sclose(data_space_id)
    ierr = H5Dclose(dset_id)

  end subroutine write_dataset_real64

# 263 "/app/src/vtkhdf_h5.F90.fypp"

  !! This private procedure extracts the part of WRITE_DATASET common to
  !! all types and ranks of the data array. It does not depend on the array
  !! directly, but only on its shape (DIMS) and its HDF5 datatype (TYPE_ID).
  !! It creates the dataset and the dataspaces for the data array and the
  !! part of the dataset where the array will be written. If ROOT is present,
  !! only that MPI rank will write, but all ranks must still participate
  !! in the subsequent collective H5Dwrite call.

  subroutine write_dataset_common(ctx, loc_id, name, type_id, dims, &
      dset_id, mem_space_id, data_space_id, root)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: loc_id, type_id
    character(*), intent(in) :: name
    integer(hsize_t), intent(in) :: dims(:)
    integer(hid_t), intent(out) :: dset_id, mem_space_id, data_space_id
    integer, intent(in), optional :: root

    integer :: ierr, rank
    integer(hsize_t), allocatable :: mem_dims(:), data_dims(:), start(:), count(:)
    integer(hsize_t) :: n
    logical :: rank_writes

    !! Convert DIMS to its equivalent C form (MEM_DIMS)
    if (size(dims) == 0) then ! scalar data
      mem_dims = [1]  ! equivalent rank-1 array
    else
      mem_dims = dims
      if (size(dims) > 1) mem_dims = dims(size(dims):1:-1) ! reverse order for the C shape
    end if

    ASSERT(all(mem_dims(2:) > 0))

    if (present(root)) then
      ASSERT((root >= 0) .and. (root < ctx%size))
      rank_writes = (root == ctx%rank)
    else
      rank_writes = .true.
    end if

    if (.not.rank_writes) mem_dims(1) = 0 ! ignore their data

    !! Create the dataset
    call ctx%global_sum(mem_dims(1), n)
    INSIST(n > 0) ! writing an empy dataset not allowed
    data_dims = mem_dims
    data_dims(1) = n
    data_space_id = H5Screate(data_dims)
    INSIST(data_space_id > 0)
    dset_id = H5Dcreate(loc_id, name, type_id, data_space_id)
    INSIST(ctx%global_all(dset_id >= 0))

    !! Starting index for the dataset hyperslab for this rank
    allocate(start, mold=mem_dims)
    start = 0
    call ctx%scan_sum(mem_dims(1), n)
    start(1) = n - mem_dims(1)

    if (mem_dims(1) > 0) then
      mem_space_id = H5Screate(mem_dims)
    else ! rank contributes nothing
      mem_space_id = H5Screate(H5S_NULL)
    end if
    INSIST(mem_space_id > 0)

    !! Create the dataspace for the part of dataset to be written by this rank
    if (rank_writes .and. mem_dims(1) > 0) then
      count = mem_dims
      ierr = H5Sselect_hyperslab(data_space_id, H5S_SELECT_SET, start, count)
    else ! rank contributes nothing
      ierr = H5Sselect_none(data_space_id)
    end if
    INSIST(ctx%global_all(ierr >= 0))

  end subroutine write_dataset_common

  !! Create an extendable dataset NAME at location LOC_ID. The dataset is
  !! extendable along the last (slowest-varying) dimension and the rank-1
  !! array MOLD specifies the fixed leading dimensions; the size of MOLD
  !! is therefore one less than the dataset rank. Note that the order of
  !! these Fortran dimensions will be reversed to conform to the C ordering
  !! used by HDF5. The type of MOLD determines the HDF5 datatype. CHUNK_SIZE
  !! specifies the chunk size in terms of the number of elements in the
  !! extendable dimension. All MPI ranks must specify the same values for
  !! MOLD and CHUNK_SIZE. The dataset's initial extent in the last dimension
  !! is 0.
# 351 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine create_unlimited_dataset_int32(ctx, loc_id, name, mold, chunk_size)
    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(int32), intent(in) :: mold(..)
    integer, intent(in) :: chunk_size
    call create_unlimited_dataset_core(ctx, loc_id, name, H5T_NATIVE_INT32, &
        shape(mold,hsize_t), chunk_size)
  end subroutine
# 351 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine create_unlimited_dataset_int64(ctx, loc_id, name, mold, chunk_size)
    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(int64), intent(in) :: mold(..)
    integer, intent(in) :: chunk_size
    call create_unlimited_dataset_core(ctx, loc_id, name, H5T_NATIVE_INT64, &
        shape(mold,hsize_t), chunk_size)
  end subroutine
# 351 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine create_unlimited_dataset_real32(ctx, loc_id, name, mold, chunk_size)
    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(real32), intent(in) :: mold(..)
    integer, intent(in) :: chunk_size
    call create_unlimited_dataset_core(ctx, loc_id, name, H5T_NATIVE_FLOAT, &
        shape(mold,hsize_t), chunk_size)
  end subroutine
# 351 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine create_unlimited_dataset_real64(ctx, loc_id, name, mold, chunk_size)
    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(real64), intent(in) :: mold(..)
    integer, intent(in) :: chunk_size
    call create_unlimited_dataset_core(ctx, loc_id, name, H5T_NATIVE_DOUBLE, &
        shape(mold,hsize_t), chunk_size)
  end subroutine
# 362 "/app/src/vtkhdf_h5.F90.fypp"


  subroutine create_unlimited_dataset_core(ctx, loc_id, name, type_id, dims, chunk_size)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(hid_t), intent(in) :: type_id
    integer(hsize_t), intent(in) :: dims(:)
    integer, intent(in) :: chunk_size

    integer :: ierr, htri
    integer(hid_t) :: space_id, dcpl_id, dset_id
    integer(hsize_t), allocatable :: c_dims(:), maxdims(:), chunk_dims(:)

    htri = H5Lexists(loc_id, name)
    INSIST(ctx%global_all(htri == 0)) ! Should not exist

    c_dims = [0_hsize_t, dims(size(dims):1:-1)] ! first dimension is unlimited

    ASSERT(all(c_dims(2:) > 0))
    ASSERT(chunk_size > 0)

    maxdims = c_dims
    maxdims(1) = H5S_UNLIMITED
    chunk_dims = c_dims
    chunk_dims(1) = chunk_size

    space_id = H5Screate(c_dims, maxdims)
    INSIST(space_id > 0)

    dcpl_id = H5Pcreate(H5P_DATASET_CREATE)
    INSIST(dcpl_id > 0)

    ierr = H5Pset_chunk(dcpl_id, size(chunk_dims), chunk_dims)
    INSIST(ierr >= 0)

    dset_id = H5Dcreate(loc_id, name, type_id, space_id, dcpl_id=dcpl_id)
    INSIST(ctx%global_all(dset_id >= 0))

    ierr = h5Sclose(space_id)
    ierr = h5Pclose(dcpl_id)
    ierr = h5Dclose(dset_id)

  end subroutine create_unlimited_dataset_core

  !! Append a distributed ARRAY to an existing extendable dataset NAME at
  !! location LOC_ID. The dataset should be one created by a previous call
  !! to CREATE_UNLIMITED_DATASET. The dataset is extended to accommodate
  !! the concatenation of the arrays from all MPI ranks. The rank of ARRAY
  !! must be the same on each MPI rank, and the extent must be identical
  !! across all ranks for all but the last (extendable) dimension. These
  !! characteristics must match those of the dataset schema.
  !!
  !! If the optional argument ROOT is present, only the specified MPI rank
  !! writes its data; the dataset is extended solely by the extent of ARRAY
  !! on ROOT. This procedure must still be called collectively.
# 420 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine append_to_dataset_int32(ctx, loc_id, name, array, root)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(int32), intent(in) :: array(..)
    integer, intent(in), optional :: root

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id, dxpl_id
    integer :: ierr

    type_id = H5T_NATIVE_INT32

    call append_to_dataset_common(ctx, loc_id, name, type_id, shape(array, hsize_t), &
        dset_id, mem_space_id, data_space_id, root)

    call ctx%get_dxpl(dxpl_id) ! collective data transfer in MPI build

    select rank (array)
    rank (0)
      ierr = H5Dwrite(dset_id, type_id, [array], mem_space_id, data_space_id, dxpl_id)
    rank (1)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl_id)
    rank (2)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl_id)
    rank default
      INSIST(.false.)
    end select
    INSIST(ctx%global_all(ierr >= 0))

    ierr = H5Pclose(dxpl_id)
    ierr = h5sclose(mem_space_id)
    ierr = h5sclose(data_space_id)
    ierr = h5dclose(dset_id)

  end subroutine append_to_dataset_int32

# 420 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine append_to_dataset_int64(ctx, loc_id, name, array, root)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(int64), intent(in) :: array(..)
    integer, intent(in), optional :: root

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id, dxpl_id
    integer :: ierr

    type_id = H5T_NATIVE_INT64

    call append_to_dataset_common(ctx, loc_id, name, type_id, shape(array, hsize_t), &
        dset_id, mem_space_id, data_space_id, root)

    call ctx%get_dxpl(dxpl_id) ! collective data transfer in MPI build

    select rank (array)
    rank (0)
      ierr = H5Dwrite(dset_id, type_id, [array], mem_space_id, data_space_id, dxpl_id)
    rank (1)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl_id)
    rank (2)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl_id)
    rank default
      INSIST(.false.)
    end select
    INSIST(ctx%global_all(ierr >= 0))

    ierr = H5Pclose(dxpl_id)
    ierr = h5sclose(mem_space_id)
    ierr = h5sclose(data_space_id)
    ierr = h5dclose(dset_id)

  end subroutine append_to_dataset_int64

# 420 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine append_to_dataset_real32(ctx, loc_id, name, array, root)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(real32), intent(in) :: array(..)
    integer, intent(in), optional :: root

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id, dxpl_id
    integer :: ierr

    type_id = H5T_NATIVE_FLOAT

    call append_to_dataset_common(ctx, loc_id, name, type_id, shape(array, hsize_t), &
        dset_id, mem_space_id, data_space_id, root)

    call ctx%get_dxpl(dxpl_id) ! collective data transfer in MPI build

    select rank (array)
    rank (0)
      ierr = H5Dwrite(dset_id, type_id, [array], mem_space_id, data_space_id, dxpl_id)
    rank (1)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl_id)
    rank (2)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl_id)
    rank default
      INSIST(.false.)
    end select
    INSIST(ctx%global_all(ierr >= 0))

    ierr = H5Pclose(dxpl_id)
    ierr = h5sclose(mem_space_id)
    ierr = h5sclose(data_space_id)
    ierr = h5dclose(dset_id)

  end subroutine append_to_dataset_real32

# 420 "/app/src/vtkhdf_h5.F90.fypp"

  subroutine append_to_dataset_real64(ctx, loc_id, name, array, root)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(..)
    integer, intent(in), optional :: root

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id, dxpl_id
    integer :: ierr

    type_id = H5T_NATIVE_DOUBLE

    call append_to_dataset_common(ctx, loc_id, name, type_id, shape(array, hsize_t), &
        dset_id, mem_space_id, data_space_id, root)

    call ctx%get_dxpl(dxpl_id) ! collective data transfer in MPI build

    select rank (array)
    rank (0)
      ierr = H5Dwrite(dset_id, type_id, [array], mem_space_id, data_space_id, dxpl_id)
    rank (1)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl_id)
    rank (2)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl_id)
    rank default
      INSIST(.false.)
    end select
    INSIST(ctx%global_all(ierr >= 0))

    ierr = H5Pclose(dxpl_id)
    ierr = h5sclose(mem_space_id)
    ierr = h5sclose(data_space_id)
    ierr = h5dclose(dset_id)

  end subroutine append_to_dataset_real64

# 459 "/app/src/vtkhdf_h5.F90.fypp"
  !! This private procedure extracts the part of APPEND_TO_DATASET common to
  !! all types and ranks of the data array. It does not depend on the array
  !! directly, but only on its shape (DIMS) and its HDF5 datatype (TYPE_ID).
  !! It opens the dataset and creates the dataspaces for the data array and the
  !! part of the dataset where the array will be written, extending the dataset
  !! dataset as needed. If ROOT is present, only that MPI rank will write, but
  !! all ranks must still participate in the subsequent collective H5Dwrite
  !! call.

  subroutine append_to_dataset_common(ctx, loc_id, name, type_id, dims, &
      dset_id, mem_space_id, data_space_id, root)

    type(vtkhdf_ctx), intent(in) :: ctx
    integer(hid_t), intent(in) :: loc_id, type_id
    character(*), intent(in) :: name
    integer(hsize_t), intent(in) :: dims(:)
    integer(hid_t), intent(out) :: dset_id, data_space_id, mem_space_id
    integer, intent(in), optional :: root

    integer :: ndims, ierr, rank
    logical :: rank_writes
    integer(hid_t) :: dset_type_id
    integer(hsize_t) :: n
    integer(hsize_t), allocatable :: mem_dims(:), data_dims(:), start(:)

    !! Convert DIMS to its equivalent C form (MEM_DIMS)
    if (size(dims) == 0) then ! scalar appended data
      mem_dims = [1]  ! compatible with a rank-1 extendable dataset
    else
      mem_dims = dims(size(dims):1:-1) ! reverse order for the C shape
    end if

    !! Open the dataset and get its current shape (DATA_DIMS)
    dset_id = H5Dopen(loc_id, name)
    INSIST(ctx%global_all(dset_id >= 0))
    data_space_id = H5Dget_space(dset_id)
    INSIST(data_space_id > 0)
    ndims = H5Sget_simple_extent_ndims(data_space_id)
    INSIST(ndims > 0)
    allocate(data_dims(ndims))
    ndims = H5Sget_simple_extent_dims(data_space_id, data_dims)
    INSIST(ndims > 0)
    ierr = H5Sclose(data_space_id)

    ASSERT(all(mem_dims(2:) > 0))

    !! Ensure the data rank and fixed extents match the dataset
    INSIST(ctx%global_all(size(mem_dims) == size(data_dims)))
    INSIST(ctx%global_all(all(mem_dims(2:) == data_dims(2:))))

    !! Ensure the data type matches the dataset
    dset_type_id = H5Dget_type(dset_id)
    INSIST(dset_type_id > 0)
    INSIST(ctx%global_all(H5Tequal(dset_type_id, type_id) > 0))
    ierr = H5Tclose(dset_type_id)

    if (present(root)) then
      ASSERT((root >= 0) .and. (root < ctx%size))
      rank_writes = (root == ctx%rank)
    else
      rank_writes = .true.
    end if

    if (.not.rank_writes) mem_dims(1) = 0 ! ignore their data

    !! Starting index for the dataset hyperslab for this rank
    allocate(start, mold=mem_dims)
    start = 0
    call ctx%scan_sum(mem_dims(1), n)
    start(1) = data_dims(1) + n - mem_dims(1)

    !! Resize the dataset to accommodate the data to be appended
    call ctx%global_sum(mem_dims(1), n)
    data_dims(1) = data_dims(1) + n
    ierr = H5Dset_extent(dset_id, data_dims)
    INSIST(ctx%global_all(ierr >= 0))

    if (mem_dims(1) > 0) then
      mem_space_id = H5Screate(mem_dims)
    else ! rank contributes nothing
      mem_space_id = H5Screate(H5S_NULL)
    end if
    INSIST(mem_space_id > 0)

    !! Create the dataspace corresponding to appended elements in dataset
    data_space_id = H5Dget_space(dset_id)
    INSIST(data_space_id > 0)
    if (rank_writes .and. mem_dims(1) > 0) then
      ierr = H5Sselect_hyperslab(data_space_id, H5S_SELECT_SET, start, mem_dims)
    else ! rank contributes nothing
      ierr = H5Sselect_none(data_space_id)
    end if
    INSIST(ctx%global_all(ierr >= 0))

  end subroutine append_to_dataset_common

end module vtkhdf_h5
