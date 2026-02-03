!!
!! VTKHDF_H5
!!
!! An application-specific layer over parallel HDF5 that provides single-call
!! procedures for writing attributes and datasets, creating and appending to
!! extendable datasets, and other utilities for private use by VTKHDF modules.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! January 2026
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! It is assumed that the HDF5 file on which these procedures operate has been
!! opened with a file access property list that enables MPI and collective
!! metadata operations and writes, as in the following:
!!

#include "f90_assert.fpp"

module vtkhdf_h5

  use,intrinsic :: iso_fortran_env
  use vtkhdf_h5_c_binding
  use mpi
  implicit none
  private

  public :: h5_write_attr, h5_create_unlimited_dataset, h5_write_dataset, h5_append_to_dataset
  public :: h5_mpi_comm

  interface h5_write_attr
    procedure write_attr_int8, write_attr_int32, write_attr_int64
    procedure write_attr_real32, write_attr_real64
    procedure write_attr_string
  end interface

  interface h5_write_dataset
    procedure write_dataset_int8, write_dataset_int32, write_dataset_int64
    procedure write_dataset_real32, write_dataset_real64
  end interface

  interface h5_create_unlimited_dataset
    procedure create_unlimited_dataset_int32, create_unlimited_dataset_int64
    procedure create_unlimited_dataset_real32, create_unlimited_dataset_real64
  end interface

  interface h5_append_to_dataset
    procedure append_to_dataset_int32, append_to_dataset_int64
    procedure append_to_dataset_real32, append_to_dataset_real64
  end interface

contains

  subroutine write_attr_int8(obj_id, name, value, stat, errmsg)

    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: name
    integer(int8), intent(in) :: value(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: type_id, attr_id
    integer :: istat

    type_id = H5T_NATIVE_UINT8

    call get_attr_id(obj_id, name, type_id, shape(value, hsize_t), attr_id, stat, errmsg)
    if (stat /= 0) return

    select rank (value)
    rank (0)
      stat = H5Awrite(attr_id, type_id, [value])
    rank (1)
      stat = H5Awrite(attr_id, type_id, value)
    rank (2)
      stat = H5Awrite(attr_id, type_id, value)
    rank default
      stat = 1
      errmsg = 'unsupported rank for attribute "' // name // '"'
      istat = H5Aclose(attr_id)
      return
    end select
    stat = min(0, stat) ! >= 0 from H5Awrite is successful

    istat = H5Aclose(attr_id)

    if (stat /= 0) then
      errmsg = 'error writing attribute "' // name // '"'
      return
    end if

  end subroutine write_attr_int8


  subroutine write_attr_int32(obj_id, name, value, stat, errmsg)

    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: name
    integer(int32), intent(in) :: value(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: type_id, attr_id
    integer :: istat

    type_id = H5T_NATIVE_INT32

    call get_attr_id(obj_id, name, type_id, shape(value, hsize_t), attr_id, stat, errmsg)
    if (stat /= 0) return

    select rank (value)
    rank (0)
      stat = H5Awrite(attr_id, type_id, [value])
    rank (1)
      stat = H5Awrite(attr_id, type_id, value)
    rank (2)
      stat = H5Awrite(attr_id, type_id, value)
    rank default
      stat = 1
      errmsg = 'unsupported rank for attribute "' // name // '"'
      istat = H5Aclose(attr_id)
      return
    end select
    stat = min(0, stat) ! >= 0 from H5Awrite is successful

    istat = H5Aclose(attr_id)

    if (stat /= 0) then
      errmsg = 'error writing attribute "' // name // '"'
      return
    end if

  end subroutine write_attr_int32


  subroutine write_attr_int64(obj_id, name, value, stat, errmsg)

    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: name
    integer(int64), intent(in) :: value(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: type_id, attr_id
    integer :: istat

    type_id = H5T_NATIVE_INT64

    call get_attr_id(obj_id, name, type_id, shape(value, hsize_t), attr_id, stat, errmsg)
    if (stat /= 0) return

    select rank (value)
    rank (0)
      stat = H5Awrite(attr_id, type_id, [value])
    rank (1)
      stat = H5Awrite(attr_id, type_id, value)
    rank (2)
      stat = H5Awrite(attr_id, type_id, value)
    rank default
      stat = 1
      errmsg = 'unsupported rank for attribute "' // name // '"'
      istat = H5Aclose(attr_id)
      return
    end select
    stat = min(0, stat) ! >= 0 from H5Awrite is successful

    istat = H5Aclose(attr_id)

    if (stat /= 0) then
      errmsg = 'error writing attribute "' // name // '"'
      return
    end if

  end subroutine write_attr_int64


  subroutine write_attr_real32(obj_id, name, value, stat, errmsg)

    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: name
    real(real32), intent(in) :: value(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: type_id, attr_id
    integer :: istat

    type_id = H5T_NATIVE_FLOAT

    call get_attr_id(obj_id, name, type_id, shape(value, hsize_t), attr_id, stat, errmsg)
    if (stat /= 0) return

    select rank (value)
    rank (0)
      stat = H5Awrite(attr_id, type_id, [value])
    rank (1)
      stat = H5Awrite(attr_id, type_id, value)
    rank (2)
      stat = H5Awrite(attr_id, type_id, value)
    rank default
      stat = 1
      errmsg = 'unsupported rank for attribute "' // name // '"'
      istat = H5Aclose(attr_id)
      return
    end select
    stat = min(0, stat) ! >= 0 from H5Awrite is successful

    istat = H5Aclose(attr_id)

    if (stat /= 0) then
      errmsg = 'error writing attribute "' // name // '"'
      return
    end if

  end subroutine write_attr_real32


  subroutine write_attr_real64(obj_id, name, value, stat, errmsg)

    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: name
    real(real64), intent(in) :: value(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: type_id, attr_id
    integer :: istat

    type_id = H5T_NATIVE_DOUBLE

    call get_attr_id(obj_id, name, type_id, shape(value, hsize_t), attr_id, stat, errmsg)
    if (stat /= 0) return

    select rank (value)
    rank (0)
      stat = H5Awrite(attr_id, type_id, [value])
    rank (1)
      stat = H5Awrite(attr_id, type_id, value)
    rank (2)
      stat = H5Awrite(attr_id, type_id, value)
    rank default
      stat = 1
      errmsg = 'unsupported rank for attribute "' // name // '"'
      istat = H5Aclose(attr_id)
      return
    end select
    stat = min(0, stat) ! >= 0 from H5Awrite is successful

    istat = H5Aclose(attr_id)

    if (stat /= 0) then
      errmsg = 'error writing attribute "' // name // '"'
      return
    end if

  end subroutine write_attr_real64


  subroutine write_attr_string(obj_id, name, value, stat, errmsg)

    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: name
#ifdef INTEL_BUG20240327
    character(*), intent(in) :: value
#else
    character(*), intent(in) :: value(..)
#endif
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: type_id, attr_id
    integer :: istat

    type_id = H5Tcopy(H5T_NATIVE_CHARACTER)
    INSIST(type_id > 0)
    istat = H5Tset_size(type_id, len(value, kind=c_size_t))
    istat = min(0, stat) ! >= 0 from H5Tset_size is successful
    INSIST(istat == 0)

    call get_attr_id(obj_id, name, type_id, shape(value, hsize_t), attr_id, stat, errmsg)
    if (stat /= 0) return

#ifdef INTEL_BUG20240327
    stat = H5Awrite(attr_id, type_id, [value])
#else
    select rank (value)
    rank (0)
      stat = H5Awrite(attr_id, type_id, [value])
    rank (1)
      stat = H5Awrite(attr_id, type_id, value)
    rank (2)
      stat = H5Awrite(attr_id, type_id, value)
    rank default
      stat = 1
      errmsg = 'unsupported rank for attribute "' // name // '"'
      istat = H5Aclose(attr_id)
      istat = H5Tclose(type_id)
      return
    end select
#endif
    stat = min(0, stat) ! >= 0 from H5Awrite is successful

    istat = H5Aclose(attr_id)
    istat = H5Tclose(type_id)

    if (stat /= 0) then
      errmsg = 'error writing attribute "' // name // '"'
      return
    end if

  end subroutine write_attr_string


  subroutine get_attr_id(obj_id, name, type_id, dims, attr_id, stat, errmsg)

    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: name
    integer(hid_t), intent(in) :: type_id
    integer(hsize_t), intent(in) :: dims(:) ! the Fortran shape
    integer(hid_t), intent(out) :: attr_id
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: space_id
    integer :: istat ! ignored status result

    stat = 0

    if (H5Aexists(obj_id, name)) then  ! open the attribute

      attr_id = H5Aopen(obj_id, name)
      if (attr_id < 0) then
        stat = 1
        errmsg = 'unable to open attribute "' // name // '"'
        return
      end if

      !NB: we assume type_id and dims match the attribute

    else ! create the attribute

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
      INSIST(space_id > 0)

      attr_id = H5Acreate(obj_id, name, type_id, space_id)
      if (attr_id < 0) then
        stat = 1
        errmsg = 'unable to create attribute "' // name // '"'
        return
      end if
      istat = H5Sclose(space_id)

    end if

  end subroutine get_attr_id

  !! Write a distributed array (DATA) to a new dataset NAME at location LOC_ID
  !! in the HDF5 file tree structure. The written dataset is the concatenation
  !! of the arrays from all the processes. The rank of DATA must be the same
  !! on each process, and in the multi-rank case, the extent must be the same
  !! in all but the last dimension, which is the distributed dimension. If the
  !! optional argument ROOT is present, the behavior is different. In this case
  !! ROOT specifies the MPI rank of the only process that writes its DATA, and
  !! the written dataset consists of that data alone. The content of DATA on
  !! all other processes is ignored, however the previous constraints on the
  !! rank and extents continues to apply. Moreover this must still be called
  !! collectively as creation of the dataset itself is a collective operation.

  subroutine write_dataset_int8(loc_id, name, array, comm, stat, errmsg, root)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(int8), intent(in) :: array(..)
    integer, intent(in) :: comm
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: root

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id, dxpl
    integer :: ierr

    type_id = H5T_NATIVE_UINT8

    call write_dataset_aux(loc_id, name, type_id, shape(array,hsize_t), comm, &
        dset_id, mem_space_id, data_space_id, dxpl, stat, errmsg, root)
    if (stat /= 0) return

    select rank (array)
    rank (0)
      ierr = H5Dwrite(dset_id, type_id, [array], mem_space_id, data_space_id, dxpl)
    rank (1)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl)
    rank (2)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl)
    rank default
      INSIST(.false.)
    end select

    if (collective_failed(ierr < 0, comm)) then
      stat = 1
      errmsg = 'error writing to dataset "' // name // '"'
    end if

    !! Cleanup
    ierr = H5Pclose(dxpl)
    ierr = H5Sclose(mem_space_id)
    ierr = H5Sclose(data_space_id)
    ierr = H5Dclose(dset_id)

  end subroutine write_dataset_int8


  subroutine write_dataset_int32(loc_id, name, array, comm, stat, errmsg, root)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(int32), intent(in) :: array(..)
    integer, intent(in) :: comm
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: root

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id, dxpl
    integer :: ierr

    type_id = H5T_NATIVE_INT32

    call write_dataset_aux(loc_id, name, type_id, shape(array,hsize_t), comm, &
        dset_id, mem_space_id, data_space_id, dxpl, stat, errmsg, root)
    if (stat /= 0) return

    select rank (array)
    rank (0)
      ierr = H5Dwrite(dset_id, type_id, [array], mem_space_id, data_space_id, dxpl)
    rank (1)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl)
    rank (2)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl)
    rank default
      INSIST(.false.)
    end select

    if (collective_failed(ierr < 0, comm)) then
      stat = 1
      errmsg = 'error writing to dataset "' // name // '"'
    end if

    !! Cleanup
    ierr = H5Pclose(dxpl)
    ierr = H5Sclose(mem_space_id)
    ierr = H5Sclose(data_space_id)
    ierr = H5Dclose(dset_id)

  end subroutine write_dataset_int32


  subroutine write_dataset_int64(loc_id, name, array, comm, stat, errmsg, root)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(int64), intent(in) :: array(..)
    integer, intent(in) :: comm
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: root

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id, dxpl
    integer :: ierr

    type_id = H5T_NATIVE_INT64

    call write_dataset_aux(loc_id, name, type_id, shape(array,hsize_t), comm, &
        dset_id, mem_space_id, data_space_id, dxpl, stat, errmsg, root)
    if (stat /= 0) return

    select rank (array)
    rank (0)
      ierr = H5Dwrite(dset_id, type_id, [array], mem_space_id, data_space_id, dxpl)
    rank (1)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl)
    rank (2)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl)
    rank default
      INSIST(.false.)
    end select

    if (collective_failed(ierr < 0, comm)) then
      stat = 1
      errmsg = 'error writing to dataset "' // name // '"'
    end if

    !! Cleanup
    ierr = H5Pclose(dxpl)
    ierr = H5Sclose(mem_space_id)
    ierr = H5Sclose(data_space_id)
    ierr = H5Dclose(dset_id)

  end subroutine write_dataset_int64


  subroutine write_dataset_real32(loc_id, name, array, comm, stat, errmsg, root)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(real32), intent(in) :: array(..)
    integer, intent(in) :: comm
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: root

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id, dxpl
    integer :: ierr

    type_id = H5T_NATIVE_FLOAT

    call write_dataset_aux(loc_id, name, type_id, shape(array,hsize_t), comm, &
        dset_id, mem_space_id, data_space_id, dxpl, stat, errmsg, root)
    if (stat /= 0) return

    select rank (array)
    rank (0)
      ierr = H5Dwrite(dset_id, type_id, [array], mem_space_id, data_space_id, dxpl)
    rank (1)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl)
    rank (2)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl)
    rank default
      INSIST(.false.)
    end select

    if (collective_failed(ierr < 0, comm)) then
      stat = 1
      errmsg = 'error writing to dataset "' // name // '"'
    end if

    !! Cleanup
    ierr = H5Pclose(dxpl)
    ierr = H5Sclose(mem_space_id)
    ierr = H5Sclose(data_space_id)
    ierr = H5Dclose(dset_id)

  end subroutine write_dataset_real32


  subroutine write_dataset_real64(loc_id, name, array, comm, stat, errmsg, root)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(real64), intent(in) :: array(..)
    integer, intent(in) :: comm
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: root

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id, dxpl
    integer :: ierr

    type_id = H5T_NATIVE_DOUBLE

    call write_dataset_aux(loc_id, name, type_id, shape(array,hsize_t), comm, &
        dset_id, mem_space_id, data_space_id, dxpl, stat, errmsg, root)
    if (stat /= 0) return

    select rank (array)
    rank (0)
      ierr = H5Dwrite(dset_id, type_id, [array], mem_space_id, data_space_id, dxpl)
    rank (1)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl)
    rank (2)
      ierr = H5Dwrite(dset_id, type_id, array, mem_space_id, data_space_id, dxpl)
    rank default
      INSIST(.false.)
    end select

    if (collective_failed(ierr < 0, comm)) then
      stat = 1
      errmsg = 'error writing to dataset "' // name // '"'
    end if

    !! Cleanup
    ierr = H5Pclose(dxpl)
    ierr = H5Sclose(mem_space_id)
    ierr = H5Sclose(data_space_id)
    ierr = H5Dclose(dset_id)

  end subroutine write_dataset_real64

  !! This auxiliary procedure extracts the part of WRITE_DATASET common to
  !! all types and ranks of the data array. It does not depend on the array
  !! directly, but only on its shape (DIMS) and its HDF5 type (TYPE_ID). It
  !! creates the dataset and the dataspaces for the data array and the part
  !! of the dataset where the array will be written.
  !!
  !! If ROOT is present, only that MPI rank will write, and the dataspace IDs
  !! are returned with value 0 (invalid) on all other processes. This procedure
  !! must still be called collectively because it creates the dataset, which
  !! is a collective operation.

  subroutine write_dataset_aux(loc_id, name, type_id, dims, comm, &
      dset_id, mem_space_id, data_space_id, dxpl, stat, errmsg, root)

    integer(hid_t), intent(in) :: loc_id, type_id
    character(*), intent(in) :: name
    integer(hsize_t), intent(in) :: dims(:)
    integer, intent(in) :: comm
    integer(hid_t), intent(out) :: dset_id, mem_space_id, data_space_id, dxpl
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: root

    integer :: ierr, rank
    logical :: local_failed
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

    !! Some parallel consistency checks on the input arguments. It is most
    !! likely that each rank calls WRITE_DATASET from the same line of code,
    !! and consequently it is generally safe to assume that the data array
    !! type, rank, and shape will be consistent across all ranks, so we
    !! limit the checks to debug mode.
    ASSERT(c_dims_is_valid(comm, mem_dims))
    ASSERT(type_id_is_valid(comm, type_id))
    ASSERT(root_is_valid(comm, root))

    if (present(root)) then
      call MPI_Comm_rank(comm, rank, ierr)
      INSIST(ierr == MPI_SUCCESS)
      rank_writes = (rank == root)
    else
      rank_writes = .true.
    end if

    if (.not.rank_writes) mem_dims(1) = 0 ! ignore their data

    !! Create the dataset (always collective!)
    call MPI_Allreduce(mem_dims(1), n, 1, MPI_INT64_T, MPI_SUM, comm, ierr)
    INSIST(ierr == MPI_SUCCESS)
    if (n == 0) then
      stat = 1
      errmsg = 'writing empty dataset "' // name // '" is not supported'
      return
    end if
    data_dims = mem_dims
    data_dims(1) = n
    data_space_id = H5Screate(data_dims)
    INSIST(data_space_id > 0)
    dset_id = H5Dcreate(loc_id, name, type_id, data_space_id)
    if (collective_failed(dset_id < 0, comm)) then
      stat = 1
      errmsg = 'error creating dataset "' // name // '"'
      ierr = H5Sclose(data_space_id)
      return
    end if

    !! Starting index for the dataset hyperslab for this rank
    allocate(start, mold=mem_dims)
    start = 0
    call MPI_Scan(mem_dims(1), n, 1, MPI_INT64_T, MPI_SUM, comm, ierr)
    INSIST(ierr == MPI_SUCCESS)
    start(1) = n - mem_dims(1)

    if (mem_dims(1) > 0) then
      mem_space_id = H5Screate(mem_dims)
    else
      mem_space_id = H5Screate(H5S_NULL)
    end if
    INSIST(mem_space_id > 0)

    !! Create the dataspace for the part of dataset to be written by this rank
    if (rank_writes .and. mem_dims(1) > 0) then
      count = mem_dims
      ierr = H5Sselect_hyperslab(data_space_id, H5S_SELECT_SET, start, count)
    else ! 0-size
      ierr = H5Sselect_none(data_space_id)
    end if

    if (collective_failed(ierr < 0,comm)) then
      stat = 1
      errmsg = 'error creating hyperslab of dataset "' // name // '"'
      ierr = H5Sclose(data_space_id)
      ierr = H5Sclose(mem_space_id)
      ierr = H5Dclose(dset_id)
      return
    end if

    dxpl = H5Pcreate(H5P_DATASET_XFER)
    INSIST(dxpl > 0)
    ierr = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE)
    INSIST(ierr == 0)

    stat = 0

  end subroutine write_dataset_aux

  logical function collective_failed(local_failed, comm) result(failed)
    logical, intent(in) :: local_failed
    integer, intent(in) :: comm
    integer :: local_stat, global_stat, ierr
    local_stat = merge(1, 0, local_failed)
    call MPI_Allreduce(local_stat, global_stat, 1, MPI_INT, MPI_MAX, comm, ierr)
    INSIST(ierr == MPI_SUCCESS)
    failed = (global_stat /= 0)
  end function

  !! Check the parallel consistency of TYPE_ID (DEBUGGING)
  logical function type_id_is_valid(comm, type_id)

    integer, intent(in) :: comm
    integer(hid_t), intent(in) :: type_id

    integer :: istat
    integer(int64) :: p_local(2), p_global(2)

    p_local = [type_id, -type_id] ! cast to a known MPI type
    call MPI_Allreduce(p_local, p_global, 2, MPI_INT64_T, MPI_MIN, comm, istat)
    INSIST(istat == MPI_SUCCESS)
    type_id_is_valid = (p_global(1) == -p_global(2))

  end function type_id_is_valid

  !! Check the parallel consistency and correctness of C_DIMS (DEBUGGING)
  logical function c_dims_is_valid(comm, c_dims)

    integer, intent(in) :: comm
    integer(hsize_t), intent(in) :: c_dims(:)

    integer :: istat
    integer :: p_local(2), p_global(2)
    integer(int64), allocatable :: q_local(:), q_global(:)

    !! Check that DIMS has the same size on all ranks
    p_local = [size(c_dims), -size(c_dims)]
    call MPI_Allreduce(p_local, p_global, 2, MPI_INT, MPI_MIN, comm, istat)
    INSIST(istat == MPI_SUCCESS)
    c_dims_is_valid = (p_global(1) == -p_global(2))

    if (c_dims_is_valid .and. size(c_dims) > 1) then
      !! Check that C_DIMS has the same extent in all dimensions but the first.
      q_local = c_dims(2:) ! cast to a known MPI type
      q_local = [q_local, -q_local(size(q_local):1:-1)]
      allocate(q_global, mold=q_local)
      call MPI_Allreduce(q_local, q_global, size(q_local), MPI_INT64_T, MPI_MIN, comm, istat)
      INSIST(istat == MPI_SUCCESS)
      c_dims_is_valid = all(q_global + q_global(size(q_global):1:-1) == 0)
      if (c_dims_is_valid) c_dims_is_valid = all(c_dims(2:) > 0)
    end if

  end function c_dims_is_valid

  !! Check the parallel consistency of CHUNK_SIZE (DEBUGGING)
  logical function chunk_is_valid(comm, chunk_size)

    integer, intent(in) :: comm
    integer, intent(in) :: chunk_size

    integer :: istat
    integer :: p_local(2), p_global(2)

    p_local = [chunk_size, -chunk_size]
    call MPI_Allreduce(p_local, p_global, 2, MPI_INT, MPI_MIN, comm, istat)
    INSIST(istat == MPI_SUCCESS)
    chunk_is_valid = (p_global(1) == -p_global(2))
    if (chunk_is_valid) chunk_is_valid = (chunk_size > 0)

  end function chunk_is_valid

  !! Check the parallel consistency of ROOT (DEBUGGING)
  logical function root_is_valid(comm, root)

    integer, intent(in) :: comm
    integer, intent(in), optional :: root

    integer :: n, p_local(2), p_global(2), istat, nproc

    !! Check that all ranks specified root or that none did.
    n = merge(1, 0, present(root))
    p_local = [n, -n]
    call MPI_Allreduce(p_local, p_global, 2, MPI_INT, MPI_MIN, comm, istat)
    INSIST(istat == MPI_SUCCESS)
    root_is_valid = (p_global(1) == -p_global(2))

    if (root_is_valid .and. present(root)) then
      !! Check that all ranks specified the same valid value.
      p_local = [root, -root]
      call MPI_Allreduce(p_local, p_global, 2, MPI_INT, MPI_MIN, comm, istat)
      INSIST(istat == MPI_SUCCESS)
      root_is_valid = (p_global(1) == -p_global(2)) ! same?
      if (root_is_valid) then! Check that the value of root is valid
        call MPI_Comm_size(comm, nproc, istat)
        INSIST(istat == MPI_SUCCESS)
        root_is_valid = (root >= 0) .and. (root < nproc)
      end if
    end if

  end function root_is_valid

  !! Create an extendible dataset NAME at location LOC_ID in the HDF5 file tree
  !! structure. The array MOLD specifies the type and rank of the dataset, and
  !! must be the same for all processes. In the multi-rank case, the extent in
  !! all dimensions but the last must be the same on all processes and sets the
  !! extent of the dataset in those dimensions. The last dimension is the
  !! extendible dimension (H5S_UNLIMITED) of the dataset, and the extent of the
  !! corresponding dimension of MOLD is ignored. Note that the DATA array passed
  !! to subsequent calls to APPEND_TO_DATASET would be suitable to use as MOLD,
  !! for example. The specified value for CHUNK_SIZE must be the same on all
  !! processes, and specifies the chunk size in terms of the number of elements
  !! in the extendible dimension. No data is written to the dataset, and its
  !! extent in the extendible dimension starts at 0.

  subroutine create_unlimited_dataset_int32(loc_id, name, mold, chunk_size, stat, errmsg)
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(int32), intent(in) :: mold(..)
    integer, intent(in) :: chunk_size
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call create_unlimited_dataset_aux(loc_id, name, H5T_NATIVE_INT32, shape(mold,hsize_t), chunk_size, stat, errmsg)
  end subroutine

  subroutine create_unlimited_dataset_int64(loc_id, name, mold, chunk_size, stat, errmsg)
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(int64), intent(in) :: mold(..)
    integer, intent(in) :: chunk_size
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call create_unlimited_dataset_aux(loc_id, name, H5T_NATIVE_INT64, shape(mold,hsize_t), chunk_size, stat, errmsg)
  end subroutine

  subroutine create_unlimited_dataset_real32(loc_id, name, mold, chunk_size, stat, errmsg)
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(real32), intent(in) :: mold(..)
    integer, intent(in) :: chunk_size
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call create_unlimited_dataset_aux(loc_id, name, H5T_NATIVE_FLOAT, shape(mold,hsize_t), chunk_size, stat, errmsg)
  end subroutine

  subroutine create_unlimited_dataset_real64(loc_id, name, mold, chunk_size, stat, errmsg)
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(real64), intent(in) :: mold(..)
    integer, intent(in) :: chunk_size
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call create_unlimited_dataset_aux(loc_id, name, H5T_NATIVE_DOUBLE, shape(mold,hsize_t), chunk_size, stat, errmsg)
  end subroutine


  subroutine create_unlimited_dataset_aux(loc_id, name, type_id, dims, chunk_size, comm, stat, errmsg)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(hid_t), intent(in) :: type_id
    integer(hsize_t), intent(in) :: dims(:)
    integer, intent(in) :: chunk_size
    integer, intent(in) :: comm
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: ierr
    integer(hid_t) :: space_id, dcpl_id, dset_id
    integer(hsize_t), allocatable :: c_dims(:), maxdims(:), chunk_dims(:)

    c_dims = dims
    if (size(dims) > 1) c_dims = dims(size(dims):1:-1)

    !! Some parallel consistency checks on the input arguments. It is most
    !! likely that each rank calls CREATE_UNLIMITED_DATASET from the same
    !! line of code, and consequently it is generally safe to assume that
    !! the data array type, rank, and shape (DIMS) will be consistent across
    !! all ranks, so we limit the checks to debug mode.
    ASSERT(c_dims_is_valid(comm, c_dims))
    ASSERT(type_id_is_valid(comm, type_id))
    ASSERT(chunk_is_valid(comm, chunk_size))

    if (size(c_dims) == 0) then
      stat = 1
      errmsg = 'invalid rank for unlimited dataset'
      return
    end if

    maxdims = c_dims
    maxdims(1) = H5S_UNLIMITED
    chunk_dims = c_dims
    chunk_dims(1) = chunk_size
    c_dims(1) = 0

    space_id = H5Screate(c_dims, maxdims)
    INSIST(space_id > 0)

    dcpl_id = H5Pcreate(H5P_DATASET_CREATE)
    INSIST(dcpl_id > 0)

    ierr = H5Pset_chunk(dcpl_id, size(chunk_dims), chunk_dims)
    INSIST(ierr >= 0)

    dset_id = H5Dcreate(loc_id, name, type_id, space_id, dcpl_id=dcpl_id)
    if (collective_failed(dset_id < 0, comm)) then
      stat = 1
      errmsg = 'error creating unlimited dataset "' // name // '"'
      ierr = H5Pclose(dcpl_id)
      ierr = H5Sclose(space_id)
      return
    end if

    ierr = h5Sclose(space_id)
    ierr = h5Pclose(dcpl_id)
    ierr = h5Dclose(dset_id)

    stat = 0

  end subroutine create_unlimited_dataset_aux

  !! Append a distributed array (DATA) to an existing extendible dataset NAME
  !! at location LOC_ID in the HDF5 file tree structure. The dataset should be
  !! one created by a previous call to CREATE_UNLIMITED_DATASET. The dataset
  !! is extended to accommodate the concatenation of the arrays from all the
  !! processes. The rank of DATA must be the same on each process, and in the
  !! multi-rank case, the extent must be the same in all but the last dimension,
  !! which is the distributed dimension. Moreover these characteristics must
  !! match those of the dataset. If the optional argument ROOT is present, the
  !! behavior is different. In this case ROOT specifies the MPI rank of the
  !! only process that writes its DATA, and the written dataset consists of
  !! that data alone. The content of DATA on all other processes is ignored,
  !! however the previous constraints on the rank and extents continues to
  !! apply. Moreover this must still be called collectively as extending the
  !! size of the dataset is itself a collective operation.

  subroutine append_to_dataset_int32(loc_id, name, data, stat, errmsg, root)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(int32), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: root

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id, dxpl
    integer :: istat

    type_id = H5T_NATIVE_INT32

    call append_to_dataset_aux(loc_id, name, type_id, shape(data, hsize_t), &
        dset_id, mem_space_id, data_space_id, stat, errmsg, root)
    if (stat /= 0) return

    if (mem_space_id == 0) then
      istat = H5Dclose(dset_id)
      return
    end if

    !! Set the transfer mode
    dxpl = H5Pcreate(H5P_DATASET_XFER)
    if (present(root)) then ! only root writes
      stat = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT)
    else  ! all write
      stat = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE)
    end if
    INSIST(stat == 0)

    select rank (data)
    rank (0)
      stat = H5Dwrite(dset_id, type_id, [data], mem_space_id, data_space_id, dxpl)
    rank (1)
      stat = H5Dwrite(dset_id, type_id, data, mem_space_id, data_space_id, dxpl)
    rank (2)
      stat = H5Dwrite(dset_id, type_id, data, mem_space_id, data_space_id, dxpl)
    rank default
      stat = 1
      errmsg = 'unsupported rank'
      return
    end select
    stat = min(0, stat) ! >= 0 from H5Dwrite is successful

    istat = H5Pclose(dxpl)
    istat = h5sclose(mem_space_id)
    istat = h5sclose(data_space_id)
    istat = h5dclose(dset_id)

    if (stat /= 0) then
      errmsg = 'error writing to dataset "' // name // '"'
      return
    end if

  end subroutine append_to_dataset_int32

  subroutine append_to_dataset_int64(loc_id, name, data, stat, errmsg, root)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(int64), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: root

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id, dxpl
    integer :: istat

    type_id = H5T_NATIVE_INT64

    call append_to_dataset_aux(loc_id, name, type_id, shape(data, hsize_t), &
        dset_id, mem_space_id, data_space_id, stat, errmsg, root)
    if (stat /= 0) return

    if (mem_space_id == 0) then
      istat = H5Dclose(dset_id)
      return
    end if

    !! Set the transfer mode
    dxpl = H5Pcreate(H5P_DATASET_XFER)
    if (present(root)) then ! only root writes
      stat = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT)
    else  ! all write
      stat = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE)
    end if
    INSIST(stat == 0)

    select rank (data)
    rank (0)
      stat = H5Dwrite(dset_id, type_id, [data], mem_space_id, data_space_id, dxpl)
    rank (1)
      stat = H5Dwrite(dset_id, type_id, data, mem_space_id, data_space_id, dxpl)
    rank (2)
      stat = H5Dwrite(dset_id, type_id, data, mem_space_id, data_space_id, dxpl)
    rank default
      stat = 1
      errmsg = 'unsupported rank'
      return
    end select
    stat = min(0, stat) ! >= 0 from H5Dwrite is successful

    istat = H5Pclose(dxpl)
    istat = h5sclose(mem_space_id)
    istat = h5sclose(data_space_id)
    istat = h5dclose(dset_id)

    if (stat /= 0) then
      errmsg = 'error writing to dataset "' // name // '"'
      return
    end if

  end subroutine append_to_dataset_int64


  subroutine append_to_dataset_real32(loc_id, name, data, stat, errmsg, root)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(real32), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: root

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id, dxpl
    integer :: istat

    type_id = H5T_NATIVE_FLOAT

    call append_to_dataset_aux(loc_id, name, type_id, shape(data, hsize_t), &
        dset_id, mem_space_id, data_space_id, stat, errmsg, root)
    if (stat /= 0) return

    if (mem_space_id == 0) then ! process doesn't write
      istat = H5Dclose(dset_id)
      return
    end if

    !! Set the transfer mode
    dxpl = H5Pcreate(H5P_DATASET_XFER)
    if (present(root)) then ! only root writes
      stat = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT)
    else  ! all write
      stat = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE)
    end if
    INSIST(stat == 0)

    select rank (data)
    rank (0)
      stat = H5Dwrite(dset_id, type_id, [data], mem_space_id, data_space_id, dxpl)
    rank (1)
      stat = H5Dwrite(dset_id, type_id, data, mem_space_id, data_space_id, dxpl)
    rank (2)
      stat = H5Dwrite(dset_id, type_id, data, mem_space_id, data_space_id, dxpl)
    rank default
      stat = 1
      errmsg = 'unsupported rank'
      return
    end select
    stat = min(0, stat) ! >= 0 from H5Dwrite is successful

    istat = H5Pclose(dxpl)
    istat = h5sclose(mem_space_id)
    istat = h5sclose(data_space_id)
    istat = h5dclose(dset_id)

    if (stat /= 0) then
      errmsg = 'error writing to dataset "' // name // '"'
      return
    end if

  end subroutine append_to_dataset_real32


  subroutine append_to_dataset_real64(loc_id, name, data, stat, errmsg, root)

    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(real64), intent(in) :: data(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: root

    integer(hid_t) :: type_id, dset_id, mem_space_id, data_space_id, dxpl
    integer :: istat

    type_id = H5T_NATIVE_DOUBLE

    call append_to_dataset_aux(loc_id, name, type_id, shape(data, hsize_t), &
        dset_id, mem_space_id, data_space_id, stat, errmsg, root)
    if (stat /= 0) return

    if (mem_space_id == 0) then ! process doesn't write
      istat = H5Dclose(dset_id)
      return
    end if

    !! Set the transfer mode
    dxpl = H5Pcreate(H5P_DATASET_XFER)
    if (present(root)) then ! only root writes
      stat = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT)
    else  ! all write
      stat = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE)
    end if
    INSIST(stat == 0)

    select rank (data)
    rank (0)
      stat = H5Dwrite(dset_id, type_id, [data], mem_space_id, data_space_id, dxpl)
    rank (1)
      stat = H5Dwrite(dset_id, type_id, data, mem_space_id, data_space_id, dxpl)
    rank (2)
      stat = H5Dwrite(dset_id, type_id, data, mem_space_id, data_space_id, dxpl)
    rank default
      stat = 1
      errmsg = 'unsupported rank'
      return
    end select
    stat = min(0, stat) ! >= 0 from H5Dwrite is successful

    istat = H5Pclose(dxpl)
    istat = h5sclose(mem_space_id)
    istat = h5sclose(data_space_id)
    istat = h5dclose(dset_id)

    if (stat /= 0) then
      errmsg = 'error writing to dataset "' // name // '"'
      return
    end if

  end subroutine append_to_dataset_real64

  !! This auxiliary procedure extracts the part of APPEND_TO_DATASET common to
  !! all types and ranks of the data array. It does not depend on the array
  !! directly, but only on its shape (DIMS) and its HDF5 type (TYPE_ID). It
  !! opens the dataset and creates the dataspaces for the data array and the
  !! part of the dataset where the array will be written, extending the dataset
  !! dataset as needed.
  !!
  !! It is most likely that each process calls APPEND_TO_DATASET from the same
  !! line of code, so that it is reasonably safe to assume that the data array
  !! type, rank, and extents in all but the last dimension are the same across
  !! processes. But since the dataset is created separately, we verify that
  !! those characteristics are consistent with the data array, and we do this
  !! independently on each process without communication.
  !!
  !! If ROOT is present, only that MPI rank will write, and the dataspace IDs
  !! are returned with value 0 (invalid) on all other processes. This procedure
  !! must still be called collectively because it extends the size of the
  !! dataset, which is a collective operation.

  subroutine append_to_dataset_aux(loc_id, name, type_id, dims, &
      dset_id, mem_space_id, data_space_id, stat, errmsg, root)

    integer(hid_t), intent(in) :: loc_id, type_id
    character(*), intent(in) :: name
    integer(hsize_t), intent(in) :: dims(:)
    integer(hid_t), intent(out) :: dset_id, data_space_id, mem_space_id
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: root

    integer :: ndims, istat, comm, rank
    integer(hid_t) :: dset_type_id
    integer(hsize_t) :: n
    integer(hsize_t), allocatable :: mem_dims(:), data_dims(:), start(:)

    !! Convert DIMS to its equivalent C form (MEM_DIMS)
    if (size(dims) == 0) then ! scalar appended data
      mem_dims = [1]  ! compatible with a rank-1 extendable dataset
    else
      mem_dims = dims
      if (size(dims) > 1) mem_dims = dims(size(dims):1:-1) ! reverse order for the C shape
    end if

    !! Open the dataset and get its current shape (DATA_DIMS)
    dset_id = H5Dopen(loc_id, name)
    if (dset_id < 0) then
      stat = 1
      errmsg = 'unable to open dataset "' // name // '"'
      return
    end if
    data_space_id = H5Dget_space(dset_id)
    INSIST(data_space_id > 0)
    ndims = H5Sget_simple_extent_ndims(data_space_id)
    INSIST(ndims > 0)
    allocate(data_dims(ndims))
    ndims = H5Sget_simple_extent_dims(data_space_id, data_dims)
    istat = H5Sclose(data_space_id)

    !! Verify data rank compatibility
    if (size(mem_dims) /= size(data_dims)) then
      stat = 1
      errmsg = 'rank of appended data is incompatible with dataset "' // name // '"'
      istat = H5Dclose(dset_id)
      return
    end if

    !! Verify data shape compatibility
    if (any(mem_dims(2:) /= data_dims(2:))) then
      stat = 1
      errmsg = 'shape of appended data is incompatible with dataset "' // name // '"'
      return
    end if

!TODO: H5Tequal doesn't recognize dset_type_id as a type ID for some reason!
!    !! Verify data type compatibility
!    dset_type_id = H5Dget_type(dset_id)
!    INSIST(dset_type_id > 0)
!    if (H5Tequal(dset_type_id, type_id) <= 0) then
!      errmsg = 'type of appended data is incompatible with dataset "' // name // '"'
!      istat = H5Tclose(dset_type_id)
!      istat = H5Dclose(dset_id)
!      return
!    end if
!    istat = H5Tclose(dset_type_id)

    comm = h5_mpi_comm(dset_id)
    call MPI_Comm_rank(comm, rank, istat)

    if (present(root)) then
      if (rank /= root) mem_dims(1) = 0
    end if

    !! Starting index for the dataset hyperslab
    allocate(start, mold=mem_dims)
    start = 0
    call MPI_Scan(mem_dims(1), n, 1, MPI_INT64_T, MPI_SUM, comm, istat)
    start(1) = data_dims(1) + n - mem_dims(1)

    !! Resize the dataset to accommodate the data to be appended
    call MPI_Allreduce(mem_dims(1), n, 1, MPI_INT64_T, MPI_SUM, comm, istat)
    call MPI_Comm_free(comm, istat)
    data_dims(1) = data_dims(1) + n
    stat = H5Dset_extent(dset_id, data_dims)
    stat = min(0, stat) ! >= 0 from H5Dset_extent is successful
    if (stat /= 0) then
      errmsg = 'error resizing dataset "' // name // '"'
      istat = H5Dclose(dset_id)
      return
    end if

    if (present(root)) then
      if (rank /= root) then ! process will not participate in the actual write
        mem_space_id = 0
        data_space_id = 0
        return
      end if
    end if

    !! Create the dataspace for the data array to be written
    mem_space_id = H5Screate(mem_dims)
    INSIST(mem_space_id > 0)

    !! Create the dataspace corresponding to appended elements in dataset
    data_space_id = H5Dget_space(dset_id)
    stat = H5Sselect_hyperslab(data_space_id, H5S_SELECT_SET, start, mem_dims)
    stat = min(0, stat) ! >= 0 from H5Sselect_hyperslab is successful
    if (stat /= 0) then
      stat = 1
      errmsg = 'error creating hyperslab of dataset "' // name // '"'
      istat = H5Sclose(data_space_id)
      istat = H5Sclose(mem_space_id)
      istat = H5Dclose(dset_id)
      return
    end if

  end subroutine append_to_dataset_aux

  !! Get a copy of the MPI communicator
  !! NB: communicator needs to be freed when finished with it
  integer function h5_mpi_comm(loc_id) result(comm)
    integer(hid_t), intent(in) :: loc_id
    integer(hid_t) :: file_id, fapl
    integer :: stat, istat
    file_id = H5Iget_file_id(loc_id)
    INSIST(file_id > 0)
    fapl = H5Fget_access_plist(file_id)
    INSIST(fapl > 0)
    stat = H5Pget_fapl_mpio(fapl, comm)
    INSIST(stat >= 0)
    istat = H5Pclose(fapl)
    istat = H5Fclose(file_id)
  end function

end module vtkhdf_h5
