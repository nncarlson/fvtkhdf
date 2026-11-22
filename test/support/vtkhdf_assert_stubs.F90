! Test-only modules; this target must not link real MPI or HDF5 libraries.
module mpi
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none
  integer, parameter :: MPI_COMM_WORLD = 123
contains
  subroutine MPI_Initialized(initialized, ierr)
    logical, intent(out) :: initialized
    integer, intent(out) :: ierr
    initialized = .true.
    ierr = 0
  end subroutine

  subroutine MPI_Comm_rank(comm, rank, ierr)
    integer, intent(in) :: comm
    integer, intent(out) :: rank, ierr
    rank = 0
    ierr = 0
  end subroutine

  subroutine MPI_Abort(comm, errorcode, ierr)
    integer, intent(in) :: comm, errorcode
    integer, intent(out) :: ierr
    write(error_unit,'(a,i0,a,i0)') 'ASSERT_MPI_ABORT_RETURNED comm=', comm, ' code=', errorcode
    flush(error_unit)
    ierr = 0
  end subroutine
end module

module vtkhdf_h5_c_binding
  use, intrinsic :: iso_c_binding, only: c_int, c_int64_t, c_ptr, c_null_ptr
  implicit none
  integer(c_int64_t), parameter :: H5E_DEFAULT = 0
contains
  function H5Eprint2(estack_id, stream) result(hdferr)
    integer(c_int64_t), value :: estack_id
    type(c_ptr), value :: stream
    integer(c_int) :: hdferr
    hdferr = 0
  end function
end module
