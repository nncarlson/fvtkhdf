subroutine f90_assert(file, line)

  use vtkhdf_h5_c_binding
  
  character(*), intent(in) :: file
  integer, intent(in) :: line

  character(256) :: string
  integer :: ierr

  write(string,fmt='(a,i4.4)') "Assertion failed at " // file // ":", line

#ifdef USE_MPI
  block
    use mpi
    use,intrinsic :: iso_fortran_env, only: error_unit
    ierr = H5Eprint2(H5E_DEFAULT, c_null_ptr)
    write(error_unit,'(a)') trim(string)
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end block
#else
  ierr = H5Eprint2(H5E_DEFAULT, c_null_ptr)
  write(string,fmt='(a,i4.4)') "Assertion failed at " // file // ":", line
  error stop trim(string)
#endif

end subroutine
