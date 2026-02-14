subroutine vtkhdf_assert(file, line)

  use vtkhdf_h5_c_binding
  use,intrinsic :: iso_fortran_env, only: error_unit

  character(*), intent(in) :: file
  integer, intent(in) :: line

  character(10) :: rank_str
  integer :: ierr

#ifdef USE_MPI
  block
    use mpi
    integer :: rank
    logical :: mpi_is_init
    call MPI_Initialized(mpi_is_init, ierr)
    if (mpi_is_init) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
      write(rank_str,'(a,i0,a)') '[R', rank, '] '
    else
      rank_str = ''
    end if
  end block
#endif

  ierr = H5Eprint2(H5E_DEFAULT, c_null_ptr)
  write(error_unit,'(/,4a,i0,a)') trim(rank_str), 'Panic! Assertion failed at ', file, ':', line

#ifdef USE_MPI
  call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
#else
  error stop
#endif

end subroutine vtkhdf_assert
