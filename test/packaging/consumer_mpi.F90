program consumer_mpi
#ifdef USE_MPI_F08
  use mpi_f08
#else
  use mpi
#endif
  use vtkhdf_ug_file_type
  implicit none
  type(vtkhdf_ug_file) :: file
  integer :: ierr
  character(:), allocatable :: errmsg
  call MPI_Init(ierr)
  call file%create('dummy.vtkhdf', MPI_COMM_WORLD, ierr, errmsg)
  call file%close
  call MPI_Finalize(ierr)
  print *, 'consumer_mpi: link OK'
end program
