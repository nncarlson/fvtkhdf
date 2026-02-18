program vtkhdf_int_comm_test
  use, intrinsic :: iso_fortran_env, only: real64
  use vtkhdf_ug_file_type
  use mpi_f08
  implicit none

  type(vtkhdf_ug_file) :: vizfile
  integer :: istat, rank, stat
  character(:), allocatable :: errmsg
  integer :: int_comm

  call MPI_Init(istat)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, istat)

  ! Extract integer handle from type(MPI_Comm)
  int_comm = MPI_COMM_WORLD%MPI_VAL

  ! Call create with the integer handle
  call vizfile%create('int_comm_test.vtkhdf', int_comm, stat, errmsg)

  if (stat /= 0) then
    if (rank == 0) print *, "Error: ", errmsg
    call MPI_Abort(MPI_COMM_WORLD, 1, istat)
  end if

  call vizfile%close()
  call MPI_Finalize(istat)

  if (rank == 0) print *, "SUCCESS: create with integer communicator worked."

end program vtkhdf_int_comm_test
