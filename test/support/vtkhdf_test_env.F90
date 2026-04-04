module vtkhdf_test_env

#ifdef USE_MPI
  use mpi
#endif
  implicit none

  private
  public :: test_env_init, test_env_finalize

contains

  subroutine test_env_init(rank, nproc)
    integer, intent(out) :: rank, nproc
    integer :: ierr

#ifdef USE_MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
#else
    rank = 0
    nproc = 1
#endif
  end subroutine test_env_init

  subroutine test_env_finalize()
    integer :: ierr

#ifdef USE_MPI
    call MPI_Finalize(ierr)
#endif
  end subroutine test_env_finalize

end module vtkhdf_test_env
