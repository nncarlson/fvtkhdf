program vtkhdf_mb_fail_late_block_test

  use, intrinsic :: iso_fortran_env, only: r8 => real64
  use vtkhdf_test_env
  use vtkhdf_mb_file_type
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  character(:), allocatable :: errmsg
  integer :: rank, nproc, stat

  type(vtkhdf_mb_file) :: vizfile
  type(vtkhdf_block_handle) :: hblk

  call test_env_init(rank, nproc)

#ifdef USE_MPI
  call vizfile%create('mb_fail_late_block_test.vtkhdf', MPI_COMM_WORLD, stat, errmsg)
#else
  call vizfile%create('mb_fail_late_block_test.vtkhdf', stat, errmsg)
#endif
  if (stat /= 0) error stop errmsg

  hblk = vizfile%add_block('Block-A', mode=UG_FIXED_MESH)

  call vizfile%start_time_step(0.0_r8)

  ! Temporal blocks must be defined before time stepping starts.
  hblk = vizfile%add_block('Late-Block', mode=UG_FIXED_MESH)
  error stop 'late temporal block add unexpectedly succeeded'

  call vizfile%close()
  call test_env_finalize()

end program vtkhdf_mb_fail_late_block_test
