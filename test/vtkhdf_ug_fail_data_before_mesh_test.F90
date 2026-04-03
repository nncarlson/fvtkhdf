program vtkhdf_ug_fail_data_before_mesh_test

  use, intrinsic :: iso_fortran_env, only: r8 => real64
  use vtkhdf_test_env
  use vtkhdf_ug_file_type
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  type(vtkhdf_ug_file) :: vizfile
  integer :: rank, nproc, stat
  character(:), allocatable :: errmsg

  call test_env_init(rank, nproc)

#ifdef USE_MPI
  call vizfile%create('ug_fail_data_before_mesh_test.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=UG_STATIC)
#else
  call vizfile%create('ug_fail_data_before_mesh_test.vtkhdf', stat, errmsg, mode=UG_STATIC)
#endif
  if (stat /= 0) error stop errmsg

  ! Must fail: static data writes require mesh geometry and topology first.
  call vizfile%write_cell_data('cell-radius', [1.0_r8], attribute='Scalars')

  error stop 'expected data-before-mesh assertion'
end program vtkhdf_ug_fail_data_before_mesh_test
