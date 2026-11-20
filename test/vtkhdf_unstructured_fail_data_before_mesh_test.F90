program vtkhdf_unstructured_fail_data_before_mesh_test

  use, intrinsic :: iso_fortran_env, only: r8 => real64
  use vtkhdf_test_env
  use vtkhdf_unstructured_file_type
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  type(vtkhdf_unstructured_file) :: vizfile
  integer :: rank, nproc, stat
  character(:), allocatable :: errmsg

  call test_env_init(rank, nproc)

#ifdef USE_MPI
  call vizfile%create('unstructured_fail_data_before_mesh_test.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=VTKHDF_STATIC)
#else
  call vizfile%create('unstructured_fail_data_before_mesh_test.vtkhdf', stat, errmsg, mode=VTKHDF_STATIC)
#endif
  if (stat /= 0) error stop errmsg

  ! Must fail: static data writes require mesh geometry and topology first.
  call vizfile%write_cell_data('cell-radius', [1.0_r8], attribute='Scalars')

  ! Return successfully if the expected assertion did not fire.
  call test_env_finalize()
end program vtkhdf_unstructured_fail_data_before_mesh_test
