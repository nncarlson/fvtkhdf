program vtkhdf_ug_fail_late_register_test

  use, intrinsic :: iso_fortran_env, only: r8 => real64, int8
  use vtkhdf_test_env
  use vtkhdf_test_data
  use vtkhdf_ug_file_type
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  type(vtkhdf_ug_file) :: vizfile
  integer :: rank, nproc, stat
  character(:), allocatable :: errmsg
  real(r8), allocatable :: points(:,:)
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  type(vtkhdf_point_data_handle) :: hpoint_radius

  call test_env_init(rank, nproc)

#ifdef USE_MPI
  call vizfile%create('ug_fail_late_register_test.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=UG_FIXED_MESH)
#else
  call vizfile%create('ug_fail_late_register_test.vtkhdf', stat, errmsg, mode=UG_FIXED_MESH)
#endif
  if (stat /= 0) error stop errmsg

  call get_tet_cube_mesh_data(points, cnode, xcnode, types)
  call shift_points(points, dx=real(rank, r8))
  call vizfile%write_mesh(points, cnode, xcnode, types)
  call vizfile%start_time_step(0.0_r8)

  ! Must fail: temporal registration is not allowed after stepping has started.
  hpoint_radius = vizfile%register_temporal_point_data('late-point-radius', 0.0_r8)

  error stop 'expected late registration assertion'
end program vtkhdf_ug_fail_late_register_test
