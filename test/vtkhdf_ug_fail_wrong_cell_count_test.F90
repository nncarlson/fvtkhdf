program vtkhdf_ug_fail_wrong_cell_count_test

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
  real(r8), allocatable :: points(:,:), scalar_cell_data(:), bad_scalar_cell_data(:)
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  type(vtkhdf_cell_data_handle) :: hcell_radius

  call test_env_init(rank, nproc)

#ifdef USE_MPI
  call vizfile%create('ug_fail_wrong_cell_count_test.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=UG_FIXED_MESH)
#else
  call vizfile%create('ug_fail_wrong_cell_count_test.vtkhdf', stat, errmsg, mode=UG_FIXED_MESH)
#endif
  if (stat /= 0) error stop errmsg

  call get_tet_cube_mesh_data(points, cnode, xcnode, types)
  call shift_points(points, dx=real(rank, r8))
  call vizfile%write_mesh(points, cnode, xcnode, types)
  hcell_radius = vizfile%register_temporal_cell_data('cell-radius', 0.0_r8)

  call get_scalar_cell_data(points, cnode, xcnode, scalar_cell_data)
  allocate(bad_scalar_cell_data(size(scalar_cell_data) - 1))
  bad_scalar_cell_data = scalar_cell_data(:size(bad_scalar_cell_data))

  call vizfile%start_time_step(0.0_r8)

  ! Must fail: temporal cell data length must match the local cell count.
  call vizfile%write_cell_data(hcell_radius, bad_scalar_cell_data)

  error stop 'expected wrong cell count assertion'
end program vtkhdf_ug_fail_wrong_cell_count_test
