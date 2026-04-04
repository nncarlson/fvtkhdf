program vtkhdf_mb_test

  use, intrinsic :: iso_fortran_env, only: r8 => real64, int8
  use vtkhdf_test_env
  use vtkhdf_test_data
  use vtkhdf_mb_file_type
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  real(r8), allocatable :: points_a(:,:), points_b(:,:)
  real(r8), allocatable :: scalar_cell_data(:), vector_point_data(:,:)
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  character(:), allocatable :: errmsg
  integer :: rank, nproc, stat

  type(vtkhdf_mb_file) :: vizfile
  type(vtkhdf_block_handle) :: hblk_a, hblk_b
  type(vtkhdf_point_data_handle) :: hpoint_velocity

  call test_env_init(rank, nproc)

#ifdef USE_MPI
  call vizfile%create('mb_test.vtkhdf', MPI_COMM_WORLD, stat, errmsg)
#else
  call vizfile%create('mb_test.vtkhdf', stat, errmsg)
#endif
  if (stat /= 0) error stop errmsg

  call get_tet_cube_mesh_data(points_a, cnode, xcnode, types)
  call shift_points(points_a, dx=real(rank, r8))
  points_b = points_a
  call shift_points(points_b, dy=1.0_r8)

  hblk_a = vizfile%add_block('Block-A', mode=UG_FIXED_MESH)
  call vizfile%write_mesh(hblk_a, points_a, cnode, xcnode, types)

  hblk_b = vizfile%add_block('Block-B', mode=UG_FIXED_MESH)
  call vizfile%write_mesh(hblk_b, points_b, cnode, xcnode, types)

  hpoint_velocity = vizfile%register_temporal_point_data(hblk_b, 'point-velocity', [real(r8) :: 0, 0, 0], &
      attribute='Vectors')

  call get_scalar_cell_data(points_a, cnode, xcnode, scalar_cell_data)
  call get_vector_point_data(points_b, vector_point_data)

  call vizfile%start_time_step(0.0_r8)
  call vizfile%write_point_data(hblk_b, hpoint_velocity, vector_point_data)
  call vizfile%finalize_time_step()

  call vizfile%start_time_step(1.0_r8)
  call vizfile%write_point_data(hblk_b, hpoint_velocity, vector_point_data + 1.0_r8)
  call vizfile%finalize_time_step()

  call vizfile%write_cell_data(hblk_a, 'static-cell-scalar', -scalar_cell_data, attribute='Scalars')

  call vizfile%close()
  call test_env_finalize()

end program vtkhdf_mb_test
