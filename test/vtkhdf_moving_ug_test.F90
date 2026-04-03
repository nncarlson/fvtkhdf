program vtkhdf_moving_ug_test

  use, intrinsic :: iso_fortran_env, only: r8 => real64, int8
  use vtkhdf_test_env
  use vtkhdf_test_data
  use vtkhdf_ug_file_type
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  real(r8), allocatable :: points(:,:), scalar_cell_data(:), vector_cell_data(:,:)
  real(r8), allocatable :: scalar_point_data(:), vector_point_data(:,:)
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  character(:), allocatable :: errmsg
  integer :: rank, nproc, stat

  type(vtkhdf_ug_file) :: vizfile
  type(vtkhdf_cell_data_handle) :: hcell_radius, hcell_velocity
  type(vtkhdf_point_data_handle) :: hpoint_radius, hpoint_velocity

  call test_env_init(rank, nproc)

#ifdef USE_MPI
  call vizfile%create('moving_ug_test.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=UG_MOVING_MESH)
#else
  call vizfile%create('moving_ug_test.vtkhdf', stat, errmsg, mode=UG_MOVING_MESH)
#endif
  if (stat /= 0) error stop errmsg

  call get_tet_cube_mesh_data(points, cnode, xcnode, types)
  call shift_points(points, dx=real(rank, r8))
  call vizfile%write_mesh_topology(cnode, xcnode, types)

  associate (scalar_mold => 0.0_r8, vector_mold => [real(r8) :: 0, 0, 0])
    hcell_radius = vizfile%register_temporal_cell_data('cell-radius', scalar_mold)
    hcell_velocity = vizfile%register_temporal_cell_data('cell-velocity', vector_mold)
    hpoint_radius = vizfile%register_temporal_point_data('point-radius', scalar_mold)
    hpoint_velocity = vizfile%register_temporal_point_data('point-velocity', vector_mold)
  end associate

  call get_scalar_cell_data(points, cnode, xcnode, scalar_cell_data)
  call get_vector_cell_data(points, cnode, xcnode, vector_cell_data)
  call get_scalar_point_data(points, scalar_point_data)
  call get_vector_point_data(points, vector_point_data)

  call vizfile%start_time_step(0.0_r8)
  call shift_points(points, dz=-0.1_r8)
  call vizfile%write_mesh_geometry(points)
  call vizfile%write_cell_data(hcell_radius, scalar_cell_data)
  call vizfile%write_cell_data(hcell_velocity, vector_cell_data)
  call vizfile%write_point_data(hpoint_radius, scalar_point_data)
  call vizfile%write_point_data(hpoint_velocity, vector_point_data)

  call vizfile%start_time_step(1.0_r8)
  call shift_points(points, dz=-0.1_r8)
  call vizfile%write_mesh_geometry(points)
  call vizfile%write_cell_data(hcell_radius, scalar_cell_data + 1.0_r8)
  call vizfile%write_cell_data(hcell_velocity, vector_cell_data + 1.0_r8)
  call vizfile%write_point_data(hpoint_radius, scalar_point_data + 1.0_r8)
  call vizfile%write_point_data(hpoint_velocity, vector_point_data + 1.0_r8)

  call vizfile%close()
  call test_env_finalize()

end program vtkhdf_moving_ug_test
