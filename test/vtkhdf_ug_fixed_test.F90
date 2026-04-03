program vtkhdf_ug_fixed_test

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
  real(r8), allocatable :: points(:,:), scalar_cell_data(:), vector_point_data(:,:)
  real(r8), allocatable :: scalar_field_data(:)
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)

  type(vtkhdf_cell_data_handle) :: hcell_radius
  type(vtkhdf_point_data_handle) :: hpoint_velocity
  type(vtkhdf_field_data_handle) :: hfield_scalar

  call test_env_init(rank, nproc)

#ifdef USE_MPI
  call vizfile%create('ug_fixed_test.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=UG_FIXED_MESH)
#else
  call vizfile%create('ug_fixed_test.vtkhdf', stat, errmsg, mode=UG_FIXED_MESH)
#endif
  if (stat /= 0) error stop errmsg

  call get_tet_cube_mesh_data(points, cnode, xcnode, types)
  call shift_points(points, dx=real(rank, r8))
  call vizfile%write_mesh(points, cnode, xcnode, types)

  hcell_radius = vizfile%register_temporal_cell_data('cell-radius', 0.0_r8, attribute='Scalars')
  hpoint_velocity = vizfile%register_temporal_point_data('point-velocity', [real(r8) :: 0, 0, 0], attribute='Vectors')
  hfield_scalar = vizfile%register_temporal_field_data('field-scalar', 0.0_r8)

  call get_scalar_cell_data(points, cnode, xcnode, scalar_cell_data)
  call get_vector_point_data(points, vector_point_data)
  call get_scalar_field_data(scalar_field_data)

  call vizfile%start_time_step(0.0_r8)
  call vizfile%write_cell_data(hcell_radius, scalar_cell_data)
  call vizfile%write_point_data(hpoint_velocity, vector_point_data)
  call vizfile%write_field_data(hfield_scalar, scalar_field_data, as_vector=.true.)
  call vizfile%finalize_time_step()

  call vizfile%start_time_step(1.0_r8)
  call vizfile%write_cell_data(hcell_radius, scalar_cell_data + 1.0_r8)
  call vizfile%write_point_data(hpoint_velocity, vector_point_data + 1.0_r8)
  ! Skip one temporal write to exercise repeated-offset behavior.
  call vizfile%finalize_time_step()

  call vizfile%close()
  call test_env_finalize()

end program vtkhdf_ug_fixed_test
