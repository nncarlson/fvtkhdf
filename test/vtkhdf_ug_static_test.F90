program vtkhdf_ug_static_test

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
  real(r8), allocatable :: points(:,:), scalar_cell_data(:), vector_cell_data(:,:)
  real(r8), allocatable :: scalar_point_data(:), vector_point_data(:,:)
  real(r8), allocatable :: scalar_field_data(:), vector_field_data(:,:)
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)

  call test_env_init(rank, nproc)

#ifdef USE_MPI
  call vizfile%create('ug_static_test.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=UG_STATIC)
#else
  call vizfile%create('ug_static_test.vtkhdf', stat, errmsg, mode=UG_STATIC)
#endif
  if (stat /= 0) error stop errmsg

  call get_tet_cube_mesh_data(points, cnode, xcnode, types)
  call shift_points(points, dx=real(rank, r8))
  call vizfile%write_mesh(points, cnode, xcnode, types)

  call get_scalar_cell_data(points, cnode, xcnode, scalar_cell_data)
  call get_vector_cell_data(points, cnode, xcnode, vector_cell_data)
  call get_scalar_point_data(points, scalar_point_data)
  call get_vector_point_data(points, vector_point_data)
  call get_scalar_field_data(scalar_field_data)
  call get_vector_field_data(vector_field_data)

  call vizfile%write_cell_data('cell-radius', scalar_cell_data, attribute='Scalars')
  call vizfile%write_cell_data('cell-velocity', vector_cell_data)
  call vizfile%write_point_data('point-radius', scalar_point_data)
  call vizfile%write_point_data('point-velocity', vector_point_data, attribute='Vectors')
  call vizfile%write_field_data('field-value', 42.0_r8)
  call vizfile%write_field_data('field-scalar', scalar_field_data, as_vector=.true.)
  call vizfile%write_field_data('field-vector', vector_field_data)

  call vizfile%close()
  call test_env_finalize()

end program vtkhdf_ug_static_test
