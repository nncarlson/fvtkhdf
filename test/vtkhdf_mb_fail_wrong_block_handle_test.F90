program vtkhdf_mb_fail_wrong_block_handle_test

  use, intrinsic :: iso_fortran_env, only: r8 => real64, int8
  use vtkhdf_test_env
  use vtkhdf_test_data
  use vtkhdf_mb_file_type
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  type(vtkhdf_mb_file) :: vizfile
  type(vtkhdf_unstructured_block), pointer :: hblk_a => null(), hblk_b => null()
  type(vtkhdf_point_data_handle) :: hpoint_radius
  integer :: rank, nproc, stat
  character(:), allocatable :: errmsg
  real(r8), allocatable :: points(:,:), scalar_point_data(:)
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)

  call test_env_init(rank, nproc)

#ifdef USE_MPI
  call vizfile%create('mb_fail_wrong_block_handle_test.vtkhdf', MPI_COMM_WORLD, stat, errmsg)
#else
  call vizfile%create('mb_fail_wrong_block_handle_test.vtkhdf', stat, errmsg)
#endif
  if (stat /= 0) error stop errmsg

  call get_tet_cube_mesh_data(points, cnode, xcnode, types)
  call shift_points(points, dx=real(rank, r8))

  call vizfile%add_block('Block-A', hblk_a, mode=VTKHDF_FIXED_MESH)
  call hblk_a%write_mesh(points, cnode, xcnode, types)

  call vizfile%add_block('Block-B', hblk_b, mode=VTKHDF_FIXED_MESH)
  call hblk_b%write_mesh(points, cnode, xcnode, types)

  hpoint_radius = hblk_a%register_temporal_point_data('point-radius', 0.0_r8)
  call get_scalar_point_data(points, scalar_point_data)

  call vizfile%start_time_step(0.0_r8)

  ! Must fail: the temporal handle belongs to Block-A, not Block-B.
  call hblk_b%write_point_data(hpoint_radius, scalar_point_data)

  ! Return successfully if the expected assertion did not fire.
  call test_env_finalize()
end program vtkhdf_mb_fail_wrong_block_handle_test
