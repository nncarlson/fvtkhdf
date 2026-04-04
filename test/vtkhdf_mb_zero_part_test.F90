program vtkhdf_mb_zero_part_test

  use, intrinsic :: iso_fortran_env, only: r8 => real64, int8
  use vtkhdf_test_env
  use vtkhdf_test_data
  use vtkhdf_mb_file_type
  use mpi
  implicit none

  real(r8), allocatable :: points_full(:,:), points_sparse(:,:), points_zero(:,:)
  real(r8), allocatable :: scalar_data(:), zero_scalar_data(:)
  integer, allocatable :: cnode(:), xcnode(:), cnode_zero(:), xcnode_zero(:)
  integer(int8), allocatable :: types(:), types_zero(:)
  character(:), allocatable :: errmsg
  integer :: rank, nproc, stat

  type(vtkhdf_mb_file) :: vizfile
  type(vtkhdf_block_handle) :: hblk_full, hblk_sparse
  type(vtkhdf_point_data_handle) :: hpoint_radius

  call test_env_init(rank, nproc)

  call vizfile%create('mb_zero_part_test.vtkhdf', MPI_COMM_WORLD, stat, errmsg)
  if (stat /= 0) error stop errmsg

  call get_tet_cube_mesh_data(points_full, cnode, xcnode, types)
  call shift_points(points_full, dx=real(rank, r8))

  hblk_full = vizfile%add_block('Block-Full', mode=UG_FIXED_MESH)
  call vizfile%write_mesh(hblk_full, points_full, cnode, xcnode, types)

  hblk_sparse = vizfile%add_block('Block-Sparse', mode=UG_FIXED_MESH)
  if (rank == 0) then
    call get_zero_mesh_data(points_zero, cnode_zero, xcnode_zero, types_zero)
    call vizfile%write_mesh(hblk_sparse, points_zero, cnode_zero, xcnode_zero, types_zero)
  else
    points_sparse = points_full
    call shift_points(points_sparse, dy=1.0_r8)
    call vizfile%write_mesh(hblk_sparse, points_sparse, cnode, xcnode, types)
  end if

  hpoint_radius = vizfile%register_temporal_point_data(hblk_sparse, 'point-radius', 0.0_r8)

  call vizfile%start_time_step(0.0_r8)
  if (rank == 0) then
    allocate(zero_scalar_data(0))
    call vizfile%write_point_data(hblk_sparse, hpoint_radius, zero_scalar_data)
  else
    call get_scalar_point_data(points_sparse, scalar_data)
    call vizfile%write_point_data(hblk_sparse, hpoint_radius, scalar_data)
  end if
  call vizfile%finalize_time_step()

  call vizfile%start_time_step(1.0_r8)
  if (rank == 0) then
    call vizfile%write_point_data(hblk_sparse, hpoint_radius, zero_scalar_data)
  else
    call vizfile%write_point_data(hblk_sparse, hpoint_radius, scalar_data + 1.0_r8)
  end if
  call vizfile%finalize_time_step()

  call vizfile%close()
  call test_env_finalize()

end program vtkhdf_mb_zero_part_test
