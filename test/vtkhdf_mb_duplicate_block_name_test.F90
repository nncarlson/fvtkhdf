program vtkhdf_mb_duplicate_block_name_test

  use, intrinsic :: iso_fortran_env, only: r8 => real64, int8
  use vtkhdf_test_env
  use vtkhdf_test_data
  use vtkhdf_mb_file_type
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  real(r8), allocatable :: points(:,:)
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  character(:), allocatable :: errmsg
  integer :: rank, nproc, stat

  type(vtkhdf_mb_file) :: vizfile
  type(vtkhdf_unstructured_block), pointer :: hblk => null()

  call test_env_init(rank, nproc)

#ifdef USE_MPI
  call vizfile%create('mb_duplicate_block_name_test.vtkhdf', MPI_COMM_WORLD, stat, errmsg)
#else
  call vizfile%create('mb_duplicate_block_name_test.vtkhdf', stat, errmsg)
#endif
  if (stat /= 0) error stop errmsg

  call get_tet_cube_mesh_data(points, cnode, xcnode, types)
  call shift_points(points, dx=real(rank, r8))

  call vizfile%add_block('Repeated Block', hblk, mode=VTKHDF_STATIC)
  call hblk%write_mesh(points, cnode, xcnode, types)

  call vizfile%add_block('Other Block', hblk, mode=VTKHDF_STATIC)
  call shift_points(points, dy=1.0_r8)
  call hblk%write_mesh(points, cnode, xcnode, types)

  call vizfile%add_block('Repeated Block', hblk, mode=VTKHDF_STATIC)
  call shift_points(points, dy=1.0_r8)
  call hblk%write_mesh(points, cnode, xcnode, types)

  call vizfile%close()
  call test_env_finalize()

end program vtkhdf_mb_duplicate_block_name_test
