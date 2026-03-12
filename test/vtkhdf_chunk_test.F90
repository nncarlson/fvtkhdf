program vtkhdf_chunk_test

  use,intrinsic :: iso_fortran_env, only: r8 => real64, int8, int64
  use vtkhdf_ug_file_type
  use vtkhdf_mb_file_type
  use vtkhdf_vtk_cell_types
#ifdef USE_MPI
  use mpi_f08
#endif
  implicit none

  real(r8), allocatable :: points(:,:)
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  character(:), allocatable :: errmsg
  integer :: stat
  integer :: rank

  type(vtkhdf_ug_file) :: ug_file
  type(vtkhdf_mb_file) :: mb_file
  type(vtkhdf_block_handle) :: block

#ifdef USE_MPI
  call MPI_Init(stat)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, stat)
#else
  rank = 0
#endif

  call get_mesh_data(points, cnode, xcnode, types)

#ifdef USE_MPI
  call ug_file%create('chunk_default.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=UG_FIXED_MESH)
#else
  call ug_file%create('chunk_default.vtkhdf', stat, errmsg, mode=UG_FIXED_MESH)
#endif
  if (stat /= 0) error stop errmsg
  call ug_file%write_mesh(points, cnode, xcnode, types)
  call ug_file%close()

#ifdef USE_MPI
  call mb_file%create('chunk_custom.vtkhdf', MPI_COMM_WORLD, stat, errmsg, data_chunk_bytes=1024_int64*1024)
#else
  call mb_file%create('chunk_custom.vtkhdf', stat, errmsg, data_chunk_bytes=1024_int64*1024)
#endif
  if (stat /= 0) error stop errmsg
  block = mb_file%add_block('Block', mode=UG_FIXED_MESH)
  call mb_file%write_mesh(block, points, cnode, xcnode, types)
  call mb_file%close()

#ifdef USE_MPI
  call MPI_Finalize(stat)
#endif

contains

  subroutine get_mesh_data(points, cnode, xcnode, types)
    real(r8), allocatable, intent(out) :: points(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable, intent(out) :: types(:)

    points = reshape([0.0_r8, 0.0_r8, 0.0_r8, &
        1.0_r8, 0.0_r8, 0.0_r8, &
        0.0_r8, 1.0_r8, 0.0_r8, &
        0.0_r8, 0.0_r8, 1.0_r8], [3,4])
    cnode = [1, 2, 3, 4]
    xcnode = [1, 5]
    types = [VTK_TETRA]
  end subroutine

end program vtkhdf_chunk_test
