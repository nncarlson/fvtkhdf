program vtkhdf_2d_test

  use,intrinsic :: iso_fortran_env, only: r8 => real64, int8
  use vtkhdf_ug_file_type
  use vtkhdf_vtk_cell_types
#ifdef USE_MPI
  use mpi_f08
#endif
  implicit none

  real(r8), allocatable :: x(:,:), scalar_cell_data(:), vector_cell_data(:,:)
  real(r8), allocatable :: scalar_point_data(:), vector_point_data(:,:)
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  character(:), allocatable :: errmsg
  integer :: stat
#ifdef USE_MPI
  integer :: istat
#endif

  type(vtkhdf_ug_file) :: vizfile

#ifdef USE_MPI
  call MPI_Init(istat)
  call vizfile%create('ug_2d_test.vtkhdf', MPI_COMM_WORLD, stat, errmsg)
#else
  call vizfile%create('ug_2d_test.vtkhdf', stat, errmsg)
#endif
  if (stat /= 0) error stop errmsg

  ! A 2D mesh (two triangles forming a square)
  x = reshape([0.0_r8,0.0_r8, 1.0_r8,0.0_r8, 1.0_r8,1.0_r8, 0.0_r8,1.0_r8], shape=[2,4])
  cnode = [1,2,3, 1,3,4]
  xcnode = [1,4,7]
  types = [VTK_TRIANGLE, VTK_TRIANGLE]

  call vizfile%write_mesh(x, cnode, xcnode, types)

  !! Add some 2D point data to make sure it still works
  allocate(scalar_point_data(4))
  scalar_point_data = [1.0_r8, 2.0_r8, 3.0_r8, 4.0_r8]
  call vizfile%write_point_data('scalar-data', scalar_point_data)

  allocate(vector_point_data(2,4))
  vector_point_data = x
  call vizfile%write_point_data('vector-data', vector_point_data)

  call vizfile%close
#ifdef USE_MPI
  call MPI_Finalize(istat)
#endif

end program
