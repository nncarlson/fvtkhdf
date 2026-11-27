program triangles_static

  use, intrinsic :: iso_fortran_env, only: int8
  use vtkhdf_unstructured_file_type
  use vtkhdf_vtk_cell_types, only: VTK_TRIANGLE
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  type(vtkhdf_unstructured_file) :: file
  real :: points(3,5), temperature(4)
  integer :: cnode(12), xcnode(5), stat
  integer(int8) :: types(4)
  character(:), allocatable :: errmsg
#ifdef USE_MPI
  integer :: nproc, ierr
#endif

#ifdef USE_MPI
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  if (nproc /= 1) then
    call MPI_Finalize(ierr)
    error stop 'Run triangles_static with exactly 1 MPI process.'
  end if
#endif

  points = reshape([ &
    0.0, 0.0, 0.0, &
    1.0, 0.0, 0.0, &
    1.0, 1.0, 0.0, &
    0.0, 1.0, 0.0, &
    0.4, 0.55, 0.0], [3,5])
  cnode = [1, 2, 5,  2, 3, 5,  3, 4, 5,  4, 1, 5]
  xcnode = [1, 4, 7, 10, 13]
  types = [VTK_TRIANGLE, VTK_TRIANGLE, VTK_TRIANGLE, VTK_TRIANGLE]

#ifdef USE_MPI
  call file%create('triangles_static.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=VTKHDF_STATIC)
#else
  call file%create('triangles_static.vtkhdf', stat, errmsg, mode=VTKHDF_STATIC)
#endif
  if (stat /= 0) then
#ifdef USE_MPI
    call MPI_Finalize(ierr)
#endif
    error stop errmsg
  end if

  call file%write_mesh(points, cnode, xcnode, types)

  temperature = [10.0, 20.0, 30.0, 40.0]
  call file%write_cell_data('temperature', temperature)
  call file%close()

#ifdef USE_MPI
  call MPI_Finalize(ierr)
#endif

end program
