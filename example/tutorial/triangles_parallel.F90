program triangles_parallel

  use, intrinsic :: iso_fortran_env, only: int8
  use vtkhdf_unstructured_file_type
  use vtkhdf_vtk_cell_types, only: VTK_TRIANGLE
  use mpi
  implicit none

  type(vtkhdf_unstructured_file) :: file
  real :: points(3,5), local_points(3,4), temperature(2)
  integer :: cnode(6), xcnode(3), cell_rank(2)
  integer(int8) :: types(2)
  integer :: rank, nproc, stat, ierr
  character(:), allocatable :: errmsg

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  if (nproc /= 2) then
    call MPI_Finalize(ierr)
    error stop 'Run triangles_parallel with exactly 2 MPI processes.'
  end if
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  points = reshape([ &
    0.0, 0.0, 0.0, &
    1.0, 0.0, 0.0, &
    1.0, 1.0, 0.0, &
    0.0, 1.0, 0.0, &
    0.4, 0.55, 0.0], [3,5])
  xcnode = [1, 4, 7]
  types = VTK_TRIANGLE

  if (rank == 0) then
    local_points = points(:, [1, 2, 3, 5])
    cnode = [1, 2, 4,  2, 3, 4]
    temperature = [10.0, 20.0]
  else
    local_points = points(:, [1, 3, 4, 5])
    cnode = [2, 3, 4,  3, 1, 4]
    temperature = [30.0, 40.0]
  end if

  call file%create('triangles_parallel.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=VTKHDF_STATIC)
  if (stat /= 0) then
    call MPI_Finalize(ierr)
    error stop errmsg
  end if

  call file%write_mesh(local_points, cnode, xcnode, types)
  cell_rank = rank
  call file%write_cell_data('mpi_rank', cell_rank)
  call file%write_cell_data('temperature', temperature)
  call file%close()
  call MPI_Finalize(ierr)

end program
