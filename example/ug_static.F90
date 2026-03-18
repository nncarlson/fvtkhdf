program ug_static_example

  use,intrinsic :: iso_fortran_env, only: int8
  use vtkhdf_ug_file_type
  use vtkhdf_vtk_cell_types, only: VTK_TRIANGLE
#ifdef USE_MPI
  use mpi_f08
#endif
  implicit none

  real, allocatable :: points(:,:)
  real, allocatable :: pressure(:), velocity(:,:)
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  character(:), allocatable :: errmsg
  integer :: stat, nproc, rank, npoints, ncells, j

  type(vtkhdf_ug_file) :: file

#ifdef USE_MPI
  call MPI_Init(stat)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, stat)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, stat)
#else
  nproc = 1
  rank = 0
#endif

  ! Create a time-independent file.
#ifdef USE_MPI
  call file%create('example.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=UG_STATIC)
  if (stat /= 0) then
    call MPI_Finalize(stat)
    error stop errmsg
  end if
#else
  call file%create('example.vtkhdf', stat, errmsg, mode=UG_STATIC)
  if (stat /= 0) error stop errmsg
#endif

  ! 2D mesh. The local mesh on each rank is a 2-triangle discretization
  ! of a unit square, suitably x-shifted so they only overlap on a boundary.
  ! The mesh topologies are identical; point coordinates differ. 3D coordinates
  ! are required; we set z=0.

  points = reshape([0,0,0, 1,0,0, 1,1,0, 0,1,0], shape=[3,4])
  points(1,:) = points(1,:) + rank
  cnode = [1,2,3,  1,3,4]
  xcnode = [1,4,7]
  types = spread(VTK_TRIANGLE, dim=1, ncopies=2)

  ! local mesh sizes
  npoints = size(points,dim=2)
  ncells  = size(xcnode) - 1

  ! Each rank writes its part of the mesh collectively.
  call file%write_mesh(points, cnode, xcnode, types)

  ! Each rank writes its part of the cell-centered pressure and
  ! point-centered velocity data.
  pressure = [(rank + j*0.5, j=0,ncells-1)]
  call file%write_cell_data('pressure', pressure)

  velocity = points
  call file%write_point_data('velocity', velocity)

  ! For field data, only the payload from rank 0 is written but the call
  ! remains collective -- all ranks must call with identical arguments.
  call file%write_field_data('cpu time', 42.0)

  call file%close

#ifdef USE_MPI
  call MPI_Finalize(stat)
#endif

end program
