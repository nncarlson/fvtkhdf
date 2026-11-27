program structured_parallel

  use vtkhdf_structured_file_type
  use mpi
  implicit none

  integer, parameter :: nx = 24, ny = 18, local_nx = nx/2
  real, parameter :: pi = acos(-1.0)
  real, parameter :: first_angle = 100.0*pi/180.0, sector_angle = pi/3.0
  real, parameter :: inner_radius = 0.65, outer_radius = 1.25
  real, parameter :: hot_radius = 0.5*(inner_radius + outer_radius)
  type(vtkhdf_structured_file) :: file
  real :: points(3,local_nx+1,ny+1,1)
  real :: temperature(local_nx+1,ny+1,1)
  integer :: point_rank(local_nx+1,ny+1,1)
  integer :: point_bounds(3,2), array_origin(3), write_bounds(3,2)
  integer :: rank, nproc, ierr, stat, i, j, global_i
  real :: radius, theta, x, y
  character(:), allocatable :: errmsg

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  if (nproc /= 2) then
    call MPI_Finalize(ierr)
    error stop 'Run structured_parallel with exactly 2 MPI processes.'
  end if
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  point_bounds(:,1) = [0,0,0]
  point_bounds(:,2) = [nx,ny,0]
  array_origin = [rank*local_nx,0,0]

  write_bounds(:,1) = array_origin
  write_bounds(:,2) = array_origin + shape(temperature) - 1
  ! Rank 0 writes the shared point column.
  if (rank == 1) write_bounds(1,1) = write_bounds(1,1) + 1

  do j = 1, ny+1
    radius = inner_radius + (outer_radius - inner_radius)*real(j-1)/ny
    do i = 1, local_nx+1
      global_i = array_origin(1) + i - 1
      theta = first_angle - sector_angle*real(global_i)/nx
      x = radius*cos(theta)
      y = radius*sin(theta)
      points(:,i,j,1) = [x,y,0.0]
      temperature(i,j,1) = exp(-40.0*(x**2 + (y - hot_radius)**2))
    end do
  end do

  call file%create('structured_parallel.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=VTKHDF_STATIC)
  if (stat /= 0) then
    call MPI_Finalize(ierr)
    error stop errmsg
  end if

  call file%write_mesh(point_bounds)
  call file%write_mesh_geometry(points, array_origin, write_bounds)
  point_rank = rank
  call file%write_point_data('mpi_rank', point_rank, array_origin, write_bounds)
  call file%write_point_data('temperature', temperature, array_origin, write_bounds)

  call file%close()
  call MPI_Finalize(ierr)

end program
