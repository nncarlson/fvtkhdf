program structured_moving

  use vtkhdf_structured_file_type
  use mpi
  implicit none

  integer, parameter :: nx = 24, ny = 18, local_nx = nx/2
  real, parameter :: pi = acos(-1.0)
  real, parameter :: first_angle = 100.0*pi/180.0, sector_angle = pi/3.0
  real, parameter :: inner_radius = 0.65, outer_radius = 1.25
  real, parameter :: hot_radius = 0.5*(inner_radius + outer_radius)
  type(vtkhdf_structured_file) :: file
  type(vtkhdf_point_data_handle) :: temperature_var
  real :: points(3,local_nx+1,ny+1,1)
  real :: temperature(local_nx+1,ny+1,1)
  integer :: point_rank(local_nx+1,ny+1,1)
  integer :: point_bounds(3,2), array_origin(3), write_bounds(3,2)
  integer :: rank, nproc, ierr, stat, step, i, j, global_i
  real :: radius, theta, x, y, time, rotation
  character(:), allocatable :: errmsg

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  if (nproc /= 2) then
    call MPI_Finalize(ierr)
    error stop 'Run structured_moving with exactly 2 MPI processes.'
  end if
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  point_bounds(:,1) = [0,0,0]
  point_bounds(:,2) = [nx,ny,0]
  array_origin = [rank*local_nx,0,0]

  write_bounds(:,1) = array_origin
  write_bounds(:,2) = array_origin + shape(temperature) - 1
  ! Rank 0 writes the shared point column.
  if (rank == 1) write_bounds(1,1) = write_bounds(1,1) + 1

  call file%create('structured_moving.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=VTKHDF_MOVING_MESH)
  if (stat /= 0) then
    call MPI_Finalize(ierr)
    error stop errmsg
  end if

  call file%write_mesh(point_bounds)
  point_rank = rank
  call file%write_point_data('mpi_rank', point_rank, array_origin, write_bounds)

  temperature_var = file%register_temporal_point_data('temperature', mold=0.0)

  do step = 0, 10
    time = real(step)/10.0
    call file%start_time_step(time)
    rotation = (40.0*pi/180.0)*time
    do j = 1, ny+1
      radius = inner_radius + (outer_radius - inner_radius)*real(j-1)/ny
      do i = 1, local_nx+1
        global_i = array_origin(1) + i - 1
        theta = first_angle - sector_angle*real(global_i)/nx + rotation
        x = radius*cos(theta)
        y = radius*sin(theta)
        points(:,i,j,1) = [x,y,0.0]
        temperature(i,j,1) = exp(-40.0*(x**2 + (y - hot_radius)**2))
      end do
    end do

    call file%write_mesh_geometry(points, array_origin, write_bounds)
    call file%write_point_data(temperature_var, temperature, array_origin, write_bounds)
    call file%finalize_time_step()
  end do

  call file%close()
  call MPI_Finalize(ierr)

end program
