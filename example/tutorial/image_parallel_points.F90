program image_parallel_points

  use vtkhdf_image_file_type
  use mpi
  implicit none

  integer, parameter :: nx = 24, ny = 18
  real, parameter :: dx = 1.0/nx, dy = 0.75/ny
  type(vtkhdf_image_file) :: file
  real, allocatable :: temperature(:,:,:)
  integer :: point_bounds(3,2), array_origin(3), write_bounds(3,2)
  integer :: i, j, stat
  real :: x, y
  character(:), allocatable :: errmsg
  type(vtkhdf_point_data_handle) :: temperature_var
  integer :: step
  real :: time, x_center
  integer :: nproc, ierr
  integer, parameter :: local_nx = nx/2
  integer :: rank
  integer :: point_rank(local_nx+1,ny+1,1)

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  if (nproc /= 2) then
    call MPI_Finalize(ierr)
    error stop 'Run image_parallel_points with exactly 2 MPI processes.'
  end if
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  point_bounds(:,1) = [0,0,0]
  point_bounds(:,2) = [nx,ny,0]
  array_origin = [rank*local_nx,0,0]
  allocate(temperature(local_nx+1,ny+1,1))

  write_bounds(:,1) = array_origin
  write_bounds(:,2) = array_origin + shape(temperature) - 1
  ! Rank 0 writes the shared point column.
  if (rank == 1) write_bounds(1,1) = write_bounds(1,1) + 1

  call file%create('image_parallel_points.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=VTKHDF_FIXED_MESH)
  if (stat /= 0) then
    call MPI_Finalize(ierr)
    error stop errmsg
  end if

  call file%write_mesh(point_bounds, origin=[0.0,0.0,0.0], spacing=[dx,dy,1.0])

  point_rank = rank
  call file%write_point_data('mpi_rank', point_rank, array_origin, write_bounds)

  temperature_var = file%register_temporal_point_data('temperature', mold=0.0)

  do step = 0, 10
    time = real(step)/10.0
    call file%start_time_step(time)
    x_center = 0.35 + 0.35*time
    do j = 1, ny+1
      y = real(j - 1)*dy
      do i = 1, local_nx+1
        x = real(array_origin(1) + i - 1)*dx
        temperature(i,j,1) = exp(-40.0*((x - x_center)**2 + (y - 0.45)**2))
      end do
    end do
    call file%write_point_data(temperature_var, temperature, array_origin, write_bounds)
    call file%finalize_time_step()
  end do

  call file%close()
  call MPI_Finalize(ierr)

end program
