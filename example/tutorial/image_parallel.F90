program image_parallel

  use vtkhdf_image_file_type
  use mpi
  implicit none

  integer, parameter :: nx = 24, ny = 18
  real, parameter :: dx = 1.0/nx, dy = 0.75/ny
  type(vtkhdf_image_file) :: file
  real, allocatable :: temperature(:,:,:)
  integer :: point_bounds(3,2), array_origin(3)
  integer :: i, j, stat
  real :: x, y
  character(:), allocatable :: errmsg
  type(vtkhdf_cell_data_handle) :: temperature_var
  integer :: step
  real :: time, x_center
  integer :: nproc, ierr
  integer, parameter :: local_nx = nx/2
  integer :: rank
  integer :: cell_rank(local_nx,ny,1)

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  if (nproc /= 2) then
    call MPI_Finalize(ierr)
    error stop 'Run image_parallel with exactly 2 MPI processes.'
  end if
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  point_bounds(:,1) = [0,0,0]
  point_bounds(:,2) = [nx,ny,0]
  array_origin = [rank*local_nx,0,0]
  allocate(temperature(local_nx,ny,1))

  call file%create('image_parallel.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=VTKHDF_FIXED_MESH)
  if (stat /= 0) then
    call MPI_Finalize(ierr)
    error stop errmsg
  end if

  call file%write_mesh(point_bounds, origin=[0.0,0.0,0.0], spacing=[dx,dy,1.0])

  cell_rank = rank
  call file%write_cell_data('mpi_rank', cell_rank, array_origin)

  temperature_var = file%register_temporal_cell_data('temperature', mold=0.0)

  do step = 0, 10
    time = real(step)/10.0
    call file%start_time_step(time)
    x_center = 0.35 + 0.35*time
    do j = 1, ny
      y = (real(j) - 0.5)*dy
      do i = 1, local_nx
        x = (array_origin(1) + real(i) - 0.5)*dx
        temperature(i,j,1) = exp(-40.0*((x - x_center)**2 + (y - 0.45)**2))
      end do
    end do
    call file%write_cell_data(temperature_var, temperature, array_origin)
    call file%finalize_time_step()
  end do

  call file%close()
  call MPI_Finalize(ierr)

end program
