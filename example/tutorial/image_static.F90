program image_static

  use vtkhdf_image_file_type
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  integer, parameter :: nx = 24, ny = 18
  real, parameter :: dx = 1.0/nx, dy = 0.75/ny
  type(vtkhdf_image_file) :: file
  real, allocatable :: temperature(:,:,:)
  integer :: point_bounds(3,2), array_origin(3)
  integer :: i, j, stat
  real :: x, y
  character(:), allocatable :: errmsg
#ifdef USE_MPI
  integer :: nproc, ierr
#endif

#ifdef USE_MPI
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  if (nproc /= 1) then
    call MPI_Finalize(ierr)
    error stop 'Run image_static with exactly 1 MPI process.'
  end if
#endif

  point_bounds(:,1) = [0,0,0]
  point_bounds(:,2) = [nx,ny,0]
  array_origin = [0,0,0]
  allocate(temperature(nx,ny,1))
  do j = 1, ny
    y = (real(j) - 0.5)*dy
    do i = 1, nx
      x = (real(i) - 0.5)*dx
      temperature(i,j,1) = exp(-40.0*((x - 0.35)**2 + (y - 0.45)**2))
    end do
  end do

#ifdef USE_MPI
  call file%create('image_static.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=VTKHDF_STATIC)
#else
  call file%create('image_static.vtkhdf', stat, errmsg, mode=VTKHDF_STATIC)
#endif
  if (stat /= 0) then
#ifdef USE_MPI
    call MPI_Finalize(ierr)
#endif
    error stop errmsg
  end if

  call file%write_mesh(point_bounds, origin=[0.0,0.0,0.0], spacing=[dx,dy,1.0])
  call file%write_cell_data('temperature', temperature, array_origin)

  call file%close()
#ifdef USE_MPI
  call MPI_Finalize(ierr)
#endif

end program
