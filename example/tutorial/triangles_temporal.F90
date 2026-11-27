program triangles_temporal

  use, intrinsic :: iso_fortran_env, only: int8
  use vtkhdf_unstructured_file_type
  use vtkhdf_vtk_cell_types, only: VTK_TRIANGLE
  use mpi
  implicit none

  integer, parameter :: nx = 24, ny = 18
  type(vtkhdf_unstructured_file) :: file
  type(vtkhdf_point_data_handle) :: temperature_var
  real, allocatable :: points(:,:), temperature(:)
  integer, allocatable :: cnode(:), xcnode(:), cell_rank(:)
  integer(int8), allocatable :: types(:)
  real :: lower(2), upper(2), time, x_center
  integer :: rank, nproc, step, stat, ierr
  character(:), allocatable :: errmsg

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  if (nproc /= 2) then
    call MPI_Finalize(ierr)
    error stop 'Run triangles_temporal with exactly 2 MPI processes.'
  end if
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  lower = [0.5*rank, 0.0]
  upper = [lower(1) + 0.5, 0.75]
  call make_triangle_mesh(lower, upper, nx/2, ny, points, cnode, xcnode, types)

  call file%create('triangles_temporal.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=VTKHDF_FIXED_MESH)
  if (stat /= 0) then
    call MPI_Finalize(ierr)
    error stop errmsg
  end if

  call file%write_mesh(points, cnode, xcnode, types)
  allocate(cell_rank(size(types)), source=rank)
  call file%write_cell_data('mpi_rank', cell_rank)

  allocate(temperature(size(points, 2)))
  temperature_var = file%register_temporal_point_data('temperature', mold=0.0)

  do step = 0, 10
    time = real(step)/10.0
    call file%start_time_step(time)
    x_center = 0.35 + 0.35*time
    temperature = exp(-40.0*((points(1,:) - x_center)**2 + (points(2,:) - 0.45)**2))
    call file%write_point_data(temperature_var, temperature)
    call file%finalize_time_step()
  end do

  call file%close()
  call MPI_Finalize(ierr)

contains

  subroutine make_triangle_mesh(lower, upper, nx, ny, points, cnode, xcnode, types)
    real, intent(in) :: lower(2), upper(2)
    integer, intent(in) :: nx, ny
    real, allocatable, intent(out) :: points(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable, intent(out) :: types(:)

    real :: spacing(2)
    integer :: nvertices, ntriangles, i, j, point, square, center, first_node
    integer :: lower_left, lower_right, upper_right, upper_left

    spacing = (upper - lower)/[nx, ny]
    nvertices = (nx + 1)*(ny + 1)
    ntriangles = 4*nx*ny
    allocate(points(3,nvertices + nx*ny))
    allocate(cnode(3*ntriangles), xcnode(ntriangles + 1), types(ntriangles))

    do j = 0, ny
      do i = 0, nx
        point = 1 + i + j*(nx + 1)
        points(:,point) = [lower + [i, j]*spacing, 0.0]
      end do
    end do

    do j = 0, ny - 1
      do i = 0, nx - 1
        square = 1 + i + j*nx
        center = nvertices + square
        lower_left = 1 + i + j*(nx + 1)
        lower_right = lower_left + 1
        upper_left = lower_left + nx + 1
        upper_right = upper_left + 1
        points(:,center) = 0.5*(points(:,lower_left) + points(:,upper_right))

        first_node = 1 + 12*(square - 1)
        cnode(first_node:first_node+11) = [ &
          lower_left,  lower_right, center, &
          lower_right, upper_right, center, &
          upper_right, upper_left,  center, &
          upper_left,  lower_left,  center]
      end do
    end do

    xcnode = [(1 + 3*i, i=0,ntriangles)]
    types = VTK_TRIANGLE
  end subroutine

end program
