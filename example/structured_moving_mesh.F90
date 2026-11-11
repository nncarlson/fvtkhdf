program structured_moving_mesh_example

  use vtkhdf_structured_file_type
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  type(vtkhdf_structured_file) :: file
  type(vtkhdf_point_data_handle) :: temp_var
  type(vtkhdf_field_data_handle) :: cpu_var
  real, allocatable :: temp(:,:,:), diff(:,:,:)
  real, allocatable :: points(:,:,:,:)
  real, parameter :: pi = acos(-1.0), dtheta = pi/12.0, dr = 0.5
  real :: theta, radius
  integer :: rank, nproc, stat, step, i, j
  integer :: point_bounds(3,2), array_origin(3), write_bounds(3,2)
  character(:), allocatable :: errmsg

#ifdef USE_MPI
  call MPI_Init(stat)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,stat)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,stat)
#else
  rank = 0
  nproc = 1
#endif

  ! Create the file with support for moving geometry and time-dependent data.
#ifdef USE_MPI
  call file%create('structured_moving_mesh.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=VTKHDF_MOVING_MESH)
  if (stat /= 0) then
    call MPI_Finalize(stat)
    error stop errmsg
  end if
#else
  call file%create('structured_moving_mesh.vtkhdf', stat, errmsg, mode=VTKHDF_MOVING_MESH)
  if (stat /= 0) error stop errmsg
#endif

  ! Define the global 2D logical grid.
  ! Each rank owns one unique 4-wide slab of cells.
  point_bounds(:,1) = [0,0,0]
  point_bounds(:,2) = [4*nproc,3,0]
  call file%write_mesh(point_bounds)
  array_origin = [4*rank,0,0] ! global index of the first local array element

  ! Map the grid to an annular strip with uniform angular and radial increments.
  ! Raise z by one per revolution so successive turns do not overlap.
  allocate(points(3,5,4,1))
  do j = 1, 4
    radius = 2.0 + dr*(j-1)
    do i = 1, 5
      theta = dtheta*(array_origin(1)+i-1)
      points(:,i,j,1) = [radius*cos(theta),radius*sin(theta),theta/(2.0*pi)]
    end do
  end do
  ! Shared points are owned by the slab to the left, for geometry and point data.
  write_bounds(:,1) = array_origin
  write_bounds(:,2) = array_origin + [4,3,0]
  write_bounds(1,1) = write_bounds(1,1) + merge(0,1,(rank==0)) ! omit unowned overlap
  call file%write_mesh_geometry(points, array_origin, write_bounds)

  ! Write time-independent (static) diffusivity variable on cells
  allocate(diff(4,3,1), source=real(rank+1))
  call file%write_cell_data('diffusivity', diff, array_origin)

  ! Register time-dependent temperature on points.
  ! Use the same shared-point ownership as for the geometry.
  allocate(temp(5,4,1))
  temp_var = file%register_temporal_point_data('temperature', mold=0.0)

  ! We'll write the cummulative cpu time as field data.
  cpu_var = file%register_temporal_field_data('cpu_time', mold=0.0)

  do step = 0, 4
    call file%start_time_step(real(step))
    if (step > 0) then
      ! Rotate about z by dtheta per step, keeping each point's height fixed.
      do j = 1, 4
        radius = 2.0 + dr*(j-1)
        do i = 1, 5
          theta = dtheta*(array_origin(1)+i-1+step)
          points(1:2,i,j,1) = [radius*cos(theta),radius*sin(theta)]
        end do
      end do
      call file%write_mesh_geometry(points, array_origin, write_bounds)
    end if
    do j = 1, 4
      do i = 1, 5
        temp(i,j,1) = array_origin(1)+i-1+10*(j-1)+step
      end do
    end do
    call file%write_point_data(temp_var, temp, array_origin, write_bounds)
    call file%write_field_data(cpu_var, 10.0*step)
    call file%finalize_time_step
  end do

  call file%close

#ifdef USE_MPI
  call MPI_Finalize(stat)
#endif

end program
