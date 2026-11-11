program rectilinear_fixed_mesh_example

  use vtkhdf_rectilinear_file_type
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  type(vtkhdf_rectilinear_file) :: file
  type(vtkhdf_point_data_handle) :: temp_var
  type(vtkhdf_field_data_handle) :: cpu_var
  real, allocatable :: temp(:,:,:), diff(:,:,:)
  real, allocatable :: x(:), y(:), z(:)
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

  ! Create the file with support for time-dependent data (mesh is fixed).
#ifdef USE_MPI
  call file%create('rectilinear_fixed_mesh.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=VTKHDF_FIXED_MESH)
  if (stat /= 0) then
    call MPI_Finalize(stat)
    error stop errmsg
  end if
#else
  call file%create('rectilinear_fixed_mesh.vtkhdf', stat, errmsg, mode=VTKHDF_FIXED_MESH)
  if (stat /= 0) error stop errmsg
#endif

  ! Define the global 2D grid with nonuniform spacing in both coordinate
  ! directions. Each rank owns one unique 4-wide slab of cells.
  point_bounds(:,1) = [0,0,0]
  point_bounds(:,2) = [4*nproc,3,0]
  ! Rank zero supplies the complete coordinate vectors; other ranks pass empty arrays.
  if (rank == 0) then
    x = [(real(i)**2,i=0,4*nproc)]
    y = [(2.0*real(j)**2,j=0,3)]
    z = [0.0]
  else
    allocate(x(0),y(0),z(0))
  end if
  call file%write_mesh(point_bounds, x, y, z)
  array_origin = [4*rank,0,0] ! global index of the first local array element

  ! Write time-independent (static) diffusivity variable on cells
  allocate(diff(4,3,1), source=real(rank+1))
  call file%write_cell_data('diffusivity', diff, array_origin)

  ! Register time-dependent temperature on points.
  ! Points on the interface between two slabs are owned by the slab to the left.
  allocate(temp(5,4,1))
  write_bounds(:,1) = array_origin
  write_bounds(:,2) = array_origin + shape(temp) - 1
  write_bounds(1,1) = write_bounds(1,1) + merge(0,1,(rank==0)) ! omit unowned overlap
  temp_var = file%register_temporal_point_data('temperature', mold=0.0)

  ! We'll write the cummulative cpu time as field data.
  cpu_var = file%register_temporal_field_data('cpu_time', mold=0.0)

  do step = 0, 4
    call file%start_time_step(real(step))
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
