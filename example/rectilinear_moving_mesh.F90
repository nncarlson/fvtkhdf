! Moving RectilinearGrid coordinates with independent temporal data reuse.
program rectilinear_moving_mesh_example
  use,intrinsic :: iso_fortran_env, only: real64
  use vtkhdf_rectilinear_file_type
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  type(vtkhdf_rectilinear_file) :: file
  type(vtkhdf_cell_data_handle) :: temperature_var
  type(vtkhdf_field_data_handle) :: energy_var
  real(real64), allocatable :: x(:), y(:), z(:)
  real(real64) :: temperature(4,3,1), time
  integer :: point_bounds(3,2), cell_origin(3), rank, nproc, stat, step, i, j
  character(:), allocatable :: errmsg
#ifdef USE_MPI
  call MPI_Init(stat)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,stat)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,stat)
#else
  rank = 0
  nproc = 1
#endif
  point_bounds(:,1) = 0
  point_bounds(:,2) = [4*nproc,3,0]
#ifdef USE_MPI
  call file%create('rectilinear_moving_mesh.vtkhdf',MPI_COMM_WORLD,stat,errmsg,mode=VTKHDF_MOVING_MESH)
#else
  call file%create('rectilinear_moving_mesh.vtkhdf',stat,errmsg,mode=VTKHDF_MOVING_MESH)
#endif
  if (stat /= 0) error stop errmsg
  call file%write_mesh(point_bounds)
  temperature_var = file%register_temporal_cell_data('temperature',0.0_real64,attribute='Scalars')
  energy_var = file%register_temporal_field_data('energy',0.0_real64)
  if (rank == 0) then
    allocate(x(4*nproc+1),y(4),z(1))
  else
    allocate(x(0),y(0),z(0))
  end if
  cell_origin = [4*rank,0,0]
  do step = 0,4
    time = real(step,real64)
    call file%start_time_step(time)
    if (modulo(step,2) == 0) then
      if (rank == 0) then
        x = [(0.1_real64*i*(1+0.1_real64*time),i=0,4*nproc)]
        y = [(0.1_real64*j,j=0,3)]
        z = [0.0_real64]
      end if
      call file%write_mesh_geometry(x,y,z)
    end if
    ! Geometry and data may reuse different preceding snapshots.
    if (step == 0 .or. step == 3) then
      do j = 1,3
        do i = 1,4
          temperature(i,j,1) = cell_origin(1)+i-1+10*(j-1)+time
        end do
      end do
      call file%write_cell_data(temperature_var,temperature,array_origin=cell_origin)
    end if
    call file%write_field_data(energy_var,time)
    call file%finalize_time_step
  end do
  call file%close
#ifdef USE_MPI
  call MPI_Finalize(stat)
#endif
end program
