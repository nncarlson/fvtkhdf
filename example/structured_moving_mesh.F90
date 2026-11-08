! Moving StructuredGrid point coordinates with fixed logical connectivity.
program structured_moving_mesh_example
  use,intrinsic :: iso_fortran_env, only: real64
  use vtkhdf_structured_file_type
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  type(vtkhdf_structured_file) :: file
  type(vtkhdf_cell_data_handle) :: temperature_var
  type(vtkhdf_field_data_handle) :: energy_var
  real(real64), allocatable :: points(:,:,:,:)
  real(real64) :: temperature(4,3,1), time
  integer :: point_bounds(3,2), cell_origin(3), first, last, rank, nproc, stat, step, i, j
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
  call file%create('structured_moving_mesh.vtkhdf',MPI_COMM_WORLD,stat,errmsg,mode=VTKHDF_MOVING_MESH)
#else
  call file%create('structured_moving_mesh.vtkhdf',stat,errmsg,mode=VTKHDF_MOVING_MESH)
#endif
  if (stat /= 0) error stop errmsg
  call file%write_mesh(point_bounds)
  temperature_var = file%register_temporal_cell_data('temperature',0.0_real64,attribute='Scalars')
  energy_var = file%register_temporal_field_data('energy',0.0_real64)
  ! Every shared slab-boundary point belongs to the rank on its right.
  first = 4*rank
  last = first+3
  if (rank == nproc-1) last = last+1
  allocate(points(3,last-first+1,4,1))
  cell_origin = [4*rank,0,0]
  do step = 0,4
    time = real(step,real64)
    call file%start_time_step(time)
    if (modulo(step,2) == 0) then
      do j = 0,3
        do i = first,last
          points(:,i-first+1,j+1,1) = [0.1_real64*i,0.1_real64*j,0.01_real64*i*j*time]
        end do
      end do
      call file%write_mesh_geometry(points,array_origin=[first,0,0])
    end if
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
