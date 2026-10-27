! Fixed ImageData geometry with one cell slab per rank and temporal reuse.
program image_fixed_mesh_example
  use,intrinsic :: iso_fortran_env, only: real64
  use vtkhdf_image_file_type
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  type(vtkhdf_image_file) :: file
  type(vtkhdf_cell_data_handle) :: temperature_var
  type(vtkhdf_field_data_handle) :: energy_var
  real(real64) :: temperature(4,3,1)
  integer :: rank, nproc, stat, step, i, j
  integer :: point_bounds(2,3), cell_origin(3), material(4,3,1)
  character(:), allocatable :: errmsg
#ifdef USE_MPI
  call MPI_Init(stat)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,stat)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,stat)
#else
  rank = 0
  nproc = 1
#endif
  point_bounds(1,:) = [0,0,0]
  point_bounds(2,:) = [4*nproc,3,0]
#ifdef USE_MPI
  call file%create('image_fixed_mesh.vtkhdf',MPI_COMM_WORLD,stat,errmsg,mode=VTKHDF_FIXED_MESH)
#else
  call file%create('image_fixed_mesh.vtkhdf',stat,errmsg,mode=VTKHDF_FIXED_MESH)
#endif
  if (stat /= 0) error stop errmsg
  call file%write_mesh(point_bounds,origin=[0.0_real64,0.0_real64,0.0_real64], &
      spacing=[0.1_real64,0.1_real64,1.0_real64])
  temperature_var = file%register_temporal_cell_data('temperature',0.0_real64,attribute='Scalars')
  energy_var = file%register_temporal_field_data('energy',0.0_real64)
  cell_origin = [4*rank,0,0]
  ! A static field needs one collective call when each rank owns one block.
  material = rank+1
  call file%write_cell_data('material',material,array_origin=cell_origin)
  do step = 0, 4
    call file%start_time_step(real(step,real64))
    if (modulo(step,2) == 0) then
      do j = 1,3
        do i = 1,4
          temperature(i,j,1) = cell_origin(1)+i-1+10*(j-1)+step
        end do
      end do
      ! The complete supplied array is this rank's disjoint global cell slab.
      call file%write_cell_data(temperature_var,temperature,array_origin=cell_origin)
    end if
    ! On odd steps every rank omits temperature: reuse the previous snapshot.
    ! FieldData is independent of decomposition; only rank 0's value is used.
    call file%write_field_data(energy_var,real(step,real64))
    call file%finalize_time_step()
  end do
  call file%close()
#ifdef USE_MPI
  call MPI_Finalize(stat)
#endif
end program
