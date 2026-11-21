! Shared Assembly API. For MB output, substitute its module and file type.
program composite_assembly
  use,intrinsic :: iso_fortran_env, only: int8
  use vtkhdf_pdc_file_type
#ifdef USE_MPI_F08
  use mpi_f08
#elif defined(USE_MPI)
  use mpi
#endif
  implicit none

  type(vtkhdf_pdc_file) :: file
  type(vtkhdf_unstructured_block), pointer :: fluid
  type(vtkhdf_image_block), pointer :: solid
  type(vtkhdf_cell_data_handle) :: pressure, temperature
  integer :: rank, nproc, stat, step, bounds(3,2), selected(3,2), regions, materials
  real :: points(3,4), values(1,1,1), time
  character(:), allocatable :: errmsg

#ifdef USE_MPI
  call MPI_Init(stat)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,stat)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,stat)
  call file%create('composite_assembly.vtkhdf',MPI_COMM_WORLD,stat,errmsg)
#else
  rank = 0
  nproc = 1
  call file%create('composite_assembly.vtkhdf',stat,errmsg)
#endif
  if (stat /= 0) error stop errmsg

  ! Each rank contributes one quad to the fluid and one cell to the solid.
  ! Use a homogeneous timeline to avoid the mixed-timeline reader limitation.
  call file%add_block('fluid',fluid,mode=VTKHDF_FIXED_MESH)
  points = reshape([0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.],[3,4])
  points(1,:) = points(1,:) + rank
  call fluid%write_mesh(points,[1,2,3,4],[1,5],[9_int8])
  pressure = fluid%register_temporal_cell_data('pressure',0.0)

  regions = file%add_node('Regions')
  call file%add_block_ref(fluid,parent=regions)

  call file%add_block('solid',solid,mode=VTKHDF_FIXED_MESH)
  bounds = reshape([0,0,0,nproc,1,0],[3,2])
  call solid%write_mesh(bounds,[0.,2.,0.],[1.,1.,1.])
  temperature = solid%register_temporal_cell_data('temperature',0.0)
  selected = reshape([rank,0,0,rank,0,0],[3,2])

  call file%add_block_ref(solid,parent=regions)
  materials = file%add_node('Materials')
  call file%add_block_ref(fluid,parent=materials)
  call file%set_assembly()
  ! Both pointers remain valid; time stepping and data output continue normally.

  do step = 0,2
    time = real(step)
    call file%start_time_step(time)
    call fluid%write_cell_data(pressure,[100.0 + rank + time])
    values = 300.0 + rank + time
    call solid%write_cell_data(temperature,values,[rank,0,0],selected)
    call file%finalize_time_step()
    call file%flush()
  end do
  call file%close()
  nullify(fluid,solid)
#ifdef USE_MPI
  call MPI_Finalize(stat)
#endif
end program composite_assembly
