! Return successfully if the expected contract check fails to reject the call.
program vtkhdf_unstructured_constant_data_fail_test
  use,intrinsic :: iso_fortran_env, only: int8, real64
  use vtkhdf_unstructured_file_type
  use vtkhdf_test_env
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  type(vtkhdf_unstructured_file) :: file
  type(vtkhdf_point_data_handle) :: point
  type(vtkhdf_cell_data_handle) :: cell
  integer :: rank, nproc, stat, mode, npoint, ncell
  character(32) :: scenario
  character(:), allocatable :: errmsg

  call test_env_init(rank,nproc)
  call get_command_argument(1,scenario)
  mode = VTKHDF_MOVING_MESH
  if (index(scenario,'dynamic') == 1) mode = VTKHDF_DYNAMIC_MESH
#ifdef USE_MPI
  call file%create('constant_data_fail_'//trim(scenario)//'.vtkhdf',MPI_COMM_WORLD,stat,errmsg,mode=mode)
#else
  call file%create('constant_data_fail_'//trim(scenario)//'.vtkhdf',stat,errmsg,mode=mode)
#endif
  if (stat /= 0) error stop errmsg
  call file%write_mesh(reshape([0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,1.],[3,4]), &
      [1,2,3,4],[1,5],[10_int8])
  npoint = 4
  ncell = 1
  select case (trim(scenario))
  case ('dynamic_point')
    call file%write_point_data('constant',[1.,2.,3.,4.])
  case ('dynamic_cell')
    call file%write_cell_data('constant',[1.])
  case ('wrong_point_count')
    if (rank == 0) npoint = 3
    call file%write_point_data('constant',spread(1.,dim=1,ncopies=npoint))
  case ('wrong_cell_count')
    if (rank == 0) ncell = 0
    call file%write_cell_data('constant',spread(1.,dim=1,ncopies=ncell))
  case ('missing_point')
    point = file%register_temporal_point_data('required',0.0_real64)
    call file%write_point_data('constant',[1.,2.,3.,4.])
    call file%start_time_step(0.)
    call file%finalize_time_step()
  case ('missing_cell')
    cell = file%register_temporal_cell_data('required',0.0_real64)
    call file%write_cell_data('constant',[1.])
    call file%start_time_step(0.)
    call file%finalize_time_step()
  case default
    error stop 'unknown scenario'
  end select
  call test_env_finalize()
end program
