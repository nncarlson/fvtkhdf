! Check pair contracts at the intended guard, rather than any failing setup.
program vtkhdf_hypertree_interface_fail_test
  use,intrinsic :: iso_fortran_env, only: real32, real64, error_unit
  use vtkhdf_hypertree_file_type
  use vtkhdf_test_env
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  type(vtkhdf_hypertree_file) :: file
  character(64) :: scenario, diagnostic_file
  character(:), allocatable :: errmsg
  integer :: rank, nproc, stat, mode
  real(real64) :: normals(3,6), intercepts(3,6)
  call test_env_init(rank,nproc)
  ! MPI_Abort can terminate the launcher before it forwards stderr. Persist
  ! each rank's assertion diagnostic for the parent checker to inspect.
  write(diagnostic_file,'(a,i0,a)') 'interface-assertion-R',rank,'.log'
  open(error_unit,file=trim(diagnostic_file),status='replace')
  call get_command_argument(1,scenario)
  normals = 1
  intercepts = 2
  mode = VTKHDF_FIXED_MESH
  if (scenario == 'static_register') mode = VTKHDF_STATIC
  if (scenario == 'dynamic_static' .or. scenario == 'changed_count' .or. &
      scenario == 'temporal_barrier') mode = VTKHDF_DYNAMIC_MESH
  if (scenario == 'static_barrier') mode = VTKHDF_MOVING_MESH
  if (scenario == 'closed_write') then
    call file%write_interface_data(normals,intercepts)
  else if (scenario == 'closed_register') then
    call file%register_temporal_interface_data(0.0_real64)
  end if
#ifdef USE_MPI
  call file%create('interface_fail_'//trim(scenario)//'.vtkhdf',MPI_COMM_WORLD,stat,errmsg,mode=mode)
#else
  call file%create('interface_fail_'//trim(scenario)//'.vtkhdf',stat,errmsg,mode=mode)
#endif
  if (stat /= 0) error stop errmsg
  if (scenario == 'static_register') then
    call file%register_temporal_interface_data(0.0_real64)
  else if (scenario == 'kind_agreement_register') then
    if (rank == 0) then
      call file%register_temporal_interface_data(0.0_real32)
    else
      call file%register_temporal_interface_data(0.0_real64)
    end if
  end if
  if (scenario == 'no_grid') call file%write_interface_data(normals,intercepts)
  call file%write_mesh([3,2,1])
  if (scenario /= 'no_geometry') call file%write_mesh_geometry([0.,1.,2.],[0.,1.],[0.])
  if (scenario /= 'no_topology') then
    if (scenario == 'missing_initial_empty' .and. rank == 0 .and. nproc > 1) then
      call file%write_mesh_topology([integer::],[1],[logical::])
    else
      call topology()
    end if
  end if
  select case (scenario)
  case ('no_geometry','no_topology','dynamic_static')
    call file%write_interface_data(normals,intercepts)
  case ('static_twice','register_after_static','static_barrier')
    if (scenario == 'static_barrier') call file%start_time_step(0.)
    call file%write_interface_data(normals,intercepts)
    if (scenario == 'static_twice') call file%write_interface_data(normals,intercepts)
    if (scenario == 'register_after_static') call file%register_temporal_interface_data(0.0_real64)
    if (scenario == 'static_barrier') call file%write_mesh_geometry([0.,1.,2.],[0.,1.],[0.])
  case ('late_register')
    call file%start_time_step(0.)
    call file%register_temporal_interface_data(0.0_real64)
  case ('normals_components')
    call file%write_interface_data(normals(:2,:),intercepts)
  case ('intercepts_components')
    call file%write_interface_data(normals,intercepts(:2,:))
  case ('normals_count')
    call file%write_interface_data(normals(:,:5),intercepts)
  case ('intercepts_count')
    call file%write_interface_data(normals,intercepts(:,:5))
  case ('kind_agreement_write')
    if (rank == 0) then
      call file%write_interface_data(real(normals,real32),real(intercepts,real32))
    else
      call file%write_interface_data(normals,intercepts)
    end if
  case ('operation_agreement')
    if (rank == 0) then
      call file%write_interface_data(normals,intercepts)
    else
      call file%register_temporal_interface_data(0.0_real64)
    end if
  case default
    if (scenario == 'kind_64_to_32') then
      call file%register_temporal_interface_data(0.0_real32)
    else
      call file%register_temporal_interface_data(0.0_real64)
    end if
    if (scenario == 'register_twice') call file%register_temporal_interface_data(0.0_real64)
    if (scenario == 'outside_step') call file%write_interface_data(normals,intercepts)
    call file%start_time_step(0.)
    select case (scenario)
    case ('missing_initial','missing_initial_empty')
      call file%finalize_time_step()
    case ('kind_64_to_32')
      call file%write_interface_data(normals,intercepts)
    case ('kind_32_to_64')
      call file%write_interface_data(real(normals,real32),real(intercepts,real32))
    case ('temporal_normals_components')
      call file%write_interface_data(normals(:2,:),intercepts)
    case ('temporal_intercepts_count')
      call file%write_interface_data(normals,intercepts(:,:5))
    case ('duplicate','changed_count','temporal_barrier','between_steps')
      call file%write_interface_data(normals,intercepts)
      if (scenario == 'duplicate') call file%write_interface_data(normals,intercepts)
      if (scenario == 'temporal_barrier') call topology()
      call file%finalize_time_step()
      if (scenario == 'between_steps') call file%write_interface_data(normals,intercepts)
      call file%start_time_step(1.)
      call file%write_mesh_topology([1,2],[1,2,3],[.false.,.false.])
      call file%finalize_time_step()
    case default
      error stop 'unknown scenario'
    end select
  end select
  ! Reaching here means the specific call was incorrectly accepted.
  call test_env_finalize()
contains
  subroutine topology()
    call file%write_mesh_topology([1,2],[1,2,7],[.false.,.true.,.false.,.false.,.false.,.false.])
  end subroutine
end program
