! Expected-failure checks: each invocation violates one public contract.
program vtkhdf_image_fail_test
  use,intrinsic :: iso_fortran_env
  use vtkhdf_image_file_type
  use vtkhdf_test_env
  use vtkhdf_ug_file_type, only: vtkhdf_ug_file
  use vtkhdf_file_modes, only: VTKHDF_MOVING_MESH, VTKHDF_DYNAMIC_MESH
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  type(vtkhdf_image_file) :: file, other
  type(vtkhdf_ug_file) :: ug
  type(vtkhdf_cell_data_handle) :: var
  type(vtkhdf_field_data_handle) :: field_var
  integer :: rank,nproc,stat,mode
  integer :: bounds(2,3), selected(2,3)
  real(real64) :: a(3,2,1), v(2,3,2,1)
  character(64) :: scenario
  character(:), allocatable :: errmsg
  call test_env_init(rank,nproc)
  call get_command_argument(1,scenario)
  mode = VTKHDF_FIXED_MESH
  select case (trim(scenario))
  case ('mode_moving')
    mode = VTKHDF_MOVING_MESH
  case ('mode_dynamic')
    mode = VTKHDF_DYNAMIC_MESH
  case ('mode_invalid')
    mode = -1
  case ('mode_mismatch')
    if (rank == 0) mode = VTKHDF_STATIC
  end select
#ifdef USE_MPI
  call file%create('image_fail_'//trim(scenario)//'.vtkhdf',MPI_COMM_WORLD,stat,errmsg,mode=mode)
#else
  call file%create('image_fail_'//trim(scenario)//'.vtkhdf',stat,errmsg,mode=mode)
#endif
  if (stat /= 0) error stop errmsg
  if (index(scenario,'mode_') == 1) then
    ! If creation accepted the invalid mode, exit successfully so WILL_FAIL
    ! reports a failure. Later lifecycle checks must not mask this violation.
    call file%close()
    call test_env_finalize()
    stop
  end if
  bounds(1,:) = 0
  bounds(2,:) = [3,2,0]
  a = 1
  v = 1
  if (scenario == 'before_mesh') then
    call file%write_cell_data('bad',a,[0,0,0])
  end if
  if (scenario == 'register_before_mesh') then
    var = file%register_temporal_cell_data('bad',0.0_real64)
  end if
  if (scenario == 'geometry_mismatch') bounds(2,1) = 3+rank
  call file%write_mesh(bounds,[0.0,0.0,0.0],[1.0,1.0,1.0])
  var = file%register_temporal_cell_data('data',0.0_real64)
  field_var = file%register_temporal_field_data('field',0.0_real64)
  if (scenario == 'outside_step') call file%write_cell_data(var,a,[0,0,0])
  if (scenario == 'field_static') call file%write_field_data('static',1.0)
  if (scenario == 'wrong_handle') then
#ifdef USE_MPI
    call other%create('image_fail_other.vtkhdf',MPI_COMM_WORLD,stat,errmsg,mode=VTKHDF_FIXED_MESH)
#else
    call other%create('image_fail_other.vtkhdf',stat,errmsg,mode=VTKHDF_FIXED_MESH)
#endif
    call other%write_mesh(bounds,[0.0,0.0,0.0],[1.0,1.0,1.0])
    var = other%register_temporal_cell_data('data',0.0_real64)
  end if
  call file%start_time_step(0.0)
  select case (trim(scenario))
  case ('bounds')
    selected(1,:) = 0
    selected(2,:) = [2,1,0]
    if (rank == 0) selected(2,1) = 3
    call file%write_cell_data(var,a,[0,0,0],selected)
  case ('source')
    call file%write_cell_data(var,a,[1,0,0],reshape([0,2,0,1,0,0],[2,3]))
  case ('type')
    call file%write_cell_data(var,real(a,real32),[0,0,0])
  case ('components')
    call file%write_cell_data(var,v,[0,0,0])
  case ('late')
    var = file%register_temporal_cell_data('late',0.0)
  case ('field_duplicate')
    call file%write_field_data(field_var,1.0_real64)
    call file%write_field_data(field_var,2.0_real64)
  case ('wrong_handle')
    call file%write_cell_data(var,a,[0,0,0])
  case ('static_type')
    call file%write_cell_data('static',a,[0,0,0],handle=var)
    call file%write_cell_data(var,real(a,real32),[0,0,0])
  case ('static_components')
    call file%write_cell_data('static',a,[0,0,0],handle=var)
    call file%write_cell_data(var,v,[0,0,0])
  case ('name_mismatch')
    block
      character(20) :: name
      write(name,'(a,i0)') 'field',rank
      call file%write_cell_data(trim(name),a,[0,0,0])
    end block
  case ('ug_duplicate')
#ifdef USE_MPI
    call ug%create('image_fail_ug_duplicate.vtkhdf',MPI_COMM_WORLD,stat,errmsg,mode=VTKHDF_FIXED_MESH)
#else
    call ug%create('image_fail_ug_duplicate.vtkhdf',stat,errmsg,mode=VTKHDF_FIXED_MESH)
#endif
    if (stat /= 0) error stop errmsg
    call ug%write_mesh(reshape([0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,1.],[3,4]), &
        [1,2,3,4],[1,5],[10_int8])
    var = ug%register_temporal_cell_data('data',0.0_real64)
    call ug%start_time_step(0.0)
    call ug%write_cell_data(var,[1.0_real64])
    call ug%write_cell_data(var,[2.0_real64])
  case ('missing')
    ! Finalization must reject a registered mesh field with no initial calls.
  end select
  if (scenario == 'missing') call file%finalize_time_step()
  ! Exit successfully if the intended assertion did not fire. Do not let
  ! missing-initial-data checks during cleanup mask an unrelated failure.
  call test_env_finalize()
end program
