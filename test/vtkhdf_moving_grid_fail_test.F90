! Each scenario violates a moving-geometry contract; normal exit is failure.
program vtkhdf_moving_grid_fail_test
  use,intrinsic :: iso_fortran_env
  use vtkhdf_rectilinear_file_type
  use vtkhdf_structured_file_type
  use vtkhdf_structured_base_class, only: vtkhdf_structured_base
  use vtkhdf_file_modes, only: VTKHDF_DYNAMIC_MESH
  use vtkhdf_test_env
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  class(vtkhdf_structured_base), allocatable :: file
  type(vtkhdf_cell_data_handle) :: cell_var
  integer :: rank, nproc, stat, mode, bounds(3,2), selected(3,2)
  real(real64) :: points(3,3,2,2), x(3), yz(2), values(2,1,1)
  character(32) :: format, scenario
  character(:), allocatable :: errmsg
  call test_env_init(rank,nproc)
  call get_command_argument(1,format)
  call get_command_argument(2,scenario)
  if (format == 'rectilinear') then
    allocate(vtkhdf_rectilinear_file :: file)
  else
    allocate(vtkhdf_structured_file :: file)
  end if
  mode = VTKHDF_MOVING_MESH
  if (scenario == 'dynamic') mode = VTKHDF_DYNAMIC_MESH
#ifdef USE_MPI
  call file%create('moving_fail_'//trim(format)//'_'//trim(scenario)//'.vtkhdf',MPI_COMM_WORLD,stat,errmsg,mode=mode)
#else
  call file%create('moving_fail_'//trim(format)//'_'//trim(scenario)//'.vtkhdf',stat,errmsg,mode=mode)
#endif
  if (stat /= 0) error stop errmsg
  if (scenario == 'dynamic') then
    call file%close
    call test_env_finalize()
    stop
  end if
  bounds(:,1) = 0
  bounds(:,2) = [2,1,1]
  points = 0
  x = [0,1,2]
  yz = [0,1]
  values = 1
  selected = bounds
  if (rank /= 0) selected(1,2) = -1
  if (scenario == 'before_grid') call geometry()
  select type (file)
  type is (vtkhdf_rectilinear_file)
    call file%write_mesh(bounds)
  type is (vtkhdf_structured_file)
    call file%write_mesh(bounds)
  end select
  select case (trim(scenario))
  case ('missing')
    call file%start_time_step(0.0)
    call file%finalize_time_step
  case ('data_before_geometry')
    cell_var = file%register_temporal_cell_data('cell',0.0_real64)
    call file%start_time_step(0.0)
    call file%write_cell_data(cell_var,values,[0,0,0])
  case ('axis_length')
    select type (file)
    type is (vtkhdf_rectilinear_file)
      call file%write_mesh_geometry(x(:1),yz,yz)
    end select
  case ('after_data')
    cell_var = file%register_temporal_cell_data('cell',0.0_real64)
    call file%start_time_step(0.0)
    call geometry()
    selected(:,2) = [1,0,0]
    if (rank /= 0) selected(1,2) = -1
    call file%write_cell_data(cell_var,values,[0,0,0],selected)
    call geometry()
  case default
    call geometry()
    select case (trim(scenario))
    case ('duplicate')
      call file%start_time_step(0.0)
      call geometry()
    case ('missing_data')
      cell_var = file%register_temporal_cell_data('cell',0.0_real64)
      call file%start_time_step(0.0)
      call file%finalize_time_step
    case ('between_steps')
      call file%start_time_step(0.0)
      call file%finalize_time_step
      call geometry()
    case ('kind')
      call file%start_time_step(0.0)
      call file%start_time_step(1.0)
      call geometry(single=.true.)
    case ('bounds_change')
      call file%start_time_step(0.0)
      select type (file)
      type is (vtkhdf_rectilinear_file)
        call file%write_mesh(bounds+1)
      type is (vtkhdf_structured_file)
        call file%write_mesh(bounds+1)
      end select
    case default
      error stop 'unknown scenario'
    end select
  end select
  call test_env_finalize()
contains
  subroutine geometry(single)
    logical, intent(in), optional :: single
    if (present(single)) then
      if (single) then
        select type (file)
        type is (vtkhdf_rectilinear_file)
          call file%write_mesh_geometry(real(x,real32),real(yz,real32),real(yz,real32))
        type is (vtkhdf_structured_file)
          call file%write_mesh_geometry(real(points,real32),[0,0,0],selected)
        end select
        return
      end if
    end if
    select type (file)
    type is (vtkhdf_rectilinear_file)
      call file%write_mesh_geometry(x,yz,yz)
    type is (vtkhdf_structured_file)
      call file%write_mesh_geometry(points,[0,0,0],selected)
    end select
  end subroutine
end program
