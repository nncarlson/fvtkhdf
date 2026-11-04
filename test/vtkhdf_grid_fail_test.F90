! Each invocation violates one geometry contract. Reaching STOP is a test failure.
program vtkhdf_grid_fail_test
  use,intrinsic :: iso_fortran_env
  use vtkhdf_rectilinear_file_type
  use vtkhdf_structured_file_type
  use vtkhdf_test_env
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  type(vtkhdf_rectilinear_file) :: rect
  type(vtkhdf_structured_file) :: grid
  real(real64) :: points(3,2,2,1), bad(2,2,2,1), axis(2)
  integer :: rank, nproc, stat, bounds(3,2), selected(3,2)
  character(64) :: scenario
  character(:), allocatable :: errmsg
  call test_env_init(rank,nproc)
  call get_command_argument(1,scenario)
  bounds = reshape([10,20,30,11,21,30],[3,2])
  points = 1
  bad = 1
  axis = [0,1]
  if (scenario == 'axis_length') then
#ifdef USE_MPI
    call rect%create('grid_fail_axis.vtkhdf',MPI_COMM_WORLD,stat,errmsg)
#else
    call rect%create('grid_fail_axis.vtkhdf',stat,errmsg)
#endif
    if (stat /= 0) error stop errmsg
    call rect%write_mesh(bounds,axis(:1),axis,axis(:1))
    ! Do not let subsequent checks mask failure to reject the bad length.
  else
#ifdef USE_MPI
    call grid%create('grid_fail_'//trim(scenario)//'.vtkhdf',MPI_COMM_WORLD,stat,errmsg,mode=VTKHDF_FIXED_MESH)
#else
    call grid%create('grid_fail_'//trim(scenario)//'.vtkhdf',stat,errmsg,mode=VTKHDF_FIXED_MESH)
#endif
    if (stat /= 0) error stop errmsg
    if (scenario == 'before_mesh') then
      call grid%write_mesh_geometry(points,bounds(:,1))
    else
      call grid%write_mesh(bounds)
      select case (trim(scenario))
      case ('missing_close')
        call grid%close()
      case ('missing_step')
        call grid%start_time_step(0.0)
      case ('components')
        call grid%write_mesh_geometry(bad,bounds(:,1))
      case ('bounds')
        call grid%write_mesh_geometry(points,bounds(:,1)+[1,0,0])
      case ('source')
        selected = bounds
        call grid%write_mesh_geometry(points,bounds(:,1)+[1,0,0],selected)
      case default
        selected = bounds
        if (rank /= 0) selected(1,2) = selected(1,1)-1
        call grid%write_mesh_geometry(points,bounds(:,1),selected)
        select case (trim(scenario))
        case ('kind')
          call grid%write_mesh_geometry(real(points,real32),bounds(:,1),selected)
        case ('late')
          call grid%start_time_step(0.0)
          call grid%write_mesh_geometry(points,bounds(:,1),selected)
        case ('duplicate')
          call grid%write_mesh(bounds)
        case default
          error stop 'unknown scenario'
        end select
      end select
    end if
  end if
  call test_env_finalize()
end program
