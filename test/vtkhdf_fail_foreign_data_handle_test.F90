! A field handle from another open file is a programming error, even when
! its field index and (for MultiBlock) block index match the target's.
program vtkhdf_fail_foreign_data_handle_test
  use,intrinsic :: iso_fortran_env, only: int8
  use vtkhdf_ug_file_type
  use vtkhdf_mb_file_type, only: vtkhdf_mb_file, vtkhdf_block_handle
  use vtkhdf_test_env
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  type(vtkhdf_ug_file) :: ug, other_ug
  type(vtkhdf_mb_file) :: mb, other_mb
  type(vtkhdf_block_handle) :: blk, other_blk
  type(vtkhdf_cell_data_handle) :: cell
  type(vtkhdf_point_data_handle) :: point
  type(vtkhdf_field_data_handle) :: field
  integer :: rank, nproc, stat
  real :: points(3,4)
  character(16) :: format, association
  character(:), allocatable :: errmsg, prefix
  call test_env_init(rank,nproc)
  call get_command_argument(1,format)
  call get_command_argument(2,association)
  prefix = 'foreign_handle_'//trim(format)//'_'//trim(association)
  points = reshape([0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,1.],[3,4])
  if (format == 'ug') then
#ifdef USE_MPI
    call ug%create(prefix//'.vtkhdf',MPI_COMM_WORLD,stat,errmsg,mode=VTKHDF_FIXED_MESH)
#else
    call ug%create(prefix//'.vtkhdf',stat,errmsg,mode=VTKHDF_FIXED_MESH)
#endif
    if (stat /= 0) error stop errmsg
#ifdef USE_MPI
    call other_ug%create(prefix//'_other.vtkhdf',MPI_COMM_WORLD,stat,errmsg,mode=VTKHDF_FIXED_MESH)
#else
    call other_ug%create(prefix//'_other.vtkhdf',stat,errmsg,mode=VTKHDF_FIXED_MESH)
#endif
    if (stat /= 0) error stop errmsg
    call ug%write_mesh(points,[1,2,3,4],[1,5],[10_int8])
    call other_ug%write_mesh(points,[1,2,3,4],[1,5],[10_int8])
    cell = ug%register_temporal_cell_data('cell',0.)
    point = ug%register_temporal_point_data('point',0.)
    field = ug%register_temporal_field_data('field',0.)
    cell = other_ug%register_temporal_cell_data('cell',0.)
    point = other_ug%register_temporal_point_data('point',0.)
    field = other_ug%register_temporal_field_data('field',0.)
    call ug%start_time_step(0.)
    select case (association)
    case ('cell')
      call ug%write_cell_data(cell,[1.])
    case ('point')
      call ug%write_point_data(point,[1.,1.,1.,1.])
    case ('field')
      call ug%write_field_data(field,1.)
    end select
  else if (format == 'mb') then
#ifdef USE_MPI
    call mb%create(prefix//'.vtkhdf',MPI_COMM_WORLD,stat,errmsg)
#else
    call mb%create(prefix//'.vtkhdf',stat,errmsg)
#endif
    if (stat /= 0) error stop errmsg
#ifdef USE_MPI
    call other_mb%create(prefix//'_other.vtkhdf',MPI_COMM_WORLD,stat,errmsg)
#else
    call other_mb%create(prefix//'_other.vtkhdf',stat,errmsg)
#endif
    if (stat /= 0) error stop errmsg
    blk = mb%add_block('block',mode=VTKHDF_FIXED_MESH)
    other_blk = other_mb%add_block('block',mode=VTKHDF_FIXED_MESH)
    call mb%write_mesh(blk,points,[1,2,3,4],[1,5],[10_int8])
    call other_mb%write_mesh(other_blk,points,[1,2,3,4],[1,5],[10_int8])
    cell = mb%register_temporal_cell_data(blk,'cell',0.)
    point = mb%register_temporal_point_data(blk,'point',0.)
    field = mb%register_temporal_field_data(blk,'field',0.)
    cell = other_mb%register_temporal_cell_data(other_blk,'cell',0.)
    point = other_mb%register_temporal_point_data(other_blk,'point',0.)
    field = other_mb%register_temporal_field_data(other_blk,'field',0.)
    call mb%start_time_step(0.)
    select case (association)
    case ('cell')
      call mb%write_cell_data(blk,cell,[1.])
    case ('point')
      call mb%write_point_data(blk,point,[1.,1.,1.,1.])
    case ('field')
      call mb%write_field_data(blk,field,1.)
    end select
  end if
  ! Success here means that the expected ownership assertion did not fire.
  ! Do not let missing-initial-write checks during close mask that failure.
  call test_env_finalize()
end program
