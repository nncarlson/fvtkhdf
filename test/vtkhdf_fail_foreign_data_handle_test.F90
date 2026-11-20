! A field handle from another open file is a programming error, even when
! its field index and (for MultiBlock) block position match the target's.
program vtkhdf_fail_foreign_data_handle_test
  use,intrinsic :: iso_fortran_env, only: int8
  use vtkhdf_unstructured_file_type
  use vtkhdf_mb_file_type, only: vtkhdf_mb_file, vtkhdf_unstructured_block
  use vtkhdf_test_env
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  type(vtkhdf_unstructured_file) :: unstructured, other_unstructured
  type(vtkhdf_mb_file) :: mb, other_mb
  type(vtkhdf_unstructured_block), pointer :: blk => null(), other_blk => null()
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
  if (format == 'unstructured') then
#ifdef USE_MPI
    call unstructured%create(prefix//'.vtkhdf',MPI_COMM_WORLD,stat,errmsg,mode=VTKHDF_FIXED_MESH)
#else
    call unstructured%create(prefix//'.vtkhdf',stat,errmsg,mode=VTKHDF_FIXED_MESH)
#endif
    if (stat /= 0) error stop errmsg
#ifdef USE_MPI
    call other_unstructured%create(prefix//'_other.vtkhdf',MPI_COMM_WORLD,stat,errmsg,mode=VTKHDF_FIXED_MESH)
#else
    call other_unstructured%create(prefix//'_other.vtkhdf',stat,errmsg,mode=VTKHDF_FIXED_MESH)
#endif
    if (stat /= 0) error stop errmsg
    call unstructured%write_mesh(points,[1,2,3,4],[1,5],[10_int8])
    call other_unstructured%write_mesh(points,[1,2,3,4],[1,5],[10_int8])
    cell = unstructured%register_temporal_cell_data('cell',0.)
    point = unstructured%register_temporal_point_data('point',0.)
    field = unstructured%register_temporal_field_data('field',0.)
    cell = other_unstructured%register_temporal_cell_data('cell',0.)
    point = other_unstructured%register_temporal_point_data('point',0.)
    field = other_unstructured%register_temporal_field_data('field',0.)
    call unstructured%start_time_step(0.)
    select case (association)
    case ('cell')
      call unstructured%write_cell_data(cell,[1.])
    case ('point')
      call unstructured%write_point_data(point,[1.,1.,1.,1.])
    case ('field')
      call unstructured%write_field_data(field,1.)
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
    call mb%add_block('block', blk,mode=VTKHDF_FIXED_MESH)
    call other_mb%add_block('block', other_blk,mode=VTKHDF_FIXED_MESH)
    call blk%write_mesh(points,[1,2,3,4],[1,5],[10_int8])
    call other_blk%write_mesh(points,[1,2,3,4],[1,5],[10_int8])
    cell = blk%register_temporal_cell_data('cell',0.)
    point = blk%register_temporal_point_data('point',0.)
    field = blk%register_temporal_field_data('field',0.)
    cell = other_blk%register_temporal_cell_data('cell',0.)
    point = other_blk%register_temporal_point_data('point',0.)
    field = other_blk%register_temporal_field_data('field',0.)
    call mb%start_time_step(0.)
    select case (association)
    case ('cell')
      call blk%write_cell_data(cell,[1.])
    case ('point')
      call blk%write_point_data(point,[1.,1.,1.,1.])
    case ('field')
      call blk%write_field_data(field,1.)
    end select
  end if
  ! Success here means that the expected ownership assertion did not fire.
  ! Do not let missing-initial-write checks during close mask that failure.
  call test_env_finalize()
end program
