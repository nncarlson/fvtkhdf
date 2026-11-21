! Return successfully if the expected programming-error assertion did not fire.
program vtkhdf_assembly_fail_test
  use,intrinsic :: iso_fortran_env, only: int8
#ifdef PDC_COLLECTION
  use vtkhdf_pdc_file_type, composite_file => vtkhdf_pdc_file
#else
  use vtkhdf_mb_file_type, composite_file => vtkhdf_mb_file
#endif
  use vtkhdf_unstructured_file_type, only: vtkhdf_unstructured_file
  use vtkhdf_test_env
#ifdef USE_MPI
  use mpi, only: MPI_COMM_WORLD
#endif
  implicit none
  type(composite_file) :: file, foreign
  type(vtkhdf_unstructured_file) :: standalone
  type(vtkhdf_unstructured_block), pointer :: a, b, outsider
  character(:), allocatable :: errmsg
  character(32) :: action, filename
  integer :: rank, nproc, stat, node, parent
  real :: points(3,3)

  call test_env_init(rank,nproc)
  call get_command_argument(1,action)
  filename = trim(action)//'.vtkhdf'
#ifdef USE_MPI
  call file%create(trim(filename),MPI_COMM_WORLD,stat,errmsg)
#else
  call file%create(trim(filename),stat,errmsg)
#endif
  if (stat /= 0) error stop errmsg
  call file%add_block('A',a)
  call file%add_block('B',b)
  points = reshape([0.,0.,0.,1.,0.,0.,0.,1.,0.],[3,3])
  call a%write_mesh(points,[1,2,3],[1,4],[5_int8])
  call b%write_mesh(points,[1,2,3],[1,4],[5_int8])
  select case (trim(action))
  case ('negative_parent')
    node = file%add_node('bad',parent=-1)
  case ('future_parent')
    call file%add_block_ref(a,parent=1)
  case ('foreign')
#ifdef USE_MPI
    call foreign%create('foreign_owner.vtkhdf',MPI_COMM_WORLD,stat,errmsg)
#else
    call foreign%create('foreign_owner.vtkhdf',stat,errmsg)
#endif
    if (stat /= 0) error stop errmsg
    call foreign%add_block('A',outsider)
    call file%add_block_ref(outsider)
  case ('standalone')
#ifdef USE_MPI
    call standalone%create('standalone_writer.vtkhdf',MPI_COMM_WORLD,stat,errmsg)
#else
    call standalone%create('standalone_writer.vtkhdf',stat,errmsg)
#endif
    if (stat /= 0) error stop errmsg
    call file%add_block_ref(standalone)
  case ('closed_node','closed_ref','closed_commit')
    call file%close()
    select case (trim(action))
    case ('closed_node')
      node = file%add_node('bad')
    case ('closed_ref')
      ! A valid foreign object avoids passing a dangling borrowed pointer.
      call file%add_block_ref(standalone)
    case ('closed_commit')
      call file%set_assembly()
    end select
  case ('mismatch_name')
    if (rank == 0) then
      node = file%add_node('first')
    else
      node = file%add_node('different name')
    end if
    call file%set_assembly()
  case ('mismatch_parent')
    node = file%add_node('first')
    parent = file%add_node('second')
    if (rank == 0) parent = node
    call file%add_block_ref(a,parent=parent)
    call file%set_assembly()
  case ('mismatch_reference')
    if (rank == 0) then
      call file%add_block_ref(a)
    else
      call file%add_block_ref(b)
    end if
    call file%set_assembly()
  case ('mismatch_custom')
    if (rank == 0) call file%add_block_ref(a)
    call file%close()
  case ('mismatch_count')
    node = file%add_node('first')
    if (rank == 0) node = file%add_node('extra')
    call file%set_assembly()
  case default
    if (trim(action) == 'duplicate_after') call file%add_block_ref(a)
    call file%set_assembly()
    select case (trim(action))
    case ('repeat_commit')
      call file%set_assembly()
    case ('block_after')
      call file%add_block('late',outsider)
    case ('node_after')
      node = file%add_node('late')
    case ('ref_after','duplicate_after')
      call file%add_block_ref(a)
    case default
      error stop 'unknown negative test'
    end select
  end select
  ! No teardown: an unrelated error must not mask a missing assertion.
  call test_env_finalize()
end program
