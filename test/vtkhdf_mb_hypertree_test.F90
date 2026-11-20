! Exercise MB's shared timeline with fixed, moving, and dynamic HTG blocks.
program vtkhdf_mb_hypertree_test
  use,intrinsic :: iso_fortran_env, only: int64, real64
  use vtkhdf_mb_file_type
  use vtkhdf_test_env
#ifdef USE_MPI_F08
  use mpi_f08, only: MPI_COMM_WORLD
#elif defined(USE_MPI)
  use mpi, only: MPI_COMM_WORLD
#endif
  implicit none
  type :: block_pointer
    type(vtkhdf_hypertree_block), pointer :: p => null()
  end type
  type(vtkhdf_mb_file) :: file
  type(block_pointer) :: blocks(3)
  type(vtkhdf_cell_data_handle) :: value(3)
  character(7), parameter :: names(3) = [character(7) :: 'fixed', 'moving', 'dynamic']
  character(:), allocatable :: errmsg
  integer :: rank, nproc, stat, mode, step, phase
  integer, allocatable :: ids(:), xtree(:)
  logical, allocatable :: refined(:), mask(:)
  real(real64), allocatable :: values(:)

  call test_env_init(rank,nproc)
#ifdef USE_MPI
  call file%create('mb_hypertree.vtkhdf',MPI_COMM_WORLD,stat,errmsg,data_chunk_bytes=64_int64)
#else
  call file%create('mb_hypertree.vtkhdf',stat,errmsg,data_chunk_bytes=64_int64)
#endif
  if (stat /= 0) error stop errmsg
  do mode = VTKHDF_FIXED_MESH,VTKHDF_DYNAMIC_MESH
    call file%add_block(trim(names(mode)),blocks(mode)%p,mode=mode)
    call blocks(mode)%p%write_mesh([3,1,1])
    call blocks(mode)%p%write_field_data('id',data=mode)
    value(mode) = blocks(mode)%p%register_temporal_cell_data('value',0.0_real64)
  end do
  do step = 0,2
    ! Starting step 2 finalizes step 1; closing MB finalizes step 2.
    call file%start_time_step(real(step,real64))
    do mode = VTKHDF_FIXED_MESH,VTKHDF_DYNAMIC_MESH
      if (.not.blocks(mode)%p%step_is_open()) error stop 'HTG step not started'
      phase = 0
      if (mode == VTKHDF_DYNAMIC_MESH .and. step > 0) phase = 1
      call local_topology(phase,mode,step)
      if (step == 0 .or. (mode == VTKHDF_DYNAMIC_MESH .and. step == 1)) then
        call blocks(mode)%p%write_mesh_topology(ids,xtree,refined,mask)
      end if
      if (step == 0 .or. (mode >= VTKHDF_MOVING_MESH .and. step == 1)) then
        phase = merge(1,0,step>0)
        call blocks(mode)%p%write_mesh_geometry([0.0_real64,1.0_real64,2.0_real64]+phase, &
            [0.0_real64],[0.0_real64])
      end if
      ! Step 2 reuses the latest cell values, including the dynamic topology.
      if (step < 2) call blocks(mode)%p%write_cell_data(value(mode),values)
    end do
    if (step == 0) then
      call file%finalize_time_step()
      do mode = VTKHDF_FIXED_MESH,VTKHDF_DYNAMIC_MESH
        if (blocks(mode)%p%step_is_open()) error stop 'HTG step not finalized'
      end do
    end if
  end do
  call file%flush()
  call file%close()
  call file%close()
  call test_env_finalize()

contains

  subroutine local_topology(phase, mode, step)
    integer, intent(in) :: phase, mode, step
    integer :: g, j, n, nactive
    ids = [integer::]
    xtree = [1]
    refined = [logical::]
    mask = [logical::]
    values = [real(real64)::]
    nactive = nproc
    if (phase == 1) nactive = max(nproc-1,1)
    do g = 0,1
      if (modulo(g+phase,nactive) /= rank) cycle
      ids = [ids,g+1]
      n = 1
      if (phase == 1 .and. g == 0) n = 3
      xtree = [xtree,xtree(size(xtree))+n]
      refined = [refined,(j==1 .and. n==3,j=1,n)]
      mask = [mask,(j==3,j=1,n)]
      values = [values,(real(1000*mode+100*g+10*step+j-1,real64),j=1,n)]
    end do
  end subroutine

end program vtkhdf_mb_hypertree_test
