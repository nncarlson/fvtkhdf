! Exercise borrowed block pointers, MB ownership, and the shared UnstructuredGrid timeline.
program vtkhdf_mb_typed_blocks_test
  use,intrinsic :: iso_fortran_env, only: int8, real32, real64
  use vtkhdf_mb_file_type
  use vtkhdf_h5_c_binding, only: hid_t, H5Iis_valid, H5Iget_file_id, H5Fclose
  use vtkhdf_test_env
#ifdef USE_MPI_F08
  use mpi_f08, only: MPI_COMM_WORLD
#elif defined(USE_MPI)
  use mpi, only: MPI_COMM_WORLD
#endif
  implicit none

  type :: block_pointer
    type(vtkhdf_unstructured_block), pointer :: p => null()
  end type
  integer, parameter :: nblocks = 24
  type(vtkhdf_mb_file) :: mb
  type(block_pointer) :: blocks(nblocks)
  type(vtkhdf_cell_data_handle) :: cell_var(3)
  type(vtkhdf_field_data_handle) :: field_var(3)
  type(vtkhdf_unstructured_block), pointer :: alias
  integer(hid_t) :: roots(nblocks), file_id
  integer :: rank, nproc, stat, i, step, ncell
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  real(real64), allocatable :: points(:,:), cell_values(:)
  character(32) :: name
  character(:), allocatable :: errmsg

  call test_env_init(rank,nproc)

  call create_mb('mb_typed_static.vtkhdf')
  if (stat /= 0) error stop errmsg
  do i = 1,nblocks
    write(name,'(a,i0)') 'Block_',i
    ! Omitting MODE exercises the default static block mode.
    call mb%add_block(trim(name),blocks(i)%p)
    if (.not.associated(blocks(i)%p)) error stop 'missing block pointer'
    roots(i) = blocks(i)%p%root()
    if (i == 1) alias => blocks(1)%p
    if (.not.associated(blocks(1)%p,alias)) error stop 'first target moved'
    if (blocks(1)%p%root() /= roots(1)) error stop 'first block changed'
  end do
  ! Write through every saved pointer after all collection reallocations.
  call make_mesh(0,.false.)
  do i = 1,nblocks
    call blocks(i)%p%write_mesh(points,cnode,xcnode,types)
    cell_values = real(100*i+rank,real64)
    call blocks(i)%p%write_cell_data('value',cell_values)
    call blocks(i)%p%write_field_data('nproc',nproc)
  end do
  call alias%write_field_data('alias',42)
  call mb%flush()
  file_id = H5Iget_file_id(roots(1))
  stat = H5Fclose(file_id) ! release the extra reference from H5Iget_file_id
  if (stat < 0) error stop 'failed to release extra file reference'
  call mb%close()
  call mb%close()
  if (H5Iis_valid(file_id) > 0) error stop 'file leaked'
  do i = 1,nblocks
    if (H5Iis_valid(roots(i)) > 0) error stop 'block root leaked'
    nullify(blocks(i)%p)
  end do
  nullify(alias)

  ! A failed create must also leave the owner usable for a fresh collection.
  call create_mb('mb_typed_static.vtkhdf/not-a-directory.vtkhdf')
  if (stat == 0 .or. .not.allocated(errmsg)) error stop 'expected creation failure'
  call mb%close()

  call create_mb('mb_typed_modes.vtkhdf')
  if (stat /= 0) error stop errmsg
  do i = 1,3
    write(name,'(a,i0)') 'Mode_',i
    call mb%add_block(trim(name),blocks(i)%p,mode=i)
    roots(i) = blocks(i)%p%root()
    call make_mesh(0,.false.)
    if (i == VTKHDF_MOVING_MESH) then
      call blocks(i)%p%write_mesh_topology(cnode,xcnode,types)
      call blocks(i)%p%write_mesh_geometry(points)
    else
      call blocks(i)%p%write_mesh(points,cnode,xcnode,types)
    end if
    cell_var(i) = blocks(i)%p%register_temporal_cell_data('value',0.0_real64,attribute='Scalars')
    field_var(i) = blocks(i)%p%register_temporal_field_data('clock',0.0_real64)
    call blocks(i)%p%write_field_data('constant',70+i)
    call blocks(i)%p%write_field_data('nproc',nproc)
  end do
  do step = 0,2
    if (step == 0) then
      call mb%start_time_step(0.0_real32)
    else
      ! Starting step 2 implicitly finalizes the preceding open step.
      call mb%start_time_step(real(step,real64))
    end if
    do i = 1,3
      if (.not.blocks(i)%p%step_is_open()) error stop 'block step not opened'
      call make_mesh(min(step,1),i == VTKHDF_DYNAMIC_MESH .and. step > 0)
      if (step == 1) then
        if (i == VTKHDF_MOVING_MESH) call blocks(i)%p%write_mesh_geometry(points)
        if (i == VTKHDF_DYNAMIC_MESH) call blocks(i)%p%write_mesh(points,cnode,xcnode,types)
      end if
      if (step == 0 .or. (step == 1 .and. i /= VTKHDF_FIXED_MESH)) then
        cell_values = real(100*i+10*step+rank,real64)
        call blocks(i)%p%write_cell_data(cell_var(i),cell_values)
      end if
      if (step /= 1) call blocks(i)%p%write_field_data(field_var(i),real(10*step+i+1000*rank,real64))
    end do
    if (step == 0) then
      call mb%finalize_time_step()
      do i = 1,3
        if (blocks(i)%p%step_is_open()) error stop 'block step not finalized'
      end do
    end if
  end do
  do i = 1,3
    call blocks(i)%p%write_field_data('late',90+i)
  end do
  call mb%flush()
  ! Closing MB finalizes step 2 on every temporal block before deallocation.
  call mb%close()
  call mb%close()
  do i = 1,3
    if (H5Iis_valid(roots(i)) > 0) error stop 'temporal block root leaked'
    nullify(blocks(i)%p)
  end do
  deallocate(points,cell_values,cnode,xcnode,types)
  call test_env_finalize()

contains

  subroutine create_mb(filename)
    character(*), intent(in) :: filename
#ifdef USE_MPI
    call mb%create(filename,MPI_COMM_WORLD,stat,errmsg)
#else
    call mb%create(filename,stat,errmsg)
#endif
  end subroutine

  subroutine make_mesh(shift, grow)
    integer, intent(in) :: shift
    logical, intent(in) :: grow
    integer :: j
    real(real64), parameter :: tetra(3,4) = reshape([0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,1.],[3,4])
    ncell = rank+1
    if (grow) ncell = ncell+1
    if (nproc > 1 .and. rank == nproc-1) ncell = 0
    if (allocated(points)) deallocate(points)
    allocate(points(3,4*ncell))
    do j = 1,ncell
      points(:,4*j-3:4*j) = tetra+spread([real(10*rank+2*(j-1),real64),0.0_real64,real(shift,real64)],2,4)
    end do
    cnode = [(j,j=1,4*ncell)]
    xcnode = [(4*j+1,j=0,ncell)]
    types = [(10_int8,j=1,ncell)]
    ! Explicit allocation gives empty ranks a zero-extent array descriptor.
    if (allocated(cell_values)) deallocate(cell_values)
    allocate(cell_values(ncell))
  end subroutine

end program vtkhdf_mb_typed_blocks_test
