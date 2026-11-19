! Exercise the shared lifecycle with heterogeneous borrowed leaves. This is
! deliberately below MB's public API, testing the dataset lifecycle directly.
program vtkhdf_simple_dataset_test
  use,intrinsic :: iso_fortran_env
  use,intrinsic :: iso_c_binding, only: c_funptr, c_ptr, c_associated
  use vtkhdf_simple_dataset_class
  use vtkhdf_structured_base_class, only: vtkhdf_structured_base
  use vtkhdf_unstructured_file_type
  use vtkhdf_unstructured_writer_type
  use vtkhdf_image_file_type
  use vtkhdf_image_writer_type
  use vtkhdf_rectilinear_file_type
  use vtkhdf_rectilinear_writer_type
  use vtkhdf_structured_file_type
  use vtkhdf_structured_writer_type
  use vtkhdf_hypertree_file_type
  use vtkhdf_hypertree_writer_type
  use vtkhdf_context_type
  use vtkhdf_h5_c_binding
  use vtkhdf_test_env
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  type :: leaf
    class(vtkhdf_simple_dataset), allocatable :: writer
  end type
  type(leaf) :: blocks(5)
  type(vtkhdf_field_data_handle) :: field_var
  type(vtkhdf_context) :: ctx
  type(c_funptr) :: saved_error_handler
  type(c_ptr) :: saved_error_data
  integer(hid_t) :: file_id, parent_id, fapl, root_ids(size(blocks))
  integer :: rank, nproc, stat, j, cycle
  integer(int32) :: bounds(3,2), selected(3,2)
  character(100) :: filename
  character(:), allocatable :: errmsg
  character(24) :: name
  real(real64) :: points(3,2,2,1)

  call test_env_init(rank,nproc)
#ifdef USE_MPI
  call ctx%init(MPI_COMM_WORLD)
#else
  call ctx%init()
#endif
  call init_hdf5
  call ctx%get_fapl(fapl)
  file_id = H5Fcreate('simple_dataset_leaves.vtkhdf',H5F_ACC_TRUNC,H5P_DEFAULT,fapl)
  if (file_id < 0) error stop 'failed to create parent file'
  stat = H5Pclose(fapl)
  parent_id = H5Gcreate(file_id,'container')
  if (parent_id < 0) error stop 'failed to create parent group'
  bounds = reshape([3,4,5,4,5,5],[3,2])
  selected = bounds
  if (rank /= 0) selected(1,2) = selected(1,1)-1
  points = reshape([0.,0.,0.,1.,0.,0.,0.,1.,0.,1.,1.,0.],[3,2,2,1])

  allocate(vtkhdf_unstructured_block :: blocks(1)%writer)
  allocate(vtkhdf_image_block :: blocks(2)%writer)
  allocate(vtkhdf_rectilinear_block :: blocks(3)%writer)
  allocate(vtkhdf_structured_block :: blocks(4)%writer)
  allocate(vtkhdf_hypertree_block :: blocks(5)%writer)
  do cycle = 1,2
    ! The second cycle reuses closed polymorphic objects and the same parent.
    do j = 1,size(blocks)
      write(name,'(a,i0,a,i0)') 'leaf_',j,'_',cycle
      select type (writer => blocks(j)%writer)
      type is (vtkhdf_unstructured_block)
        call writer%init(parent_id,trim(name),ctx,mode=VTKHDF_FIXED_MESH)
        call register_fields(writer)
        if (cycle == 2) call writer%start_time_step(0.0_real32)
        if (rank == 0) then
          call writer%write_mesh(reshape(points,[3,4]),[1,2,4,3],[1,5],[9_int8])
        else
          call writer%write_mesh(reshape([real(real64)::],[3,0]),[integer::],[1],[integer(int8)::])
        end if
      type is (vtkhdf_image_block)
        call writer%init(parent_id,trim(name),ctx,mode=VTKHDF_FIXED_MESH)
        call register_fields(writer)
        if (cycle == 2) call writer%start_time_step(0.0_real32)
        call writer%write_mesh(bounds,[0.,0.,0.],[1.,1.,1.])
      type is (vtkhdf_rectilinear_block)
        call writer%init(parent_id,trim(name),ctx,mode=merge(VTKHDF_MOVING_MESH,VTKHDF_FIXED_MESH,cycle==2))
        call register_fields(writer)
        if (cycle == 2) call writer%start_time_step(0.0_real32)
        call writer%write_mesh(bounds,[0.,1.],[0.,1.],[0.])
      type is (vtkhdf_structured_block)
        call writer%init(parent_id,trim(name),ctx,mode=merge(VTKHDF_MOVING_MESH,VTKHDF_FIXED_MESH,cycle==2))
        call register_fields(writer)
        if (cycle == 2) call writer%start_time_step(0.0_real32)
        call writer%write_mesh(bounds)
        call writer%write_mesh_geometry(points,bounds(:,1),selected)
      type is (vtkhdf_hypertree_block)
        call writer%init(parent_id,trim(name),ctx,mode=merge(VTKHDF_MOVING_MESH,VTKHDF_FIXED_MESH,cycle==2))
        call register_fields(writer)
        if (cycle == 2) call writer%start_time_step(0.0_real32)
        call write_hypertree_mesh(writer)
      class default
        error stop 'unexpected leaf type'
      end select
      ! Lifecycle operations dispatch through the abstract ancestor.
      root_ids(j) = blocks(j)%writer%root()
      if (root_ids(j) < 0) error stop 'missing root'
      if (.not.blocks(j)%writer%temporal_enabled()) error stop 'missing temporal mode'
      if (cycle == 1) then
        if (blocks(j)%writer%step_is_open()) error stop 'step already open'
        call blocks(j)%writer%start_time_step(0.0_real32)
      end if
      if (.not.blocks(j)%writer%step_is_open()) error stop 'step did not open'
      call blocks(j)%writer%start_time_step(1.0_real64)
      select type (writer => blocks(j)%writer)
      class is (vtkhdf_unstructured_writer)
        call writer%write_field_data(field_var,real(cycle,real64))
      class is (vtkhdf_structured_base)
        call writer%write_field_data(field_var,real(cycle,real64))
      class is (vtkhdf_hypertree_writer)
        if (cycle == 2) call writer%write_mesh_geometry([1.0,2.0,3.0],[0.0],[0.0])
        call writer%write_field_data(field_var,real(cycle,real64))
      end select
      if (cycle == 1) then
        call blocks(j)%writer%finalize_time_step()
        if (blocks(j)%writer%step_is_open()) error stop 'step did not close'
      end if
    end do
    do j = 1,size(blocks)
      call blocks(j)%writer%close()
      call blocks(j)%writer%close()
      if (blocks(j)%writer%root() >= 0) error stop 'root not reset'
      if (blocks(j)%writer%temporal_enabled()) error stop 'mode not reset'
      if (H5Iis_valid(root_ids(j)) > 0) error stop 'root leaked'
    end do
    if (H5Iis_valid(file_id) <= 0) error stop 'borrowed file closed'
    if (H5Iis_valid(parent_id) <= 0) error stop 'borrowed group closed'
    if (.not.ctx%global_all(.true.)) error stop 'borrowed communicator invalid'
  end do
  ! File and block writers are siblings; replace the borrowed block objects
  ! with standalone file objects for the owned-resource checks below.
  do j = 1,size(blocks)
    deallocate(blocks(j)%writer)
  end do
  allocate(vtkhdf_unstructured_file :: blocks(1)%writer)
  allocate(vtkhdf_image_file :: blocks(2)%writer)
  allocate(vtkhdf_rectilinear_file :: blocks(3)%writer)
  allocate(vtkhdf_structured_file :: blocks(4)%writer)
  allocate(vtkhdf_hypertree_file :: blocks(5)%writer)
  ! Failed standalone creation must release its duplicated communicator.
  ! The same polymorphic object can then be closed and reopened normally.
  stat = H5Eget_auto2(H5E_DEFAULT,saved_error_handler,saved_error_data)
  if (stat < 0) error stop 'failed to query HDF5 error handler'
  do j = 1,size(blocks)
    do cycle = 1,3
      filename = 'simple_dataset_leaves.vtkhdf/not-a-directory.vtkhdf'
      if (cycle > 1) write(filename,'(a,i0,a)') 'simple_dataset_owned_',j,'.vtkhdf'
      select type (writer => blocks(j)%writer)
      type is (vtkhdf_unstructured_file)
#ifdef USE_MPI
        call writer%create(trim(filename),MPI_COMM_WORLD,stat,errmsg)
#else
        call writer%create(trim(filename),stat,errmsg)
#endif
        call check_error_handler()
        if (stat == 0) call writer%flush()
      type is (vtkhdf_image_file)
#ifdef USE_MPI
        call writer%create(trim(filename),MPI_COMM_WORLD,stat,errmsg)
#else
        call writer%create(trim(filename),stat,errmsg)
#endif
        call check_error_handler()
        if (stat == 0) call writer%flush()
      type is (vtkhdf_rectilinear_file)
#ifdef USE_MPI
        call writer%create(trim(filename),MPI_COMM_WORLD,stat,errmsg)
#else
        call writer%create(trim(filename),stat,errmsg)
#endif
        call check_error_handler()
        if (stat == 0) call writer%flush()
      type is (vtkhdf_structured_file)
#ifdef USE_MPI
        call writer%create(trim(filename),MPI_COMM_WORLD,stat,errmsg)
#else
        call writer%create(trim(filename),stat,errmsg)
#endif
        call check_error_handler()
        if (stat == 0) then
          call writer%flush()
          ! Cleanup must allow a logical grid whose coordinates were omitted.
          call writer%write_mesh(bounds)
        end if
      type is (vtkhdf_hypertree_file)
#ifdef USE_MPI
        call writer%create(trim(filename),MPI_COMM_WORLD,stat,errmsg)
#else
        call writer%create(trim(filename),stat,errmsg)
#endif
        call check_error_handler()
        if (stat == 0) then
          call writer%flush()
          ! Reopening must reset the inherited root grid and topology state.
          call write_hypertree_mesh(writer)
        end if
      class default
        error stop 'unexpected leaf type'
      end select
      call check_error_handler()
      if (cycle == 1) then
        if (stat == 0 .or. .not.allocated(errmsg)) error stop 'expected creation failure'
        if (blocks(j)%writer%root() >= 0) error stop 'failed creation left an open root'
      else
        if (stat /= 0) error stop errmsg
        if (blocks(j)%writer%root() < 0) error stop 'standalone creation failed'
        if (blocks(j)%writer%temporal_enabled()) error stop 'wrong standalone mode'
      end if
      call blocks(j)%writer%close()
      call check_error_handler()
    end do
  end do
  stat = H5Gclose(parent_id)
  if (stat < 0) error stop 'parent group close failed'
  stat = H5Fclose(file_id)
  if (stat < 0) error stop 'parent file close failed'
  call ctx%free()
  call test_env_finalize()

contains

  subroutine register_fields(dataset)
    class(vtkhdf_simple_dataset), intent(inout) :: dataset
    select type (dataset)
    class is (vtkhdf_unstructured_writer)
      call dataset%write_field_data('cycle',cycle)
      field_var = dataset%register_temporal_field_data('time_value',0.0_real64)
    class is (vtkhdf_hypertree_writer)
      call dataset%write_field_data('cycle',cycle)
      field_var = dataset%register_temporal_field_data('time_value',0.0_real64)
    class is (vtkhdf_structured_base)
      call dataset%write_field_data('cycle',cycle)
      field_var = dataset%register_temporal_field_data('time_value',0.0_real64)
    end select
  end subroutine

  subroutine write_hypertree_mesh(writer)
    class(vtkhdf_hypertree_writer), intent(inout) :: writer
    integer :: i, branch_factor
    ! Change the grid and branch factor on reuse; only rank zero owns a tree.
    branch_factor = 2+modulo(cycle+1,2)
    call writer%write_mesh([cycle+1,1,1],branch_factor=branch_factor)
    call writer%write_mesh_geometry([(real(i,real64),i=0,cycle)],[0.0_real64],[0.0_real64])
    if (rank == 0) then
      call writer%write_mesh_topology([1],[1,branch_factor+2],[.true.,(.false.,i=1,branch_factor)])
      call writer%write_cell_data('value',[(real(i,real64),i=1,branch_factor+1)])
    else
      call writer%write_mesh_topology([integer::],[1],[logical::])
#ifdef GNU_PR123899
      block
        real(real64) :: values(0)
        call writer%write_cell_data('value',values)
      end block
#else
      call writer%write_cell_data('value',[real(real64)::])
#endif
    end if
  end subroutine

  subroutine check_error_handler()
    type(c_funptr) :: handler
    type(c_ptr) :: data
    integer :: ierr
    ierr = H5Eget_auto2(H5E_DEFAULT,handler,data)
    if (ierr < 0) error stop 'failed to query HDF5 error handler'
    if (c_associated(handler) .neqv. c_associated(saved_error_handler)) error stop 'HDF5 error handler changed'
    if (c_associated(handler)) then
      if (.not.c_associated(handler,saved_error_handler)) error stop 'HDF5 error handler changed'
    end if
    if (c_associated(data) .neqv. c_associated(saved_error_data)) error stop 'HDF5 error handler data changed'
    if (c_associated(data)) then
      if (.not.c_associated(data,saved_error_data)) error stop 'HDF5 error handler data changed'
    end if
  end subroutine

end program
