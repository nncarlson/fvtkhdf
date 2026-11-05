! Exercise the shared lifecycle with heterogeneous borrowed leaves. This is
! deliberately below MB's public API, which still creates only UG blocks.
program vtkhdf_simple_dataset_test
  use,intrinsic :: iso_fortran_env
  use vtkhdf_simple_dataset_class
  use vtkhdf_ug_file_type
  use vtkhdf_image_file_type
  use vtkhdf_rectilinear_file_type
  use vtkhdf_structured_file_type
  use vtkhdf_ctx_type
  use vtkhdf_h5_c_binding
  use vtkhdf_test_env
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  type :: leaf
    class(vtkhdf_simple_dataset), allocatable :: writer
  end type
  type(leaf) :: blocks(4)
  type(vtkhdf_ctx) :: ctx
  integer(hid_t) :: file_id, parent_id, fapl, root_ids(4)
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

  allocate(vtkhdf_ug_file :: blocks(1)%writer)
  allocate(vtkhdf_image_file :: blocks(2)%writer)
  allocate(vtkhdf_rectilinear_file :: blocks(3)%writer)
  allocate(vtkhdf_structured_file :: blocks(4)%writer)
  do cycle = 1,2
    ! The second cycle reuses closed polymorphic objects and the same parent.
    do j = 1,4
      write(name,'(a,i0,a,i0)') 'leaf_',j,'_',cycle
      select type (writer => blocks(j)%writer)
      type is (vtkhdf_ug_file)
        call writer%init(parent_id,trim(name),ctx,mode=VTKHDF_FIXED_MESH)
        if (rank == 0) then
          call writer%write_mesh(reshape(points,[3,4]),[1,2,4,3],[1,5],[9_int8])
        else
          call writer%write_mesh(reshape([real(real64)::],[3,0]),[integer::],[1],[integer(int8)::])
        end if
      type is (vtkhdf_image_file)
        call writer%init(parent_id,trim(name),ctx,mode=VTKHDF_FIXED_MESH)
        call writer%write_mesh(bounds,[0.,0.,0.],[1.,1.,1.])
      type is (vtkhdf_rectilinear_file)
        call writer%init(parent_id,trim(name),ctx,mode=VTKHDF_FIXED_MESH)
        call writer%write_mesh(bounds,[0.,1.],[0.,1.],[0.])
      type is (vtkhdf_structured_file)
        call writer%init(parent_id,trim(name),ctx,mode=VTKHDF_FIXED_MESH)
        call writer%write_mesh(bounds)
        call writer%write_mesh_geometry(points,bounds(:,1),selected)
      class default
        error stop 'unexpected leaf type'
      end select
      ! All subsequent operations dispatch through the abstract ancestor.
      root_ids(j) = blocks(j)%writer%root()
      if (root_ids(j) < 0) error stop 'missing root'
      if (.not.blocks(j)%writer%temporal_enabled()) error stop 'missing temporal mode'
      if (blocks(j)%writer%step_is_open()) error stop 'step already open'
      call blocks(j)%writer%start_time_step(0.0_real32)
      if (.not.blocks(j)%writer%step_is_open()) error stop 'step did not open'
      call blocks(j)%writer%start_time_step(1.0_real64)
      call blocks(j)%writer%finalize_time_step()
      if (blocks(j)%writer%step_is_open()) error stop 'step did not close'
    end do
    do j = 1,4
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
  ! Failed standalone creation must release its duplicated communicator.
  ! The same polymorphic object can then be closed and reopened normally.
  do j = 1,4
    do cycle = 1,3
      filename = 'simple_dataset_leaves.vtkhdf/not-a-directory.vtkhdf'
      if (cycle > 1) write(filename,'(a,i0,a)') 'simple_dataset_owned_',j,'.vtkhdf'
      select type (writer => blocks(j)%writer)
      type is (vtkhdf_ug_file)
#ifdef USE_MPI
        call writer%create(trim(filename),MPI_COMM_WORLD,stat,errmsg)
#else
        call writer%create(trim(filename),stat,errmsg)
#endif
      type is (vtkhdf_image_file)
#ifdef USE_MPI
        call writer%create(trim(filename),MPI_COMM_WORLD,stat,errmsg)
#else
        call writer%create(trim(filename),stat,errmsg)
#endif
      type is (vtkhdf_rectilinear_file)
#ifdef USE_MPI
        call writer%create(trim(filename),MPI_COMM_WORLD,stat,errmsg)
#else
        call writer%create(trim(filename),stat,errmsg)
#endif
      type is (vtkhdf_structured_file)
#ifdef USE_MPI
        call writer%create(trim(filename),MPI_COMM_WORLD,stat,errmsg)
#else
        call writer%create(trim(filename),stat,errmsg)
#endif
      class default
        error stop 'unexpected leaf type'
      end select
      if (cycle == 1) then
        if (stat == 0 .or. .not.allocated(errmsg)) error stop 'expected creation failure'
        if (blocks(j)%writer%root() >= 0) error stop 'failed creation left an open root'
      else
        if (stat /= 0) error stop errmsg
        if (blocks(j)%writer%root() < 0) error stop 'standalone creation failed'
        if (blocks(j)%writer%temporal_enabled()) error stop 'wrong standalone mode'
      end if
      call blocks(j)%writer%close()
    end do
  end do
  stat = H5Gclose(parent_id)
  if (stat < 0) error stop 'parent group close failed'
  stat = H5Fclose(file_id)
  if (stat < 0) error stop 'parent file close failed'
  call ctx%free()
  call test_env_finalize()
end program
