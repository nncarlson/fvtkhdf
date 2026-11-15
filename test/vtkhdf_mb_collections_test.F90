! Exercise MB ownership and orchestration of heterogeneous block collections.
program vtkhdf_mb_collections_test
  use,intrinsic :: iso_fortran_env, only: int8, real64
  use vtkhdf_mb_file_type
  use vtkhdf_simple_dataset_class, only: vtkhdf_simple_dataset
  use vtkhdf_ug_writer_type, only: vtkhdf_ug_writer
  use vtkhdf_structured_base_class, only: vtkhdf_structured_base
  use vtkhdf_h5_c_binding, only: hid_t, H5Iis_valid, H5Iget_file_id, H5Fclose
  use vtkhdf_test_env
#ifdef USE_MPI_F08
  use mpi_f08, only: MPI_COMM_WORLD
#elif defined(USE_MPI)
  use mpi, only: MPI_COMM_WORLD
#endif
  implicit none

  type :: block_pointer
    class(vtkhdf_simple_dataset), pointer :: p => null()
  end type
  type(vtkhdf_mb_file) :: mb
  type(vtkhdf_ug_block), pointer :: ug
  type(vtkhdf_image_block), pointer :: image
  type(vtkhdf_rectilinear_block), pointer :: rect
  type(vtkhdf_structured_block), pointer :: struct
  type(block_pointer) :: blocks(4)
  type(vtkhdf_field_data_handle) :: clock(4)
  integer(hid_t) :: roots(4), file_id
  integer :: rank, nproc, stat, cycle, mode, image_mode, j, step
  integer :: bounds(3,2), selected(3,2)
  real(real64) :: points(3,2,2,1)
  character(:), allocatable :: errmsg
  character(12), parameter :: cases(3) = [character(12) :: 'static', 'temporal', 'mixed']

  call test_env_init(rank,nproc)
  bounds = reshape([0,0,0,1,1,0],[3,2])
  selected = bounds
  if (rank /= 0) selected(1,2) = -1
  points = reshape([0.,0.,0.,1.,0.,0.,0.,1.,0.,1.,1.,0.],[3,2,2,1])

  do cycle = 1,3
    ! Reuse the owner after closing each collection. Keep mixed static/temporal
    ! blocks in a separate file: VTK readers have trouble with that mixture.
#ifdef USE_MPI
    call mb%create('mb_collection_'//trim(cases(cycle))//'.vtkhdf',MPI_COMM_WORLD,stat,errmsg)
#else
    call mb%create('mb_collection_'//trim(cases(cycle))//'.vtkhdf',stat,errmsg)
#endif
    if (stat /= 0) error stop errmsg
    mode = VTKHDF_STATIC
    if (cycle > 1) mode = VTKHDF_FIXED_MESH
    image_mode = mode
    if (cycle == 3) image_mode = VTKHDF_STATIC

    ! Every insertion grows the heterogeneous collection. Repeated names must
    ! be resolved in one namespace, independent of the concrete block type.
    if (cycle == 1) then
      call mb%add_block('Repeated Block',image)
    else
      call mb%add_block('Repeated Block',image,mode=image_mode)
    end if
    blocks(1)%p => image
    roots(1) = image%root()
    call mb%add_block('Repeated Block',rect,mode=mode)
    blocks(2)%p => rect
    roots(2) = rect%root()
    call mb%add_block('Repeated Block',struct,mode=mode)
    blocks(3)%p => struct
    roots(3) = struct%root()
    call mb%add_block('Repeated Block',ug,mode=mode)
    blocks(4)%p => ug
    roots(4) = ug%root()

    ! Use all saved pointers after collection growth, both directly and through
    ! their common base. MB must own exactly the same targets.
    if (rank == 0) then
      call ug%write_mesh(reshape(points,[3,4]),[1,2,4,3],[1,5],[9_int8])
    else
      call ug%write_mesh(reshape([real(real64)::],[3,0]),[integer::],[1],[integer(int8)::])
    end if
    call image%write_mesh(bounds,[0.,0.,0.],[1.,1.,1.])
    call rect%write_mesh(bounds,[0.,1.],[0.,1.],[0.])
    call struct%write_mesh(bounds)
    call struct%write_mesh_geometry(points,[0,0,0],selected)
    do j = 1,4
      if (blocks(j)%p%root() /= roots(j)) error stop 'block target changed'
      select type (dataset => blocks(j)%p)
      class is (vtkhdf_ug_writer)
        call dataset%write_field_data('id',j)
        if (dataset%temporal_enabled()) clock(j) = dataset%register_temporal_field_data('clock',0.0_real64)
      class is (vtkhdf_structured_base)
        call dataset%write_field_data('id',j)
        if (dataset%temporal_enabled()) clock(j) = dataset%register_temporal_field_data('clock',0.0_real64)
      end select
    end do

    if (cycle > 1) then
      do step = 0,2
        ! Step 2 implicitly finalizes step 1; CLOSE will finalize step 2.
        call mb%start_time_step(real(step,real64))
        do j = 1,4
          if (blocks(j)%p%step_is_open() .neqv. blocks(j)%p%temporal_enabled()) &
              error stop 'incorrect block step state'
          if (.not.blocks(j)%p%temporal_enabled()) cycle
          select type (dataset => blocks(j)%p)
          class is (vtkhdf_ug_writer)
            call dataset%write_field_data(clock(j),real(10*j+step,real64))
          class is (vtkhdf_structured_base)
            call dataset%write_field_data(clock(j),real(10*j+step,real64))
          end select
        end do
        if (step == 0) then
          call mb%finalize_time_step()
          do j = 1,4
            if (blocks(j)%p%step_is_open()) error stop 'block step not finalized'
          end do
        end if
      end do
    end if

    call mb%flush()
    file_id = H5Iget_file_id(roots(1))
    stat = H5Fclose(file_id) ! release the extra reference from H5Iget_file_id
    if (stat < 0) error stop 'failed to release extra file reference'
    call mb%close()
    call mb%close()
    if (H5Iis_valid(file_id) > 0) error stop 'collection file leaked'
    do j = 1,4
      if (H5Iis_valid(roots(j)) > 0) error stop 'block root leaked'
      nullify(blocks(j)%p)
    end do
    nullify(ug,image,rect,struct)
  end do
  call test_env_finalize()
end program vtkhdf_mb_collections_test
