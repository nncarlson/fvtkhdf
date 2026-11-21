!! Shared implementation for MultiBlockDataSet and PartitionedDataSetCollection.
!! Composite-owned typed blocks expose mesh and data writes through borrowed
!! pointers. This module is internal; use the concrete file modules.
!! SPDX-License-Identifier: BSD-2-Clause

#include "vtkhdf_assert.inc"

module vtkhdf_composite_file_type

  use,intrinsic :: iso_fortran_env
  use vtkhdf_h5_c_binding
  use vtkhdf_h5
  use vtkhdf_data_support, only: named_record, resolve_record_name
  use vtkhdf_data_handle_types
  use vtkhdf_simple_dataset_class
  use vtkhdf_unstructured_writer_type, only: vtkhdf_unstructured_block
  use vtkhdf_image_writer_type, only: vtkhdf_image_block
  use vtkhdf_rectilinear_writer_type, only: vtkhdf_rectilinear_block
  use vtkhdf_structured_writer_type, only: vtkhdf_structured_block
  use vtkhdf_hypertree_writer_type, only: vtkhdf_hypertree_block
  use vtkhdf_file_modes, only: VTKHDF_STATIC, VTKHDF_FIXED_MESH, VTKHDF_MOVING_MESH, VTKHDF_DYNAMIC_MESH
  use vtkhdf_file_context_type
  use vtkhdf_h5e_context_type
  use vtkhdf_version_param
  implicit none
  private

  public :: vtkhdf_version, vtkhdf_unstructured_block
  public :: vtkhdf_image_block, vtkhdf_rectilinear_block, vtkhdf_structured_block
  public :: vtkhdf_hypertree_block
  public :: vtkhdf_cell_data_handle, vtkhdf_point_data_handle, vtkhdf_field_data_handle
  public :: VTKHDF_STATIC, VTKHDF_FIXED_MESH, VTKHDF_MOVING_MESH, VTKHDF_DYNAMIC_MESH

  type, extends(named_record) :: block_item
    class(vtkhdf_simple_dataset), pointer :: b => null()
    character(:), allocatable :: storage_name
  end type

  type :: block_registry
    type(block_item), allocatable :: blocks(:) ! a 1-indexed array
  contains
    final :: delete_registry
  end type

  integer, parameter :: ASSEMBLY_NODE_CHILD = 1, ASSEMBLY_BLOCK_CHILD = 2

  type, extends(named_record) :: assembly_child
    ! INTERNAL_NAME remains unallocated until commitment resolves the names.
    integer :: kind = 0, index = 0
  end type

  type :: assembly_node
    character(:), allocatable :: requested_name
    type(assembly_child), allocatable :: children(:)
    integer :: nchildren = 0
  end type

  ! Entirely local, allocatable Assembly state. Root is node 0; positive
  ! node IDs and one-based block registry slots occupy separate namespaces.
  type :: assembly_description
    type(assembly_node), allocatable :: nodes(:)
    integer :: nnodes = 0
    logical :: committed = .false.
  end type

  type, abstract, public :: vtkhdf_composite_file
    private
    type(vtkhdf_file_context) :: file
    integer(hid_t) :: vtk_id=H5I_INVALID_HID, ass_id=H5I_INVALID_HID
    integer(int64) :: data_chunk_bytes = 128_int64*1024
    integer :: next_block_index = 0 ! zero-based index in the file
    type(block_registry) :: registry
    type(assembly_description) :: assembly
    logical :: is_temporal = .false., start_time_step_called = .false.
  contains
    procedure :: internal_create ! internal entry point for the concrete CREATE wrappers
    procedure :: close
    procedure :: flush
    generic :: add_block => add_unstructured_block, add_image_block, add_rectilinear_block, &
        add_structured_block, add_hypertree_block
    procedure, private :: add_unstructured_block, add_image_block, add_rectilinear_block, &
        add_structured_block, add_hypertree_block
    procedure, private :: add_block_impl
    generic :: start_time_step => start_time_step_real32, start_time_step_real64
    procedure, private :: start_time_step_real32, start_time_step_real64
    procedure :: finalize_time_step
    procedure :: add_node, add_block_ref, set_assembly
  end type

contains

  ! The concrete registry component can be finalized even though its owner is
  ! abstract. Implicit finalization releases heap targets only, never MPI/HDF5
  ! resources; collective cleanup requires an explicit CLOSE.
  impure elemental subroutine delete_registry(this)
    type(block_registry), intent(inout) :: this
    integer :: j
    if (allocated(this%blocks)) then
      do j = 1, size(this%blocks)
        if (associated(this%blocks(j)%b)) deallocate(this%blocks(j)%b)
      end do
    end if
  end subroutine

  !! Cleanly closes H5 identifiers and the file, and default initializes.
  subroutine close(this)
    class(vtkhdf_composite_file), intent(inout) :: this
    integer :: ierr, j
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    if (this%file%id >= 0) then
#ifndef NDEBUG
      call check_assembly_agreement(this)
#endif
      if (allocated(this%assembly%nodes) .and. .not.this%assembly%committed) call this%set_assembly()
    end if
    if (allocated(this%registry%blocks)) then
      do j = 1, size(this%registry%blocks)
        if (associated(this%registry%blocks(j)%b)) then
          call this%registry%blocks(j)%b%close
          deallocate(this%registry%blocks(j)%b)
        end if
      end do
      deallocate(this%registry%blocks)
    end if
    if (H5Iis_valid(this%ass_id) > 0) ierr = H5Gclose(this%ass_id)
    if (H5Iis_valid(this%vtk_id) > 0) ierr = H5Gclose(this%vtk_id)
    call this%file%close
    call err_ctx%restore(0)
    call finalize(this) ! free local memory
  contains
    subroutine finalize(this)
      class(vtkhdf_composite_file), intent(out) :: this
    end subroutine
  end subroutine

  subroutine flush(this)
    class(vtkhdf_composite_file), intent(in) :: this
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%file%flush
    call err_ctx%restore(0)
  end subroutine

#ifdef USE_MPI
  subroutine internal_create(this, filename, comm, stat, errmsg, data_chunk_bytes, composite_type)
#else
  subroutine internal_create(this, filename, stat, errmsg, data_chunk_bytes, composite_type)
#endif

    class(vtkhdf_composite_file), intent(inout) :: this
    character(*), intent(in) :: filename
#ifdef USE_MPI
    integer, intent(in) :: comm
#endif
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer(int64), intent(in), optional :: data_chunk_bytes
    character(*), intent(in) :: composite_type

    integer(hid_t) :: gcpl
    integer(c_int) :: flag
    integer :: ierr
    type(vtkhdf_h5e_context) :: err_ctx

    INSIST(composite_type == "MultiBlockDataSet" .or. composite_type == "PartitionedDataSetCollection")

    call init_hdf5
    call err_ctx%capture

    this%data_chunk_bytes = 128_int64*1024
    if (present(data_chunk_bytes)) this%data_chunk_bytes = data_chunk_bytes
    INSIST(this%data_chunk_bytes > 0)

#ifdef USE_MPI
    call this%file%open(filename, comm, stat, errmsg)
#else
    call this%file%open(filename, stat, errmsg)
#endif
    if (stat /= 0) then
      errmsg = errmsg // ' "' // filename // '"'
      call err_ctx%restore(stat)
      return
    end if

    !! Group creation properties for PDC datasets (no harm to MB)
    !! Applies to VTKHDF and Assembly groups
    flag = ior(H5P_CRT_ORDER_TRACKED, H5P_CRT_ORDER_INDEXED)
    gcpl = H5Pcreate(H5P_GROUP_CREATE)
    stat = H5Pset_link_creation_order(gcpl, flag)
    INSIST(stat == 0)

    !! Create the root VTKHDF group and write its attributes
    this%vtk_id = H5Gcreate(this%file%id, 'VTKHDF', gcpl_id=gcpl)
    INSIST(this%file%ctx%global_all(this%vtk_id >= 0))

    call h5_write_attr(this%file%ctx, this%vtk_id, 'Version', vtkhdf_version)

    call h5_write_attr(this%file%ctx, this%vtk_id, 'Type', composite_type)

    !! Create the Assembly group
    this%ass_id = H5Gcreate(this%vtk_id, 'Assembly', gcpl_id=gcpl)
    ierr = H5Pclose(gcpl)
    INSIST(this%file%ctx%global_all(this%ass_id >= 0))

    stat = 0
    call err_ctx%restore(stat)

  end subroutine internal_create

  ! The returned pointer borrows an individually allocated target owned by the composite.
  ! Growing BLOCKS preserves that target; CLOSE invalidates borrowed pointers.
  subroutine add_unstructured_block(this, name, dataset, mode)
    class(vtkhdf_composite_file), intent(inout) :: this
    character(*), intent(in) :: name
    type(vtkhdf_unstructured_block), pointer, intent(out) :: dataset
    integer, intent(in), optional :: mode
    allocate(dataset)
    call this%add_block_impl(name, dataset, mode)
  end subroutine

  subroutine add_image_block(this, name, dataset, mode)
    class(vtkhdf_composite_file), intent(inout) :: this
    character(*), intent(in) :: name
    type(vtkhdf_image_block), pointer, intent(out) :: dataset
    integer, intent(in), optional :: mode
    allocate(dataset)
    call this%add_block_impl(name, dataset, mode)
  end subroutine

  subroutine add_rectilinear_block(this, name, dataset, mode)
    class(vtkhdf_composite_file), intent(inout) :: this
    character(*), intent(in) :: name
    type(vtkhdf_rectilinear_block), pointer, intent(out) :: dataset
    integer, intent(in), optional :: mode
    allocate(dataset)
    call this%add_block_impl(name, dataset, mode)
  end subroutine

  subroutine add_structured_block(this, name, dataset, mode)
    class(vtkhdf_composite_file), intent(inout) :: this
    character(*), intent(in) :: name
    type(vtkhdf_structured_block), pointer, intent(out) :: dataset
    integer, intent(in), optional :: mode
    allocate(dataset)
    call this%add_block_impl(name, dataset, mode)
  end subroutine

  subroutine add_hypertree_block(this, name, dataset, mode)
    class(vtkhdf_composite_file), intent(inout) :: this
    character(*), intent(in) :: name
    type(vtkhdf_hypertree_block), pointer, intent(out) :: dataset
    integer, intent(in), optional :: mode
    allocate(dataset)
    call this%add_block_impl(name, dataset, mode)
  end subroutine

  ! DATASET is an individually allocated pointer target supplied by an overload.
  ! The composite retains that target after this procedure returns.
  subroutine add_block_impl(this, name, dataset, mode)

    class(vtkhdf_composite_file), intent(inout) :: this
    character(*), intent(in) :: name
    class(vtkhdf_simple_dataset), target, intent(inout) :: dataset
    integer, intent(in), optional :: mode

    integer :: ierr, mode_
    type(vtkhdf_h5e_context) :: err_ctx
    character(16) :: group_name
    character(:), allocatable :: external_name, internal_name

    INSIST(this%file%id >= 0)
    INSIST(.not.this%assembly%committed)

    mode_ = VTKHDF_STATIC
    if (present(mode)) mode_ = mode
    ! Temporal blocks must exist before the common timeline starts.
    if (mode_ > VTKHDF_STATIC) then
      INSIST(.not.this%start_time_step_called)
    end if

    call resolve_record_name(name, 'Unamed_Block', external_name, internal_name, this%registry%blocks)

    if (.not.allocated(this%registry%blocks)) then
      allocate(this%registry%blocks(1))
    else
      block
        integer :: j, n
        type(block_item), allocatable :: tmp(:)
        n = size(this%registry%blocks)
        allocate(tmp(n+1))
        do j = 1, n
          tmp(j)%b => this%registry%blocks(j)%b
          call move_alloc(this%registry%blocks(j)%internal_name, tmp(j)%internal_name)
          call move_alloc(this%registry%blocks(j)%storage_name, tmp(j)%storage_name)
        end do
        call move_alloc(tmp, this%registry%blocks)
      end block
    end if

    this%registry%blocks(this%next_block_index+1)%b => dataset
    this%registry%blocks(this%next_block_index+1)%internal_name = internal_name

    call err_ctx%capture
    write(group_name,'(a,i0)') 'Group_', this%next_block_index
    this%registry%blocks(this%next_block_index+1)%storage_name = trim(group_name)
    select type (dataset)
    type is (vtkhdf_unstructured_block)
      call dataset%init(this%vtk_id, trim(group_name), this%file%ctx, mode_, this%data_chunk_bytes)
      call dataset%set_block_metadata(this%next_block_index, external_name)
    type is (vtkhdf_image_block)
      call dataset%init(this%vtk_id, trim(group_name), this%file%ctx, mode_, this%data_chunk_bytes)
      call dataset%set_block_metadata(this%next_block_index, external_name)
    type is (vtkhdf_rectilinear_block)
      call dataset%init(this%vtk_id, trim(group_name), this%file%ctx, mode_, this%data_chunk_bytes)
      call dataset%set_block_metadata(this%next_block_index, external_name)
    type is (vtkhdf_structured_block)
      call dataset%init(this%vtk_id, trim(group_name), this%file%ctx, mode_, this%data_chunk_bytes)
      call dataset%set_block_metadata(this%next_block_index, external_name)
    type is (vtkhdf_hypertree_block)
      call dataset%init(this%vtk_id, trim(group_name), this%file%ctx, mode_, this%data_chunk_bytes)
      call dataset%set_block_metadata(this%next_block_index, external_name)
    class default
      INSIST(.false.)
    end select
    this%is_temporal = this%is_temporal .or. dataset%temporal_enabled()
    this%next_block_index = this%next_block_index + 1

    ierr = H5Lcreate_soft('/VTKHDF/' // trim(group_name), this%ass_id, internal_name)
    INSIST(this%file%ctx%global_all(ierr >= 0))
    call err_ctx%restore(0)

  end subroutine add_block_impl

  !! Mark the start of a new time step with time value TIME. Subsequent output
  !! of time-dependent datasets will be associated with this time step.

  subroutine start_time_step_real32(this, time)
    class(vtkhdf_composite_file), intent(inout) :: this
    real(real32), intent(in) :: time
    call this%start_time_step(real(time, real64))
  end subroutine

  subroutine start_time_step_real64(this, time)
    class(vtkhdf_composite_file), intent(inout) :: this
    real(real64), intent(in) :: time
    type(vtkhdf_h5e_context) :: err_ctx
    integer :: j
    call err_ctx%capture
    INSIST(this%is_temporal)
    this%start_time_step_called = .true.
    if (allocated(this%registry%blocks)) then
      do j = 1, size(this%registry%blocks)
        associate (b => this%registry%blocks(j)%b)
          if (b%temporal_enabled()) call b%start_time_step(time)
        end associate
      end do
    end if
    call err_ctx%restore(0)
  end subroutine

  subroutine finalize_time_step(this)
    class(vtkhdf_composite_file), intent(inout) :: this
    type(vtkhdf_h5e_context) :: err_ctx
    integer :: j
    call err_ctx%capture
    INSIST(this%is_temporal)
    if (allocated(this%registry%blocks)) then
      do j = 1, size(this%registry%blocks)
        associate (b => this%registry%blocks(j)%b)
          if (b%temporal_enabled() .and. b%step_is_open()) call b%finalize_time_step()
        end associate
      end do
    end if
    call err_ctx%restore(0)
  end subroutine

  ! Assembly construction, name resolution, and serialization. Keep these procedures
  ! together so the local description can later be extracted into a helper type.

  integer function add_node(this, name, parent) result(node)

    class(vtkhdf_composite_file), intent(inout) :: this
    character(*), intent(in) :: name
    integer, intent(in), optional :: parent

    type(assembly_node), allocatable :: grown(:)
    integer :: p, j, capacity

    INSIST(this%file%id >= 0)
    INSIST(.not.this%assembly%committed)

    p = 0
    if (present(parent)) p = parent
    INSIST(p >= 0 .and. p <= this%assembly%nnodes)
    if (.not.allocated(this%assembly%nodes)) allocate(this%assembly%nodes(0:7))

    node = this%assembly%nnodes + 1
    if (node > ubound(this%assembly%nodes,1)) then
      capacity = 2*size(this%assembly%nodes)
      allocate(grown(0:capacity-1))
      do j = 0, this%assembly%nnodes
        call move_alloc(this%assembly%nodes(j)%requested_name, grown(j)%requested_name)
        call move_alloc(this%assembly%nodes(j)%children, grown(j)%children)
        grown(j)%nchildren = this%assembly%nodes(j)%nchildren
      end do
      call move_alloc(grown, this%assembly%nodes)
    end if
    this%assembly%nnodes = node
    this%assembly%nodes(node)%requested_name = name
    call append_assembly_child(this%assembly%nodes(p), assembly_child(kind=ASSEMBLY_NODE_CHILD,index=node))

  end function add_node


  subroutine add_block_ref(this, dataset, parent)

    class(vtkhdf_composite_file), intent(inout) :: this
    class(vtkhdf_simple_dataset), target, intent(in) :: dataset
    integer, intent(in), optional :: parent

    integer :: p, j, slot

    INSIST(this%file%id >= 0)
    INSIST(.not.this%assembly%committed)

    ! Allocating NODES also marks the start of custom Assembly construction.
    if (.not.allocated(this%assembly%nodes)) allocate(this%assembly%nodes(0:7))

    p = 0
    if (present(parent)) p = parent
    INSIST(p >= 0 .and. p <= this%assembly%nnodes)

    ! Find DATASET's registry slot by pointer identity. A match establishes
    ! that DATASET belongs to this file.
    slot = 0
    if (allocated(this%registry%blocks)) then
      do j = 1, size(this%registry%blocks)
        if (associated(this%registry%blocks(j)%b, dataset)) then
          slot = j
          exit
        end if
      end do
    end if
    INSIST(slot > 0)

    ! Append a block reference under node P unless it is already present.
    ! Duplicate references retain their original position.
    associate (node => this%assembly%nodes(p))
      do j = 1, node%nchildren
        if (node%children(j)%kind /= ASSEMBLY_BLOCK_CHILD) cycle
        if (node%children(j)%index == slot) return
      end do
    end associate
    call append_assembly_child(this%assembly%nodes(p), assembly_child(kind=ASSEMBLY_BLOCK_CHILD,index=slot))

  end subroutine add_block_ref

  subroutine append_assembly_child(node, child)
    type(assembly_node), intent(inout) :: node
    type(assembly_child), intent(in) :: child
    type(assembly_child), allocatable :: grown(:)
    if (.not.allocated(node%children)) allocate(node%children(4))
    if (node%nchildren >= size(node%children)) then
      allocate(grown(2*size(node%children)))
      grown(:node%nchildren) = node%children(:node%nchildren)
      call move_alloc(grown, node%children)
    end if
    node%nchildren = node%nchildren + 1
    node%children(node%nchildren) = child
  end subroutine

  subroutine set_assembly(this)

    class(vtkhdf_composite_file), intent(inout) :: this

    type(vtkhdf_h5e_context) :: err_ctx
    integer(hid_t) :: gcpl
    integer :: j, ierr

    INSIST(this%file%id >= 0)
#ifndef NDEBUG
    call check_assembly_agreement(this)
#endif
    INSIST(.not.this%assembly%committed)
    if (allocated(this%assembly%nodes)) then
      call resolve_assembly(this)
      call err_ctx%capture
      ! The registry knows every existing flat link; no HDF5 traversal is needed.
      do j = 1, this%next_block_index
        ierr = H5Ldelete(this%ass_id, this%registry%blocks(j)%internal_name)
        INSIST(this%file%ctx%global_all(ierr >= 0))
      end do
      gcpl = H5Pcreate(H5P_GROUP_CREATE)
      INSIST(this%file%ctx%global_all(gcpl >= 0))
      ierr = H5Pset_link_creation_order(gcpl, int(ior(H5P_CRT_ORDER_TRACKED,H5P_CRT_ORDER_INDEXED),c_int))
      INSIST(this%file%ctx%global_all(ierr >= 0))
      call write_assembly(this, 0, this%ass_id, gcpl)
      ierr = H5Pclose(gcpl)
      INSIST(this%file%ctx%global_all(ierr >= 0))
      call err_ctx%restore(0)
    end if
    ! With no custom construction the flat hierarchy is already correct.
    this%assembly%committed = .true.

  end subroutine set_assembly

  ! Add implicit root references and resolve names in place before writing
  ! the HDF5 links.

  subroutine resolve_assembly(this)

    class(vtkhdf_composite_file), intent(inout) :: this

    type(named_record), allocatable :: names(:)
    logical, allocatable :: placed(:)
    character(:), allocatable :: external_name, internal_name
    integer :: p, j, n, nused

    ! Append implicit root references for blocks with no explicit placement.
    allocate(placed(this%next_block_index), source=.false.)
    do p = 0, this%assembly%nnodes
      associate (node => this%assembly%nodes(p))
        do j = 1, node%nchildren
          if (node%children(j)%kind == ASSEMBLY_BLOCK_CHILD) then
            placed(node%children(j)%index) = .true.
          end if
        end do
      end associate
    end do

    do j = 1, size(placed)
      if (placed(j)) cycle
      call append_assembly_child(this%assembly%nodes(0), assembly_child(kind=ASSEMBLY_BLOCK_CHILD,index=j))
    end do

    ! Resolve names independently for each parent. Reserve the established
    ! names of its block references first, including implicit root references.
    ! Normalize and disambiguate ordinary node names in insertion order
    ! against block names and previously resolved node names.

    do p = 0, this%assembly%nnodes
      n = this%assembly%nodes(p)%nchildren
      allocate(names(n))
      ! Reserve only the block names actually present under this parent.
      nused = 0
      do j = 1, n
        associate (child => this%assembly%nodes(p)%children(j))
          if (child%kind /= ASSEMBLY_BLOCK_CHILD) cycle
          child%internal_name = this%registry%blocks(child%index)%internal_name
          nused = nused + 1
          names(nused)%internal_name = child%internal_name
        end associate
      end do
      ! Resolve node names against block names and earlier siblings.
      do j = 1, n
        associate (child => this%assembly%nodes(p)%children(j))
          if (child%kind /= ASSEMBLY_NODE_CHILD) cycle
          call resolve_record_name(this%assembly%nodes(child%index)%requested_name, 'Unnamed_Node', &
                                   external_name, internal_name, names(:nused))
          child%internal_name = internal_name
          nused = nused + 1
          names(nused)%internal_name = internal_name
        end associate
      end do
      deallocate(names)
    end do

  end subroutine resolve_assembly

  recursive subroutine write_assembly(this, node, parent_id, gcpl)

    class(vtkhdf_composite_file), intent(in) :: this
    integer, intent(in) :: node
    integer(hid_t), intent(in) :: parent_id, gcpl

    integer(hid_t) :: id
    integer :: j, ierr

    do j = 1, this%assembly%nodes(node)%nchildren
      associate (child => this%assembly%nodes(node)%children(j))
        select case (child%kind)
        case (ASSEMBLY_NODE_CHILD)
          id = H5Gcreate(parent_id, child%internal_name, gcpl_id=gcpl)
          INSIST(this%file%ctx%global_all(id >= 0))
          call write_assembly(this, child%index, id, gcpl)
          ierr = H5Gclose(id)
          INSIST(this%file%ctx%global_all(ierr >= 0))
        case (ASSEMBLY_BLOCK_CHILD)
          ierr = H5Lcreate_soft('/VTKHDF/'//this%registry%blocks(child%index)%storage_name, &
                               parent_id, child%internal_name)
          INSIST(this%file%ctx%global_all(ierr >= 0))
        end select
      end associate
    end do

  end subroutine write_assembly

  ! Compare counts before traversing variable-length buffers. CLOSE calls this
  ! before branching on implicit commitment, so a mismatch cannot strand ranks
  ! in different collective paths. These follow the existing debug-check policy.

  subroutine check_assembly_agreement(this)

    class(vtkhdf_composite_file), intent(in) :: this

    integer :: p, j, state(4)

    state = [merge(1,0,this%assembly%committed), merge(1,0,allocated(this%assembly%nodes)), &
             this%assembly%nnodes, this%next_block_index]
    ASSERT(this%file%ctx%identical(state))
    do j = 1, this%next_block_index
      ASSERT(this%file%ctx%identical(this%registry%blocks(j)%internal_name))
      ASSERT(this%file%ctx%identical(this%registry%blocks(j)%storage_name))
    end do
    if (.not.allocated(this%assembly%nodes)) return
    do p = 0, this%assembly%nnodes
      associate (node => this%assembly%nodes(p))
        if (p /= 0) then
          ASSERT(this%file%ctx%identical(node%requested_name))
        end if
        ASSERT(this%file%ctx%identical([node%nchildren]))
        do j = 1, node%nchildren
          ASSERT(this%file%ctx%identical([node%children(j)%kind, node%children(j)%index]))
        end do
      end associate
    end do

  end subroutine check_assembly_agreement

end module vtkhdf_composite_file_type
