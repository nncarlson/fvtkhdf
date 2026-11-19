!!
!! VTKHDF_MB_FILE_TYPE
!!
!! This module defines a derived type for exporting mesh-based solution data
!! to a VTKHDF file readable by the ParaView visualization tool. The type
!! produces a MultiBlockDataSet (MB) dataset, and supports static and
!! time-dependent cell, point, and field data. Composite-owned typed blocks
!! expose their mesh and data operations directly through borrowed pointers.
!!
!! The output is limited to a flat collection of simple dataset blocks;
!! hierarchical nesting of MultiBlockDataSets is not supported.
!!
!! SPDX-License-Identifier: BSD-2-Clause
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! The latest specification for the VTKHDF file format is available at:
!! https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html#
!!

#include "vtkhdf_assert.inc"

module vtkhdf_mb_file_type

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

  type, public :: vtkhdf_mb_file
    private
    type(vtkhdf_file_context) :: file
    integer(hid_t) :: vtk_id=H5I_INVALID_HID, ass_id=H5I_INVALID_HID
    integer(int64) :: data_chunk_bytes = 128_int64*1024
    integer :: next_block_index = 0 ! zero-based index in the file
    type(block_item), allocatable :: blocks(:) ! a 1-indexed array
    logical :: is_temporal = .false., start_time_step_called = .false.
  contains
#ifdef USE_MPI_F08
    generic :: create => create_f08, create_int
    procedure, private :: create_f08, create_int
#elif defined(USE_MPI)
    generic :: create => create_int
    procedure, private :: create_int
#else
    procedure :: create
#endif
    procedure, private :: create_impl
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
    final :: vtkhdf_mb_file_delete
  end type

  type, extends(named_record) :: block_item
    class(vtkhdf_simple_dataset), pointer :: b => null()
  end type

contains

  !! Finalizer for vtkhdf_mb_file objects. We free heap memory we own but avoid
  !! doing things that may require syncronization with other ranks (MPI/PHDF5)
  !! because where implicit finalization occurs it is not guaranteed to be
  !! collective or ordered with respect to other ranks. This can leak HDF5
  !! IDs and the MPI communicator, but that is unavoidable. Users should
  !! always use CLOSE to do a proper collective cleanup and close of the file.

  impure elemental subroutine vtkhdf_mb_file_delete(this)
    type(vtkhdf_mb_file), intent(inout) :: this
    integer :: j
    if (allocated(this%blocks)) then
      do j = 1, size(this%blocks)
        if (associated(this%blocks(j)%b)) deallocate(this%blocks(j)%b)
      end do
    end if
    this%ass_id  = H5I_INVALID_HID
    this%vtk_id  = H5I_INVALID_HID
    this%file%id = H5I_INVALID_HID
    this%next_block_index = 0
  end subroutine

  !! Cleanly closes H5 identifiers and the file, and default initializes.
  subroutine close(this)
    class(vtkhdf_mb_file), intent(inout) :: this
    integer :: ierr, j
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    if (allocated(this%blocks)) then
      do j = 1, size(this%blocks)
        if (associated(this%blocks(j)%b)) then
          call this%blocks(j)%b%close
          deallocate(this%blocks(j)%b)
        end if
      end do
      deallocate(this%blocks)
    end if
    if (H5Iis_valid(this%ass_id) > 0) ierr = H5Gclose(this%ass_id)
    if (H5Iis_valid(this%vtk_id) > 0) ierr = H5Gclose(this%vtk_id)
    call this%file%close
    call err_ctx%restore(0)
    call finalize(this) ! free local memory
  contains
    subroutine finalize(this)
      class(vtkhdf_mb_file), intent(out) :: this
    end subroutine
  end subroutine

  subroutine flush(this)
    class(vtkhdf_mb_file), intent(in) :: this
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%file%flush
    call err_ctx%restore(0)
  end subroutine

#ifdef USE_MPI_F08
  subroutine create_f08(this, filename, comm, stat, errmsg, data_chunk_bytes)
    use mpi_f08, only: MPI_Comm
    class(vtkhdf_mb_file), intent(out) :: this
    character(*), intent(in) :: filename
    type(MPI_Comm), intent(in) :: comm
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer(int64), intent(in), optional :: data_chunk_bytes
    call this%create(filename, comm%mpi_val, stat, errmsg, data_chunk_bytes)
  end subroutine create_f08

#endif
#ifdef USE_MPI
  subroutine create_int(this, filename, comm, stat, errmsg, data_chunk_bytes)
    class(vtkhdf_mb_file), intent(out) :: this
    character(*), intent(in) :: filename
    integer, intent(in) :: comm
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer(int64), intent(in), optional :: data_chunk_bytes
    call this%create_impl(filename, comm, stat, errmsg, data_chunk_bytes)
  end subroutine create_int
#else
  subroutine create(this, filename, stat, errmsg, data_chunk_bytes)
    class(vtkhdf_mb_file), intent(out) :: this
    character(*), intent(in) :: filename
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer(int64), intent(in), optional :: data_chunk_bytes
    call this%create_impl(filename, stat, errmsg, data_chunk_bytes)
  end subroutine create
#endif

#ifdef USE_MPI
  subroutine create_impl(this, filename, comm, stat, errmsg, data_chunk_bytes)
#else
  subroutine create_impl(this, filename, stat, errmsg, data_chunk_bytes)
#endif

    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: filename
#ifdef USE_MPI
    integer, intent(in) :: comm
#endif
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer(int64), intent(in), optional :: data_chunk_bytes

    integer(hid_t) :: gcpl
    integer(c_int) :: flag
    integer :: ierr
    type(vtkhdf_h5e_context) :: err_ctx

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

    !NB: We stick with the older MB type due to an issue with the modern PDC
    !type; see https://gitlab.kitware.com/vtk/vtk/-/issues/19902
    !call h5_write_attr(this%file%ctx, this%vtk_id, 'Type', 'PartitionedDataSetCollection', stat, errmsg)
    call h5_write_attr(this%file%ctx, this%vtk_id, 'Type', 'MultiBlockDataSet')

    !! Create the Assembly group
    this%ass_id = H5Gcreate(this%vtk_id, 'Assembly', gcpl_id=gcpl)
    ierr = H5Pclose(gcpl)
    INSIST(this%file%ctx%global_all(this%ass_id >= 0))

    stat = 0
    call err_ctx%restore(stat)

  end subroutine create_impl

  ! The returned pointer borrows an individually allocated target owned by MB.
  ! Growing BLOCKS preserves that target; CLOSE invalidates borrowed pointers.
  subroutine add_unstructured_block(this, name, dataset, mode)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: name
    type(vtkhdf_unstructured_block), pointer, intent(out) :: dataset
    integer, intent(in), optional :: mode
    allocate(dataset)
    call this%add_block_impl(name, dataset, mode)
  end subroutine add_unstructured_block

  subroutine add_image_block(this, name, dataset, mode)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: name
    type(vtkhdf_image_block), pointer, intent(out) :: dataset
    integer, intent(in), optional :: mode
    allocate(dataset)
    call this%add_block_impl(name, dataset, mode)
  end subroutine add_image_block

  subroutine add_rectilinear_block(this, name, dataset, mode)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: name
    type(vtkhdf_rectilinear_block), pointer, intent(out) :: dataset
    integer, intent(in), optional :: mode
    allocate(dataset)
    call this%add_block_impl(name, dataset, mode)
  end subroutine add_rectilinear_block

  subroutine add_structured_block(this, name, dataset, mode)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: name
    type(vtkhdf_structured_block), pointer, intent(out) :: dataset
    integer, intent(in), optional :: mode
    allocate(dataset)
    call this%add_block_impl(name, dataset, mode)
  end subroutine add_structured_block

  subroutine add_hypertree_block(this, name, dataset, mode)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: name
    type(vtkhdf_hypertree_block), pointer, intent(out) :: dataset
    integer, intent(in), optional :: mode
    allocate(dataset)
    call this%add_block_impl(name, dataset, mode)
  end subroutine add_hypertree_block

  ! DATASET is an individually allocated pointer target supplied by an overload.
  ! MB retains that target after this procedure returns.
  subroutine add_block_impl(this, name, dataset, mode)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: name
    class(vtkhdf_simple_dataset), target, intent(inout) :: dataset
    integer, intent(in), optional :: mode

    integer :: ierr, mode_
    type(vtkhdf_h5e_context) :: err_ctx
    character(16) :: group_name
    character(:), allocatable :: external_name, internal_name

    INSIST(this%file%id >= 0)

    mode_ = VTKHDF_STATIC
    if (present(mode)) mode_ = mode
    ! Temporal blocks must exist before the common timeline starts.
    if (mode_ > VTKHDF_STATIC) then
      INSIST(.not.this%start_time_step_called)
    end if

    call resolve_record_name(name, 'Unamed_Block', external_name, internal_name, this%blocks)

    if (.not.allocated(this%blocks)) then
      allocate(this%blocks(1))
    else
      block
        integer :: j, n
        type(block_item), allocatable :: tmp(:)
        n = size(this%blocks)
        allocate(tmp(n+1))
        do j = 1, n
          tmp(j)%b => this%blocks(j)%b
          call move_alloc(this%blocks(j)%internal_name, tmp(j)%internal_name)
        end do
        call move_alloc(tmp, this%blocks)
      end block
    end if

    this%blocks(this%next_block_index+1)%b => dataset
    this%blocks(this%next_block_index+1)%internal_name = internal_name

    call err_ctx%capture
    write(group_name,'(a,i0)') 'Group_', this%next_block_index
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
    class(vtkhdf_mb_file), intent(inout) :: this
    real(real32), intent(in) :: time
    call this%start_time_step(real(time, real64))
  end subroutine start_time_step_real32

  subroutine start_time_step_real64(this, time)
    class(vtkhdf_mb_file), intent(inout) :: this
    real(real64), intent(in) :: time
    type(vtkhdf_h5e_context) :: err_ctx
    integer :: j
    call err_ctx%capture
    INSIST(this%is_temporal)
    this%start_time_step_called = .true.
    if (allocated(this%blocks)) then
      do j = 1, size(this%blocks)
        associate (b => this%blocks(j)%b)
          if (b%temporal_enabled()) call b%start_time_step(time)
        end associate
      end do
    end if
    call err_ctx%restore(0)
  end subroutine start_time_step_real64

  subroutine finalize_time_step(this)
    class(vtkhdf_mb_file), intent(inout) :: this
    type(vtkhdf_h5e_context) :: err_ctx
    integer :: j
    call err_ctx%capture
    INSIST(this%is_temporal)
    if (allocated(this%blocks)) then
      do j = 1, size(this%blocks)
        associate (b => this%blocks(j)%b)
          if (b%temporal_enabled() .and. b%step_is_open()) call b%finalize_time_step()
        end associate
      end do
    end if
    call err_ctx%restore(0)
  end subroutine finalize_time_step

end module vtkhdf_mb_file_type
