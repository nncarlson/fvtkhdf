# 1 "/app/src/vtkhdf_mb_file_type.F90.fypp"
!!
!! VTKHDF_MB_FILE_TYPE
!!
!! This module defines a derived type for exporting mesh-based solution data
!! to a VTKHDF file readable by the ParaView visualization tool. The type
!! produces a MultiBlockDataSet (MB) dataset, and supports static and
!! time-dependent cell and point data, assuming a static mesh for each block.
!!
!! The output is limited to a flat collection of UnstructuredGrid blocks;
!! hierarchical nesting of MultiBlockDataSets is not supported.
!!
!! Copyright (c) 2026 Neil Carlson <neil.n.carlson@gmail.com>
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

# 27 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 29 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 30 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 31 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 32 "/app/src/vtkhdf_mb_file_type.F90.fypp"

module vtkhdf_mb_file_type

  use,intrinsic :: iso_fortran_env
  use vtkhdf_h5_c_binding
  use vtkhdf_h5
  use vtkhdf_ug_type
  use vtkhdf_ctx_type
  use vtkhdf_h5e_context_type
  use vtkhdf_version_param
  implicit none
  private

  public :: vtkhdf_version

  type, public :: vtkhdf_mb_file
    private
    type(vtkhdf_ctx) :: ctx
    integer(hid_t) :: file_id=H5I_INVALID_HID, vtk_id=H5I_INVALID_HID, ass_id=H5I_INVALID_HID
    integer :: next_bid = 0
    type(pdc_block), pointer :: blocks => null()
    type(pdc_block), pointer :: last_block => null()
    logical :: is_temporal = .false., write_time_step_called = .false.
  contains
    procedure :: create
    procedure :: close
    procedure :: add_block
# 60 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 61 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_block_mesh => write_block_mesh_real32_int32
    procedure, private :: write_block_mesh_real32_int32
# 61 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_block_mesh => write_block_mesh_real32_int64
    procedure, private :: write_block_mesh_real32_int64
# 64 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 60 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 61 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_block_mesh => write_block_mesh_real64_int32
    procedure, private :: write_block_mesh_real64_int32
# 61 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_block_mesh => write_block_mesh_real64_int64
    procedure, private :: write_block_mesh_real64_int64
# 64 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 65 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    procedure :: write_time_step
# 67 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 68 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 69 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_int32_scalar
    procedure, private :: write_cell_data_int32_scalar
    generic :: register_temporal_cell_data => register_temporal_cell_data_int32_scalar
    procedure, private :: register_temporal_cell_data_int32_scalar
    generic :: write_temporal_cell_data => write_temporal_cell_data_int32_scalar
    procedure, private :: write_temporal_cell_data_int32_scalar
# 69 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_int32_vector
    procedure, private :: write_cell_data_int32_vector
    generic :: register_temporal_cell_data => register_temporal_cell_data_int32_vector
    procedure, private :: register_temporal_cell_data_int32_vector
    generic :: write_temporal_cell_data => write_temporal_cell_data_int32_vector
    procedure, private :: write_temporal_cell_data_int32_vector
# 76 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 68 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 69 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_point_data => write_point_data_int32_scalar
    procedure, private :: write_point_data_int32_scalar
    generic :: register_temporal_point_data => register_temporal_point_data_int32_scalar
    procedure, private :: register_temporal_point_data_int32_scalar
    generic :: write_temporal_point_data => write_temporal_point_data_int32_scalar
    procedure, private :: write_temporal_point_data_int32_scalar
# 69 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_point_data => write_point_data_int32_vector
    procedure, private :: write_point_data_int32_vector
    generic :: register_temporal_point_data => register_temporal_point_data_int32_vector
    procedure, private :: register_temporal_point_data_int32_vector
    generic :: write_temporal_point_data => write_temporal_point_data_int32_vector
    procedure, private :: write_temporal_point_data_int32_vector
# 76 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 77 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 67 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 68 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 69 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_int64_scalar
    procedure, private :: write_cell_data_int64_scalar
    generic :: register_temporal_cell_data => register_temporal_cell_data_int64_scalar
    procedure, private :: register_temporal_cell_data_int64_scalar
    generic :: write_temporal_cell_data => write_temporal_cell_data_int64_scalar
    procedure, private :: write_temporal_cell_data_int64_scalar
# 69 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_int64_vector
    procedure, private :: write_cell_data_int64_vector
    generic :: register_temporal_cell_data => register_temporal_cell_data_int64_vector
    procedure, private :: register_temporal_cell_data_int64_vector
    generic :: write_temporal_cell_data => write_temporal_cell_data_int64_vector
    procedure, private :: write_temporal_cell_data_int64_vector
# 76 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 68 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 69 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_point_data => write_point_data_int64_scalar
    procedure, private :: write_point_data_int64_scalar
    generic :: register_temporal_point_data => register_temporal_point_data_int64_scalar
    procedure, private :: register_temporal_point_data_int64_scalar
    generic :: write_temporal_point_data => write_temporal_point_data_int64_scalar
    procedure, private :: write_temporal_point_data_int64_scalar
# 69 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_point_data => write_point_data_int64_vector
    procedure, private :: write_point_data_int64_vector
    generic :: register_temporal_point_data => register_temporal_point_data_int64_vector
    procedure, private :: register_temporal_point_data_int64_vector
    generic :: write_temporal_point_data => write_temporal_point_data_int64_vector
    procedure, private :: write_temporal_point_data_int64_vector
# 76 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 77 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 67 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 68 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 69 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_real32_scalar
    procedure, private :: write_cell_data_real32_scalar
    generic :: register_temporal_cell_data => register_temporal_cell_data_real32_scalar
    procedure, private :: register_temporal_cell_data_real32_scalar
    generic :: write_temporal_cell_data => write_temporal_cell_data_real32_scalar
    procedure, private :: write_temporal_cell_data_real32_scalar
# 69 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_real32_vector
    procedure, private :: write_cell_data_real32_vector
    generic :: register_temporal_cell_data => register_temporal_cell_data_real32_vector
    procedure, private :: register_temporal_cell_data_real32_vector
    generic :: write_temporal_cell_data => write_temporal_cell_data_real32_vector
    procedure, private :: write_temporal_cell_data_real32_vector
# 76 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 68 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 69 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_point_data => write_point_data_real32_scalar
    procedure, private :: write_point_data_real32_scalar
    generic :: register_temporal_point_data => register_temporal_point_data_real32_scalar
    procedure, private :: register_temporal_point_data_real32_scalar
    generic :: write_temporal_point_data => write_temporal_point_data_real32_scalar
    procedure, private :: write_temporal_point_data_real32_scalar
# 69 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_point_data => write_point_data_real32_vector
    procedure, private :: write_point_data_real32_vector
    generic :: register_temporal_point_data => register_temporal_point_data_real32_vector
    procedure, private :: register_temporal_point_data_real32_vector
    generic :: write_temporal_point_data => write_temporal_point_data_real32_vector
    procedure, private :: write_temporal_point_data_real32_vector
# 76 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 77 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 67 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 68 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 69 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_real64_scalar
    procedure, private :: write_cell_data_real64_scalar
    generic :: register_temporal_cell_data => register_temporal_cell_data_real64_scalar
    procedure, private :: register_temporal_cell_data_real64_scalar
    generic :: write_temporal_cell_data => write_temporal_cell_data_real64_scalar
    procedure, private :: write_temporal_cell_data_real64_scalar
# 69 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_cell_data => write_cell_data_real64_vector
    procedure, private :: write_cell_data_real64_vector
    generic :: register_temporal_cell_data => register_temporal_cell_data_real64_vector
    procedure, private :: register_temporal_cell_data_real64_vector
    generic :: write_temporal_cell_data => write_temporal_cell_data_real64_vector
    procedure, private :: write_temporal_cell_data_real64_vector
# 76 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 68 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 69 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_point_data => write_point_data_real64_scalar
    procedure, private :: write_point_data_real64_scalar
    generic :: register_temporal_point_data => register_temporal_point_data_real64_scalar
    procedure, private :: register_temporal_point_data_real64_scalar
    generic :: write_temporal_point_data => write_temporal_point_data_real64_scalar
    procedure, private :: write_temporal_point_data_real64_scalar
# 69 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    generic :: write_point_data => write_point_data_real64_vector
    procedure, private :: write_point_data_real64_vector
    generic :: register_temporal_point_data => register_temporal_point_data_real64_vector
    procedure, private :: register_temporal_point_data_real64_vector
    generic :: write_temporal_point_data => write_temporal_point_data_real64_vector
    procedure, private :: write_temporal_point_data_real64_vector
# 76 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 77 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 78 "/app/src/vtkhdf_mb_file_type.F90.fypp"
    procedure, private :: get_block_ptr
    final :: vtkhdf_mb_file_delete
  end type

  type :: pdc_block
    character(:), allocatable :: name
    type(vtkhdf_ug) :: b
    type(pdc_block), pointer :: next => null()
  contains
    final :: pdc_block_delete
  end type

contains

  !! Finalizer for vtkhdf_mb_file objects. We free heap memory we own but avoid
  !! doing things that may require syncronization with other ranks (MPI/PHDF5)
  !! because where implicit finalization occurs it is not guaranteed to be
  !! collective or ordered with respect to other ranks. This can leak HDF5
  !! IDs and the MPI communicator, but that is unavoidable. Users should
  !! always use CLOSE to do a proper collective cleanup and close of the file.

  subroutine vtkhdf_mb_file_delete(this)
    type(vtkhdf_mb_file), intent(inout) :: this
    if (associated(this%blocks)) deallocate(this%blocks)
    this%ass_id  = H5I_INVALID_HID
    this%vtk_id  = H5I_INVALID_HID
    this%file_id = H5I_INVALID_HID
    this%next_bid = 0
    this%last_block => null()
  end subroutine

  recursive subroutine pdc_block_delete(this)
    type(pdc_block), intent(inout) :: this
    if (associated(this%next)) deallocate(this%next)
  end subroutine

  !! Cleanly closes H5 identifiers and the file, and default initializes.
  subroutine close(this)
    class(vtkhdf_mb_file), intent(inout) :: this
    integer :: ierr
    type(pdc_block), pointer :: p
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    p => this%blocks
    do while (associated(p))
      call p%b%close
      p => p%next
    end do
    if (H5Iis_valid(this%ass_id) > 0) ierr = H5Gclose(this%ass_id)
    if (H5Iis_valid(this%vtk_id) > 0) ierr = H5Gclose(this%vtk_id)
    if (H5Iis_valid(this%file_id) > 0) ierr = H5Fclose(this%file_id)
    call err_ctx%restore(0)
    call this%ctx%free  ! frees duplicated comm in MPI build
    call finalize(this) ! free local memory
  contains
    subroutine finalize(this)
      class(vtkhdf_mb_file), intent(out) :: this
    end subroutine
  end subroutine

#ifdef USE_MPI
  subroutine create(this, filename, comm, stat, errmsg)
#else
  subroutine create(this, filename, stat, errmsg)
#endif

    class(vtkhdf_mb_file), intent(out) :: this
    character(*), intent(in) :: filename
#ifdef USE_MPI
    integer, intent(in) :: comm
#endif
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: fapl, gcpl
    integer(c_int) :: flag
    integer :: ierr
    type(vtkhdf_h5e_context) :: err_ctx

    ! MPI build: duplicates comm, caches rank/size
    ! Serial build: no comm, sets rank=0, size=1
#ifdef USE_MPI
    call this%ctx%init(comm)
#else
    call this%ctx%init()
#endif

    call init_hdf5
    call err_ctx%capture

    !! Open the file (read/write), overwrite any existing file
    call this%ctx%get_fapl(fapl) ! configures for parallel I/O in MPI build
    this%file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl)
    ierr = H5Pclose(fapl)
    if (this%ctx%global_any(this%file_id < 0)) then
      stat = 1
      errmsg = 'failed to open file "' // filename // '"'
      if (this%file_id > 0) ierr = H5Fclose(this%file_id)
      this%file_id = H5I_INVALID_HID
      call this%ctx%free
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
    this%vtk_id = H5Gcreate(this%file_id, 'VTKHDF', gcpl_id=gcpl)
    INSIST(this%ctx%global_all(this%vtk_id >= 0))

    call h5_write_attr(this%ctx, this%vtk_id, 'Version', vtkhdf_version)

    !NB: We stick with the older MB type due to an issue with the modern PDC
    !type; see https://gitlab.kitware.com/vtk/vtk/-/issues/19902
    !call h5_write_attr(this%ctx, this%vtk_id, 'Type', 'PartitionedDataSetCollection', stat, errmsg)
    call h5_write_attr(this%ctx, this%vtk_id, 'Type', 'MultiBlockDataSet')

    !! Create the Assembly group
    this%ass_id = H5Gcreate(this%vtk_id, 'Assembly', gcpl_id=gcpl)
    ierr = H5Pclose(gcpl)
    INSIST(this%ctx%global_all(this%ass_id >= 0))

    stat = 0
    call err_ctx%restore(stat)

  end subroutine create


  subroutine add_block(this, name, stat, errmsg, is_temporal)

    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: name
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    logical, intent(in), optional :: is_temporal

    integer :: n, ierr
    type(pdc_block), pointer :: b, new
    type(vtkhdf_h5e_context) :: err_ctx
    logical :: temporal_flag

    INSIST(this%file_id >= 0)

    !! Ensure the specified NAME is valid
    if (name == '') then
      stat = 1
      errmsg = 'invalid empty block name'
      return
    end if

    n = scan(name, './')
    if (n /= 0) then
      stat = 1
      errmsg = 'invalid character "' // name(n:n) // '" in block name'
      return
    end if

    if (name == 'Assembly') then
      stat = 1
      errmsg = 'invalid block name "Assembly"'
      return
    end if

    b => this%blocks
    do while (associated(b))
      if (b%name == name) then
        stat = 1
        errmsg = 'block "' // name // '" already defined'
        return
      end if
      b => b%next
    end do

    ! If a block is temporal, it must be added before the
    ! timeline starts to ensure step-count consistency.
    temporal_flag = .false.
    if (present(is_temporal)) temporal_flag = is_temporal
    if (temporal_flag) then
      INSIST(.not.this%write_time_step_called)
      this%is_temporal = .true.
    end if

    !! Create the UG block group and its HDF5 hierarchy
    allocate(new)
    new%name = name
    new%next => this%blocks
    this%blocks => new

    call err_ctx%capture ! begin HDF5 calls

    call new%b%init(this%vtk_id, name, this%ctx, is_temporal)

    !NB: Unused for MultiBlockDataSet, but required for PartitionedDataSetCollection.
    call h5_write_attr(this%ctx, new%b%root_id, 'Index', this%next_bid)
    this%next_bid = this%next_bid + 1

    !! Create a softlink in Assembly group to the block group.
    ierr = H5Lcreate_soft('/VTKHDF/'//name, this%ass_id, name)
    INSIST(this%ctx%global_all(ierr >= 0))

    stat = 0
    call err_ctx%restore(stat)

  end subroutine add_block

  !! Write the UnstructuredGrid data for the specified block. The unstructured
  !! mesh is described in the conventional manner by the X, CNODE, and XCNODE
  !! arrays. The additional array TYPES specifies the VTK cell types. This
  !! procedure must be called for each of the blocks before any of the
  !! following procedures.
# 294 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 295 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_block_mesh_real32_int32(this, block_name, x, cnode, xcnode, types)

    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name
    real(real32), intent(in) :: x(:,:)
    integer(int32), intent(in) :: cnode(:), xcnode(:)
    integer(int8), intent(in) :: types(:)

    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx

    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_mesh(x, cnode, xcnode, types)
    call err_ctx%restore(0)

  end subroutine write_block_mesh_real32_int32

# 295 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_block_mesh_real32_int64(this, block_name, x, cnode, xcnode, types)

    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name
    real(real32), intent(in) :: x(:,:)
    integer(int64), intent(in) :: cnode(:), xcnode(:)
    integer(int8), intent(in) :: types(:)

    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx

    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_mesh(x, cnode, xcnode, types)
    call err_ctx%restore(0)

  end subroutine write_block_mesh_real32_int64

# 315 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 294 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 295 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_block_mesh_real64_int32(this, block_name, x, cnode, xcnode, types)

    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name
    real(real64), intent(in) :: x(:,:)
    integer(int32), intent(in) :: cnode(:), xcnode(:)
    integer(int8), intent(in) :: types(:)

    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx

    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_mesh(x, cnode, xcnode, types)
    call err_ctx%restore(0)

  end subroutine write_block_mesh_real64_int32

# 295 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_block_mesh_real64_int64(this, block_name, x, cnode, xcnode, types)

    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name
    real(real64), intent(in) :: x(:,:)
    integer(int64), intent(in) :: cnode(:), xcnode(:)
    integer(int8), intent(in) :: types(:)

    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx

    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_mesh(x, cnode, xcnode, types)
    call err_ctx%restore(0)

  end subroutine write_block_mesh_real64_int64

# 315 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 316 "/app/src/vtkhdf_mb_file_type.F90.fypp"
  !! Writes the cell/point-based data ARRAY to a new named cell/point dataset
  !! for the specified mesh block. Scalar, vector, and tensor data are
  !! supported. In the case of a temporal block supporting time-dependent
  !! datasets, this dataset is static and not associated with any time step.
# 321 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 322 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 323 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_cell_data_int32_scalar(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int32), intent(in) :: array(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 323 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_cell_data_int32_vector(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int32), intent(in) :: array(:,:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 336 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 322 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 323 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_cell_data_int64_scalar(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int64), intent(in) :: array(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 323 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_cell_data_int64_vector(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int64), intent(in) :: array(:,:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 336 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 322 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 323 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_cell_data_real32_scalar(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real32), intent(in) :: array(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 323 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_cell_data_real32_vector(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real32), intent(in) :: array(:,:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 336 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 322 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 323 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_cell_data_real64_scalar(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real64), intent(in) :: array(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 323 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_cell_data_real64_vector(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real64), intent(in) :: array(:,:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 336 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 337 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 321 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 322 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 323 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_point_data_int32_scalar(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int32), intent(in) :: array(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 323 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_point_data_int32_vector(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int32), intent(in) :: array(:,:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 336 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 322 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 323 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_point_data_int64_scalar(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int64), intent(in) :: array(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 323 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_point_data_int64_vector(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int64), intent(in) :: array(:,:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 336 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 322 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 323 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_point_data_real32_scalar(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real32), intent(in) :: array(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 323 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_point_data_real32_vector(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real32), intent(in) :: array(:,:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 336 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 322 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 323 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_point_data_real64_scalar(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real64), intent(in) :: array(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 323 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_point_data_real64_vector(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real64), intent(in) :: array(:,:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 336 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 337 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 338 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  !! Register the specified NAME as a time-dependent cell/point dataset for
  !! the specified mesh block. This writes no data, but only configures some
  !! necessary internal metadata. The MOLD array argument shall have the same
  !! type, kind, and rank as the actual dataset, and the same extent in all
  !! but the last dimension, whose extent is ignored, but the array values
  !! themselves are not accessed. Scalar, vector, and tensor-valued mesh
  !! data are supported (rank-1, 2, and 3 MOLD).
# 347 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 348 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 349 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine register_temporal_cell_data_int32_scalar(this, block_name, name, mold)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int32), intent(in) :: mold
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%register_temporal_cell_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 349 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine register_temporal_cell_data_int32_vector(this, block_name, name, mold)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int32), intent(in) :: mold(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%register_temporal_cell_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 362 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 348 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 349 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine register_temporal_cell_data_int64_scalar(this, block_name, name, mold)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int64), intent(in) :: mold
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%register_temporal_cell_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 349 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine register_temporal_cell_data_int64_vector(this, block_name, name, mold)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int64), intent(in) :: mold(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%register_temporal_cell_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 362 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 348 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 349 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine register_temporal_cell_data_real32_scalar(this, block_name, name, mold)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real32), intent(in) :: mold
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%register_temporal_cell_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 349 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine register_temporal_cell_data_real32_vector(this, block_name, name, mold)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real32), intent(in) :: mold(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%register_temporal_cell_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 362 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 348 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 349 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine register_temporal_cell_data_real64_scalar(this, block_name, name, mold)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real64), intent(in) :: mold
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%register_temporal_cell_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 349 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine register_temporal_cell_data_real64_vector(this, block_name, name, mold)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real64), intent(in) :: mold(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%register_temporal_cell_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 362 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 363 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 347 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 348 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 349 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine register_temporal_point_data_int32_scalar(this, block_name, name, mold)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int32), intent(in) :: mold
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%register_temporal_point_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 349 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine register_temporal_point_data_int32_vector(this, block_name, name, mold)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int32), intent(in) :: mold(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%register_temporal_point_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 362 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 348 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 349 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine register_temporal_point_data_int64_scalar(this, block_name, name, mold)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int64), intent(in) :: mold
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%register_temporal_point_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 349 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine register_temporal_point_data_int64_vector(this, block_name, name, mold)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int64), intent(in) :: mold(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%register_temporal_point_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 362 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 348 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 349 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine register_temporal_point_data_real32_scalar(this, block_name, name, mold)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real32), intent(in) :: mold
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%register_temporal_point_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 349 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine register_temporal_point_data_real32_vector(this, block_name, name, mold)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real32), intent(in) :: mold(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%register_temporal_point_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 362 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 348 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 349 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine register_temporal_point_data_real64_scalar(this, block_name, name, mold)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real64), intent(in) :: mold
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%register_temporal_point_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 349 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine register_temporal_point_data_real64_vector(this, block_name, name, mold)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real64), intent(in) :: mold(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%register_temporal_point_data(name, mold)
    call err_ctx%restore(0)
  end subroutine
# 362 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 363 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 364 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  !! Mark the start of a new time step with time value TIME. Subsequent output
  !! of time-dependent datasets will be associated with this time step.

  subroutine write_time_step(this, time)
    class(vtkhdf_mb_file), intent(inout) :: this
    real(real64), intent(in) :: time
    type(pdc_block), pointer :: b
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    INSIST(this%is_temporal)
    this%write_time_step_called = .true.
    b => this%blocks
    do while (associated(b))
      if (b%b%nsteps >= 0) call b%b%write_time_step(time)
      b => b%next
    end do
    call err_ctx%restore(0)
  end subroutine

  !! Writes the cell/point-based data ARRAY to the named time-dependent
  !! dataset for the specified mesh block. The data is associated with the
  !! current time step.
# 388 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 389 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 390 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_temporal_cell_data_int32_scalar(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int32), intent(in) :: array(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_temporal_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 390 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_temporal_cell_data_int32_vector(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int32), intent(in) :: array(:,:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_temporal_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 403 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 389 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 390 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_temporal_cell_data_int64_scalar(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int64), intent(in) :: array(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_temporal_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 390 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_temporal_cell_data_int64_vector(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int64), intent(in) :: array(:,:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_temporal_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 403 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 389 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 390 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_temporal_cell_data_real32_scalar(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real32), intent(in) :: array(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_temporal_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 390 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_temporal_cell_data_real32_vector(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real32), intent(in) :: array(:,:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_temporal_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 403 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 389 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 390 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_temporal_cell_data_real64_scalar(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real64), intent(in) :: array(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_temporal_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 390 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_temporal_cell_data_real64_vector(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real64), intent(in) :: array(:,:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_temporal_cell_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 403 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 404 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 388 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 389 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 390 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_temporal_point_data_int32_scalar(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int32), intent(in) :: array(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_temporal_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 390 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_temporal_point_data_int32_vector(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int32), intent(in) :: array(:,:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_temporal_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 403 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 389 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 390 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_temporal_point_data_int64_scalar(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int64), intent(in) :: array(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_temporal_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 390 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_temporal_point_data_int64_vector(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    integer(int64), intent(in) :: array(:,:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_temporal_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 403 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 389 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 390 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_temporal_point_data_real32_scalar(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real32), intent(in) :: array(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_temporal_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 390 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_temporal_point_data_real32_vector(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real32), intent(in) :: array(:,:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_temporal_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 403 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 389 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 390 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_temporal_point_data_real64_scalar(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real64), intent(in) :: array(:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_temporal_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 390 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  subroutine write_temporal_point_data_real64_vector(this, block_name, name, array)
    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real64), intent(in) :: array(:,:)
    type(vtkhdf_ug), pointer :: bptr
    type(vtkhdf_h5e_context) :: err_ctx
    call err_ctx%capture
    call this%get_block_ptr(block_name, bptr)
    call bptr%write_temporal_point_data(name, array)
    call err_ctx%restore(0)
  end subroutine
# 403 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 404 "/app/src/vtkhdf_mb_file_type.F90.fypp"
# 405 "/app/src/vtkhdf_mb_file_type.F90.fypp"

  !! This private procedure returns a pointer to the named block.
  subroutine get_block_ptr(this, name, bptr)

    class(vtkhdf_mb_file), intent(inout) :: this
    character(*), intent(in) :: name
    type(vtkhdf_ug), pointer, intent(out) :: bptr

    type(pdc_block), pointer :: b

    if (associated(this%last_block)) then
      if (this%last_block%name == name) then
        bptr => this%last_block%b
        return
      end if
    end if

    b => this%blocks
    do while (associated(b))
      if (b%name == name) exit
      b => b%next
    end do

    INSIST(this%ctx%global_all(associated(b)))
    this%last_block => b
    bptr => b%b

  end subroutine get_block_ptr

end module vtkhdf_mb_file_type
