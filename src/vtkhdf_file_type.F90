!!
!! VTKHDF_FILE_TYPE
!!
!! This module defines a derived type for exporting mesh-based solution data
!! to a VTKHDF format file that can be read by the ParaView visualizaton tool.
!! The format uses HDF5 for on-disk storage.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! January 2026
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! The up-to-date specification for VTKHDF is at
!! https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html
!!
!! This module was written for version 2.5 of the format.
!!
!! This module uses the "MultiBlockDataSet" type of format. It supports static
!! and time dependent data (point and cell), but assumes a single static mesh.
!! The "Assembly" dataset hierarchy is limited to a flat collection
!! of UnstructuredGrid blocks.
!!

#include "f90_assert.fpp"

module vtkhdf_file_type

  use,intrinsic :: iso_fortran_env
  use vtkhdf_h5_c_binding
  use vtkhdf_hl_h5
  use vtkhdf_ug_type
  implicit none
  private

  type, public :: vtkhdf_file
    private
    integer(hid_t) :: file_id = -1, vtk_id=-1, ass_id=-1
    integer :: next_bid = 0
    type(pdc_block), pointer :: blocks => null()
  contains
    procedure :: create, create_block
    procedure :: write_block_mesh
    procedure :: write_time_step
    generic :: write_cell_data  => write_cell_data_real64, write_cell_data_int32
    generic :: write_point_data => write_point_data_real64
    generic :: register_temporal_cell_data  => register_temporal_cell_data_real64
    generic :: register_temporal_point_data => register_temporal_point_data_real64
    generic :: write_temporal_cell_data  => write_temporal_cell_data_real64
    generic :: write_temporal_point_data => write_temporal_point_data_real64
    procedure, private :: write_cell_data_real64, write_cell_data_int32
    procedure, private :: write_point_data_real64
    procedure, private :: register_temporal_cell_data_real64
    procedure, private :: register_temporal_point_data_real64
    procedure, private :: write_temporal_cell_data_real64
    procedure, private :: write_temporal_point_data_real64
    procedure :: close
    procedure, private :: get_block_ptr
    !final :: vtkhdf_file_delete
  end type

  type :: pdc_block
    character(:), allocatable :: name
    type(vtkhdf_ug) :: b
    type(pdc_block), pointer :: next => null()
  contains
    final :: pdc_block_delete
  end type

  !! The VTK cell types that are relevant to Truchas.
  !! NB: The local node ordering for a VTK wedge cell may differ from Truchas
  integer(int8), parameter, public :: VTK_TETRA = 10
  integer(int8), parameter, public :: VTK_HEXAHEDRON = 12
  integer(int8), parameter, public :: VTK_WEDGE = 13
  integer(int8), parameter, public :: VTK_PYRAMID = 14

  integer, parameter :: vtkhdf_version(*) = [2,5]

contains

  subroutine vtkhdf_file_delete(this)
    type(vtkhdf_file), intent(inout) :: this
    if (associated(this%blocks)) deallocate(this%blocks)
  end subroutine

  recursive subroutine pdc_block_delete(this)
    type(pdc_block), intent(inout) :: this
    call this%b%close
    if (associated(this%next)) deallocate(this%next)
  end subroutine

  subroutine create(this, filename, comm, stat, errmsg)

    use,intrinsic :: iso_c_binding, only: c_bool

    class(vtkhdf_file), intent(out) :: this
    character(*), intent(in) :: filename
    integer, intent(in) :: comm
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer(hid_t) :: fapl, crt_prop
    integer(c_int) :: flag
    integer :: istat ! ignored status result

    call init_hdf5

    fapl = H5Pcreate(H5P_FILE_ACCESS)
    stat = H5Pset_fapl_mpio(fapl, comm)
    stat = H5Pset_all_coll_metadata_ops(fapl, is_collective=.true._c_bool)
    stat = H5Pset_coll_metadata_write(fapl, is_collective=.true._c_bool)
    this%file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl)
    if (this%file_id < 0) then
      stat = 1
      errmsg = 'h5fcreate error'  !TODO: refine msg
      return
    end if
    istat = H5Pclose(fapl)
    INSIST(istat == 0)

    crt_prop = H5Pcreate(H5P_GROUP_CREATE)
    flag = ior(H5P_CRT_ORDER_TRACKED, H5P_CRT_ORDER_INDEXED)
    stat = H5Pset_link_creation_order(crt_prop, flag)
    INSIST(stat == 0)

    this%vtk_id = H5Gcreate(this%file_id, 'VTKHDF', gcpl_id=crt_prop)
    INSIST(this%vtk_id > 0)

    call h5_write_attr(this%vtk_id, 'Version', vtkhdf_version, stat, errmsg)
    INSIST(stat == 0)
    !NB: We stick with the older MB type due to an issue with the modern PDC
    !type; see https://gitlab.kitware.com/vtk/vtk/-/issues/19902
    !call h5_write_attr(this%vtk_id, 'Type', 'PartitionedDataSetCollection', stat, errmsg)
    call h5_write_attr(this%vtk_id, 'Type', 'MultiBlockDataSet', stat, errmsg)
    INSIST(stat == 0)

    this%ass_id = H5Gcreate(this%vtk_id, 'Assembly', gcpl_id=crt_prop)
    INSIST(this%ass_id > 0)
    istat = H5Pclose(crt_prop)

  end subroutine create

  subroutine create_block(this, name, stat, errmsg, temporal)

    class(vtkhdf_file), intent(inout) :: this
    character(*), intent(in) :: name
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    logical, intent(in), optional :: temporal

    integer :: n
    type(pdc_block), pointer :: b, new

    if (name == '') then
      stat = -1
      errmsg = 'invalid empty block name'
      return
    end if

    n = scan(name, './')
    if (n /= 0) then
      stat = -1
      errmsg = 'invalid character "' // name(n:n) // '" in block name'
      return
    end if

    if (name == 'Assembly') then
      stat = -1
      errmsg = 'invalid block name "Assembly"'
      return
    end if

    !! Check that the named block has not already been defined.
    b => this%blocks
    do while (associated(b))
      if (b%name == name) then
        stat = -1
        errmsg = 'block "' // name // '" already defined'
        return
      end if
      b => b%next
    end do

    allocate(new)
    new%name = name
    new%next => this%blocks
    this%blocks => new

    call new%b%init(this%vtk_id, name, stat, errmsg, temporal)
    if (stat /= 0) return
    call h5_write_attr(new%b%root_id, 'Index', this%next_bid, stat, errmsg)
    INSIST(stat == 0)
    this%next_bid = this%next_bid + 1

    !! Create an Assembly group softlink to the block.
    stat = H5Lcreate_soft('/VTKHDF/'//name, this%ass_id, name)
    INSIST(stat == 0)

  end subroutine

  subroutine close(this)
    class(vtkhdf_file), intent(inout) :: this
    integer :: istat
    if (associated(this%blocks)) deallocate(this%blocks)
    if (this%ass_id > 0) istat = H5Gclose(this%ass_id)
    if (this%vtk_id > 0) istat = H5Gclose(this%vtk_id)
    if (this%file_id > 0) istat = H5Fclose(this%file_id)
    !call default_initialize(this)
  contains
    subroutine default_initialize(this)
      class(vtkhdf_file), intent(out) :: this
    end subroutine
  end subroutine

  !! Write the UnstructuredGrid data for the specified block. The unstructured
  !! mesh is described in the conventional manner by the X, CNODE, and XCNODE
  !! arrays. The additional array TYPES provides an unambiguous specification
  !! of the cell types that one might otherwise infer from the XCNODE array.
  !! This procedure must be called for each of the blocks before any of the
  !! following procedures.

  subroutine write_block_mesh(this, block_name, x, cnode, xcnode, types, stat, errmsg)

    class(vtkhdf_file), intent(inout) :: this
    character(*), intent(in) :: block_name
    real(real64), intent(in) :: x(:,:)
    integer, intent(in) :: cnode(:), xcnode(:)
    integer(int8), intent(in) :: types(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(vtkhdf_ug), pointer :: bptr

    call this%get_block_ptr(block_name, bptr, stat, errmsg)
    if (stat /= 0) return
    call bptr%write_mesh(x, cnode, xcnode, types, stat, errmsg)

  end subroutine write_block_mesh

  !! Writes the cell-based data ARRAY to a new named cell dataset for the
  !! specified mesh block. Scalar, vector, and tensor cell-based data are
  !! supported. In the case of a temporal block supporting time-dependent
  !! datasets, this dataset is static and not associated with any time step.

  subroutine write_cell_data_real64(this, block_name, name, array, stat, errmsg)
    class(vtkhdf_file), intent(in) :: this
    character(*), intent(in) :: block_name, name
    real(real64), intent(in) :: array(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    type(vtkhdf_ug), pointer :: bptr
    call this%get_block_ptr(block_name, bptr, stat, errmsg)
    if (stat /= 0) return
    call bptr%write_cell_data_real64(name, array, stat, errmsg)
  end subroutine

  subroutine write_cell_data_int32(this, block_name, name, array, stat, errmsg)
    class(vtkhdf_file), intent(in) :: this
    character(*), intent(in) :: block_name, name
    integer(int32), intent(in) :: array(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    type(vtkhdf_ug), pointer :: bptr
    call this%get_block_ptr(block_name, bptr, stat, errmsg)
    if (stat /= 0) return
    call bptr%write_cell_data_int32(name, array, stat, errmsg)
  end subroutine

  !! Writes the point-based data ARRAY to a new named cell dataset for the
  !! specified mesh block. Scalar, vector, and tensor cell-based data are
  !! supported. In the case of a temporal block supporting time-dependent
  !! datasets, this dataset is static and not associated with any time step.

  subroutine write_point_data_real64(this, block_name, name, array, stat, errmsg)
    class(vtkhdf_file), intent(in) :: this
    character(*), intent(in) :: block_name, name
    real(real64), intent(in) :: array(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    type(vtkhdf_ug), pointer :: bptr
    call this%get_block_ptr(block_name, bptr, stat, errmsg)
    if (stat /= 0) return
    call bptr%write_point_data_real64(name, array, stat, errmsg)
  end subroutine

  !! Register the specified NAME as a time-dependent point dataset for the
  !! specified mesh block. This writes no data, but only configures some
  !! necessary internal metadata. The MOLD array argument shall have the same
  !! type, kind, and rank as the actual dataset, and the same extent in all
  !! but the last dimension, whose extent is ignored, but the array values
  !! themselves are not accessed. Scalar, vector, and tensor-valued mesh
  !! data are supported (rank-1, 2, and 3 MOLD).

  subroutine register_temporal_point_data_real64(this, block_name, name, mold, stat, errmsg)
    class(vtkhdf_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real64), intent(in) :: mold(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    type(vtkhdf_ug), pointer :: bptr
    call this%get_block_ptr(block_name, bptr, stat, errmsg)
    if (stat /= 0) return
    call bptr%register_temporal_point_data_real64(name, mold, stat, errmsg)
  end subroutine

  !! Register the specified NAME as a time-dependent cell dataset for the
  !! specified mesh block. This writes no data, but only configures some
  !! necessary internal metadata. The MOLD array argument shall have the same
  !! type, kind, and rank as the actual dataset, and the same extent in all
  !! but the last dimension, whose extent is ignored, but the array values
  !! themselves are not accessed. Scalar, vector, and tensor-valued mesh
  !! data are supported (rank-1, 2, and 3 MOLD).

  subroutine register_temporal_cell_data_real64(this, block_name, name, mold, stat, errmsg)
    class(vtkhdf_file), intent(inout) :: this
    character(*), intent(in) :: block_name, name
    real(real64), intent(in) :: mold(..)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    type(vtkhdf_ug), pointer :: bptr
    call this%get_block_ptr(block_name, bptr, stat, errmsg)
    if (stat /= 0) return
    call bptr%register_temporal_cell_data_real64(name, mold, stat, errmsg)
  end subroutine

  !! Mark the start of a new time step with time value TIME. Subsequent output
  !! of time-dependent datasets will be associated with this time step.

  subroutine write_time_step(this, time)
    class(vtkhdf_file), intent(inout) :: this
    real(real64), intent(in) :: time
    type(pdc_block), pointer :: b
    b => this%blocks
    do while (associated(b))
      if (b%b%nsteps >= 0) call b%b%write_time_step(time)
      b => b%next
    end do
  end subroutine

  !! Write the point-based data ARRAY to the named time-dependent cell dataset
  !! for the specified mesh block. The data is associated with the current
  !! time step.

  subroutine write_temporal_point_data_real64(this, block_name, name, array, stat, errmsg)
    class(vtkhdf_file), intent(in) :: this
    character(*), intent(in) :: block_name, name
    real(real64), intent(in) :: array(..)
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg
    type(vtkhdf_ug), pointer :: bptr
    call this%get_block_ptr(block_name, bptr, stat, errmsg)
    if (stat /= 0) return
    call bptr%write_temporal_point_data_real64(name, array, stat, errmsg)
  end subroutine

  !! Write the cell-based data ARRAY to the named time-dependent cell dataset
  !! for the specified mesh block. The data is associated with the current
  !! time step.

  subroutine write_temporal_cell_data_real64(this, block_name, name, array, stat, errmsg)
    class(vtkhdf_file), intent(in) :: this
    character(*), intent(in) :: block_name, name
    real(real64), intent(in) :: array(..)
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg
    type(vtkhdf_ug), pointer :: bptr
    call this%get_block_ptr(block_name, bptr, stat, errmsg)
    if (stat /= 0) return
    call bptr%write_temporal_cell_data_real64(name, array, stat, errmsg)
  end subroutine

  !! This auxiliary procedure returns a pointer to the named block, or a null
  !! pointer if the block does not exist. If successful (block exists) STAT
  !! returns 0; otherwise it returns a non-zero value and ERRMSG returns an
  !! informative error message.

  subroutine get_block_ptr(this, name, bptr, stat, errmsg)

    class(vtkhdf_file), intent(in) :: this
    character(*), intent(in) :: name
    type(vtkhdf_ug), pointer, intent(out) :: bptr
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(pdc_block), pointer :: b

    b => this%blocks
    do while (associated(b))
      if (b%name == name) exit
      b => b%next
    end do

    stat = merge(0, -1, associated(b))
    if (stat == 0) then
      bptr => b%b
    else
      bptr => null()
      errmsg = 'no such block "' // name // '"'
    end if

  end subroutine get_block_ptr

end module vtkhdf_file_type
