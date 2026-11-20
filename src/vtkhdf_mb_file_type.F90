!! MultiBlockDataSet file writer. Shared behavior belongs to the composite base.
!! SPDX-License-Identifier: BSD-2-Clause

module vtkhdf_mb_file_type

  use,intrinsic :: iso_fortran_env, only: int64
  use vtkhdf_composite_file_type
  implicit none
  private

  public :: vtkhdf_version, vtkhdf_unstructured_block
  public :: vtkhdf_image_block, vtkhdf_rectilinear_block, vtkhdf_structured_block
  public :: vtkhdf_hypertree_block
  public :: vtkhdf_cell_data_handle, vtkhdf_point_data_handle, vtkhdf_field_data_handle
  public :: VTKHDF_STATIC, VTKHDF_FIXED_MESH, VTKHDF_MOVING_MESH, VTKHDF_DYNAMIC_MESH

  type, extends(vtkhdf_composite_file), public :: vtkhdf_mb_file
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
  end type

contains

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
  end subroutine

#endif
#ifdef USE_MPI
  subroutine create_int(this, filename, comm, stat, errmsg, data_chunk_bytes)
    class(vtkhdf_mb_file), intent(out) :: this
    character(*), intent(in) :: filename
    integer, intent(in) :: comm
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer(int64), intent(in), optional :: data_chunk_bytes
    call this%internal_create(filename, comm, stat, errmsg, data_chunk_bytes, "MultiBlockDataSet")
  end subroutine
#else
  subroutine create(this, filename, stat, errmsg, data_chunk_bytes)
    class(vtkhdf_mb_file), intent(out) :: this
    character(*), intent(in) :: filename
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer(int64), intent(in), optional :: data_chunk_bytes
    call this%internal_create(filename, stat, errmsg, data_chunk_bytes, "MultiBlockDataSet")
  end subroutine
#endif

end module vtkhdf_mb_file_type
