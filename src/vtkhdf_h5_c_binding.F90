!!
!! VTKHDF_H5_C_BINDING
!!
!! Private bindings to a relatively small subset of HDF5 C functions that
!! are used by the VTKHDF modules.
!!
!! Copyright (c) 2026 Neil Carlson <neil.n.carlson@gmail.com>
!! SPDX-License-Identifier: BSD-2-Clause
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! This module provides temporary C bindings to a subset of the HDF5 API.
!! Ideally it will be replaced by the HDF5 Fortran interface once current
!! build issues with the NAG compiler (HDF5 1.14.3) are resolved.
!!
!! The HDF5_INIT subroutine must be called before any other procedure in this
!! module. It initializes several read-only module variables (type and property
!! identifiers) whose values can only be determined at runtime as copies of
!! corresponding identifiers in the HDF5 C library.
!!
!! Creating these copies has the side effect of initializing the HDF5 C library.
!! Unlike the Fortran interface, which performs this initialization explicitly
!! via h5open_f, the C library initializes itself on the first HDF5 call.
!! HDF5_INIT ensures that this initialization occurs before any other bindings
!! are used.
!!

module vtkhdf_h5_c_binding

  use iso_c_binding
  implicit none
  public

  !! Header file constants from from H5Ipublic.h
  integer, parameter :: hid_t = c_int64_t
  integer(c_int), parameter :: H5I_INVALID_HID = -1

  !! Header file constants from H5public.h
  integer, parameter :: hsize_t = c_int64_t ! really uint64_t

  !! Header file constants from H5Spublic.h
  integer(c_int),   parameter :: H5S_ALL = 0
  integer(hsize_t), parameter :: H5S_UNLIMITED = -1 ! UINT64_MAX

  !! Header file constants from H5Fpublic.h
  integer(c_int), parameter :: H5F_ACC_RDONLY = int(Z'0000', c_int)
  integer(c_int), parameter :: H5F_ACC_RDWR   = int(Z'0001', c_int)
  integer(c_int), parameter :: H5F_ACC_TRUNC  = int(Z'0002', c_int)
  integer(c_int), parameter :: H5F_ACC_EXCL   = int(Z'0004', c_int)

  !! Header file constants from H5Ppublic.h
  integer(hid_t), parameter :: H5P_DEFAULT = 0

  !! Header file constants from H5Epublic.h
  integer(hid_t), parameter :: H5E_DEFAULT = 0

  !! Object IDs that are run-time *copies* of objects on the C side.
  !! These need to be initialized by a call to init_hdf5.
  integer(hid_t), protected :: H5T_NATIVE_UINT8
  integer(hid_t), protected :: H5T_NATIVE_INT32
  integer(hid_t), protected :: H5T_NATIVE_INT64
  integer(hid_t), protected :: H5T_NATIVE_FLOAT
  integer(hid_t), protected :: H5T_NATIVE_DOUBLE
  integer(hid_t), protected :: H5T_NATIVE_CHARACTER
  integer(hid_t), protected :: H5P_DATASET_CREATE
  integer(hid_t), protected :: H5P_GROUP_CREATE
  integer(hid_t), protected :: H5P_CRT_ORDER_TRACKED
  integer(hid_t), protected :: H5P_CRT_ORDER_INDEXED
  integer(hid_t), protected :: H5P_DATASET_XFER
  integer(hid_t), protected :: H5P_FILE_ACCESS

  enum, bind(c) ! H5S_class_t from H5Spublic.h
    enumerator :: H5S_SCALAR = 0, H5S_SIMPLE, H5S_NULL
  end enum

  enum, bind(c) ! H5S_seloper_t from H5Spublic.h
    enumerator :: H5S_SELECT_SET = 0
  end enum

  enum, bind(c) ! H5T_str_t from H5Tpublic.h
    enumerator :: H5T_STR_NULLTERM = 0, H5T_STR_NULLPAD, H5T_STR_SPACEPAD
  end enum

  enum, bind(c) ! H5FD_mpio_xfer_t from H5FDmpi.h
    enumerator :: H5FD_MPIO_INDEPENDENT = 0 ! Use independent I/O access
    enumerator :: H5FD_MPIO_COLLECTIVE      ! Use collective I/O access
  end enum

  !!!! H5F functions that can be used as-is !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface
    function H5Fclose(file_id) result(h5err) bind(c,name='H5Fclose')
      import :: hid_t, c_int
      integer(hid_t), value :: file_id
      integer(c_int) :: h5err
    end function

    function H5Fget_access_plist(file_id) &
        result(plist_id) bind(c,name='H5Fget_access_plist')
      import :: hid_t
      integer(hid_t), value :: file_id
      integer(hid_t) :: plist_id
    end function
  end interface

  public :: H5Fclose, H5Fget_access_plist

  !!!! H5G functions that can be used as-is !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface
    function H5Gclose(grp_id) result(h5err) bind(c,name='H5Gclose')
      import :: hid_t, c_int
      integer(hid_t), value :: grp_id
      integer(c_int) :: h5err
    end function
  end interface

  public :: H5Gclose
  public :: H5Gcreate ! module procedure

  !!!! H5I functions that can be used as-is !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface
    function H5Iget_file_id(id) result(file_id) bind(c,name='H5Iget_file_id')
      import :: hid_t
      integer(hid_t), value :: id
      integer(hid_t) :: file_id
    end function

    function H5Iis_valid(id) result(htri) bind(c,name='H5Iis_valid')
      import :: hid_t, c_int
      integer(hid_t), value :: id
      integer(c_int) :: htri ! >0, true; ==0, false; <0, failure
    end function
  end interface

  public :: H5Iget_file_id, H5Iis_valid

  !!!! H5A functions that can be used as-is !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface
    function H5Aclose(attr_id) result(h5err) bind(c,name='H5Aclose')
      import :: hid_t, c_int
      integer(hid_t), value :: attr_id
      integer(c_int) :: h5err
    end function

    function H5Awrite(attr_id, type_id, buf) result(h5err) bind(c,name='H5Awrite')
      import :: hid_t, c_int
      integer(hid_t), value :: attr_id, type_id
      type(*), intent(in) :: buf(*)
      integer(c_int) :: h5err
    end function

    function H5Aget_type(attr_id) result(type_id) bind(c,name='H5Aget_type')
      import :: hid_t
      integer(hid_t), value :: attr_id
      integer(hid_t) :: type_id
    end function

    function H5Aget_space(attr_id) result(space_id) bind(c,name='H5Aget_space')
      import :: hid_t
      integer(hid_t), value :: attr_id
      integer(hid_t) :: space_id
    end function
  end interface

  public :: H5Aclose, H5Awrite, H5Aget_type, H5Aget_space
  public :: H5Acreate, H5Aopen, H5Aexists

  interface H5Acreate
    module procedure H5Acreate, H5Acreate_by_name
  end interface

  interface H5Aopen
    module procedure H5Aopen, H5Aopen_by_name
  end interface

  !!!! H5S functions that can be used as-is !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface
!    function H5Screate(type) result(space_id) bind(c,name='H5Screate')
!      import :: c_int, hid_t
!      integer(c_int), value :: type
!      integer(hid_t) :: space_id
!    end function

    function H5Sget_simple_extent_ndims(space_id) &
        result(ndims) bind(c,name='H5Sget_simple_extent_ndims')
      import :: hid_t, c_int
      integer(hid_t), value :: space_id
      integer(c_int) :: ndims
    end function

    function H5Sget_simple_extent_dims(space_id, dims, maxdims) &
        result(ndim) bind(c,name='H5Sget_simple_extent_dims')
      import :: hid_t, hsize_t, c_int
      integer(hid_t), value :: space_id
      integer(hsize_t), optional :: dims(*), maxdims(*)
      integer(c_int) :: ndim
    end function

    function H5Sselect_elements(space_id, op, num_elem, coord) &
        result(h5err) bind(c,name='H5Sselect_elements')
      import :: hid_t, c_int, c_size_t, hsize_t
      integer(hid_t), value :: space_id
      integer(c_int), value :: op
      integer(c_size_t), value :: num_elem
      integer(hsize_t), intent(in) :: coord(*)
      integer(c_int) :: h5err
    end function

    function H5Sselect_none(space_id) result(h5err) bind(c,name='H5Sselect_none')
      import :: hid_t, c_int
      integer(hid_t), value :: space_id
      integer(c_int) :: h5err
    end function

    function H5Sclose(space_id) result(h5err) bind(c,name='H5Sclose')
      import :: hid_t, c_int
      integer(hid_t), value :: space_id
      integer(c_int) :: h5err
    end function
  end interface

  public :: H5Sget_simple_extent_ndims, H5Sget_simple_extent_dims, H5Sselect_elements, H5Sclose
  public :: H5Screate, H5Sselect_hyperslab ! module procedures

  interface H5Screate
    module procedure H5Screate, H5Screate_array
  end interface

  !!!! H5T functions that can be used as-is !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface
    function H5Tcopy(type_id) result(copy_id) bind(c,name='H5Tcopy')
      import :: hid_t
      integer(hid_t), value :: type_id
      integer(hid_t) :: copy_id
    end function

    function H5Tset_size(type_id, size) result(h5err) bind(c,name='H5Tset_size')
      import :: hid_t, c_size_t, c_int
      integer(hid_t), value :: type_id
      integer(c_size_t), value :: size
      integer(c_int) :: h5err
    end function

    function H5Tset_strpad(type_id, strpad) result(h5err) bind(c,name='H5Tset_strpad')
      import :: hid_t, c_int
      integer(hid_t), value :: type_id
      integer(c_int), value :: strpad
      integer(c_int) :: h5err
    end function

    function H5Tclose(type_id) result(h5err) bind(c,name='H5Tclose')
      import :: hid_t, c_int
      integer(hid_t), value :: type_id
      integer(c_int) :: h5err
    end function

    function H5Tequal(type1_id, type2_id) result(htri) bind(c,name='H5Tequal')
      import :: hid_t, c_int
      integer(hid_t), value :: type1_id, type2_id
      integer(c_int) :: htri ! >0, true; ==0, false; <0, failure
    end function
  end interface

  public :: H5Tcopy, H5Tset_size, H5Tset_strpad, H5Tclose, H5Tequal

  !!!! H5D functions that can be used as-is !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface
    function H5Dget_space(dset_id) result(space_id) bind(c,name='H5Dget_space')
      import :: hid_t
      integer(hid_t), value :: dset_id
      integer(hid_t) :: space_id
    end function

    function H5Dset_extent(dset_id, size) result(h5err) bind(c,name='H5Dset_extent')
      import :: hid_t, hsize_t, c_int
      integer(hid_t), value :: dset_id
      integer(hsize_t), intent(in) :: size(*)
      integer(c_int) :: h5err
    end function

    function H5Dclose(dset_id) result(h5err) bind(c,name='H5Dclose')
      import :: hid_t, c_int
      integer(hid_t), value :: dset_id
      integer(c_int) :: h5err
    end function

    function H5Dget_type(dset_id) result(type_id) bind(c,name='H5Dget_type')
      import :: hid_t
      integer(hid_t), value :: dset_id
      integer(hid_t) :: type_id
    end function
  end interface

  public :: H5Dget_space, H5Dset_extent, H5Dclose, H5Dget_type
  public :: H5Dopen, H5Dcreate, H5Dwrite  ! module procedures

  !!!! H5P functions that can be used as-is !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface
    function H5Pcreate(cls_id) result(prop_id) bind(c,name='H5Pcreate')
      import :: hid_t
      integer(hid_t), value :: cls_id
      integer(hid_t) :: prop_id
    end function

    function H5Pset_chunk(plist_id, ndims, dim) result(hdferr) bind(c,name='H5Pset_chunk')
      import :: hid_t, hsize_t, c_int
      integer(hid_t), value :: plist_id
      integer(c_int), value :: ndims
      integer(hsize_t), intent(in) :: dim(*)
      integer(c_int) :: hdferr
    end function

    function H5Pset_link_creation_order(plist_id, crt_order_flags) &
        result(hdferr) bind(c,name='H5Pset_link_creation_order')
      import :: hid_t, c_int
      integer(hid_t), value :: plist_id
      integer(c_int), value :: crt_order_flags ! really unsigned
      integer(c_int) :: hdferr
    end function

    function H5Pset_all_coll_metadata_ops(plist_id, is_collective) &
        result(hdferr) bind(c,name='H5Pset_all_coll_metadata_ops')
      import :: hid_t, c_int, c_bool
      integer(hid_t), value :: plist_id
      logical(c_bool), value :: is_collective
      integer(c_int) :: hdferr
    end function

    function H5Pset_coll_metadata_write(plist_id, is_collective) &
        result(hdferr) bind(c,name='H5Pset_coll_metadata_write')
      import :: hid_t, c_int, c_bool
      integer(hid_t), value :: plist_id
      logical(c_bool), value :: is_collective
      integer(c_int) :: hdferr
    end function

    function H5Pset_dxpl_mpio(dxpl_id, xfer_mode) &
        result(hdferr) bind(c,name='H5Pset_dxpl_mpio')
      import hid_t, c_int
      integer(hid_t), value :: dxpl_id
      integer(c_int), value :: xfer_mode
      integer(c_int) :: hdferr
    end function

    function H5Pclose(prp_id) result(h5err) bind(c,name='H5Pclose')
      import :: hid_t, c_int
      integer(hid_t), value :: prp_id
      integer(c_int) :: h5err
    end function
  end interface

  public :: H5Pcreate, H5Pclose, H5Pset_chunk, H5Pset_link_creation_order, &
            H5Pset_all_coll_metadata_ops, H5Pset_coll_metadata_write, H5Pset_dxpl_mpio

  interface ! to wrapper functions that take Fortran comm instead of C comm
    function H5Pset_fapl_mpio(fapl_id, comm) &
        result(hdferr) bind(c,name='H5Pset_fapl_mpio_Fcomm')
      import :: hid_t, c_int
      integer(hid_t), value :: fapl_id
      integer, value :: comm
      integer(c_int) :: hdferr
    end function

    function H5Pget_fapl_mpio(fapl, comm) &
        result(hdferr) bind(c,name='H5Pget_fapl_mpio_Fcomm')
      import :: hid_t, c_int
      integer(hid_t), value :: fapl
      integer :: comm
      integer(c_int) :: hdferr
    end function
  end interface

  public :: H5Pset_fapl_mpio, H5Pget_fapl_mpio

  !!!! H5L functions

  public :: H5Lexists, H5Lcreate_soft, H5Lcreate_hard ! module procedures

  !!!! H5E functions

  interface
    function H5Eset_auto2(estack_id, func, client_data) &
        result(hdferr) bind(c,name='H5Eset_auto2')
      import :: hid_t, c_int, c_ptr, c_funptr
      integer(hid_t), value :: estack_id
      type(c_funptr), value :: func
      type(c_ptr), value :: client_data
      integer(c_int) :: hdferr
    end function

    function H5Eget_auto2(estack_id, func, client_data) &
        result(hdferr) bind(c,name='H5Eget_auto2')
      import :: hid_t, c_int, c_ptr, c_funptr
      integer(hid_t), value :: estack_id
      type(c_funptr) :: func           ! Received by reference
      type(c_ptr)    :: client_data    ! Received by reference
      integer(c_int) :: hdferr
    end function

    function H5Eclear2(estack_id) result(hdferr) bind(c,name='H5Eclear2')
      import :: hid_t, c_int
      integer(hid_t), value :: estack_id
      integer(c_int) :: hdferr
    end function

    function H5Eprint2(estack_id, stream) result(hdferr) bind(c, name="H5Eprint2")
      import :: hid_t, c_int, c_ptr
      integer(hid_t), value :: estack_id
      type(c_ptr), value :: stream
      integer(c_int) :: hdferr
    end function
  end interface

  public :: H5Eset_auto2, H5Eget_auto2, H5Eclear2, H5Eprint2

contains

  subroutine init_hdf5
    logical, save :: initialized = .false.

    interface
      function H5P_DATASET_CREATE_value() result(type_id) &
          bind(c,name='H5P_DATASET_CREATE_value')
        import :: hid_t
        integer(hid_t) :: type_id
      end function
      function H5P_GROUP_CREATE_value() result(type_id) &
          bind(c,name='H5P_GROUP_CREATE_value')
        import :: hid_t
        integer(hid_t) :: type_id
      end function
      function H5P_CRT_ORDER_TRACKED_value() result(type_id) &
          bind(c,name='H5P_CRT_ORDER_TRACKED_value')
        import :: hid_t
        integer(hid_t) :: type_id
      end function
      function H5P_CRT_ORDER_INDEXED_value() result(type_id) &
          bind(c,name='H5P_CRT_ORDER_INDEXED_value')
        import :: hid_t
        integer(hid_t) :: type_id
      end function
      function H5P_DATASET_XFER_value() result(type_id) &
          bind(c,name='H5P_DATASET_XFER_value')
        import :: hid_t
        integer(hid_t) :: type_id
      end function
      function H5P_FILE_ACCESS_value() result(type_id) &
          bind(c,name='H5P_FILE_ACCESS_value')
        import :: hid_t
        integer(hid_t) :: type_id
      end function
      function H5T_NATIVE_UINT8_value() result(type_id) &
          bind(c,name='H5T_NATIVE_UINT8_value')
        import :: hid_t
        integer(hid_t) :: type_id
      end function
      function H5T_NATIVE_INT32_value() result(type_id) &
          bind(c,name='H5T_NATIVE_INT32_value')
        import :: hid_t
        integer(hid_t) :: type_id
      end function
      function H5T_NATIVE_INT64_value() result(type_id) &
          bind(c,name='H5T_NATIVE_INT64_value')
        import :: hid_t
        integer(hid_t) :: type_id
      end function
      function H5T_NATIVE_FLOAT_value() result(type_id) &
          bind(c,name='H5T_NATIVE_FLOAT_value')
        import :: hid_t
        integer(hid_t) :: type_id
      end function
      function H5T_NATIVE_DOUBLE_value() result(type_id) &
          bind(c,name='H5T_NATIVE_DOUBLE_value')
        import :: hid_t
        integer(hid_t) :: type_id
      end function
      function H5T_NATIVE_CHARACTER_value() result(type_id) &
          bind(c,name='H5T_NATIVE_CHARACTER_value')
        import :: hid_t
        integer(hid_t) :: type_id
      end function
    end interface

    if (initialized) return

    H5P_DATASET_CREATE = H5P_DATASET_CREATE_value()
    H5P_GROUP_CREATE = H5P_GROUP_CREATE_value()
    H5P_CRT_ORDER_TRACKED = H5P_CRT_ORDER_TRACKED_value()
    H5P_CRT_ORDER_INDEXED = H5P_CRT_ORDER_INDEXED_value()
    H5P_DATASET_XFER = H5P_DATASET_XFER_value()
    H5P_FILE_ACCESS = H5P_FILE_ACCESS_value()
    H5T_NATIVE_UINT8 = H5T_NATIVE_UINT8_value()
    H5T_NATIVE_INT32 = H5T_NATIVE_INT32_value()
    H5T_NATIVE_INT64 = H5T_NATIVE_INT64_value()
    H5T_NATIVE_FLOAT = H5T_NATIVE_FLOAT_value()
    H5T_NATIVE_DOUBLE = H5T_NATIVE_DOUBLE_value()
    H5T_NATIVE_CHARACTER = H5T_NATIVE_CHARACTER_value()

    initialized = .true.
  end subroutine

  function H5Fcreate(filename, flags, fcpl_id, fapl_id) result(file_id)
    character(*), intent(in) :: filename
    integer(c_int), intent(in) :: flags
    integer(hid_t), intent(in), optional :: fcpl_id, fapl_id
    integer(hid_t) :: file_id
    interface
      function H5Fcreate_c(filename, flags, fcpl_id, fapl_id) &
          result(file_id) bind(c,name="H5Fcreate")
        import :: hid_t, c_int
        character, intent(in) :: filename(*)
        integer(c_int), value :: flags ! really unsigned
        integer(hid_t), value :: fcpl_id, fapl_id
        integer(hid_t) :: file_id
      end function
    end interface
    integer(hid_t) :: fcpl_id_, fapl_id_
    fcpl_id_ = H5P_DEFAULT; if (present(fcpl_id)) fcpl_id_ = fcpl_id
    fapl_id_ = H5P_DEFAULT; if (present(fapl_id)) fapl_id_ = fapl_id
    file_id = H5Fcreate_c(filename//c_null_char, flags, fcpl_id_, fapl_id_)
  end function

  function H5Gcreate(loc_id, name, lcpl_id, gcpl_id, gapl_id) result(grp_id)
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(hid_t), intent(in), optional :: lcpl_id, gcpl_id, gapl_id
    integer(hid_t) :: grp_id
    interface
      function H5Gcreate2_c(loc_id, name, lcpl_id, gcpl_id, gapl_id) &
          result(grp_id) bind(c,name="H5Gcreate2")
        import :: hid_t
        integer(hid_t), value :: loc_id
        character, intent(in) :: name(*)
        integer(hid_t), value :: lcpl_id, gcpl_id, gapl_id
        integer(hid_t) :: grp_id
      end function
    end interface
    integer(hid_t) :: lcpl_id_, gcpl_id_, gapl_id_
    lcpl_id_ = H5P_DEFAULT; if (present(lcpl_id)) lcpl_id_ = lcpl_id
    gcpl_id_ = H5P_DEFAULT; if (present(gcpl_id)) gcpl_id_ = gcpl_id
    gapl_id_ = H5P_DEFAULT; if (present(gapl_id)) gapl_id_ = gapl_id
    grp_id = H5Gcreate2_c(loc_id, name//c_null_char, lcpl_id_, gcpl_id_, gapl_id_)
  end function

  function H5Acreate(loc_id, attr_name, type_id, space_id, acpl_id, aapl_id) result(attr_id)
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: attr_name
    integer(hid_t), intent(in) :: type_id, space_id
    integer(hid_t), intent(in), optional :: acpl_id, aapl_id
    integer(hid_t) :: attr_id
    interface
      function H5Acreate2_c(loc_id, attr_name, type_id, space_id, acpl_id, aapl_id) &
          result(attr_id) bind(c,name='H5Acreate2')
        import :: hid_t
        integer(hid_t), value :: loc_id
        character, intent(in) :: attr_name(*)
        integer(hid_t), value :: type_id, space_id, acpl_id, aapl_id
        integer(hid_t) :: attr_id
      end function
    end interface
    integer(hid_t) :: acpl_id_, aapl_id_
    acpl_id_ = H5P_DEFAULT; if (present(acpl_id)) acpl_id_ = acpl_id
    aapl_id_ = H5P_DEFAULT; if (present(aapl_id)) aapl_id_ = aapl_id
    attr_id = H5Acreate2_c(loc_id, attr_name//c_null_char, type_id, space_id, acpl_id_, aapl_id_)
  end function

  function H5Acreate_by_name(loc_id, obj_name, attr_name, type_id, space_id, acpl_id, aapl_id, lapl_id) result(attr_id)
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: obj_name, attr_name
    integer(hid_t), intent(in) :: type_id, space_id
    integer(hid_t), intent(in), optional :: acpl_id, aapl_id, lapl_id
    integer(hid_t) :: attr_id
    interface
      function H5Acreate_by_name_c(loc_id, obj_name, attr_name, type_id, space_id, acpl_id, aapl_id, lapd_id) &
          result(attr_id) bind(c,name='H5Acreate_by_name')
        import :: hid_t
        integer(hid_t), value :: loc_id
        character, intent(in) :: obj_name(*), attr_name(*)
        integer(hid_t), value :: type_id, space_id, acpl_id, aapl_id, lapd_id
        integer(hid_t) :: attr_id
      end function
    end interface
    integer(hid_t) :: acpl_id_, aapl_id_, lapl_id_
    acpl_id_ = H5P_DEFAULT; if (present(acpl_id)) acpl_id_ = acpl_id
    aapl_id_ = H5P_DEFAULT; if (present(aapl_id)) aapl_id_ = aapl_id
    lapl_id_ = H5P_DEFAULT; if (present(lapl_id)) lapl_id_ = lapl_id
    attr_id = H5Acreate_by_name_c(loc_id, obj_name//c_null_char, attr_name//c_null_char, type_id, space_id, &
        acpl_id_, aapl_id_, lapl_id_)
  end function

  function H5Aopen(obj_id, attr_name) result(attr_id)
    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: attr_name
    integer(hid_t) :: attr_id
    interface
      function H5Aopen_c(obj_id, attr_name, aapl_id) &
          result(attr_id) bind(c,name='H5Aopen')
        import :: hid_t
        integer(hid_t), value :: obj_id
        character, intent(in) :: attr_name(*)
        integer(hid_t), value :: aapl_id
        integer(hid_t) :: attr_id
      end function
    end interface
    attr_id = H5Aopen_c(obj_id, attr_name//c_null_char, H5P_DEFAULT)
  end function

  function H5Aopen_by_name(loc_id, obj_name, attr_name, lapl_id) result(attr_id)
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: obj_name, attr_name
    integer(hid_t), intent(in), optional :: lapl_id
    integer(hid_t) :: attr_id
    interface
      function H5Aopen_by_name_c(loc_id, obj_name, attr_name, aapl_id, lapl_id) &
          result(attr_id) bind(c,name='H5Aopen_by_name')
        import :: hid_t
        integer(hid_t), value :: loc_id
        character, intent(in) :: obj_name(*), attr_name(*)
        integer(hid_t), value :: aapl_id, lapl_id
        integer(hid_t) :: attr_id
      end function
    end interface
    integer(hid_t) :: lapl_id_
    lapl_id_ = H5P_DEFAULT; if (present(lapl_id)) lapl_id_ = lapl_id
    attr_id = H5Aopen_by_name_c(loc_id, obj_name//c_null_char, attr_name//c_null_char, H5P_DEFAULT, lapl_id_)
  end function

  integer function H5Aexists(obj_id, attr_name) result(res)
    integer(hid_t), intent(in) :: obj_id
    character(*), intent(in) :: attr_name
    interface
      function H5Aexists_c(obj_id, attr_name) result(res) bind(c,name='H5Aexists')
        import :: hid_t, c_int
        integer(hid_t), value :: obj_id
        character, intent(in) :: attr_name(*)
        integer(c_int) :: res
      end function
    end interface
    res = H5Aexists_c(obj_id, attr_name//c_null_char) ! >0 yes, 0 no, <0 error
  end function

  function H5Screate(type) result(space_id)
    integer(c_int), intent(in) :: type
    integer(hid_t) :: space_id
    interface
      function H5Screate_c(type) result(space_id) bind(c,name='H5Screate')
        import :: c_int, hid_t
        integer(c_int), value :: type
        integer(hid_t) :: space_id
      end function
    end interface
    space_id = H5Screate_c(type)
  end function

  function H5Screate_array(dims, maxdims) result(space_id)
    integer(hsize_t), intent(in) :: dims(:)
    integer(hsize_t), intent(in), optional :: maxdims(:)
    integer(hid_t) :: space_id
    interface
      function H5Screate_simple_c(rank, dims, maxdims) result(space_id) bind(c,name='H5Screate_simple')
        import :: c_int, hsize_t, hid_t
        integer(c_int), value :: rank
        integer(hsize_t), intent(in) :: dims(*)
        integer(hsize_t), intent(in), optional :: maxdims(*)
        integer(hid_t) :: space_id
      end function
    end interface
    space_id = H5Screate_simple_c(size(dims), dims, maxdims)
  end function

  function H5Sselect_hyperslab(space_id, op, start, count, stride, block) result(h5err)
    integer(hid_t), intent(in) :: space_id
    integer(c_int), intent(in) :: op
    integer(hsize_t), intent(in) :: start(:), count(:)
    integer(hsize_t), intent(in), optional :: stride(:), block(:)
    integer(c_int) :: h5err
    interface
      function H5Sselect_hyperslab_c(space_id, op, start, stride, count, block) &
          result(h5err) bind(c,name='H5Sselect_hyperslab')
        import :: hid_t, c_int, hsize_t
        integer(hid_t), value :: space_id
        integer(c_int), value :: op
        integer(hsize_t), intent(in) :: start(*), count(*)
        integer(hsize_t), intent(in), optional :: stride(*), block(*)
        integer(c_int) :: h5err
      end function
    end interface
    h5err = H5Sselect_hyperslab_c(space_id, op, start, stride, count, block)
  end function

  function H5Dopen(loc_id, name, dapl_id) result(dset_id)
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(hid_t), intent(in), optional :: dapl_id
    integer(hid_t) :: dset_id
    interface
      function H5Dopen2_c(loc_id, name, dapl) result(dset_id) bind(c,name='H5Dopen2')
        import :: hid_t
        integer(hid_t), value :: loc_id, dapl
        character, intent(in) :: name(*)
        integer(hid_t) :: dset_id
      end function
    end interface
    integer(hid_t) :: dapl_id_
    dapl_id_ = H5P_DEFAULT; if (present(dapl_id)) dapl_id_ = dapl_id
    dset_id = H5Dopen2_c(loc_id, name//c_null_char, dapl_id_)
  end function

  function H5Dwrite(dset_id, mem_type_id, buf, mem_space_id, file_space_id, dxpl_id) result(h5err)
    integer(hid_t), intent(in) :: dset_id, mem_type_id
    type(*), intent(in) :: buf(*)
    integer(hid_t), intent(in), optional :: mem_space_id, file_space_id, dxpl_id
    integer(c_int) :: h5err
    interface
      function H5Dwrite_c(dset_id, mem_type_id, mem_space_id, file_space_id, dxpl_id, buf) &
          result(h5err) bind(c,name='H5Dwrite')
        import :: hid_t, c_int
        integer(hid_t), value :: dset_id, mem_type_id, mem_space_id, file_space_id, dxpl_id
        type(*), intent(in) :: buf(*)
        integer(c_int) :: h5err
      end function
    end interface
    integer(hid_t) :: mem_space_id_, file_space_id_, dxpl_id_
    mem_space_id_ = H5S_ALL; if (present(mem_space_id)) mem_space_id_ = mem_space_id
    file_space_id_ = H5S_ALL; if (present(file_space_id)) file_space_id_ = file_space_id
    dxpl_id_ = H5P_DEFAULT; if (present(dxpl_id)) dxpl_id_ = dxpl_id
    h5err = H5Dwrite_c(dset_id, mem_type_id, mem_space_id_, file_space_id_, dxpl_id_, buf)
  end function

  function H5Dcreate(loc_id, name, type_id, space_id, lcpl_id, dcpl_id, dapl_id) result(dset_id)
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(hid_t), intent(in) :: type_id, space_id
    integer(hid_t), intent(in), optional :: lcpl_id, dcpl_id, dapl_id
    integer(hid_t) :: dset_id
    interface
      function H5Dcreate2_c(loc_id, name, type_id, space_id, lcpl_id, dcpl_id, dapl_id) &
          result(dset_id) bind(c,name='H5Dcreate2')
        import :: hid_t
        integer(hid_t), value :: loc_id
        character, intent(in) :: name(*)
        integer(hid_t), value :: type_id, space_id, lcpl_id, dcpl_id, dapl_id
        integer(hid_t) :: dset_id
      end function
    end interface
    integer(hid_t) :: lcpl_id_, dcpl_id_, dapl_id_
    lcpl_id_ = H5P_DEFAULT; if (present(lcpl_id)) lcpl_id_ = lcpl_id
    dcpl_id_ = H5P_DEFAULT; if (present(dcpl_id)) dcpl_id_ = dcpl_id
    dapl_id_ = H5P_DEFAULT; if (present(dapl_id)) dapl_id_ = dapl_id
    dset_id = H5Dcreate2_c(loc_id, name//c_null_char, type_id, space_id, lcpl_id_, dcpl_id_, dapl_id_)
  end function

  function H5Lexists(loc_id, name, lapl_id) result(exists)
    integer(hid_t), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer(hid_t), intent(in), optional :: lapl_id
    integer(c_int) :: exists
    interface
      function H5Lexists_c(loc_id, name, lapl_id) result(htri) bind(c,name='H5Lexists')
        import :: hid_t, c_int
        integer(hid_t), value :: loc_id
        character, intent(in) :: name(*)
        integer(hid_t), value :: lapl_id
        integer(c_int) :: htri
      end function
    end interface
    integer(hid_t) :: lapl_id_
    lapl_id_ = H5P_DEFAULT; if (present(lapl_id)) lapl_id_ = lapl_id
    exists = H5Lexists_c(loc_id, name//c_null_char, lapl_id_)
  end function

  function H5Lcreate_soft(link_target, link_loc_id, link_name, lcpl_id, lapl_id) result(h5err)
    character(*), intent(in) :: link_target, link_name
    integer(hid_t), intent(in) :: link_loc_id
    integer(hid_t), intent(in), optional :: lcpl_id, lapl_id
    integer(c_int) :: h5err
    interface
      function H5Lcreate_soft_c(link_target, link_loc_id, link_name, lcpl_id, lapl_id) &
          result(h5err) bind(c,name='H5Lcreate_soft')
        import :: hid_t, c_int
        character, intent(in) :: link_target(*), link_name(*)
        integer(hid_t), value :: link_loc_id, lcpl_id, lapl_id
        integer(c_int) :: h5err
      end function
    end interface
    integer(hid_t) :: lcpl_id_, lapl_id_
    lcpl_id_ = H5P_DEFAULT; if (present(lcpl_id)) lcpl_id_ = lcpl_id
    lapl_id_ = H5P_DEFAULT; if (present(lapl_id)) lapl_id_ = lapl_id
    h5err = H5Lcreate_soft_c(link_target//c_null_char, link_loc_id, link_name//c_null_char, &
                             lcpl_id_, lapl_id_)
  end function

  function H5Lcreate_hard(cur_loc, cur_name, dst_loc, dst_name, lcpl_id, lapl_id) result(h5err)
    integer(hid_t), intent(in) :: cur_loc, dst_loc
    character(*), intent(in) :: cur_name, dst_name
    integer(hid_t), intent(in), optional :: lcpl_id, lapl_id
    integer(c_int) :: h5err
    interface
      function H5Lcreate_hard_c(cur_loc, cur_name, dst_loc, dst_name, lcpl_id, lapl_id) &
          result(h5err) bind(c,name='H5Lcreate_hard')
        import :: hid_t, c_int
        integer(hid_t), value :: cur_loc, dst_loc
        character, intent(in) :: cur_name(*), dst_name(*)
        integer(hid_t), value :: lcpl_id, lapl_id
        integer(c_int) :: h5err
      end function
    end interface
    integer(hid_t) :: lcpl_id_, lapl_id_
    lcpl_id_ = H5P_DEFAULT; if (present(lcpl_id)) lcpl_id_ = lcpl_id
    lapl_id_ = H5P_DEFAULT; if (present(lapl_id)) lapl_id_ = lapl_id
    h5err = H5Lcreate_hard_c(cur_loc, cur_name//c_null_char, dst_loc, dst_name//c_null_char, &
                             lcpl_id_, lapl_id_)
  end function

end module vtkhdf_h5_c_binding
