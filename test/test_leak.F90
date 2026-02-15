program test_leak
  use vtkhdf_h5_c_binding
  implicit none
  integer(hid_t) :: val1, val2

  ! First call to init_hdf5
  call init_hdf5
  val1 = H5T_NATIVE_UINT8
  print *, "First call value: ", val1

  ! Second call to init_hdf5
  call init_hdf5
  val2 = H5T_NATIVE_UINT8
  print *, "Second call value: ", val2

  ! Check if values are the same
  if (val1 == val2) then
    print *, "SUCCESS: Handles are the same. Leak prevented."
  else
    print *, "FAILURE: Handles are different. Leak detected."
    stop 1
  end if

end program test_leak
