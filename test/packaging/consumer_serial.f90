program consumer_serial
  use vtkhdf_ug_file_type
  implicit none
  integer :: ierr
  character(:), allocatable :: errmsg
  type(vtkhdf_ug_file) :: file
  call file%create('dummy.vtkhdf', ierr, errmsg)
  call file%close
  print *, 'consumer_serial: link OK'
end program
