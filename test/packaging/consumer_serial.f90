program consumer_serial
  use vtkhdf_ug_file_type
  use vtkhdf_image_file_type
  use vtkhdf_file_modes, only: VTKHDF_FIXED_MESH
  implicit none
  integer :: ierr
  character(:), allocatable :: errmsg
  type(vtkhdf_ug_file) :: file
  type(vtkhdf_image_file) :: image
  type(vtkhdf_cell_data_handle) :: variable
  integer :: selected(2,3), rank
  real :: values(1,1,1)
  call file%create('dummy.vtkhdf', ierr, errmsg)
  call file%close
  selected(1,:) = 0
  selected(2,:) = [0,0,0]
  values = 2
  rank = 0
  call image%create('consumer_image.vtkhdf',ierr,errmsg,mode=VTKHDF_FIXED_MESH)
  if (ierr /= 0) error stop errmsg
  call image%write_mesh(reshape([0,1,0,1,0,0],[2,3]),[0.,0.,0.],[1.,1.,1.])
  variable = image%register_temporal_cell_data('value',0.0)
  call image%start_time_step(0.0)
  call image%write_cell_data(variable,values,[0,0,0],selected)
  call image%close()
  print *, 'consumer_serial: link OK'
end program
