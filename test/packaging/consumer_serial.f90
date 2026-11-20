program consumer_serial
  use,intrinsic :: iso_fortran_env, only: int8
  use vtkhdf_pdc_file_type, only: vtkhdf_pdc_file
  use vtkhdf_unstructured_file_type
  use vtkhdf_mb_file_type, only: vtkhdf_mb_file, vtkhdf_unstructured_block, vtkhdf_image_block, &
      vtkhdf_rectilinear_block, vtkhdf_structured_block, vtkhdf_hypertree_block
  use vtkhdf_rectilinear_file_type
  use vtkhdf_structured_file_type
  use vtkhdf_image_file_type
  use vtkhdf_file_modes, only: VTKHDF_FIXED_MESH
  implicit none
  integer :: ierr
  character(:), allocatable :: errmsg
  type(vtkhdf_unstructured_file) :: file
  type(vtkhdf_image_file) :: image
  type(vtkhdf_rectilinear_file) :: rect
  type(vtkhdf_structured_file) :: grid
  type(vtkhdf_mb_file) :: mb
  type(vtkhdf_pdc_file) :: pdc
  type(vtkhdf_unstructured_block), pointer :: block
  type(vtkhdf_image_block), pointer :: image_block
  type(vtkhdf_rectilinear_block), pointer :: rect_block
  type(vtkhdf_structured_block), pointer :: struct_block
  type(vtkhdf_hypertree_block), pointer :: hypertree_block
  real :: points(3,2,2,1)
  type(vtkhdf_cell_data_handle) :: variable
  integer :: selected(3,2), rank
  real :: values(1,1,1)
  call file%create('dummy.vtkhdf', ierr, errmsg)
  call file%close
  selected(:,1) = 0
  selected(:,2) = [0,0,0]
  values = 2
  rank = 0
  call image%create('consumer_image.vtkhdf',ierr,errmsg,mode=VTKHDF_FIXED_MESH)
  if (ierr /= 0) error stop errmsg
  call image%write_mesh(reshape([0,0,0,1,1,0],[3,2]),[0.,0.,0.],[1.,1.,1.])
  variable = image%register_temporal_cell_data('value',0.0)
  call image%start_time_step(0.0)
  call image%write_cell_data(variable,values,[0,0,0],selected)
  call image%close()
  call rect%create('consumer_rect.vtkhdf',ierr,errmsg)
  if (ierr /= 0) error stop errmsg
  call rect%write_mesh(reshape([0,0,0,1,1,0],[3,2]),[0.,1.],[0.,1.],[0.])
  call rect%close()
  call grid%create('consumer_grid.vtkhdf',ierr,errmsg)
  if (ierr /= 0) error stop errmsg
  call grid%write_mesh(reshape([0,0,0,1,1,0],[3,2]))
  points = reshape([0.,0.,0.,1.,0.,0.,0.,1.,0.,1.,1.,0.],[3,2,2,1])
  selected = reshape([0,0,0,1,1,0],[3,2])
  if (rank /= 0) selected(1,2) = -1
  call grid%write_mesh_geometry(points,[0,0,0],selected)
  call grid%close()

  call mb%create('consumer_mb.vtkhdf',ierr,errmsg)
  if (ierr /= 0) error stop errmsg
  call mb%add_block('mesh',block)
  call block%write_mesh(reshape(points,[3,4]),[1,2,4,3],[1,5],[9_int8])
  call block%write_cell_data('value',[2.])
  call mb%add_block('image',image_block)
  call image_block%write_mesh(reshape([0,0,0,1,1,0],[3,2]),[0.,0.,0.],[1.,1.,1.])
  call mb%add_block('rect',rect_block)
  call rect_block%write_mesh(reshape([0,0,0,1,1,0],[3,2]),[0.,1.],[0.,1.],[0.])
  call mb%add_block('struct',struct_block)
  call struct_block%write_mesh(reshape([0,0,0,1,1,0],[3,2]))
  call struct_block%write_mesh_geometry(points,[0,0,0],selected)
  call mb%add_block('hypertree',hypertree_block)
  call hypertree_block%write_mesh([2,1,1])
  call hypertree_block%write_mesh_geometry([0.,1.],[0.],[0.])
  if (rank == 0) then
    call hypertree_block%write_mesh_topology([1],[1,2],[.false.])
  else
    call hypertree_block%write_mesh_topology([integer::],[1],[logical::])
  end if
  call mb%close()

  call pdc%create('consumer_pdc.vtkhdf',ierr,errmsg)
  if (ierr /= 0) error stop errmsg
  call pdc%add_block('mesh',block)
  call block%write_mesh(reshape(points,[3,4]),[1,2,4,3],[1,5],[9_int8])
  call block%write_cell_data('value',[2.])
  call pdc%flush()
  call pdc%close()

  print *, 'consumer_serial: link OK'
end program
