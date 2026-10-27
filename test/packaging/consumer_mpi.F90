program consumer_mpi
#ifdef USE_MPI_F08
  use mpi_f08
#else
  use mpi
#endif
  use vtkhdf_ug_file_type
  use vtkhdf_image_file_type
  use vtkhdf_file_modes, only: VTKHDF_FIXED_MESH
  implicit none
  type(vtkhdf_ug_file) :: file
  type(vtkhdf_image_file) :: image
  type(vtkhdf_cell_data_handle) :: variable
  integer :: selected(2,3), rank
  real :: values(1,1,1)
  integer :: ierr
  character(:), allocatable :: errmsg
  call MPI_Init(ierr)
  call file%create('dummy.vtkhdf', MPI_COMM_WORLD, ierr, errmsg)
  call file%close
  selected(1,:) = 0
  selected(2,:) = [0,0,0]
  values = 2
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  if (rank /= 0) selected(2,1) = -1
  call image%create('consumer_image.vtkhdf',MPI_COMM_WORLD,ierr,errmsg,mode=VTKHDF_FIXED_MESH)
  if (ierr /= 0) error stop errmsg
  call image%write_mesh(reshape([0,1,0,1,0,0],[2,3]),[0.,0.,0.],[1.,1.,1.])
  variable = image%register_temporal_cell_data('value',0.0)
  call image%start_time_step(0.0)
  call image%write_cell_data(variable,values,[0,0,0],selected)
  call image%close()
  call MPI_Finalize(ierr)
  print *, 'consumer_mpi: link OK'
end program
