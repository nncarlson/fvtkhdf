program vtkhdf_demo

  use,intrinsic :: iso_fortran_env, only: r8 => real64, int8
  use vtkhdf_file_type
  use mpi
  implicit none

  type(vtkhdf_file) :: vizfile
  integer :: istat, nproc, rank, stat, j
  character(:), allocatable :: errmsg
  character(7), allocatable :: name(:)
  !! 0-sized mesh
  real(r8) :: x0(3,0)
  integer  :: cnode0(0), xcnode0(1) = 1
  integer(int8) :: types0(0)
  !! Basic mesh unit
  real(r8), allocatable :: x(:,:), y(:,:)
  integer,  allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  real(r8), allocatable :: s(:), v(:,:) ! scalar and vector data arrays
  real(r8) :: s0(0), v0(3,0) ! 0-sized scalar and vector arrays; also molds for same

  call MPI_Init(istat)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, istat)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, istat)

  call vizfile%create('demo.vtkhdf', MPI_COMM_WORLD, stat, errmsg)
  if (stat /= 0) error stop errmsg

  !! The unstructured mesh data for a basic mesh unit. The full mesh will
  !! be a collection of non-overlapping shifts of this basic unit. Each
  !! rank has one of these, which is right-shifted proportional to the rank.
  call get_mesh1_data(x, cnode, xcnode, types)
  x(1,:) = x(1,:) + rank ! shift right

  !! Create the mesh blocks and export the mesh for each. Each mesh block
  !! consists of 2 or 3 basic mesh units from select ranks; other ranks
  !! contribute an empty mesh. Not all blocks (or any) need to be temporal.
  !!
  !! NB: all VTKHDF_FILE methods are collective and ranks not contributing
  !! any data need to be called with appropriate 0-sized data.

  allocate(name(0:nproc-1))
  do j = 0, nproc-1! Block names
    write(name(j),'("block",i2.2)') j
  end do

  y = x ! initial node coordinates
  do j = 0, nproc-1 ! block loop
    call vizfile%create_block(name(j), stat, errmsg, temporal=.true.)
    if (stat /= 0) error stop errmsg
    if (abs(rank-j) <= 1) then
      call vizfile%write_block_mesh(name(j), y, cnode, xcnode, types, stat, errmsg)
    else ! pass a 0-sized mesh
      call vizfile%write_block_mesh(name(j), x0, cnode0, xcnode0, types0, stat, errmsg)
    endif
    if (stat /= 0) error stop errmsg
    y(2,:) = y(2,:) + 1 ! everyone shifts up
  end do

  !!!! Register the data arrays that evolve with time.

  do j = 0, nproc-1
    call vizfile%register_temporal_cell_data(name(j), 'cell-scalar', s0, stat, errmsg)
    if (stat /= 0) error stop errmsg
    call vizfile%register_temporal_cell_data(name(j), 'cell-vector', v0, stat, errmsg)
    if (stat /= 0) error stop errmsg
    call vizfile%register_temporal_point_data(name(j), 'point-scalar', s0, stat, errmsg)
    if (stat /= 0) error stop errmsg
    call vizfile%register_temporal_point_data(name(j), 'point-vector', v0, stat, errmsg)
    if (stat /= 0) error stop errmsg
  end do

  !!!! Write the datasets for the first time step !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call vizfile%write_time_step(0.0_r8)

  y = x ! restore initial node coordinates
  do j = 0, nproc-1
    if (abs(rank-j) <= 1) then
      call get_scalar_cell_data(y, cnode, xcnode, s)
      call vizfile%write_temporal_cell_data(name(j), 'cell-scalar', s, stat, errmsg)
      call get_vector_cell_data(y, cnode, xcnode, v)
      call vizfile%write_temporal_cell_data(name(j), 'cell-vector', v, stat, errmsg)
      call get_scalar_point_data(y, s)
      call vizfile%write_temporal_point_data(name(j), 'point-scalar', s, stat, errmsg)
      call get_vector_point_data(y, v)
      call vizfile%write_temporal_point_data(name(j), 'point-vector', v, stat, errmsg)
    else ! pass 0-sized data
      call vizfile%write_temporal_cell_data(name(j), 'cell-scalar', s0, stat, errmsg)
      call vizfile%write_temporal_cell_data(name(j), 'cell-vector',  v0, stat, errmsg)
      call vizfile%write_temporal_point_data(name(j), 'point-scalar', s0, stat, errmsg)
      call vizfile%write_temporal_point_data(name(j), 'point-vector',  v0, stat, errmsg)
    end if
    y(2,:) = y(2,:) + 2 ! everyone shifts up
  end do

  !!!! Write the datasets for the second time step !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call vizfile%write_time_step(1.0_r8)

  y = x ! restore initial node coordinates
  do j = 0, nproc-1
    if (abs(rank-j) <= 1) then
      call get_scalar_cell_data(y, cnode, xcnode, s)
      call vizfile%write_temporal_cell_data(name(j), 'cell-scalar', s+1, stat, errmsg)
      call get_vector_cell_data(y, cnode, xcnode, v)
      call vizfile%write_temporal_cell_data(name(j), 'cell-vector', v+1, stat, errmsg)
      call get_scalar_point_data(y, s)
      call vizfile%write_temporal_point_data(name(j), 'point-scalar', s+1, stat, errmsg)
      call get_vector_point_data(y, v)
      call vizfile%write_temporal_point_data(name(j), 'point-vector', v+1, stat, errmsg)
    else ! pass 0-sized data
      call vizfile%write_temporal_cell_data(name(j), 'cell-scalar', s0, stat, errmsg)
      call vizfile%write_temporal_cell_data(name(j), 'cell-vector',  v0, stat, errmsg)
      call vizfile%write_temporal_point_data(name(j), 'point-scalar', s0, stat, errmsg)
      call vizfile%write_temporal_point_data(name(j), 'point-vector',  v0, stat, errmsg)
    end if
    y(2,:) = y(2,:) + 2 ! everyone shifts up
  end do

  !! At any point you can write a dataset that isn't time dependent, but its name must
  !! be unique from any other dataset, temporal or not, of the same type (cell or point).

  y = x ! restore initial node coordinates
  do j = 0, nproc-1
    if (abs(rank-j) <= 1) then
      call get_scalar_cell_data(y, cnode, xcnode, s)
      call vizfile%write_cell_data(name(j), 'static-cell-scalar', s, stat, errmsg)
      call get_vector_cell_data(y, cnode, xcnode, v)
      call vizfile%write_cell_data(name(j), 'static-cell-vector', v, stat, errmsg)
      call get_scalar_point_data(y, s)
      call vizfile%write_point_data(name(j), 'static-point-scalar', s, stat, errmsg)
      call get_vector_point_data(y, v)
      call vizfile%write_point_data(name(j), 'static-point-vector', v, stat, errmsg)
    else ! pass 0-sized data
      call vizfile%write_cell_data(name(j), 'static-cell-scalar', s0, stat, errmsg)
      call vizfile%write_cell_data(name(j), 'static-cell-vector',  v0, stat, errmsg)
      call vizfile%write_point_data(name(j), 'static-point-scalar', s0, stat, errmsg)
      call vizfile%write_point_data(name(j), 'static-point-vector',  v0, stat, errmsg)
    end if
    y(2,:) = y(2,:) + 2 ! everyone shifts up
  end do

  call vizfile%close
  call MPI_Finalize(istat)

contains

  ! A 5-tet subdivision of the squished unit cube.
  subroutine get_mesh1_data(x, cnode, xcnode, types)
    real(r8), allocatable, intent(out) :: x(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable :: types(:)
    x = reshape([0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1], shape=[3,8])
    ! distort to catch C/Fortran index ordering errors
    x(1,:) = 0.9_r8*x(1,:)
    x(2,:) = 0.7_r8*x(2,:)
    x(3,:) = 0.5_r8*x(3,:)
    cnode = [1,2,4,5, 2,3,4,7, 2,5,6,7, 4,5,7,8, 2,4,5,7]
    xcnode = [1,5,9,13,17,21]
    types = spread(VTK_TETRA, dim=1, ncopies=5)
  end subroutine

  ! A single squished unit hex cell.
  subroutine get_mesh2_data(x, cnode, xcnode, types)
    real(r8), allocatable, intent(out) :: x(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable :: types(:)
    x = reshape([0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1], shape=[3,8])
    ! distort to catch C/Fortran index ordering errors
    x(1,:) = 0.9_r8*x(1,:)
    x(2,:) = 0.7_r8*x(2,:)
    x(3,:) = 0.5_r8*x(3,:)
    cnode = [1, 2, 3, 4, 5, 6, 7, 8]
    xcnode = [1,9]
    types = [VTK_HEXAHEDRON]
  end subroutine

  ! norm of the node coordinates
  subroutine get_scalar_point_data(x, pdata)
    real(r8), intent(in) :: x(:,:)
    real(r8), allocatable, intent(out) :: pdata(:)
    integer :: j
    allocate(pdata(size(x,dim=2)))
    do j = 1, size(x,dim=2)
      pdata(j) = norm2(x(:,j))
    end do
  end subroutine

  ! node coordinates
  subroutine get_vector_point_data(x, pdata)
    real(r8), intent(in) :: x(:,:)
    real(r8), allocatable, intent(out) :: pdata(:,:)
    integer :: j
    pdata = x
  end subroutine

  ! norm of the cell centroid
  subroutine get_scalar_cell_data(x, cnode, xcnode, cdata)
    real(r8), intent(in) :: x(:,:)
    integer, intent(in) :: cnode(:), xcnode(:)
    real(r8), allocatable, intent(out) :: cdata(:)
    integer :: j
    allocate(cdata(size(xcnode)-1))
    do j = 1, size(cdata)
      associate(pid => cnode(xcnode(j):xcnode(j+1)-1))
        cdata(j) = norm2(sum(x(:,pid),dim=2)/size(pid))
      end associate
    end do
  end subroutine

  ! cell centroids
  subroutine get_vector_cell_data(x, cnode, xcnode, cdata)
    real(r8), intent(in) :: x(:,:)
    integer, intent(in) :: cnode(:), xcnode(:)
    real(r8), allocatable, intent(out) :: cdata(:,:)
    integer :: j
    allocate(cdata(size(x,dim=1),size(xcnode)-1))
    do j = 1, size(cdata,dim=2)
      associate(pid => cnode(xcnode(j):xcnode(j+1)-1))
        cdata(:,j) = sum(x(:,pid),dim=2)/size(pid)
      end associate
    end do
  end subroutine

end program
