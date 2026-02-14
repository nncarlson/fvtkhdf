program vtkhdf_mb_test

  use,intrinsic :: iso_fortran_env, only: r8 => real64, int8
  use vtkhdf_mb_file_type
  implicit none

  real(r8), allocatable :: x(:,:), scalar_cell_data(:), vector_cell_data(:,:)
  real(r8), allocatable :: scalar_point_data(:), vector_point_data(:,:)
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  character(:), allocatable :: errmsg
  integer :: stat

  type(vtkhdf_mb_file) :: vizfile

  call vizfile%create('mb_test.vtkhdf', stat, errmsg)
  if (stat /= 0) error stop errmsg

  !! The unstructured mesh data for a basic mesh unit. The full mesh will
  !! consist of two non-overlapping shifts of this basic unit.
  call get_mesh_data(x, cnode, xcnode, types)

  call vizfile%add_block('Block-A', stat, errmsg, is_temporal=.true.)
  if (stat /= 0) error stop errmsg

  call vizfile%write_block_mesh('Block-A', x, cnode, xcnode, types)

  call vizfile%add_block('Block-B', stat, errmsg, is_temporal=.true.)
  if (stat /= 0) error stop errmsg

  call vizfile%write_block_mesh('Block-B', x+1, cnode, xcnode, types)

  !! Register the datasets that evolve with time. At this stage the data arrays
  !! are only used to glean their types and shapes.

  associate (scalar_mold => 0.0_r8, vector_mold => [real(r8) :: 0, 0, 0])
    ! Block-A has time-dependent cell data
    call vizfile%register_temporal_cell_data('Block-A', 'cell-radius', scalar_mold)
    call vizfile%register_temporal_cell_data('Block-A', 'cell-velocity', vector_mold)
    ! Block-B has time-dependent point data
    call vizfile%register_temporal_point_data('Block-B', 'point-radius', scalar_mold)
    call vizfile%register_temporal_point_data('Block-B', 'point-velocity', vector_mold)
  end associate

  !! Generate some cell and point data to use for output
  call get_scalar_cell_data(x, cnode, xcnode, scalar_cell_data)
  call get_vector_cell_data(x, cnode, xcnode, vector_cell_data)

  call get_scalar_point_data(x+1, scalar_point_data)
  call get_vector_point_data(x+1, vector_point_data)

  !!!! Write the data for the first time step !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call vizfile%write_time_step(0.0_r8)

  call vizfile%write_temporal_cell_data('Block-A', 'cell-radius', scalar_cell_data)
  call vizfile%write_temporal_cell_data('Block-A', 'cell-velocity', vector_cell_data)
  call vizfile%write_temporal_point_data('Block-B', 'point-radius', scalar_point_data)
  call vizfile%write_temporal_point_data('Block-B', 'point-velocity', vector_point_data)

  !!!! Write the data for the second time step !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call vizfile%write_time_step(1.0_r8)

  call vizfile%write_temporal_cell_data('Block-A', 'cell-radius', scalar_cell_data+1)
  call vizfile%write_temporal_cell_data('Block-A', 'cell-velocity', vector_cell_data+1)
  call vizfile%write_temporal_point_data('Block-B', 'point-radius', scalar_point_data+1)
  call vizfile%write_temporal_point_data('Block-B', 'point-velocity', vector_point_data+1)

  !! At any point you can write a data that isn't time dependent, but its name must
  !! be unique from any other data temporal or not of the same type (cell or point).

  call vizfile%write_cell_data('Block-A', 'static-cell-scalar', -scalar_cell_data)
  call vizfile%write_cell_data('Block-A', 'static-cell-vector', -vector_cell_data)
  call vizfile%write_point_data('Block-B', 'static-point-scalar', -scalar_point_data)
  call vizfile%write_point_data('Block-B', 'static-point-vector', -vector_point_data)

  call vizfile%close

contains

  ! A 5-tet subdivision of a squished unit cube.
  subroutine get_mesh_data(x, cnode, xcnode, types)
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

  ! Point scalar is the magnitude of the node coordinate
  subroutine get_scalar_point_data(x, pdata)
    real(r8), intent(in) :: x(:,:)
    real(r8), allocatable, intent(out) :: pdata(:)
    integer :: j
    allocate(pdata(size(x,dim=2)))
    do j = 1, size(x,dim=2)
      pdata(j) = norm2(x(:,j))
    end do
  end subroutine

  ! Point vector is the node coordinate itself
  subroutine get_vector_point_data(x, pdata)
    real(r8), intent(in) :: x(:,:)
    real(r8), allocatable, intent(out) :: pdata(:,:)
    pdata = x
  end subroutine

  ! Cell scalar is the magnitude of the cell centroid
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

  ! Cell vector is the cell centroid
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
