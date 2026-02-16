program vtkhdf_ug_test

  use,intrinsic :: iso_fortran_env, only: r8 => real64, r4 => real32, int8
  use vtkhdf_ug_file_type
  use vtkhdf_vtk_cell_types
  use mpi
  implicit none

  type(vtkhdf_ug_file) :: vizfile
  integer :: istat, nproc, rank, stat
  character(:), allocatable :: errmsg
  real(r8), allocatable :: x(:,:)
  integer,  allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  real(r8), allocatable :: s(:), v(:,:) ! scalar and vector data arrays

  call MPI_Init(istat)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, istat)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, istat)

  call vizfile%create('ug_test.vtkhdf', MPI_COMM_WORLD, stat, errmsg, is_temporal=.true.)
  if (stat /= 0) error stop errmsg

  !! The unstructured mesh data for a basic mesh unit. The full mesh will
  !! be a collection of non-overlapping shifts of this basic unit. Each
  !! rank has one of these, which is right-shifted proportional to the rank.
  call get_mesh_data(x, cnode, xcnode, types)
  x(1,:) = x(1,:) + rank ! shift right

  call vizfile%write_mesh(x, cnode, xcnode, types)

  !!!! Register the data arrays that evolve with time.

  associate (scalar_mold => 0.0_r8, vector_mold => [real(r8) :: 0, 0, 0])
    call vizfile%register_temporal_cell_data('cell-scalar', scalar_mold)
    call vizfile%register_temporal_cell_data('cell-vector', vector_mold)
    call vizfile%register_temporal_point_data('point-scalar', scalar_mold)
    call vizfile%register_temporal_point_data('point-vector', vector_mold)
  end associate

  !!!! Write the datasets for the first time step !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call vizfile%write_time_step(0.0_r4)

  call get_scalar_cell_data(x, cnode, xcnode, s)
  call vizfile%write_temporal_cell_data('cell-scalar', s)

  call get_vector_cell_data(x, cnode, xcnode, v)
  call vizfile%write_temporal_cell_data('cell-vector', v)

  call get_scalar_point_data(x, s)
  call vizfile%write_temporal_point_data('point-scalar', s)

  call get_vector_point_data(x, v)
  call vizfile%write_temporal_point_data('point-vector', v)

  call vizfile%flush()

  !!!! Write the datasets for the second time step !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call vizfile%write_time_step(1.0_r4)

  call get_scalar_cell_data(x, cnode, xcnode, s)
  call vizfile%write_temporal_cell_data('cell-scalar', s+1)

  call get_vector_cell_data(x, cnode, xcnode, v)
  call vizfile%write_temporal_cell_data('cell-vector', v+1)

  call get_scalar_point_data(x, s)
  call vizfile%write_temporal_point_data('point-scalar', s+1)

  call get_vector_point_data(x, v)
  call vizfile%write_temporal_point_data('point-vector', v+1)

  !! Some time-independent cell and point data

  call get_scalar_cell_data(x, cnode, xcnode, s)
  call vizfile%write_cell_data('static-cell-scalar', -s)

  call get_vector_cell_data(x, cnode, xcnode, v)
  call vizfile%write_cell_data('static-cell-vector', -v)

  call get_scalar_point_data(x, s)
  call vizfile%write_point_data('static-point-scalar', -s)

  call get_vector_point_data(x, v)
  call vizfile%write_point_data('static-point-vector', -v)

  call vizfile%close
  call MPI_Finalize(istat)

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
