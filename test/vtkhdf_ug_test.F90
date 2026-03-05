program vtkhdf_ug_test

  use,intrinsic :: iso_fortran_env, only: r8 => real64, r4 => real32, int8
  use vtkhdf_ug_file_type
  use vtkhdf_vtk_cell_types
  use mpi_f08
  implicit none

  type(vtkhdf_ug_file) :: vizfile
  integer :: istat, nproc, rank, stat
  character(:), allocatable :: errmsg
  real(r8), allocatable :: points(:,:)
  integer,  allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  real(r8), allocatable :: s(:), v(:,:) ! scalar and vector data arrays
  real(r8), allocatable :: fs(:), fv(:,:) ! scalar and vector field arrays
  type(vtkhdf_cell_data_handle) :: hcell_scalar, hcell_vector
  type(vtkhdf_point_data_handle) :: hpoint_scalar, hpoint_vector
  type(vtkhdf_field_data_handle) :: hfield_value, hfield_scalar, hfield_vector

  call MPI_Init(istat)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, istat)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, istat)

  call vizfile%create('ug_test.vtkhdf', MPI_COMM_WORLD, stat, errmsg, is_temporal=.true.)
  if (stat /= 0) error stop errmsg

  !! The unstructured mesh data for a basic mesh unit. The full mesh will
  !! be a collection of non-overlapping shifts of this basic unit. Each
  !! rank has one of these, which is right-shifted proportional to the rank.
  call get_mesh_data(points, cnode, xcnode, types)
  points(1,:) = points(1,:) + rank ! shift right

  call vizfile%write_mesh(points, cnode, xcnode, types)

  !!!! Register the data arrays that evolve with time.

  associate (scalar_mold => 0.0_r8, vector_mold => [real(r8) :: 0, 0, 0])
    hcell_scalar = vizfile%register_temporal_cell_data('cell-scalar', scalar_mold)
    hcell_vector = vizfile%register_temporal_cell_data('cell-vector', vector_mold)
    hpoint_scalar = vizfile%register_temporal_point_data('point-scalar', scalar_mold)
    hpoint_vector = vizfile%register_temporal_point_data('point-vector', vector_mold)
    hfield_value = vizfile%register_temporal_field_data('field-value', scalar_mold)
    hfield_scalar = vizfile%register_temporal_field_data('field-scalar', scalar_mold)
    hfield_vector = vizfile%register_temporal_field_data('field-vector', scalar_mold)
  end associate

  !!!! Write the datasets for the first time step !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call vizfile%start_time_step(0.0_r4)

  call get_scalar_cell_data(points, cnode, xcnode, s)
  call vizfile%write_temporal_cell_data(hcell_scalar, s)

  call get_vector_cell_data(points, cnode, xcnode, v)
  call vizfile%write_temporal_cell_data(hcell_vector, v)

  call get_scalar_point_data(points, s)
  call vizfile%write_temporal_point_data(hpoint_scalar, s)

  call get_vector_point_data(points, v)
  call vizfile%write_temporal_point_data(hpoint_vector, v)

  fs = [1.0_r8, 2.0_r8]
  fv = reshape([1.0_r8, 2.0_r8, 3.0_r8, 11.0_r8, 12.0_r8, 13.0_r8], [3,2])
  call vizfile%write_temporal_field_data(hfield_value, 42.0_r8)
  call vizfile%write_temporal_field_data(hfield_scalar, fs, as_vector=.true.)
  call vizfile%write_temporal_field_data(hfield_vector, fv)
  call vizfile%finalize_time_step()

  call vizfile%flush()

  !!!! Write the datasets for the second time step !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call vizfile%start_time_step(1.0_r4)

  call get_scalar_cell_data(points, cnode, xcnode, s)
  call vizfile%write_temporal_cell_data(hcell_scalar, s+1)

  call get_vector_cell_data(points, cnode, xcnode, v)
  call vizfile%write_temporal_cell_data(hcell_vector, v+1)

  call get_scalar_point_data(points, s)
  call vizfile%write_temporal_point_data(hpoint_scalar, s+1)

  fs = [10.0_r8, 20.0_r8, 30.0_r8]
  call vizfile%write_temporal_field_data(hfield_scalar, fs)

  !! Skip one temporal dataset write; offset should repeat the last written value.
  call vizfile%finalize_time_step()

  !! Some time-independent cell and point data

  call get_scalar_cell_data(points, cnode, xcnode, s)
  call vizfile%write_cell_data('static-cell-scalar', -s)

  call get_vector_cell_data(points, cnode, xcnode, v)
  call vizfile%write_cell_data('static-cell-vector', -v)

  call get_scalar_point_data(points, s)
  call vizfile%write_point_data('static-point-scalar', -s)

  call get_vector_point_data(points, v)
  call vizfile%write_point_data('static-point-vector', -v)

  fs = [-1.0_r8, -2.0_r8, -3.0_r8, -4.0_r8]
  fv = reshape([-1.0_r8, -2.0_r8, -3.0_r8, -11.0_r8, -12.0_r8, -13.0_r8], [3,2])
  call vizfile%write_field_data('static-field-scalar', fs, as_vector=.true.)
  call vizfile%write_field_data('static-field-vector', fv)
  call vizfile%write_field_data('static-field-value', -9.0_r8)

  call vizfile%close
  call MPI_Finalize(istat)

contains

  ! A 5-tet subdivision of a squished unit cube.
  subroutine get_mesh_data(points, cnode, xcnode, types)
    real(r8), allocatable, intent(out) :: points(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable :: types(:)
    points = reshape([0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1], shape=[3,8])
    ! distort to catch C/Fortran index ordering errors
    points(1,:) = 0.9_r8*points(1,:)
    points(2,:) = 0.7_r8*points(2,:)
    points(3,:) = 0.5_r8*points(3,:)
    cnode = [1,2,4,5, 2,3,4,7, 2,5,6,7, 4,5,7,8, 2,4,5,7]
    xcnode = [1,5,9,13,17,21]
    types = spread(VTK_TETRA, dim=1, ncopies=5)
  end subroutine

  ! Point scalar is the magnitude of the node coordinate
  subroutine get_scalar_point_data(points, pdata)
    real(r8), intent(in) :: points(:,:)
    real(r8), allocatable, intent(out) :: pdata(:)
    integer :: j
    allocate(pdata(size(points,dim=2)))
    do j = 1, size(points,dim=2)
      pdata(j) = norm2(points(:,j))
    end do
  end subroutine

  ! Point vector is the node coordinate itself
  subroutine get_vector_point_data(points, pdata)
    real(r8), intent(in) :: points(:,:)
    real(r8), allocatable, intent(out) :: pdata(:,:)
    pdata = points
  end subroutine

  ! Cell scalar is the magnitude of the cell centroid
  subroutine get_scalar_cell_data(points, cnode, xcnode, cdata)
    real(r8), intent(in) :: points(:,:)
    integer, intent(in) :: cnode(:), xcnode(:)
    real(r8), allocatable, intent(out) :: cdata(:)
    integer :: j
    allocate(cdata(size(xcnode)-1))
    do j = 1, size(cdata)
      associate(pid => cnode(xcnode(j):xcnode(j+1)-1))
        cdata(j) = norm2(sum(points(:,pid),dim=2)/size(pid))
      end associate
    end do
  end subroutine

  ! Cell vector is the cell centroid
  subroutine get_vector_cell_data(points, cnode, xcnode, cdata)
    real(r8), intent(in) :: points(:,:)
    integer, intent(in) :: cnode(:), xcnode(:)
    real(r8), allocatable, intent(out) :: cdata(:,:)
    integer :: j
    allocate(cdata(size(points,dim=1),size(xcnode)-1))
    do j = 1, size(cdata,dim=2)
      associate(pid => cnode(xcnode(j):xcnode(j+1)-1))
        cdata(:,j) = sum(points(:,pid),dim=2)/size(pid)
      end associate
    end do
  end subroutine

end program
