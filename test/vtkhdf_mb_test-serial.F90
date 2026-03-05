program vtkhdf_mb_test

  use,intrinsic :: iso_fortran_env, only: r8 => real64, int8
  use vtkhdf_mb_file_type
  use vtkhdf_vtk_cell_types
  implicit none

  real(r8), allocatable :: points(:,:), scalar_cell_data(:), vector_cell_data(:,:)
  real(r8), allocatable :: scalar_point_data(:), vector_point_data(:,:)
  real(r8), allocatable :: scalar_field_data(:), vector_field_data(:,:)
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  character(:), allocatable :: errmsg
  integer :: stat

  type(vtkhdf_mb_file) :: vizfile
  type(vtkhdf_block_handle) :: hblk_a, hblk_b
  type(vtkhdf_cell_data_handle) :: hcell_radius, hcell_velocity
  type(vtkhdf_point_data_handle) :: hpoint_radius, hpoint_velocity
  type(vtkhdf_field_data_handle) :: hfield_value, hfield_scalar, hfield_vector

  call vizfile%create('mb_test.vtkhdf', stat, errmsg)
  if (stat /= 0) error stop errmsg

  !! The unstructured mesh data for a basic mesh unit. The full mesh will
  !! consist of two non-overlapping shifts of this basic unit.
  call get_mesh_data(points, cnode, xcnode, types)

  hblk_a = vizfile%add_block('Block-A', is_temporal=.true.)
  call vizfile%write_mesh(hblk_a, points, cnode, xcnode, types)

  hblk_b = vizfile%add_block('Block-B', is_temporal=.true.)
  call vizfile%write_mesh(hblk_b, points+1, cnode, xcnode, types)

  !! Register the datasets that evolve with time. At this stage the data arrays
  !! are only used to glean their types and shapes.

  associate (scalar_mold => 0.0_r8, vector_mold => [real(r8) :: 0, 0, 0])
    ! Block-A has time-dependent cell data
    hcell_radius = vizfile%register_temporal_cell_data(hblk_a, 'cell-radius', scalar_mold)
    hcell_velocity = vizfile%register_temporal_cell_data(hblk_a, 'cell-velocity', vector_mold)
    ! Block-B has time-dependent point data
    hpoint_radius = vizfile%register_temporal_point_data(hblk_b, 'point-radius', scalar_mold)
    hpoint_velocity = vizfile%register_temporal_point_data(hblk_b, 'point-velocity', vector_mold)
    hfield_value = vizfile%register_temporal_field_data(hblk_a, 'field-value', scalar_mold)
    hfield_scalar = vizfile%register_temporal_field_data(hblk_a, 'field-scalar', scalar_mold)
    hfield_vector = vizfile%register_temporal_field_data(hblk_a, 'field-vector', scalar_mold)
  end associate

  !! Generate some cell and point data to use for output
  call get_scalar_cell_data(points, cnode, xcnode, scalar_cell_data)
  call get_vector_cell_data(points, cnode, xcnode, vector_cell_data)

  call get_scalar_point_data(points+1, scalar_point_data)
  call get_vector_point_data(points+1, vector_point_data)
  scalar_field_data = [1.0_r8, 2.0_r8]
  vector_field_data = reshape([1.0_r8, 2.0_r8, 3.0_r8, 11.0_r8, 12.0_r8, 13.0_r8], [3,2])

  !!!! Write the data for the first time step !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call vizfile%start_time_step(0.0_r8)

  call vizfile%write_temporal_cell_data(hblk_a, hcell_radius, scalar_cell_data)
  call vizfile%write_temporal_cell_data(hblk_a, hcell_velocity, vector_cell_data)
  call vizfile%write_temporal_point_data(hblk_b, hpoint_radius, scalar_point_data)
  call vizfile%write_temporal_point_data(hblk_b, hpoint_velocity, vector_point_data)
  call vizfile%write_temporal_field_data(hblk_a, hfield_value, 42.0_r8)
  call vizfile%write_temporal_field_data(hblk_a, hfield_scalar, scalar_field_data, as_vector=.true.)
  call vizfile%write_temporal_field_data(hblk_a, hfield_vector, vector_field_data)
  call vizfile%finalize_time_step()

  !!!! Write the data for the second time step !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call vizfile%start_time_step(1.0_r8)

  call vizfile%write_temporal_cell_data(hblk_a, hcell_radius, scalar_cell_data+1)
  call vizfile%write_temporal_cell_data(hblk_a, hcell_velocity, vector_cell_data+1)
  call vizfile%write_temporal_point_data(hblk_b, hpoint_radius, scalar_point_data+1)
  call vizfile%write_temporal_point_data(hblk_b, hpoint_velocity, vector_point_data+1)
  call vizfile%write_temporal_field_data(hblk_a, hfield_scalar, [10.0_r8, 20.0_r8, 30.0_r8])
  call vizfile%finalize_time_step()

  !! At any point you can write a data that isn't time dependent, but its name must
  !! be unique from any other data temporal or not of the same type (cell or point).

  call vizfile%write_cell_data(hblk_a, 'static-cell-scalar', -scalar_cell_data)
  call vizfile%write_cell_data(hblk_a, 'static-cell-vector', -vector_cell_data)
  call vizfile%write_point_data(hblk_b, 'static-point-scalar', -scalar_point_data)
  call vizfile%write_point_data(hblk_b, 'static-point-vector', -vector_point_data)
  call vizfile%write_field_data(hblk_a, 'static-field-scalar', [-1.0_r8, -2.0_r8, -3.0_r8, -4.0_r8], as_vector=.true.)
  call vizfile%write_field_data(hblk_a, 'static-field-vector', -vector_field_data)
  call vizfile%write_field_data(hblk_a, 'static-field-value', -9.0_r8)

  call vizfile%close

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
