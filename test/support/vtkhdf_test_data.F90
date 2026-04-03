module vtkhdf_test_data

  use, intrinsic :: iso_fortran_env, only: r8 => real64, int8
  use vtkhdf_vtk_cell_types, only: VTK_TETRA, VTK_TRIANGLE, VTK_QUAD
  implicit none

  private
  public :: get_zero_mesh_data
  public :: get_tet_cube_mesh_data
  public :: get_quad_mesh_data
  public :: get_tri_remesh_data
  public :: shift_points
  public :: get_scalar_point_data
  public :: get_vector_point_data
  public :: get_scalar_cell_data
  public :: get_vector_cell_data
  public :: get_scalar_field_data
  public :: get_vector_field_data

contains

  subroutine get_zero_mesh_data(points, cnode, xcnode, types)
    real(r8), allocatable, intent(out) :: points(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable, intent(out) :: types(:)

    allocate(points(3,0), cnode(0), xcnode(1), types(0))
    xcnode = [1]
  end subroutine get_zero_mesh_data

  ! A 5-tet subdivision of a squished unit cube.
  subroutine get_tet_cube_mesh_data(points, cnode, xcnode, types)
    real(r8), allocatable, intent(out) :: points(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable, intent(out) :: types(:)

    points = reshape([0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1], shape=[3,8])
    ! Distort to catch C/Fortran index ordering errors.
    points(1,:) = 0.9_r8*points(1,:)
    points(2,:) = 0.7_r8*points(2,:)
    points(3,:) = 0.5_r8*points(3,:)
    cnode = [1,2,4,5, 2,3,4,7, 2,5,6,7, 4,5,7,8, 2,4,5,7]
    xcnode = [1,5,9,13,17,21]
    types = spread(VTK_TETRA, dim=1, ncopies=5)
  end subroutine get_tet_cube_mesh_data

  ! Single quad cell mesh.
  subroutine get_quad_mesh_data(points, cnode, xcnode, types)
    real(r8), allocatable, intent(out) :: points(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable, intent(out) :: types(:)

    points = 0.8_r8 * reshape([0,0,0, 1,0,0, 1,1,0, 0,1,0], shape=[3,4])
    cnode = [1,2,3,4]
    xcnode = [1,5]
    types = [VTK_QUAD]
  end subroutine get_quad_mesh_data

  ! Four-triangle remeshing of get_quad_mesh_data.
  subroutine get_tri_remesh_data(points, cnode, xcnode, types)
    real(r8), allocatable, intent(out) :: points(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable, intent(out) :: types(:)

    points = 0.4_r8 * reshape([0,0,0, 2,0,0, 2,2,0, 0,2,0, 1,1,0], shape=[3,5])
    cnode = [1,2,5, 2,3,5, 3,4,5, 4,1,5]
    xcnode = [1,4,7,10,13]
    types = spread(VTK_TRIANGLE, dim=1, ncopies=4)
  end subroutine get_tri_remesh_data

  subroutine shift_points(points, dx, dy, dz)
    real(r8), intent(inout) :: points(:,:)
    real(r8), intent(in), optional :: dx, dy, dz

    if (present(dx)) points(1,:) = points(1,:) + dx
    if (present(dy)) points(2,:) = points(2,:) + dy
    if (present(dz)) points(3,:) = points(3,:) + dz
  end subroutine shift_points

  ! Point scalar is the magnitude of the node coordinate.
  subroutine get_scalar_point_data(points, pdata)
    real(r8), intent(in) :: points(:,:)
    real(r8), allocatable, intent(out) :: pdata(:)
    integer :: j

    allocate(pdata(size(points,dim=2)))
    do j = 1, size(points,dim=2)
      pdata(j) = norm2(points(:,j))
    end do
  end subroutine get_scalar_point_data

  ! Point vector is the node coordinate itself.
  subroutine get_vector_point_data(points, pdata)
    real(r8), intent(in) :: points(:,:)
    real(r8), allocatable, intent(out) :: pdata(:,:)

    pdata = points
  end subroutine get_vector_point_data

  ! Cell scalar is the magnitude of the cell centroid.
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
  end subroutine get_scalar_cell_data

  ! Cell vector is the cell centroid.
  subroutine get_vector_cell_data(points, cnode, xcnode, cdata)
    real(r8), intent(in) :: points(:,:)
    integer, intent(in) :: cnode(:), xcnode(:)
    real(r8), allocatable, intent(out) :: cdata(:,:)
    integer :: j

    allocate(cdata(size(points,dim=1), size(xcnode)-1))
    do j = 1, size(cdata,dim=2)
      associate(pid => cnode(xcnode(j):xcnode(j+1)-1))
        cdata(:,j) = sum(points(:,pid),dim=2)/size(pid)
      end associate
    end do
  end subroutine get_vector_cell_data

  subroutine get_scalar_field_data(data)
    real(r8), allocatable, intent(out) :: data(:)

    data = [1.0_r8, 2.0_r8]
  end subroutine get_scalar_field_data

  subroutine get_vector_field_data(data)
    real(r8), allocatable, intent(out) :: data(:,:)

    data = reshape([1.0_r8, 2.0_r8, 3.0_r8, 11.0_r8, 12.0_r8, 13.0_r8], [3,2])
  end subroutine get_vector_field_data

end module vtkhdf_test_data
