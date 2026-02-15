module vtkhdf_vtk_cell_types
  use, intrinsic :: iso_fortran_env, only: int8
  implicit none
  public

  !! Standard VTK Cell Types
  !! See https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html#
  !! and Common/DataModel/vtkCellType.h in the VTK source code.

  integer(int8), parameter :: VTK_EMPTY_CELL = 0
  integer(int8), parameter :: VTK_VERTEX = 1
  integer(int8), parameter :: VTK_POLY_VERTEX = 2
  integer(int8), parameter :: VTK_LINE = 3
  integer(int8), parameter :: VTK_POLY_LINE = 4
  integer(int8), parameter :: VTK_TRIANGLE = 5
  integer(int8), parameter :: VTK_TRIANGLE_STRIP = 6
  integer(int8), parameter :: VTK_POLYGON = 7
  integer(int8), parameter :: VTK_PIXEL = 8
  integer(int8), parameter :: VTK_QUAD = 9
  integer(int8), parameter :: VTK_TETRA = 10
  integer(int8), parameter :: VTK_VOXEL = 11
  integer(int8), parameter :: VTK_HEXAHEDRON = 12
  integer(int8), parameter :: VTK_WEDGE = 13
  integer(int8), parameter :: VTK_PYRAMID = 14
  integer(int8), parameter :: VTK_PENTAGONAL_PRISM = 15
  integer(int8), parameter :: VTK_HEXAGONAL_PRISM = 16

  !! Quadratic, isoparametric cells
  integer(int8), parameter :: VTK_QUADRATIC_EDGE = 21
  integer(int8), parameter :: VTK_QUADRATIC_TRIANGLE = 22
  integer(int8), parameter :: VTK_QUADRATIC_QUAD = 23
  integer(int8), parameter :: VTK_QUADRATIC_TETRA = 24
  integer(int8), parameter :: VTK_QUADRATIC_HEXAHEDRON = 25
  integer(int8), parameter :: VTK_QUADRATIC_WEDGE = 26
  integer(int8), parameter :: VTK_QUADRATIC_PYRAMID = 27
  integer(int8), parameter :: VTK_BIQUADRATIC_QUAD = 28
  integer(int8), parameter :: VTK_TRIQUADRATIC_HEXAHEDRON = 29
  integer(int8), parameter :: VTK_QUADRATIC_LINEAR_QUAD = 30
  integer(int8), parameter :: VTK_QUADRATIC_LINEAR_WEDGE = 31
  integer(int8), parameter :: VTK_BIQUADRATIC_QUADRATIC_WEDGE = 32
  integer(int8), parameter :: VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON = 33
  integer(int8), parameter :: VTK_BIQUADRATIC_TRI_QUADRATIC_HEXAHEDRON = 34
  integer(int8), parameter :: VTK_QUADRATIC_POLYGON = 36
  integer(int8), parameter :: VTK_TRIQUADRATIC_PYRAMID = 37

end module vtkhdf_vtk_cell_types
