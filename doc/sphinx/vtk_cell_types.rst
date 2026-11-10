VTK Cell Types
==============

The ``vtkhdf_vtk_cell_types`` module defines named ``int8`` parameters
corresponding to VTK cell type identifiers.

These parameters are intended for use in the ``types`` argument of
the UnstructuredGrid mesh-writing methods.

The set of named parameters provided here is not intended to be
exhaustive. In general, any VTK cell type code that is valid for
``UnstructuredGrid`` may be used in the ``types`` array; fVTKHDF
does not attempt to restrict or reinterpret these codes.

The integer values are defined by VTK and are reproduced here for
convenience. The authoritative definitions may be found in:

* The `VTKHDF file format specification
  <https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html>`_.
* ``Common/DataModel/vtkCellType.h`` in the VTK source tree.

Standard Cell Types
-------------------

.. list-table::
   :header-rows: 1
   :widths: 40 20

   * - Name
     - Value
   * - ``VTK_EMPTY_CELL``
     - 0
   * - ``VTK_VERTEX``
     - 1
   * - ``VTK_POLY_VERTEX``
     - 2
   * - ``VTK_LINE``
     - 3
   * - ``VTK_POLY_LINE``
     - 4
   * - ``VTK_TRIANGLE``
     - 5
   * - ``VTK_TRIANGLE_STRIP``
     - 6
   * - ``VTK_POLYGON``
     - 7
   * - ``VTK_PIXEL``
     - 8
   * - ``VTK_QUAD``
     - 9
   * - ``VTK_TETRA``
     - 10
   * - ``VTK_VOXEL``
     - 11
   * - ``VTK_HEXAHEDRON``
     - 12
   * - ``VTK_WEDGE``
     - 13
   * - ``VTK_PYRAMID``
     - 14
   * - ``VTK_PENTAGONAL_PRISM``
     - 15
   * - ``VTK_HEXAGONAL_PRISM``
     - 16

Quadratic and Higher-Order Cell Types
-------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 20

   * - Name
     - Value
   * - ``VTK_QUADRATIC_EDGE``
     - 21
   * - ``VTK_QUADRATIC_TRIANGLE``
     - 22
   * - ``VTK_QUADRATIC_QUAD``
     - 23
   * - ``VTK_QUADRATIC_TETRA``
     - 24
   * - ``VTK_QUADRATIC_HEXAHEDRON``
     - 25
   * - ``VTK_QUADRATIC_WEDGE``
     - 26
   * - ``VTK_QUADRATIC_PYRAMID``
     - 27
   * - ``VTK_BIQUADRATIC_QUAD``
     - 28
   * - ``VTK_TRIQUADRATIC_HEXAHEDRON``
     - 29
   * - ``VTK_QUADRATIC_LINEAR_QUAD``
     - 30
   * - ``VTK_QUADRATIC_LINEAR_WEDGE``
     - 31
   * - ``VTK_BIQUADRATIC_QUADRATIC_WEDGE``
     - 32
   * - ``VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON``
     - 33
   * - ``VTK_BIQUADRATIC_TRI_QUADRATIC_HEXAHEDRON``
     - 34
   * - ``VTK_QUADRATIC_POLYGON``
     - 36
   * - ``VTK_TRIQUADRATIC_PYRAMID``
     - 37
