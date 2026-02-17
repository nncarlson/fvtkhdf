.. toctree::
   :hidden:

.. role:: f(code)
   :language: fortran

===================
The fVTKHDF Library
===================

Parallel Concerns
=================


The vtkhdf_ug_file_type module
==============================

This module defines the `vtkhdf_ug_file` derived type for exporting mesh-based
solution data to a VTKHDF file readable by ParaView. It produces an
UnstructuredGrid (UG) dataset and supports both static and time-dependent cell
and point data on a single static mesh.

.. code-block:: fortran

   use vtkhdf_ug_file_type
   type(vtkhdf_ug_file) :: file

File procedures
---------------

:f:`call file%create(filename, comm, stat, errmsg [,is_temporal])`
   Initializes a new VTKHDF unstructured grid file.

   - `filename` (character(*), intent(in)): The name of the file to create.
   - `comm` (integer, intent(in)): MPI communicator.
   - `stat` (integer, intent(out)): Status code (0 for success).
   - `errmsg` (character(:), allocatable, intent(out)): Error message if `stat` is non-zero.
   - `is_temporal` (logical, intent(in), optional): Set to `.true.` if the file will contain time-dependent data.

``call file%close()``
    Close the file.

``call file%flush()``
    Flush the internal HDF5 buffers associated with the file and request
    the operating system to flush the system buffers for the open file.

Mesh data
-----------
This subroutine writes the unstructured mesh geometry and topology to the
file. It **must** be called before calling any of the procedures that
follow.

``call file%write_mesh(x, cnode, xcnode, types)``
  * ``x`` is a rank-2 `real32` or `real64` array storing the mesh node
    coordinates. [which dimension is the cell/point index]
  * ``cnode`` and ``xcnode`` stores the mesh topology in the common CSR
    format: ``cnode`` is a flat array storing the concatenated cell node
    connectivity data and the array ``xcnode`` stores the starting index
    in ``cnode`` of the connectivity list for each cell. The size of
    ``xcnode`` is 1 more than the number of cells [and the value of the
    last element is ...]
    ``int32`` and  ``int64`` arrays are supported.
  * ``types`` is a ``int8`` array storing the VTK types of the cells.
    The module ``vtkhdf_vtk_cell_types`` provides parameters that may
    be used to initialize this array; e.g.
    ``VTK_TETRA`` for tetrahedral cells.

Static mesh data
----------------
These subroutines write mesh-based data to the file. These can be called
at any time after writing the mesh. For temporal files, this data is
static and not associated with any time step.

.. glossary::

   ``call file%write_cell_data(name, array)``
   ``call file%write_point_data(name, array)``
      Write the cell/point-based data ``array`` to a new cell/point dataset
      ``name``. Scalar-valued data (rank-1 ``array``) and vector-valued data
      (rank-2 ``array``) are supported. [which dimension is the mesh index]

Time-dependent mesh data
------------------------
Temporal files support time-dependent mesh data. Such data must be registered
beforehand. Once the first time step is started, no further registrations are
allowed.

.. glossary::

   ``call file%register_temporal_cell_data(name, mold)``
   ``call file%register_temporal_point_data(name, mold)``
      Register ``name`` as a time-dependent cell/point dataset.
      The type and kind of `mold` determines the dataset type.
      For scalar-valued cell/point data, `mold` shall be a scalar,
      and for vector-valued data, `mold` shall be a rank-1 array
      whose size is the number of vector components. The value of
      ``mold`` is never referenced.


   ``call file%write_time_step(time)``
      Start a new time step with time value ``time``.


   ``call file%write_temporal_cell_data(name, array)``
   ``call file%write_temporal_point_data(name, array)``
      Write the cell/point-based data ``array`` to the time-dependent
      dataset ``name``. The data is associated with the current time
      step. [characteristics of array must agree with the mold specified
      when name was registered.  size of data must equal number of cells
      or points].

The vtkhdf_mb_file_type module
==============================

The vtkhdf_vtk_cell_types module
================================
This module defines a collection of parameter ...

[TODO: (1) put the values into a compact table; (2) keep a reference to the source of the data]
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

