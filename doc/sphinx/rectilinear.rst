The vtkhdf_rectilinear_file_type module
=======================================

This module defines ``vtkhdf_rectilinear_file`` for writing VTKHDF
RectilinearGrid files. The geometry is the Cartesian product of three
coordinate vectors, allowing nonuniform spacing along each axis. One-, two-,
and three-dimensional grids are supported.

File management, point and cell data, and time stepping use :doc:`structured`.
Field data uses :doc:`field_data`. The module exports ``VTKHDF_STATIC``
(default) and ``VTKHDF_FIXED_MESH``. Geometry is fixed for the lifetime of
the file; the latter mode allows temporal point, cell, and field data.

.. code-block:: fortran

   use vtkhdf_rectilinear_file_type
   type(vtkhdf_rectilinear_file) :: file

Mesh Data
---------

.. method-doc::

   call file%write_mesh(point_bounds, x, y, z)

   Define the global grid and write its coordinate vectors. Call exactly
   once, before writing or registering point/cell data or starting time stepping.

   * ``point_bounds(3,2)``: inclusive global logical point bounds of kind
     ``int32``. Column 1 contains lower bounds and column 2 upper bounds,
     in x, y, z order. All ranks supply identical bounds. Upper bounds must
     be at least their lower bounds, with at least one nonsingleton direction.
     Each point dimension must fit in a signed 32-bit integer, and the total
     point count must fit in a signed 64-bit integer.
   * ``x(:)``, ``y(:)``, ``z(:)``: full physical coordinate vectors in logical
     index order, supplied by rank zero. Their lengths equal the corresponding
     global point dimensions, including one coordinate for a singleton axis.
     Values and lengths on other ranks are ignored; zero-sized arrays suffice.
     All three vectors, on all ranks, must have the same kind, either
     ``real32`` or ``real64``. Storage preserves that kind.

   This is a collective call, including on ranks supplying no coordinates.
   The first coordinate in each vector corresponds to the lower global
   logical point bound on that axis.

Logical indices used by fVTKHDF may have arbitrary lower bounds. VTKHDF stores
only the dimensions for this grid type, so the reader exposes zero-based
indices. Physical coordinates and data placement are preserved.

Reader Compatibility
--------------------

RectilinearGrid requires ParaView 6.2+ (VTK 9.7+). Readers before VTK 9.7.2
have a shared structured-data defect for nontrailing singleton axes. With an
affected reader, use ``(nx,ny,1)`` for 2D and ``(nx,1,1)`` for 1D. Arbitrary
logical lower bounds remain supported because the file uses zero-based storage.

Example
-------

``example/rectilinear_fixed_mesh.F90`` demonstrates rank-zero coordinate
vectors, distributed cell data, and temporal reuse in serial and MPI builds.

.. literalinclude:: ../../example/rectilinear_fixed_mesh.F90
   :language: fortran
