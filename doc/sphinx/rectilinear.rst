RectilinearGrid
===============

The ``vtkhdf_rectilinear_file_type`` module defines the ``vtkhdf_rectilinear_file``
derived type for writing VTKHDF RectilinearGrid files. The geometry is the
Cartesian product of three coordinate vectors, allowing nonuniform spacing
along each axis. One-, two-, and three-dimensional grids are supported.
Reading these files requires ParaView 6.2 or later (VTK 9.7 or later).

File management, point and cell data, and time stepping use the
:doc:`shared interface for structured VTK datasets <structured>`.
Field data uses the :doc:`common field data interface <field_data>`.

The RectilinearGrid writer supports ``VTKHDF_STATIC``
(default), ``VTKHDF_FIXED_MESH``, and ``VTKHDF_MOVING_MESH``. The logical grid
is fixed for the lifetime of the file. Moving-mesh mode permits changing
coordinate vectors and static or temporal point, cell, and field data.

.. code-block:: fortran

   use vtkhdf_rectilinear_file_type
   type(vtkhdf_rectilinear_file) :: file

Mesh Data
---------

.. method-doc::

   call file%write_mesh(point_bounds)

   Define the global logical grid. Call exactly once, before writing geometry
   or writing/registering point/cell data.

   * ``point_bounds(3,2)``: inclusive global logical point bounds of kind
     ``int32``. Column 1 contains lower bounds and column 2 upper bounds,
     in x, y, z order. All ranks supply identical bounds. Upper bounds must
     be at least their lower bounds, with at least one nonsingleton direction.
     Each point dimension must fit in a signed 32-bit integer, and the total
     point count must fit in a signed 64-bit integer.

.. method-doc::

   call file%write_mesh_geometry(x, y, z)

   Write the complete geometry as three coordinate vectors. In static and
   fixed-mesh modes, call once. In moving-mesh mode, call for the initial
   geometry and whenever it changes, at most once per time step. Supply all
   three vectors together, before writing point/cell data for that step.

   * ``x(:)``, ``y(:)``, ``z(:)``: full physical coordinate vectors in logical
     index order, supplied by rank zero. Their lengths equal the corresponding
     global point dimensions, including one coordinate for a singleton axis.
     Values and lengths on other ranks are ignored; zero-sized arrays suffice.
     All three vectors, on all ranks, must have the same kind, either
     ``real32`` or ``real64``. Storage preserves that kind, which must remain
     the same for all geometry writes.

   This is a collective call, including on ranks supplying no coordinates.
   The first coordinate in each vector corresponds to the lower global
   logical point bound on that axis.

.. method-doc::

   call file%write_mesh(point_bounds, x, y, z)

   Define the logical grid and write its initial coordinate vectors. This is
   equivalent to calling the two methods above in sequence. Call exactly once;
   subsequent moving-mesh updates use ``write_mesh_geometry``.

Omitting ``write_mesh_geometry`` during a step reuses the most recently
written coordinates. Geometry and point/cell data may be omitted or replaced
independently. The global point dimensions and logical bounds remain fixed.
Initial geometry may be supplied before or during the first time step.

Logical indices used by fVTKHDF may have arbitrary lower bounds. VTKHDF stores
only the dimensions for this grid type, so the reader exposes zero-based
indices. Physical coordinates and data placement are preserved.

.. warning::

   Readers earlier than VTK 9.7.2 (ParaView 6.2.1) require singleton axes to be
   trailing: use ``(nx,ny,1)`` for 2D grids and ``(nx,1,1)`` for 1D grids. Other
   arrangements may produce incorrect results.

Examples
--------

:repo:`example/rectilinear_fixed_mesh.F90` follows the ImageData example, replacing
its geometry definition with nonuniform coordinate vectors supplied by rank
zero. It demonstrates static cell data, temporal point and field data, and
shared-point ownership using ``write_bounds`` in serial and MPI builds.

:repo:`example/rectilinear_moving_mesh.F90` follows the fixed-mesh example,
shrinking the x coordinates and expanding the y coordinates at each successive
step. Static diffusivity and temporal point and field data are written as in
the fixed-mesh example. It supports serial and MPI builds.
