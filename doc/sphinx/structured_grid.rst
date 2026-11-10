StructuredGrid
==============

The ``vtkhdf_structured_file_type`` module defines the ``vtkhdf_structured_file``
derived type for writing VTKHDF StructuredGrid files. Every grid point has
explicit physical coordinates, allowing warped or curved grids with implicit
structured connectivity. One-, two-, and three-dimensional grids are supported.
Reading these files requires ParaView 6.2 or later (VTK 9.7 or later).

File management, point and cell data, and time stepping use the
:doc:`shared interface for structured VTK datasets <structured>`.
Field data uses the :doc:`common field data interface <field_data>`.

The StructuredGrid writer supports ``VTKHDF_STATIC``
(default), ``VTKHDF_FIXED_MESH``, and ``VTKHDF_MOVING_MESH``. The logical grid
is fixed for the lifetime of the file. Moving-mesh mode permits changing
point coordinates and static or temporal point, cell, and field data.

.. code-block:: fortran

   use vtkhdf_structured_file_type
   type(vtkhdf_structured_file) :: file

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

   call file%write_mesh_geometry(points, array_origin [, write_bounds])

   Write a rectangular region of physical point coordinates. Repeated calls
   assemble the global geometry. Complete these writes before writing point or
   cell data and before finalizing the corresponding time step. Initial
   geometry may be supplied before or during the first step. Static and
   fixed-mesh modes permit only this initial geometry; moving-mesh mode permits
   a replacement during each subsequent step.

   * ``points(3,nx,ny,nz)``: physical x, y, z coordinates, with all three
     spatial dimensions present even for 1D or 2D. Accepted kinds are
     ``real32`` and ``real64``. The first call creates the coordinate dataset
     with the supplied kind, including when some contributions are empty.
     All ranks and all subsequent calls must use that same kind.
   * ``array_origin(3)``: global logical point indices corresponding to the
     first spatial element of ``points``; kind ``int32``.
   * ``write_bounds(3,2)``: optional inclusive logical point bounds selecting
     the region to write; kind ``int32``. The selection must lie within both
     the global point bounds and the mapped footprint of ``points``.

   Geometry uses the same region selection, halo exclusion, and empty
   contribution conventions as :doc:`structured`. Visible points must have
   coordinates supplied exactly once. Unwritten coordinates default to quiet
   NaN. A zero spatial extent with ``write_bounds`` omitted supplies an empty
   contribution; the component extent must still be three. ``array_origin``
   remains required, and ``[0,0,0]`` suffices for an empty contribution.

In moving-mesh mode, the first geometry call in a step begins a replacement;
subsequent calls in that step supply further regions of the same replacement.
An all-empty collective call also begins a replacement. Unwritten regions use
fill values, rather than inheriting preceding coordinates. Complete all regions
before writing point/cell data or finalizing the step.

Omitting all geometry calls during a step reuses the most recently written
geometry. A stored geometry snapshot is therefore distinct from a time step:
seven steps with geometry written at steps 0, 2, and 4 require three snapshots.
Geometry and point/cell data may be reused independently because the logical
grid dimensions remain fixed.

Logical indices used by fVTKHDF may have arbitrary lower bounds. VTKHDF stores
only the dimensions for this grid type, so the reader exposes zero-based
indices. Physical coordinates and data placement are preserved.

.. warning::

   Readers earlier than VTK 9.7.2 (ParaView 6.2.1) require singleton axes to be
   trailing: use ``(nx,ny,1)`` for 2D grids and ``(nx,1,1)`` for 1D grids. Other
   arrangements may produce incorrect results.

Examples
--------

:repo:`example/structured_fixed_mesh.F90` follows the fixed-mesh examples, mapping
the logical grid to an annular strip whose height increases with azimuth.
It demonstrates distributed point coordinates, static cell data, temporal point
and field data, and shared-point ownership in serial and MPI builds.

:repo:`example/structured_moving_mesh.F90` follows the fixed-mesh example, rotating
the strip about z by one angular increment at each successive step. Static cell
data and temporal point and field data are written as in the fixed-mesh example.
It supports serial and MPI builds.
