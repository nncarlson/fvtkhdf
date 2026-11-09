The vtkhdf_structured_file_type module
======================================

This module defines ``vtkhdf_structured_file`` for writing VTKHDF StructuredGrid
files. Every grid point has explicit physical coordinates, allowing warped or
curved grids with implicit structured connectivity. One-, two-, and
three-dimensional grids are supported.

File management, point and cell data, and time stepping use :doc:`structured`.
Field data uses :doc:`field_data`. The module exports ``VTKHDF_STATIC``
(default), ``VTKHDF_FIXED_MESH``, and ``VTKHDF_MOVING_MESH``. The logical grid
is fixed for the lifetime of the file. Moving-mesh mode permits changing
point coordinates and requires temporal point and cell data. Field data may
be static or temporal.

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

Reader Compatibility
--------------------

StructuredGrid requires ParaView 6.2+ (VTK 9.7+). Readers before VTK 9.7.2
have a shared structured-data defect for nontrailing singleton axes. With an
affected reader, use ``(nx,ny,1)`` for 2D and ``(nx,1,1)`` for 1D. Arbitrary
logical lower bounds remain supported because the file uses zero-based storage.

Example
-------

``example/structured_fixed_mesh.F90`` demonstrates distributed point
coordinates and temporal cell data in serial and MPI builds.

.. literalinclude:: ../../example/structured_fixed_mesh.F90
   :language: fortran

``example/structured_moving_mesh.F90`` demonstrates distributed moving
coordinates with independent geometry and data reuse.

.. literalinclude:: ../../example/structured_moving_mesh.F90
   :language: fortran
