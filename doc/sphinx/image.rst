ImageData
=========

The ``vtkhdf_image_file_type`` module defines the ``vtkhdf_image_file``
derived type for writing VTKHDF ImageData files. ImageData represents a regular
structured grid. One-, two-, and three-dimensional grids are supported. The
writer supports ``VTKHDF_STATIC`` (default) and ``VTKHDF_FIXED_MESH``. The mesh
geometry cannot change during the lifetime of the file.

File management, point and cell data, and time stepping use the
:doc:`shared interface for structured VTK datasets <structured>`.
Field data uses the :doc:`common field data interface <field_data>`.

.. code-block:: fortran

   use vtkhdf_image_file_type
   type(vtkhdf_image_file) :: file

Mesh Data
---------

.. method-doc::

   call file%write_mesh(point_bounds, origin, spacing [, direction])

   Define the global grid geometry. Call exactly once, before writing or
   registering point/cell data or starting time stepping. All ranks must
   supply identical arguments.

   * ``point_bounds(3,2)``: inclusive global logical point bounds of kind
     ``int32``. Column 1 contains lower bounds and column 2 contains
     upper bounds, in x, y, z order. Each upper bound must be at least its
     lower bound, and at least one direction must contain more than one point.
   * ``origin(3)``: physical position of logical index zero, including when
     zero lies outside the global bounds.
   * ``spacing(3)``: physical spacing along the logical axes.
   * ``direction(3,3)``: optional direction matrix, defaulting to identity.
     Its columns describe the logical axes. Column lengths are preserved,
     without normalization or an orthogonality requirement.

   ``origin``, ``spacing``, and ``direction`` accept ``real32`` or ``real64``;
   all supplied real arguments must have the same kind. Each point dimension
   must fit in a signed 32-bit integer, and the total point count must fit in
   a signed 64-bit integer.

The physical position at logical index ``(i,j,k)`` is

.. code-block:: text

   origin + i*spacing(1)*direction(:,1)
          + j*spacing(2)*direction(:,2)
          + k*spacing(3)*direction(:,3)

A singleton axis has equal lower and upper point bounds. A two-dimensional
XY grid, for example, may have point bounds ``[0,nx-1]``, ``[0,ny-1]``, and
``[0,0]``. Its point arrays have spatial shape ``(nx,ny,1)`` and its cell arrays
have spatial shape ``(nx-1,ny-1,1)``. Other choices of singleton axes are
supported by the writer, subject to the reader constraints below. A single
isolated point is not supported.

.. _image-reader-compatibility:

.. warning::

   Readers earlier than VTK 9.7.2 (ParaView 6.2.1) have the following limitations:

   * Use zero-based global logical point lower bounds on every axis. Positive
     global lower bounds, including one-based indexing, can cause out-of-range
     data reads. **This restriction applies only to ImageData**; RectilinearGrid
     and StructuredGrid support arbitrary global logical lower bounds.
   * Keep singleton axes trailing: use ``(nx,ny,1)`` for 2D grids and
     ``(nx,1,1)`` for 1D grids. Other arrangements may produce incorrect results.

   These constraints apply to the global logical grid. A nonzero physical
   ``origin`` and nonzero local ``array_origin`` or ``write_bounds`` remain
   supported with a zero-based global grid.

Examples
--------

The following serial example defines an XY grid with eight by six points:

.. code-block:: fortran

   use vtkhdf_image_file_type
   type(vtkhdf_image_file) :: file
   integer :: point_bounds(3,2)
   integer :: stat
   character(:), allocatable :: errmsg

   call file%create('image.vtkhdf',stat,errmsg)
   if (stat /= 0) error stop errmsg
   point_bounds(:,1) = [0,0,0]
   point_bounds(:,2) = [7,5,0]
   call file%write_mesh(point_bounds, origin=[0.0,0.0,0.0], spacing=[0.1,0.2,1.0])
   ! Write point, cell, or field data here.
   call file%close()

:repo:`example/image_fixed_mesh.F90` demonstrates static cell data, temporal point
and field data, and shared-point ownership using ``write_bounds``. It supports
serial and MPI builds.
