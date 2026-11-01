The vtkhdf_image_file_type module
=================================

This module defines the ``vtkhdf_image_file`` derived type for writing
VTKHDF ImageData files. ImageData represents a regular structured grid with
fixed geometry. One-, two-, and three-dimensional grids are supported.

File management, point and cell data, and time stepping are described in
:doc:`structured`. Field data is described in :doc:`field_data`.

.. code-block:: fortran

   use vtkhdf_image_file_type
   type(vtkhdf_image_file) :: file

Writer Modes
------------

The module exports the following mode constants:

* ``VTKHDF_STATIC`` (default): fixed geometry and static point, cell, and
  field data; no time stepping.
* ``VTKHDF_FIXED_MESH``: fixed geometry, static or temporal point and cell
  data, and temporal field data.

The mesh geometry cannot change during the lifetime of the file.

Mesh Data
---------

``call file%write_mesh(point_bounds, origin, spacing [, direction])``
   Define the global grid geometry. Call exactly once, before writing or
   registering point/cell data or starting time stepping. All ranks must
   supply identical arguments.

   * ``point_bounds(2,3)``: inclusive global logical point bounds of type
     ``integer(int32)``. Row 1 contains lower bounds and row 2 contains upper
     bounds, in x, y, z order. Each upper bound must be at least its lower
     bound, and at least one direction must contain more than one point.
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

Reader Compatibility
--------------------

VTK 9.7.0 and 9.7.1 have ImageData reader defects. For these versions and
ParaView builds using the affected reader:

* Use zero-based global logical point lower bounds on every axis. Positive
  global lower bounds, including one-based indexing, can cause out-of-range
  data reads.
* Keep singleton axes trailing: use ``(nx,ny,1)`` for two-dimensional grids
  and ``(nx,1,1)`` for one-dimensional grids. Other singleton-axis choices
  can fail or produce incorrect tuple counts.

These constraints apply to the global logical grid. A nonzero physical
``origin`` and nonzero local ``array_origin`` or ``write_bounds`` remain
supported with a zero-based global grid. fVTKHDF preserves the supplied
bounds and does not automatically rebase indices.

Both reader defects are fixed by upstream VTK
`commit 32b8631901b9 <https://gitlab.kitware.com/vtk/vtk/-/commit/32b8631901b90fac46c470a312c30221c452bd54>`_
(`merge request 13507 <https://gitlab.kitware.com/vtk/vtk/-/merge_requests/13507>`_).
The fixes have been backported to the VTK 9.7 release branch for VTK 9.7.2.
A VTK or ParaView build containing the fixes removes these constraints.

Example
-------

The following serial example defines an XY grid with eight by six points:

.. code-block:: fortran

   use, intrinsic :: iso_fortran_env, only: int32, real64
   use vtkhdf_image_file_type
   type(vtkhdf_image_file) :: file
   integer(int32) :: point_bounds(2,3)
   integer :: stat
   character(:), allocatable :: errmsg

   call file%create('image.vtkhdf',stat,errmsg)
   if (stat /= 0) error stop errmsg
   point_bounds(:,1) = [0,7]
   point_bounds(:,2) = [0,5]
   point_bounds(:,3) = [0,0]
   call file%write_mesh(point_bounds, &
       origin=[0.0_real64,0.0_real64,0.0_real64], &
       spacing=[0.1_real64,0.2_real64,1.0_real64])
   ! Write point, cell, or field data here.
   call file%close()

Complete Temporal Example
-------------------------

``example/image_fixed_mesh.F90`` supports serial and MPI builds and demonstrates
geometry definition, temporal registration, rank-local regions, omitted-write
reuse, field data, and explicit close.

.. literalinclude:: ../../example/image_fixed_mesh.F90
   :language: fortran
