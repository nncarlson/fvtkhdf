=========
ImageData
=========

``vtkhdf_image_file_type`` provides ``vtkhdf_image_file`` for a single regular
structured grid. Geometry is fixed for the lifetime of the file. Both serial
and MPI builds support static and temporal point, cell, and field data.
RectilinearGrid, StructuredGrid, and ImageData blocks in MultiBlock files are
not implemented by this interface.

Creation and geometry
=====================

``create`` follows the UnstructuredGrid convention, including ``stat``,
``errmsg``, MPI communicator overloads, and optional ``data_chunk_bytes``.
The supported writer modes, exported by ``vtkhdf_image_file_type`` and defined
in ``vtkhdf_file_modes``, are:

* ``VTKHDF_STATIC`` (default): geometry and static data.
* ``VTKHDF_FIXED_MESH``: fixed geometry, static mesh-associated data, and
  registered temporal data.

.. important::

   For VTK 9.7.0/9.7.1 or a ParaView build using the affected reader, use
   global point lower bounds ``[0,0,0]``. One-based whole extents such as
   ``[1,nx,1,ny,1,nz]`` fail in those readers. The example below uses zero-based
   bounds. See :ref:`image-reader-compatibility` for the current constraints
   and the verified upstream fixes.

For example, in a serial build:

.. code-block:: fortran

   use,intrinsic :: iso_fortran_env, only: real64
   use vtkhdf_image_file_type
   type(vtkhdf_image_file) :: file
   integer :: point_bounds(2,3), stat
   character(:), allocatable :: errmsg

   call file%create('image.vtkhdf',stat,errmsg,mode=VTKHDF_FIXED_MESH)
   if (stat /= 0) error stop errmsg
   point_bounds(:,1) = [0,7]
   point_bounds(:,2) = [0,5]
   point_bounds(:,3) = [0,0]
   call file%write_mesh(point_bounds, &
       origin=[0.0_real64,0.0_real64,0.0_real64], &
       spacing=[0.1_real64,0.2_real64,1.0_real64])

In an MPI build, pass the communicator immediately after the filename.
All ranks provide identical geometry. Define it exactly once, before writing
or registering mesh-associated fields or beginning time steps.

``point_bounds(2,3)`` contains inclusive global logical point indices: lower
bounds in row 1 and upper bounds in row 2. Logical indices and the placement
arguments below are ``integer(int32)``. Internal dimension arithmetic, HDF5
selections, and temporal offsets use wider integers. ``WholeExtent`` is stored
as six 32-bit integers, ordered ``[ilo,ihi,jlo,jhi,klo,khi]``.

``origin(3)`` and ``spacing(3)`` accept either ``real32`` or ``real64``; use the
same kind for the optional ``direction(3,3)``. Geometry attributes are stored
as double precision. Direction defaults to identity. Its columns describe
the logical axes, so the physical position at logical index ``(i,j,k)`` is

.. code-block:: text

   origin + i*spacing(1)*direction(:,1)
          + j*spacing(2)*direction(:,2)
          + k*spacing(3)*direction(:,3)

Origin is the location of **logical index zero**, including when zero lies
outside the global bounds. Direction is preserved without normalization or
orthogonality requirements. Non-unit column lengths scale physical spacing.
The nine entries are serialized in row-major order for VTK.

Cells use the index of their lower-index corner. In each nonsingleton
direction, the cell upper bound is the point upper bound minus one. A
singleton direction retains its point bound. Thus the example above has
point dimensions ``(8,6,1)`` and cell dimensions ``(7,5,1)``. Images with one,
two, or three active axes are supported, with any choice of singleton axes.
A single isolated point is outside this interface's scope. Each axis dimension
must fit in VTK's 32-bit dimension representation, and the total point count
must fit in signed 64-bit arithmetic.

Point and cell arrays
=====================

Supported numeric kinds are ``int8``, ``int32``, ``int64``, ``real32``, and
``real64``. As in UG, ``int8`` values are stored using the unsigned-byte HDF5
type, suitable for VTK visibility and ghost arrays.

* Scalars have shape ``(nx,ny,nz)``.
* Multicomponent arrays have shape ``(ncomp,nx,ny,nz)``.

Keep all three spatial dimensions in public arrays, including singleton
directions. Only the file representation omits singleton spatial axes;
components remain the fastest-varying dimension. An active point direction
with two points retains its spatial dimension in cell datasets even though
that direction has only one cell.

Every write supplies ``array_origin(3)``, the global logical indices of the
array's first spatial element. Optional ``write_bounds(2,3)`` specifies an
inclusive rectangular subset in the applicable cell or point index space.
Without it, the entire supplied array is written.

.. code-block:: fortran

   call file%write_cell_data('Pressure',pressure, &
       array_origin=cell_origin,write_bounds=cell_bounds,attribute='Scalars')

Each name-based static point/cell write creates a new field, with one block
contributed by each process. Repeated names and names that collide after
sanitization receive distinct internal HDF5 names, following UG's best-effort
naming policy. Static and temporal fields may share an external name.

To contribute additional blocks to a static field, request its handle on the
first write and use that handle on subsequent writes:

.. code-block:: fortran

   type(vtkhdf_cell_data_handle) :: pressure_var
   call file%write_cell_data('Pressure',pressure_a,origin_a,bounds_a, &
       attribute='Scalars',handle=pressure_var)
   call file%write_cell_data(pressure_var,pressure_b,origin_b,bounds_b)

The optional ``attribute`` designation follows UG's conventions and is set
when creating the field. All blocks of a field must have the same numeric
kind, scalar/component interpretation, and component count. The handle types
are shared with temporal fields; each handle identifies its field independently
of the user-supplied name. The optional ``handle`` output is available for
both scalar and multicomponent cell/point writes.

Actual Fortran subscripts are independent of placement. These two arrays
represent the same global region:

.. code-block:: fortran

   real(real64) :: local(4,3,1), global(10:13,-2:0,5:5)
   ! Fill LOCAL and GLOBAL with the same values.
   call file%write_cell_data('local_copy',local,array_origin=[10,-2,5])
   call file%write_cell_data('global_copy',global,array_origin=[10,-2,5])

For arrays that include halos, supply their global logical origin and select
only the owned interior:

.. code-block:: fortran

   real(real64) :: pressure(0:5,0:4,1)
   integer :: owned(2,3)
   ! Interior subscripts (1:4,1:3,1) correspond to global (10:13,-2:0,0).
   owned(:,1) = [10,13]
   owned(:,2) = [-2,0]
   owned(:,3) = [0,0]
   call file%write_cell_data('Pressure',pressure, &
       array_origin=[9,-3,0],write_bounds=owned)

The HDF5 memory selection excludes halos without explicitly packing an
interior buffer. Array dummies are contiguous; noncontiguous Fortran actual
arguments may cause compiler-generated temporaries. Passing a contiguous
allocation containing halos avoids that additional packing.

Collective blocks and ownership
===============================

All ranks call compatible operations in the same order for the same field,
with identical numeric/component characteristics. Each call may select a
different rectangular region on each rank. There is no requirement for a
Cartesian processor arrangement, and a rank can contribute several blocks.

The caller owns the decomposition: every required point or cell value must
be supplied exactly once. Assign shared boundary points to one owner, omit
computational halos, and write any desired blanking/visibility arrays
explicitly. The library checks local selections and lifecycle state; it does
not detect holes, overlap, or complete global coverage.

An explicit bound with upper less than lower in any direction is an empty
contribution. Its array may have zero spatial size. It still needs a valid
numeric kind and component count. Use empty calls to pad unequal local block
counts:

.. code-block:: fortran

   type(vtkhdf_point_data_handle) :: temperature_var
   do b = 1, max_blocks_across_ranks
     if (b == 1) then
       if (b <= local_blocks) then
         call file%write_point_data('Temperature',block(b)%values, &
             block(b)%array_origin,block(b)%owned_bounds,handle=temperature_var)
       else
         ! EMPTY has shape (0,0,0); EMPTY_BOUNDS(:,1) = [0,-1].
         call file%write_point_data('Temperature',empty,[0,0,0],empty_bounds, &
             handle=temperature_var)
       end if
     else if (b <= local_blocks) then
       call file%write_point_data(temperature_var,block(b)%values, &
           block(b)%array_origin,block(b)%owned_bounds)
     else
       call file%write_point_data(temperature_var,empty,[0,0,0],empty_bounds)
     end if
   end do

The caller determines the padded call count. There is no block object or
special empty-write method in the library. An empty first contribution still
creates the field and returns a valid handle on that rank.

Temporal snapshots
==================

Register fields after geometry and before starting the first step:

.. code-block:: fortran

   type(vtkhdf_cell_data_handle) :: pressure_var
   pressure_var = file%register_temporal_cell_data('Pressure',0.0_real64)

A scalar mold registers scalar data; a rank-one mold registers its number
of components. Mold values are ignored. Registration retains no placement
arguments. Handles identify fields and can be used for multiple blocks:

.. code-block:: fortran

   call file%start_time_step(time)
   call file%write_cell_data(pressure_var,pressure_a,origin_a,bounds_a)
   call file%write_cell_data(pressure_var,pressure_b,origin_b,bounds_b)
   call file%finalize_time_step()

For each temporal point/cell field:

* No calls during a step reuse the preceding snapshot and its stored offset.
* One or more calls assemble a **complete replacement** snapshot. Allocation
  occurs once on the first call, and subsequent calls write into that slice.
  Unwritten regions do not inherit preceding values.

Every registered point/cell field must begin an initial snapshot in the first
step. Once a replacement begins, the caller must supply all its global values
exactly once before finalizing. Blocks and rank ownership may change between
snapshots. An all-empty collective call still begins a replacement: omit the
field's calls collectively to request reuse.

Static point/cell fields in temporal files have one stored slice and zero
offsets at every step. They can be written after mesh definition, including
after time stepping begins. Starting a new step or closing the file finalizes
an open step, as in UG. Call ``close`` explicitly for collective cleanup;
implicit Fortran finalization performs no collective operations.

FieldData
=========

FieldData follows the existing UG API and behavior independently of the block
rules above:

* Scalar, rank-one, and rank-two values accept ``int32``, ``int64``, ``real32``,
  and ``real64``. ``as_vector`` gives a rank-one value multiple components;
  otherwise it represents scalar tuples. Rank-two arrays are component first.
* Calls are collective, and only rank 0 contributes values.
* Each static write supplies the complete field and creates a unique dataset,
  including repeated names. There is no returned handle or incremental write.
  Static FieldData is available only in ``VTKHDF_STATIC`` mode.
* Temporal registration uses a scalar mold. Each field may be written at most
  once per step. Omission reuses the previous value; initial omission is an
  empty value. Shapes may change between writes, matching UG.
* FieldData has no ``array_origin`` or ``write_bounds`` arguments and no
  full-first-snapshot requirement.

Temporal FieldData uses the same tuple-by-component storage as UG. Both
dataset dimensions are extendible, so component and tuple counts may change
between steps without declaring a maximum at registration. Component storage
grows to the largest count encountered; unused entries are padding.
``Steps/FieldDataSizes`` records each step's actual component and tuple counts,
and ``Steps/FieldDataOffsets`` measures offsets in tuple rows. Readers select
only the valid rectangle and exclude padding. Omitted writes repeat the
previous offset and sizes without appending values.

.. _image-reader-compatibility:

Storage and reader compatibility
================================

Point/cell datasets use disk ordering ``[z,y,x,(components)]`` with singleton
spatial axes omitted. Temporal files prepend a stored-snapshot dimension,
including for static fields. Repeated offsets represent reuse without copying.

Datasets are chunked deterministically. The policy fills chunks along the
fastest disk dimensions first, bounded by ``data_chunk_bytes`` (128 KiB by
default, or at least one numeric value). A temporal chunk spans one snapshot.
This is a general default, with no public chunk-shape interface.

Current constraints for affected readers
----------------------------------------

The installed VTK releases 9.7.0 and 9.7.1 have geometry reader defects.
For these readers, including ParaView builds that use the affected code:

* Use zero-based global point bounds on every axis, for example
  ``WholeExtent = [0,nx-1,0,ny-1,0,nz-1]``. Positive global lower bounds,
  including the common one-based convention, cause out-of-range field reads.
  A negative-bound full-domain control passed, but this does not establish
  correct subextent reads; use zero-based bounds for affected readers.
* Keep singleton axes trailing: ``(nx,ny,1)`` for 2D and ``(nx,1,1)`` for 1D.
  Singleton x/y cases can fail or silently produce incorrect tuple counts.

These are reader constraints on the global logical extent. A nonzero physical
``origin`` and positive ``array_origin`` or ``write_bounds`` for local blocks
are supported with a zero-based whole extent.

fVTKHDF always maps the supplied ``point_bounds`` directly to ``WholeExtent``
and preserves the supplied physical origin. There is no automatic index
rebasing or compatibility option. The library continues to support nonzero
global bounds and arbitrary singleton-axis choices for readers containing the
upstream fixes.

Verified upstream fixes
-----------------------

Both geometry reader defects are fixed in the unmodified VTK development
snapshot reporting ``9.7.20261022``
(`commit 9f0588e2fe0d <https://github.com/Kitware/VTK/commit/9f0588e2fe0d10573555834a66ac3dab95be4114>`_).
All ImageData geometry and mesh-associated data cases passed with a local
build of that source. The tested snapshot's
`ReadStructuredData implementation <https://github.com/Kitware/VTK/blob/9f0588e2fe0d10573555834a66ac3dab95be4114/IO/HDF/vtkHDFReader.cxx#L272-L320>`_
subtracts global lower bounds when forming file offsets and handles each
singleton axis independently.

A VTK or ParaView build containing those fixes removes the geometry constraints
above. The first released versions containing them have not been verified;
no specific future release is promised here. No ParaView release was tested
for this milestone.

The failures were reproduced with independent small HDF5 fixtures. The test
suite probes reader capability and skips affected semantic checks while
retaining direct HDF5 layout checks. A standalone reproducer, measured results,
and source analysis are recorded in
``vtk-9.7.1-vtkhdf-reader-failures.md`` at the repository root.

A complete serial/MPI temporal example is provided in
``example/image_fixed_mesh.F90``. It shows file creation, global mesh
definition, registration, disjoint rank slabs, omitted-write reuse, FieldData,
and explicit close.

Complete example
================

.. literalinclude:: ../../example/image_fixed_mesh.F90
   :language: fortran
