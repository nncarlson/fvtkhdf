Structured VTK Dataset Types
============================

Structured VTK datasets use a logical grid indexed by three integer coordinates.
Their point and cell data are written as rectangular regions of that grid.
The interface below is provided by the ImageData writer. RectilinearGrid and
StructuredGrid writers are not yet available.

See :doc:`image` for geometry definition, supported modes, and reader
requirements.

File Creation and Management
----------------------------

The declaration of ``file`` and its geometry-writing method depend on the
VTK dataset type. For ImageData, use ``vtkhdf_image_file`` from
``vtkhdf_image_file_type``.

``call file%create(filename, [comm,] stat, errmsg [, mode] [, data_chunk_bytes])``
   Create a new file.

   * ``filename``: path to the file to create. The recommended extension is
     ``.vtkhdf``.
   * ``comm``: MPI communicator, omitted in serial builds. With
     ``USE_MPI_F08=ON`` (default), either ``integer`` or ``type(MPI_Comm)`` is
     accepted. With ``USE_MPI_F08=OFF``, only ``integer`` is accepted.
   * ``stat``: output status code; zero indicates success.
   * ``errmsg``: allocated error message when ``stat /= 0``.
   * ``mode``: optional supported writer mode; defaults to ``VTKHDF_STATIC``.
   * ``data_chunk_bytes``: optional positive ``integer(int64)`` target chunk
     size in bytes for point and cell data; defaults to ``128*1024`` (128 KiB).

   In MPI builds, all input arguments are identical on all ranks, and all
   ranks return identical ``stat`` and ``errmsg`` values.

``call file%flush()``
   Flush the file's HDF5 buffers.

``call file%close()``
   Finalize any open time step, close the file, and release resources. Call
   this explicitly; automatic Fortran finalization cannot perform collective
   cleanup.

Grid Indices and Data Layout
----------------------------

Global point bounds specify inclusive logical indices in each of three
directions. Point and cell indices, ``array_origin``, and ``write_bounds``
use ``integer(int32)``. For a direction with more than one point, cells are
indexed by their lower-index corner: the cell lower bound equals the point
lower bound, and the cell upper bound is one less than the point upper bound.
A singleton direction retains its point bound for cell indexing.

Point and cell arrays use the following shapes, including all three spatial
dimensions even for one- or two-dimensional grids:

.. list-table::
   :header-rows: 1

   * - Data
     - Fortran shape
   * - Scalar
     - ``(nx,ny,nz)``
   * - Multicomponent
     - ``(ncomp,nx,ny,nz)``

The supported numeric kinds are ``int8``, ``int32``, ``int64``, ``real32``,
and ``real64`` from ``iso_fortran_env``. All writes to one field must have
the same type, kind, scalar/multicomponent interpretation, and component
count. The component count must be positive.

.. warning::

   ``int8`` data uses VTK/ParaView-compatible unsigned-byte storage. Use only
   nonnegative values; negative ``int8`` values are not a portable interface.

Static Point and Cell Data
--------------------------

.. code-block:: fortran

   type(vtkhdf_cell_data_handle) :: cell_var
   type(vtkhdf_point_data_handle) :: point_var

.. code-block:: text

   call file%write_cell_data(name, array, array_origin [, write_bounds] [, attribute] [, handle])
   call file%write_point_data(name, array, array_origin [, write_bounds] [, attribute] [, handle])

Create a static field and write one rectangular region per calling rank.

* ``name``: character name; see :doc:`common` for naming rules.
* ``array``: rank-3 scalar or rank-4 multicomponent data, as above.
* ``array_origin(3)``: global logical indices of the first spatial
  element of ``array``, independently of its Fortran lower bounds.
* ``write_bounds(2,3)``: optional inclusive global bounds selecting the
  region to write. Row 1 gives lower bounds and row 2 gives upper
  bounds. The selection must lie within both the supplied array and
  the global bounds for the association. If omitted, the entire array
  is written.
* ``attribute``: optional character VTK role; see :doc:`common`.
* ``handle``: optional output handle of the corresponding cell or point
  data type, for subsequent writes to this field.

Each name-based call creates a new field. Repeating the name creates
another field with a disambiguated name.

.. code-block:: text

   call file%write_cell_data(cell_var, array, array_origin [, write_bounds])
   call file%write_point_data(point_var, array, array_origin [, write_bounds])

Write another region to the field identified by a handle. Placement and
array arguments have the same meaning as in the name-based calls.

Static point and cell data are supported in ``VTKHDF_STATIC`` and
``VTKHDF_FIXED_MESH`` modes. Define the mesh before creating or writing these
fields. In a temporal file, static fields may also be written after time
stepping begins; they remain independent of time.

Region Ownership and Collective Writes
--------------------------------------

Each rank may contribute different rectangular regions, and may contribute
multiple regions to one field. No Cartesian arrangement of ranks is required.
All ranks must make the same sequence of calls for the same fields, with
matching numeric and component characteristics.

The caller must supply every required point or cell exactly once. Shared
boundary points must have one owner, and computational halos must be excluded
unless intentionally part of the output. The library does not detect global
holes, overlaps, or incomplete coverage.

A selection with an upper bound less than its lower bound in any direction
is an empty contribution. Arrays with zero spatial size are also allowed.
Empty contributions retain the field's numeric and component characteristics.
Ranks with fewer regions must make empty calls so that every rank makes the
same number of collective calls. An empty first call still creates the field
and returns a valid handle when requested.

For an array containing halos, pass the full contiguous allocation and use
``write_bounds`` to select its owned interior. Passing a noncontiguous array
section may cause a compiler-generated copy.

Temporal Point and Cell Data
----------------------------

.. code-block:: text

   cell_var = file%register_temporal_cell_data(name, mold [, attribute])
   point_var = file%register_temporal_point_data(name, mold [, attribute])

Register a time-dependent field and return its opaque handle.

* ``name``: character name; see :doc:`common`.
* ``mold``: scalar for scalar data, or rank-1 array whose size is the
  number of components for multicomponent data. Its type, kind, and
  scalar/component interpretation determine those of later writes;
  its values are ignored.
* ``attribute``: optional character VTK role, set at registration.

Define the mesh before registration. Complete all registration before
the first time step starts. Placement is supplied on each write.

``call file%start_time_step(time)``
   Start a new time step with the scalar ``real32`` or ``real64`` time value
   ``time``. In MPI builds, all ranks supply the same value.

.. code-block:: text

   call file%write_cell_data(cell_var, array, array_origin [, write_bounds])
   call file%write_point_data(point_var, array, array_origin [, write_bounds])

Write a region of the temporal field for the current time step. The
array and placement arguments have the same meaning as for static
fields, and must match the registered numeric and component
characteristics. Multiple calls may assemble one complete field.

``call file%finalize_time_step()``
   Complete the current time step after all its temporal writes.

Every registered point or cell field must be written in the first time step.
In each subsequent step, the caller must either omit all writes to the field,
reusing its previous complete value, or supply a complete replacement. An
unwritten region in a replacement does not inherit values from the previous
step. All-empty collective writes still begin a replacement; reuse requires
omitting the calls on all ranks.

Region shapes and rank ownership may change between replacements. Every
replacement must cover the entire field exactly once before the step is
finalized.

Field Data
----------

See :doc:`field_data` for the shared interface. Field data writes supply a
complete value; they have no region-placement arguments and do not use the
incremental point/cell writing rules.

Example: Static Data with Halos
-------------------------------

After creating the file and defining a sufficiently large mesh, the following
calls create a cell field and contribute two regions. The first array contains
halos excluded by ``owned``:

.. code-block:: fortran

   real(real64) :: pressure(0:5,0:4,1)
   integer(int32) :: owned(2,3)
   type(vtkhdf_cell_data_handle) :: pressure_var

   ! Interior subscripts (1:4,1:3,1) map to global indices (10:13,-2:0,0).
   owned(:,1) = [10,13]
   owned(:,2) = [-2,0]
   owned(:,3) = [0,0]
   call file%write_cell_data('Pressure',pressure, &
       array_origin=[9,-3,0],write_bounds=owned,handle=pressure_var)
   call file%write_cell_data(pressure_var,other_pressure, &
       array_origin=other_origin,write_bounds=other_owned)
