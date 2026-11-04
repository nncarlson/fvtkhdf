Structured VTK Dataset Types
============================

Structured VTK datasets use a logical grid indexed by three integer coordinates.
Their point and cell data are written as rectangular regions of that grid.
The interface below is shared by the ImageData, RectilinearGrid, and
StructuredGrid writers.

See :doc:`image`, :doc:`rectilinear`, and :doc:`structured_grid` for geometry
definition, supported modes, and reader requirements.

File Creation and Management
----------------------------

The declaration of ``file`` and its geometry-writing method depend on the
VTK dataset type: ``vtkhdf_image_file``, ``vtkhdf_rectilinear_file``, or
``vtkhdf_structured_file``, from the corresponding ``*_file_type`` module.

.. method-doc::

   call file%create(filename, [comm,] stat, errmsg [, mode] [, data_chunk_bytes])

   Create a new file.

   * ``filename``: path to the file to create. The recommended extension is
     ``.vtkhdf``.
   * ``comm``: MPI communicator, omitted in serial builds. With
     ``USE_MPI_F08=ON`` (default), either ``integer`` or ``type(MPI_Comm)`` is
     accepted. With ``USE_MPI_F08=OFF``, only ``integer`` is accepted.
   * ``stat``: output status code; zero indicates success.
   * ``errmsg``: allocated error message when ``stat /= 0``.
   * ``mode``: optional supported writer mode; defaults to ``VTKHDF_STATIC``.
   * ``data_chunk_bytes``: optional positive ``int64`` target chunk
     size in bytes for point and cell data; defaults to ``128*1024`` (128 KiB).

   In MPI builds, all input arguments are identical on all ranks, and all
   ranks return identical ``stat`` and ``errmsg`` values.

.. method-doc::

   call file%flush()

   Flush the file's HDF5 buffers.

.. method-doc::

   call file%close()

   Finalize any open time step, close the file, and release resources. Users
   must call this explicitly; automatic Fortran finalization cannot perform
   collective cleanup.

Grid Indices and Data Layout
----------------------------

Global point bounds specify inclusive logical indices in each of three
directions. Point and cell indices, ``array_origin``, and ``write_bounds``
use kind ``int32``. Cells are indexed by their lower-index corner.
For example, point indices ``3:8`` in one direction give cell indices
``3:7`` in that direction. In a singleton direction, points and cells use
the same single index.

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
and ``real64`` from ``iso_fortran_env``. All writes to one dataset must have
the same type, kind, scalar/multicomponent interpretation, and component
count. The component count must be positive.

Unwritten point and cell data values default to quiet NaN for real kinds
and zero for integer kinds. NaN makes missing real data
distinguishable from a valid zero; it does not mark a location hidden
(:ref:`see below <structured-hidden-locations>`).

.. warning::

   ``int8`` data uses VTK/ParaView-compatible unsigned-byte storage. Use only
   nonnegative values; negative ``int8`` values are not a portable interface.

Static Point and Cell Data
--------------------------

.. code-block:: fortran

   type(vtkhdf_cell_data_handle) :: cell_var
   type(vtkhdf_point_data_handle) :: point_var

.. method-doc::

   call file%write_cell_data(name, array, array_origin [, write_bounds] [, attribute] [, handle])
   call file%write_point_data(name, array, array_origin [, write_bounds] [, attribute] [, handle])

   Create a static dataset and write one rectangular region per calling rank.

   * ``name``: character name; see :doc:`common` for naming rules.
   * ``array``: rank-3 scalar or rank-4 multicomponent data, as above.
   * ``array_origin(3)``: global logical indices of the first spatial
     element of ``array``, independently of its Fortran lower bounds.
   * ``write_bounds(3,2)``: optional inclusive global bounds selecting the
     region to write. Column 1 gives lower bounds and column 2 gives upper
     bounds. The selection must lie within the global bounds for the association
     and within the global region represented by the array, from ``array_origin``
     through ``array_origin + spatial_shape - 1``, where ``spatial_shape`` gives
     the array's three spatial extents. If omitted, the entire array is written.
   * ``attribute``: optional character VTK role; see :doc:`common`.
   * ``handle``: optional output handle of the corresponding cell or point
     data type, for subsequent writes to this dataset.

   For example, with an actual array declared ``a(0:4,1,1)``,
   ``array_origin=[10,0,0]`` places ``a(0,1,1)`` at global index ``(10,0,0)``.
   A ``write_bounds`` selection from ``[11,0,0]`` to ``[13,0,0]`` writes
   ``a(1:3,1,1)``.

   Each name-based call creates a new dataset. Repeating the name creates
   another dataset with a disambiguated name.

.. method-doc::

   call file%write_cell_data(cell_var, array, array_origin [, write_bounds])
   call file%write_point_data(point_var, array, array_origin [, write_bounds])

   Write another region to the dataset identified by a handle. Placement and
   array arguments have the same meaning as in the name-based calls.

   Static point and cell data are supported in ``VTKHDF_STATIC`` and
   ``VTKHDF_FIXED_MESH`` modes. Define the mesh before creating or writing these
   datasets. In a temporal file, static datasets may also be written after time
   stepping begins; they remain independent of time.

Region Ownership and Collective Writes
--------------------------------------

Each rank may contribute different rectangular regions, and may contribute
multiple regions to one dataset. No Cartesian arrangement of ranks is required.
All ranks must make the same sequence of calls for the same datasets, with
matching numeric and component characteristics. Ranks with fewer regions
must therefore supply :ref:`empty contributions <structured-empty-regions>`
to complete that sequence.

For each dataset, the caller must supply data for every visible global point
or cell index exactly once. Hidden locations may be left unwritten. Overlapping
writes are forbidden; shared boundary points must have one writer. This rule
also applies to boundary-condition regions included in the global grid.
The library does not detect overlaps or incomplete coverage.

For an array containing halos, pass the full contiguous allocation and use
``write_bounds`` to select its owned interior. Passing a noncontiguous array
section may cause a compiler-generated copy.

.. _structured-empty-regions:

Empty Regions
~~~~~~~~~~~~~

To contribute an empty region, pass an array with zero extent in any spatial
dimension and omit ``write_bounds``. Alternatively, supply ``write_bounds``
with an upper bound less than its lower bound in any direction; the array
may then be nonempty.

``array_origin`` is still required, but its values are not checked against
the global grid for an empty region. Use ``[0,0,0]``, even if zero lies outside
the grid. The array must still have the correct numeric type, kind, and
scalar/multicomponent interpretation, including the component count.
An empty first call still creates the dataset and returns a valid handle
when requested.

.. _structured-hidden-locations:

Hidden Locations
~~~~~~~~~~~~~~~~

The global grid is rectangular, but parts of it may be hidden to represent
holes or exclude boundary-condition regions from display. Write scalar
``int8`` point and/or cell data arrays named ``vtkGhostType``, as needed,
using the usual write methods. Each array specifies hidden locations for
its own association. No ``attribute`` argument is needed. Use these values:

.. list-table::
   :header-rows: 1

   * - Association
     - Visible
     - Hidden
   * - Point
     - ``0``
     - ``2`` (``HIDDENPOINT``)
   * - Cell
     - ``0``
     - ``32`` (``HIDDENCELL``)

For example, a cell mask for a region starting at global index ``(0,0,0)``
could be written as follows:

.. code-block:: fortran

   integer(int8) :: ghost(8,6,1)

   ghost = 0_int8
   ghost(3:4,2:3,1) = 32_int8
   call file%write_cell_data('vtkGhostType',ghost,array_origin=[0,0,0])

Write the mask for every location whose visibility must be specified.
Unwritten mask entries default to zero (visible). Other datasets may be left
unwritten at hidden locations; HDF5 supplies their fill values. NaN values
alone do not mark locations hidden.

Temporal Point and Cell Data
----------------------------

.. code-block:: fortran

   type(vtkhdf_cell_data_handle) :: cell_var
   type(vtkhdf_point_data_handle) :: point_var

.. method-doc::

   cell_var = file%register_temporal_cell_data(name, mold [, attribute])
   point_var = file%register_temporal_point_data(name, mold [, attribute])

   Register a time-dependent dataset and return its opaque handle.

   * ``name``: character name; see :doc:`common`.
   * ``mold``: scalar for scalar data, or rank-1 array whose size is the
     number of components for multicomponent data. Its type, kind, and
     scalar/component interpretation determine those of later writes;
     its values are ignored.
   * ``attribute``: optional character VTK role, set at registration.

Define the mesh before registration. Complete all registration before
the first time step starts. Placement is supplied on each write.

.. method-doc::

   call file%start_time_step(time)

   Start a new time step with the scalar ``real32`` or ``real64`` time value
   ``time``. In MPI builds, all ranks supply the same value.

.. method-doc::

   call file%write_cell_data(cell_var, array, array_origin [, write_bounds])
   call file%write_point_data(point_var, array, array_origin [, write_bounds])

   Write a region of the temporal dataset for the current time step. The
   array and placement arguments have the same meaning as for static
   datasets, and must match the registered numeric and component
   characteristics. Multiple calls may assemble one complete dataset.

.. method-doc::

   call file%finalize_time_step()

   Complete the current time step after all its temporal writes.

Every registered point or cell dataset must be written in the first time step.
In each subsequent step, the caller must either omit all writes to the dataset,
reusing its previous complete value, or supply a complete replacement. An
unwritten region in a replacement does not inherit values from the previous
step. All-empty collective writes still begin a replacement; reuse requires
omitting the calls on all ranks.

Region shapes and rank ownership may change between replacements. Every
replacement must cover the entire dataset exactly once before the step is
finalized.

Field Data
----------

See :doc:`field_data` for the shared interface. Field data writes supply a
complete value; they have no region-placement arguments and do not use the
incremental point/cell writing rules.

Example: Static Data with Halos
-------------------------------

After creating the file and defining a sufficiently large mesh, the following
calls create a cell dataset and contribute two regions. The first array contains
halos excluded by ``owned``:

.. code-block:: fortran

   real(real64) :: pressure(0:5,0:4,1)
   integer(int32) :: owned(3,2)
   type(vtkhdf_cell_data_handle) :: pressure_var

   ! Interior subscripts (1:4,1:3,1) map to global indices (10:13,-2:0,0).
   owned(:,1) = [10,-2,0]
   owned(:,2) = [13,0,0]
   call file%write_cell_data('Pressure',pressure, &
       array_origin=[9,-3,0],write_bounds=owned,handle=pressure_var)
   call file%write_cell_data(pressure_var,other_pressure, &
       array_origin=other_origin,write_bounds=other_owned)
