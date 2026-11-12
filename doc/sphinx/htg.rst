HyperTreeGrid
=============

The ``vtkhdf_htg_file_type`` module defines the ``vtkhdf_htg_file``
derived type for writing VTKHDF HyperTreeGrid files, VTK's representation
of tree-based adaptive mesh refinement (AMR).

A HyperTreeGrid is a rectilinear grid of *root cells*. Each root cell is the
root of a *tree*: a cell is either a leaf or is *refined* into children by
dividing it into ``branch_factor`` equal parts along every nonsingleton axis.
A refined cell thus has ``branch_factor**d`` children in a ``d``-dimensional
grid; the default ``branch_factor`` of 2 gives binary trees, quadtrees, and
octrees. Every cell of every tree, refined or leaf, carries cell data.
There is no point data.

Reading these files requires ParaView 6.2 or later (VTK 9.7 or later).

Optional VTK attribute roles, specified with the ``attribute`` argument, are
described in :doc:`common`.

In the MPI build, all type-bound procedures are collective over the
communicator passed to ``create``. Every rank must call the same method
in the same order and supply identical values for all non-distributed
arguments.

File Creation and Management
----------------------------

.. code-block:: fortran

   use vtkhdf_htg_file_type
   type(vtkhdf_htg_file) :: file

.. method-doc::

   call file%create(filename, [comm,] stat, errmsg [, mode] [, data_chunk_bytes])

   Create a new VTKHDF HyperTreeGrid file. The arguments are those of the
   :doc:`UnstructuredGrid <ug>` ``create`` method.

.. method-doc::

   call file%close()

   Finalize any open time step, close the file, and release resources. Call
   this explicitly; automatic Fortran finalization cannot perform collective
   cleanup.

.. method-doc::

   call file%flush()

   Flush the file's HDF5 buffers.

Writer Modes
~~~~~~~~~~~~
The writer supports all four modes:

* ``VTKHDF_STATIC`` (default): fixed mesh and static cell and field data;
  no time stepping.
* ``VTKHDF_FIXED_MESH``: fixed mesh and static or temporal cell and field
  data.
* ``VTKHDF_MOVING_MESH``: fixed trees with temporal root grid coordinates;
  cell data may be static or temporal.
* ``VTKHDF_DYNAMIC_MESH``: temporal trees and root grid coordinates; cell
  datasets must be temporal. This is the mode for adaptive refinement.

Field data may be static or temporal in all temporal modes. The dimensions
of the root grid, the branch factor, and the root indexing are fixed for the
lifetime of the file in every mode.

Mesh Data
---------

The mesh has three parts: the definition of the root grid, its geometry, and
its topology, the trees.

.. code-block:: fortran

   use vtkhdf_htg_file_type
   type(vtkhdf_htg_file) :: file

.. method-doc::

   call file%write_mesh(dims [, branch_factor] [, transposed_root_indexing])

   Define the grid of root cells. Call exactly once, before writing geometry
   or topology. All ranks supply identical values.

   * ``dims(3)``: the number of grid *points* along the x, y, and z axes, as
     for a RectilinearGrid; there are ``dims-1`` root cells along a
     nonsingleton axis. Use 1 for an unused axis of a one- or two-dimensional
     grid. At least one axis must have more than one point.
   * ``branch_factor``: the number of parts, 2 (default) or 3, into which a
     refined cell is divided along each nonsingleton axis.
   * ``transposed_root_indexing``: root cells are numbered with x varying
     fastest by default, and with z varying fastest when this is ``.true.``.

.. method-doc::

   call file%write_mesh_geometry(x, y, z)

   Write the coordinate vectors of the root grid.

   * ``x(:)``, ``y(:)``, ``z(:)``: point coordinates along each axis,
     supplied by rank zero, with lengths ``dims``. A singleton axis has one
     coordinate. Values and lengths on other ranks are ignored; zero-sized
     arrays suffice. All three vectors have the same kind, ``real32`` or
     ``real64``. Coordinates are stored as ``real64``.

.. method-doc::

   call file%write_mesh_topology(tree_ids, xtree, refined [, mask])

   Write the trees provided by the calling MPI rank. Each tree of the grid
   is supplied, complete, by exactly one rank; a rank may supply none.

   * ``tree_ids``: the root cell of each of the rank's trees, as ``int32`` or
     ``int64``. Root cells are numbered from 1, following Fortran conventions:
     root cell ``(i,j,k)`` of a grid with ``nx`` by ``ny`` by ``nz`` root
     cells has number ``i + (j-1)*nx + (k-1)*nx*ny`` under the default root
     indexing. Conversion to the 0-based VTKHDF representation is handled
     internally. A root cell without a tree is absent from the dataset.

   * ``xtree``: the starting index of each tree's cells in the rank's cell
     list, of the same kind as ``tree_ids``. Tree ``n`` consists of cells
     ``xtree(n):xtree(n+1)-1``. ``size(xtree)`` must equal
     ``size(tree_ids) + 1``, its first element is 1, and its final element
     equals the rank's number of cells plus 1.

   * ``refined``: a default ``logical`` array with one element per cell;
     ``.true.`` marks a refined cell and ``.false.`` a leaf.

   * ``mask``: an optional default ``logical`` array with one element per
     cell; ``.true.`` hides a leaf from view. A refined cell cannot be masked.
     All ranks supply ``mask`` or none does.

   The cells of a tree are listed in *breadth-first* order: the root cell,
   then all cells of the first level of refinement, then all cells of the
   second, and so on. Within a level, the children of each refined cell of
   the preceding level are contiguous, and follow the order of their parents.
   The children of a cell are ordered with x varying fastest, then y, then z,
   over the nonsingleton axes. For example, a quadtree whose root is refined,
   and whose second child is refined again, has 9 cells::

      refined = [.true., .false., .true., .false., .false., &
                 .false., .false., .false., .false.]

   The same order applies to cell data: the rank's cell data arrays hold one
   value per cell, tree by tree, in the order of ``tree_ids`` and ``refined``.
   The writer checks that ``refined`` describes complete trees.

Temporal Mesh Data
~~~~~~~~~~~~~~~~~~
The geometry and topology must both be written before any cell data, in
either order. In ``VTKHDF_STATIC`` and ``VTKHDF_FIXED_MESH`` modes each is
written once. ``VTKHDF_MOVING_MESH`` additionally allows
``write_mesh_geometry``, and ``VTKHDF_DYNAMIC_MESH`` both
``write_mesh_geometry`` and ``write_mesh_topology``, to be called again in
later time steps, subject to the same rules as for
:doc:`UnstructuredGrid <ug>`:

* Any mesh component written before the first call to ``start_time_step``
  is regarded as being written for the first time step.

* It is an error to write a mesh component more than once for a single
  time step.

* If a mesh component is not written for a given time step, its most
  recently written value is used.

* Once cell data is written for a time step, no mesh component may be
  written until the next step.

Trees may be supplied by different ranks at different time steps.

Cell Data
---------
The ``write_cell_data`` method supports two interfaces: writes by name create
static cell datasets, and writes by handle supply temporal data. In both
cases, write the mesh before writing the data. See :doc:`common` for naming
and handle conventions.

.. code-block:: fortran

   use vtkhdf_htg_file_type
   type(vtkhdf_htg_file) :: file
   type(vtkhdf_cell_data_handle) :: cell_var

.. method-doc::

   call file%write_cell_data(name, array [, attribute])

   Write the data ``array`` to a new static cell dataset ``name``. This is
   supported in all modes except ``VTKHDF_DYNAMIC_MESH``. Scalar data uses a
   rank-1 ``array``; multicomponent data uses a rank-2 ``array`` whose first
   dimension indexes components. The last dimension indexes the rank's cells
   in the order described above. Supported types are ``real32``, ``real64``,
   ``int8``, ``int32`` and ``int64``.

.. method-doc::

   cell_var = file%register_temporal_cell_data(name, mold [, attribute])

   Register ``name`` as a time-dependent cell dataset, before the first call
   to ``start_time_step``. Pass a scalar ``mold`` for scalar data, or a
   rank-1 ``mold`` whose size equals the number of components.

.. method-doc::

   call file%start_time_step(time)
   call file%write_cell_data(cell_var, array)
   call file%finalize_time_step()

   These behave as for :doc:`UnstructuredGrid <ug>`. A temporal dataset that
   is not written during a step reuses its most recently written values,
   which requires that the rank's number of cells is unchanged; after a
   change of topology the dataset must normally be rewritten.

.. warning::

   ``int8`` data is written using VTK/ParaView-compatible 8-bit
   unsigned storage. Use only nonnegative values.

.. note::

   The VTK reader labels HyperTreeGrid cell arrays by their HDF5 dataset
   names. A ``name`` containing ``/``, ``.``, or a space therefore appears
   with those characters replaced by underscores; see :doc:`common`.

Field Data
----------
See :doc:`field_data` for the shared static and temporal field data interface,
supported values, array interpretation, and MPI requirements.

Parallel Storage
----------------
Each MPI rank's trees are stored as a separate partition, and read as one
HyperTreeGrid of a partitioned dataset, each covering the whole root grid.
A rank that supplies no trees contributes an empty partition, so the number
of partitions always equals the number of ranks. No tree may be divided
between ranks.

Examples
--------

:repo:`example/htg_dynamic_mesh.F90` refines a two-dimensional grid of
quadtrees around a circle that moves across the grid, writing a new topology
and temporal cell and field data at every step. It shows how to generate the
breadth-first cell order, and supports serial and MPI builds.
