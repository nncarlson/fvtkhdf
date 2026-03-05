The vtkhdf_mb_file_type module
==============================

This module defines the ``vtkhdf_mb_file`` derived type for writing
VTKHDF MultiBlockDataSet files.

A MultiBlockDataSet file contains a flat collection of named
UnstructuredGrid blocks; hierarchical nesting of
MultiBlockDataSet objects is not supported.

Each block behaves semantically like a ``vtkhdf_ug_file`` mesh: it has
its own static mesh and associated point, cell, and field datasets.

As in ``vtkhdf_ug_file``, partitions correspond 1:1 with MPI ranks.
Each MPI rank contributes one VTKHDF partition for every block.

In the MPI build, all type-bound procedures are collective over the
communicator passed to ``create`` and must be called in the same order
on all ranks.

Most methods of ``vtkhdf_mb_file`` correspond directly to those of
``vtkhdf_ug_file``, but operate on an opaque ``vtkhdf_block_handle``
returned by ``add_block``.

.. code-block:: fortran

   use vtkhdf_mb_file_type
   type(vtkhdf_mb_file) :: file
   type(vtkhdf_block_handle) :: block
   type(vtkhdf_cell_data_handle) :: cell_var
   type(vtkhdf_point_data_handle) :: point_var
   type(vtkhdf_field_data_handle) :: field_var

File Creation and Management
----------------------------

``call file%create(filename, [comm,] stat, errmsg)``
   Create a new VTKHDF MultiBlockDataSet file.

   Arguments are identical to those of ``vtkhdf_ug_file%create``,
   except that there is no ``is_temporal`` argument.

``call file%flush()``
   Collectively flush the file's HDF5 buffers and request the operating
   system to flush file buffers.

``call file%close()``
   Close the file and release internal resources. Users should always
   call this to "finalize" the object; automatic finalization
   cannot perform a proper collective cleanup.

Block definition
----------------

``block = file%add_block(name [, is_temporal])``
   Define a new UnstructuredGrid block and return its handle.

   ``name`` is normalized and disambiguated before it is used as the block
   name seen by the VTKHDF reader and ParaView. Input is trimmed, empty input
   is replaced by a default name, the characters ``/``, ``.``, and space are
   replaced by ``_``, and duplicate internal names are made unique by
   appending a suffix.

   The original user-facing input name is still written to the block group's
   ``Name`` attribute for informational purposes.

   If ``is_temporal`` is present and ``.true.``, the block supports
   time-dependent datasets. Temporal blocks must be defined before the first
   call to ``start_time_step``. Non-temporal blocks may be defined at any time.

   The returned ``vtkhdf_block_handle`` is opaque. Its components are private,
   so user code cannot inspect or construct handles directly; only store the
   returned value and pass it to later block-scoped operations.

The file is considered temporal if at least one block is temporal.
All temporal blocks share a common timeline defined by calls to
``start_time_step``.

Mesh Data
---------

``call file%write_mesh(block, points, cnode, xcnode, types)``
   Write the mesh geometry and topology for the block identified by ``block``.

   Mesh arguments and semantics are identical to those of
   ``vtkhdf_ug_file``. The mesh for each block must be written
   exactly once and before writing any data for that block.

Static mesh-centered data
-------------------------

.. glossary::

   ``call file%write_cell_data(block, name, array)``
   ``call file%write_point_data(block, name, array)``
      Write static cell or point datasets for the block identified by
      ``block``.

      Dataset semantics are identical to those of ``vtkhdf_ug_file``.
      ``array`` must conform to the same type and shape requirements
      described for ``vtkhdf_ug_file``.

      Within a block, cell datasets share a ``CellData`` namespace and point
      datasets share a ``PointData`` namespace. Input names are trimmed, empty
      input is replaced by a default name, the characters ``/``, ``.``, and
      space are replaced by ``_``, and duplicate internal names are made
      unique by appending a suffix. The sanitized internal dataset name is
      what the VTKHDF reader and ParaView use. The original user-facing name
      is written to the dataset ``Name`` attribute for informational purposes.

Static field data
-----------------

``call file%write_field_data(block, name, array [, as_vector])``
   Write static field data for ``block``.

   Semantics are identical to ``vtkhdf_ug_file%write_field_data``:
   scalar, rank-1, and rank-2 arrays are supported. For rank-1 arrays,
   ``as_vector=.true.`` marks the data as one tuple with ``n`` components;
   otherwise rank-1 is interpreted as ``n`` tuples of one component.
   In MPI builds, calls are collective and only rank 0 contributes the payload;
   all ranks must still pass matching array type/rank/shape.

Time-dependent mesh-centered data
---------------------------------

MultiBlockDataSet files support time-dependent point and cell datasets on
a per-block basis. A block is temporal if it was defined with
``is_temporal = .true.`` in ``add_block``.

Temporal blocks must be defined before the first call to ``start_time_step``.
After the first time step is started, no further temporal block definitions
or temporal dataset registrations are allowed.

.. glossary::

   ``cell_var = file%register_temporal_cell_data(block, name, mold)``
   ``point_var = file%register_temporal_point_data(block, name, mold)``
      Register ``name`` as a time-dependent dataset on the block identified by
      ``block``. Registration semantics are identical to those of
      ``vtkhdf_ug_file``. The returned handle is opaque; user code should
      store it and pass it to later temporal writes for the same block.

``call file%start_time_step(time)``
   Start a new time step with time value ``time``. The timeline is shared by
   all temporal blocks in the file. A call to ``start_time_step`` is required
   before writing any temporal dataset in any block.

   For the first time step, each temporal dataset registered on each temporal
   block must be written once for that block.

.. glossary::

   ``call file%write_temporal_cell_data(block, cell_var, array)``
   ``call file%write_temporal_point_data(block, point_var, array)``
      Write ``array`` to the temporal dataset identified by ``cell_var`` or
      ``point_var`` for the specified ``block``, associating it with the
      current time step.
      Write semantics are identical to those of ``vtkhdf_ug_file``.
      After the first time step, writes may be skipped on later time steps; if
      omitted, the most recently written value is used.

      The data handle and block handle must match the same block registration.

``call file%finalize_time_step()``
   Finalize the current shared time step for all temporal blocks. Until this
   call, the file is in an in-progress state for that step. Best practice is
   to call ``finalize_time_step`` immediately after all temporal writes are
   complete.

   ``finalize_time_step`` is called implicitly when ``start_time_step`` begins
   a new step while one is still open, and when ``close`` is called.

Time-dependent field data
-------------------------

.. glossary::

   ``field_var = file%register_temporal_field_data(block, name, mold)``
      Register time-dependent field data for ``block``. Semantics match
      ``vtkhdf_ug_file%register_temporal_field_data``.

   ``call file%write_temporal_field_data(block, field_var, array [, as_vector])``
      Write time-dependent field data for the current shared time step.
      Semantics match ``vtkhdf_ug_file%write_temporal_field_data``.
      In MPI builds, calls are collective and only rank 0 contributes the
      payload; all ranks must still pass matching array type/rank/shape.

      The data handle and block handle must match the same block registration.
