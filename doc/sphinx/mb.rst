The vtkhdf_mb_file_type module
==============================

This module defines the ``vtkhdf_mb_file`` derived type for writing
VTKHDF MultiBlockDataSet files.

A MultiBlockDataSet file contains a flat collection of named
UnstructuredGrid blocks; hierarchical nesting of
MultiBlockDataSet objects is not supported.

Each block behaves semantically like a ``vtkhdf_ug_file`` mesh: it has
its own static mesh and associated point and cell datasets.

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

   ``name`` is the user-facing block name stored in the file metadata.
   Names that are non-empty, unique, and avoid ``/``, ``.``, and spaces are
   still recommended. However, ``add_block`` no longer fails for empty,
   duplicate, or otherwise awkward names. It sanitizes invalid characters,
   substitutes a default name for empty input, and appends a suffix when
   needed so the file remains valid.

   If ``is_temporal`` is present and ``.true.``, the block supports
   time-dependent datasets. Temporal blocks must be defined before the first
   call to ``write_time_step``. Non-temporal blocks may be defined at any time.

   The returned ``vtkhdf_block_handle`` is opaque. Its components are private,
   so user code cannot inspect or construct handles directly; only store the
   returned value and pass it to later block-scoped operations.

The file is considered temporal if at least one block is temporal.
All temporal blocks share a common timeline defined by calls to
``write_time_step``.

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

Time-dependent mesh-centered data
---------------------------------

MultiBlockDataSet files support time-dependent point and cell datasets on
a per-block basis. A block is temporal if it was defined with
``is_temporal = .true.`` in ``add_block``.

Temporal blocks must be defined before the first call to ``write_time_step``.
After the first time step is started, no further temporal block definitions
or temporal dataset registrations are allowed.

.. glossary::

   ``call file%register_temporal_cell_data(block, name, mold)``
   ``call file%register_temporal_point_data(block, name, mold)``
      Register ``name`` as a time-dependent dataset on the block identified by
      ``block``. Registration semantics are identical to those of
      ``vtkhdf_ug_file``.

``call file%write_time_step(time)``
   Start a new time step with time value ``time``. The timeline is shared by
   all temporal blocks in the file. A call to ``write_time_step`` is required
   before writing any temporal dataset in any block.

.. glossary::

   ``call file%write_temporal_cell_data(block, name, array)``
   ``call file%write_temporal_point_data(block, name, array)``
      Write ``array`` to the temporal dataset ``name`` on the block identified
      by ``block``, associating it with the current time step. Write semantics
      are identical to those of ``vtkhdf_ug_file``.
