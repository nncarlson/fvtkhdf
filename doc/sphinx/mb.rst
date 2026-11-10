MultiBlockDataSet
=================

The ``vtkhdf_mb_file_type`` module defines the ``vtkhdf_mb_file``
derived type for writing VTKHDF MultiBlockDataSet files.

Optional VTK attribute roles, specified with the ``attribute`` argument, are
recognized by VTK 9.7.0 / ParaView 6.2.0 and later; older readers ignore these
designations while still reading the data. See :doc:`common`.

MultiBlockDataSet is a composite VTK dataset: a container for simple VTK
datasets. This interface supports a flat collection of named UnstructuredGrid
blocks; hierarchical nesting and other simple VTK dataset types are not
supported.

Each block has the same semantics as a VTK dataset written by
``vtkhdf_ug_file``, with its own selected mode, mesh, and associated cell,
point, and field datasets. Each block has its own mode and data; all temporal
blocks share a common timeline.

As in ``vtkhdf_ug_file``, partitions correspond 1:1 with MPI ranks.
Each MPI rank contributes one VTKHDF partition for every block.

In the MPI build, all type-bound procedures are collective over the
communicator passed to ``create`` and must be called in the same order
on all ranks.

Most methods of ``vtkhdf_mb_file`` correspond directly to those of
``vtkhdf_ug_file``, but operate on the block identified by an opaque
handle returned by ``add_block``.

File Creation and Management
----------------------------

.. code-block:: fortran

   use vtkhdf_mb_file_type
   type(vtkhdf_mb_file) :: file

.. method-doc::

   call file%create(filename, [comm,] stat, errmsg [, data_chunk_bytes])

   Create a new VTKHDF MultiBlockDataSet file.

   Arguments are identical to those of ``vtkhdf_ug_file%create``,
   except that there is no file-level ``mode`` argument.

.. method-doc::

   call file%flush()

   Flush the file's HDF5 buffers.

.. method-doc::

   call file%close()

   Finalize any open time step, close the file, and release resources. Call
   this explicitly; automatic Fortran finalization cannot perform collective
   cleanup.

Block Definition
----------------

The shared ``VTKHDF_*`` mode constants are exported by
``vtkhdf_mb_file_type``. The original ``UG_*`` names remain available as
compatibility aliases.

.. code-block:: fortran

   use vtkhdf_mb_file_type
   type(vtkhdf_mb_file) :: file
   type(vtkhdf_block_handle) :: block

.. method-doc::

   block = file%add_block(name [, mode])

   Define a new UnstructuredGrid block and return its handle.

   ``mode`` is one of ``VTKHDF_STATIC``, ``VTKHDF_FIXED_MESH``,
   ``VTKHDF_MOVING_MESH``, or ``VTKHDF_DYNAMIC_MESH`` and defaults to
   ``VTKHDF_STATIC``.

   Block names follow the naming conventions in :doc:`common`.
   The original name is preserved as the value of the block group
   ``Name`` attribute for user reference.

   The returned ``vtkhdf_block_handle`` is opaque. Store it and reuse it for
   later calls for that block.

The file is temporal if at least one block uses a temporal UG mode
(``VTKHDF_FIXED_MESH``, ``VTKHDF_MOVING_MESH``, or ``VTKHDF_DYNAMIC_MESH``). All
temporal blocks share calls to ``start_time_step`` and ``finalize_time_step``.

Mesh Data
---------

.. code-block:: fortran

   use vtkhdf_mb_file_type
   type(vtkhdf_mb_file) :: file
   type(vtkhdf_block_handle) :: block

.. method-doc::

   call file%write_mesh(block, points, cnode, xcnode, types)

   Write the mesh geometry and topology for the block identified by ``block``.

.. method-doc::

   call file%write_mesh_topology(block, cnode, xcnode, types)

   Write only mesh topology for the block identified by ``block``. This is
   mainly used with ``VTKHDF_MOVING_MESH`` or ``VTKHDF_DYNAMIC_MESH`` when geometry
   is written separately.

.. method-doc::

   call file%write_mesh_geometry(block, points)

   Write only mesh geometry for the block identified by ``block``.

Mesh arguments and semantics are identical to those of ``vtkhdf_ug_file``.
The mesh for each block must be written before writing its point or cell data.


Point, Cell, and Field Data
---------------------------

Except for the added block handle, all arguments and semantics are the
same as for ``vtkhdf_ug_file``; see its description for details omitted
here. Field data values, array interpretation, and reuse rules are described
in :doc:`field_data`.

For field data, ``file`` is a ``vtkhdf_mb_file`` object, and each registration
or write takes the block handle returned by ``file%add_block`` as its first
argument, as shown below.

Within each block, cell data, point data, and field data have their
own namespaces. Names are normalized and deduplicated exactly as in
``vtkhdf_ug_file``.

If multiple blocks contain datasets with the same name and
association (cell or point), ParaView treats them as a single variable across
the composite VTK dataset. Selecting that variable uses the dataset from each
block that contains it. Datasets with the same name and association should
therefore represent the same quantity and have compatible structure (e.g.,
same type and number of components).

.. code-block:: fortran

   use vtkhdf_mb_file_type
   type(vtkhdf_mb_file) :: file
   type(vtkhdf_block_handle) :: block
   type(vtkhdf_cell_data_handle) :: cell_var
   type(vtkhdf_point_data_handle) :: point_var
   type(vtkhdf_field_data_handle) :: field_var

Static Data
~~~~~~~~~~~
Static point and cell writes are supported in ``VTKHDF_STATIC`` and
``VTKHDF_FIXED_MESH`` block modes. Static field data is supported in every
block mode and may be written at any point after block creation and before
closing the file.

.. method-doc::

   call file%write_cell_data(block, name, array [, attribute])
   call file%write_point_data(block, name, array [, attribute])
   call file%write_field_data(block, name, array [, as_vector])

   Write ``array`` to a new cell, point, or field dataset ``name``
   for the block identified by ``block``.

   For point and cell data, the optional ``attribute`` argument designates
   a VTK role; see :doc:`common` for attribute conventions.

Temporal Data
~~~~~~~~~~~~~
All temporal blocks share a common timeline created by ``start_time_step``.

* Temporal blocks must be added before the first time step starts.
* Temporal data registration is block-scoped, and registrations for all
  temporal blocks must occur before the first time step starts.
* Temporal variable handles are tied to the block for which they were
  registered and must not be used with a different block.

.. method-doc::

   cell_var  = file%register_temporal_cell_data(block, name, mold [, attribute])
   point_var = file%register_temporal_point_data(block, name, mold [, attribute])
   field_var = file%register_temporal_field_data(block, name, mold)

   Register ``name`` as a time-dependent cell, point, or field dataset
   for the block identified by ``block``, and return an opaque handle
   to it.

   For point and cell data, the optional ``attribute`` argument designates
   a VTK role at registration; see :doc:`common` for attribute conventions.

.. method-doc::

   call file%start_time_step(time)

   Start a new shared time step with the scalar ``real32`` or ``real64``
   time value ``time`` for every temporal block in the file. In MPI builds,
   only rank zero's value is used; the call is collective.

.. method-doc::

   call file%write_cell_data(block, cell_var, array)
   call file%write_point_data(block, point_var, array)
   call file%write_field_data(block, field_var, array [, as_vector])

   Write ``array`` to the temporal dataset identified by the variable
   handle for the block identified by ``block``, associating it with
   the current time step.

.. method-doc::

   call file%finalize_time_step()

   Complete the current shared time step after all its temporal writes.
   Starting another step or closing the file also finalizes an open step.

Static blocks may coexist with temporal blocks and are unaffected by
``start_time_step``.
