The vtkhdf_mb_file_type module
==============================

This module defines the ``vtkhdf_mb_file`` derived type for writing
VTKHDF MultiBlockDataSet files.

A MultiBlockDataSet file contains a flat collection of named
UnstructuredGrid blocks; hierarchical nesting is not supported.

Each block behaves semantically like a ``vtkhdf_ug_file`` dataset; each with
its own selected mode, mesh and associated cell, point and field datasets.
However, all temporal blocks share a common timeline.

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

``call file%create(filename, [comm,] stat, errmsg [, data_chunk_bytes])``
   Create a new VTKHDF MultiBlockDataSet file.

   Arguments are identical to those of ``vtkhdf_ug_file%create``,
   except that there is no file-level ``mode`` argument.

``call file%flush()``
   Collectively flush the file's HDF5 buffers and request the operating
   system to flush file buffers.

``call file%close()``
   Close the file and release internal resources. Users should always
   call this to "finalize" the object; automatic finalization
   cannot perform a proper collective cleanup.

Block Definition
----------------
.. code-block:: fortran

   use vtkhdf_mb_file_type
   type(vtkhdf_mb_file) :: file
   type(vtkhdf_block_handle) :: block

``block = file%add_block(name [, mode])``
   Define a new UnstructuredGrid block and return its handle.

   ``mode`` is one of ``UG_STATIC``, ``UG_FIXED_MESH``,
   ``UG_MOVING_MESH``, or ``UG_DYNAMIC_MESH`` and defaults to
   ``UG_STATIC``.

   The block name is normalized in the same way as ``vtkhdf_ug_file`` dataset
   names. The original name is preserved as the value of the block group
   ``Name`` attribute for user reference.

   The returned ``vtkhdf_block_handle`` is opaque. Store it and reuse it for
   later block-scoped calls.

The file is temporal if at least one block uses a temporal UG mode
(``UG_FIXED_MESH``, ``UG_MOVING_MESH``, or ``UG_DYNAMIC_MESH``). All
temporal blocks share calls to ``start_time_step`` and ``finalize_time_step``.

Mesh Data
---------
.. code-block:: fortran

   use vtkhdf_mb_file_type
   type(vtkhdf_mb_file) :: file
   type(vtkhdf_block_handle) :: block

``call file%write_mesh(block, points, cnode, xcnode, types)``
   Write the mesh geometry and topology for the block identified by ``block``.

``call file%write_mesh_topology(block, cnode, xcnode, types)``
   Write only mesh topology for the block identified by ``block``. This is
   mainly used with ``UG_MOVING_MESH`` or ``UG_DYNAMIC_MESH`` when geometry
   is written separately.

``call file%write_mesh_geometry(block, points)``
   Write only mesh geometry for the block identified by ``block``.

Mesh arguments and semantics are identical to those of ``vtkhdf_ug_file``.
The mesh for each block must be written before writing any associated
mesh-centered data for that block.


Cell, Point, and Field Data
---------------------------
Except for the added block handle, all arguments and semantics are the
same as for ``vtkhdf_ug_file``; see its description for details omitted
here.

Within each block, cell data, point data, and field data have their
own namespaces. Names are normalized and deduplicated exactly as in
``vtkhdf_ug_file``.

If multiple blocks contain mesh-centered datasets with the same name and
association (cell or point), ParaView treats them as a single variable across
the composite dataset. Selecting that variable uses the dataset from each
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
These by-name writes are allowed only where the selected ``mode`` for the
block permits static data of that category (mesh-centered or field).

.. glossary::

   ``call file%write_cell_data(block, name, array [, attribute])``
   ``call file%write_point_data(block, name, array [, attribute])``
   ``call file%write_field_data(block, name, array [, as_vector])``
      Write ``array`` to a new cell, point, or field dataset ``name``
      for the block identified by ``block``.

      For cell and point data, the optional ``attribute`` argument writes
      the VTKHDF dataset ``Attribute`` HDF5 attribute for arrays with a specific
      VTK role, for example ``Tensors`` or ``GlobalIds``. For the full
      set of allowed values, see the VTKHDF specification section
      `Attribute Data <https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/vtkhdf_specifications.html#attribute-data>`_.

Temporal Data
~~~~~~~~~~~~~
All temporal blocks share a common timeline created by ``start_time_step``.

* Temporal blocks must be added before the first time step starts.
* Temporal data registration is block-scoped, and registrations for all
  temporal blocks must occur before the first time step starts.
* Temporal variable handles are tied to the block for which they were
  registered and must not be used with a different block.

.. glossary::

   ``cell_var  = file%register_temporal_cell_data(block, name, mold [, attribute])``
   ``point_var = file%register_temporal_point_data(block, name, mold [, attribute])``
   ``field_var = file%register_temporal_field_data(block, name, mold)``
      Register ``name`` as a time-dependent cell, point, or field dataset
      for the block identified by ``block``, and return an opaque handle
      to it.

      For cell and point data, the optional ``attribute`` argument writes
      the dataset ``Attribute`` HDF5 attribute at registration time, for example
      ``Tensors`` or ``GlobalIds``. For the full set of allowed values,
      see the VTKHDF specification section `Attribute Data
      <https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/vtkhdf_specifications.html#attribute-data>`_.

``call file%start_time_step(time)``
   Start a new shared time step with time value ``time`` for every
   temporal block in the file.

.. glossary::

   ``call file%write_cell_data(block, cell_var, array)``
   ``call file%write_point_data(block, point_var, array)``
   ``call file%write_field_data(block, field_var, array [, as_vector])``
      Write ``array`` to the temporal dataset identified by the variable
      handle for the block identified by ``block``, associating it with
      the current time step.

``call file%finalize_time_step()``
   Finalize the current shared time step. Until this call, the file
   is in an in-progress state for that step. Best practice is to call
   ``finalize_time_step`` immediately after all temporal writes for
   the step are complete.

Static blocks may coexist with temporal blocks and are unaffected by
``start_time_step``.
