UnstructuredGrid
================

The ``vtkhdf_ug_file_type`` module defines the ``vtkhdf_ug_file``
derived type for writing VTKHDF UnstructuredGrid files.

Optional VTK attribute roles, specified with the ``attribute`` argument, are
recognized by VTK 9.7.0 / ParaView 6.2.0 and later; older readers ignore these
designations while still reading the data. See :doc:`common`.

In the MPI build, all type-bound procedures are collective over the
communicator passed to ``create``. Every rank must call the same method
in the same order and supply identical values for all non-distributed
arguments.

File Creation and Management
----------------------------

.. code-block:: fortran

   use vtkhdf_ug_file_type
   type(vtkhdf_ug_file) :: file

.. method-doc::

   call file%create(filename, [comm,] stat, errmsg [, mode] [, data_chunk_bytes])

   Create a new VTKHDF UnstructuredGrid file.

   * ``filename``: path to the file to create. The recommended extension is
     ``.vtkhdf``.
   * ``comm``: MPI communicator. In MPI builds with ``USE_MPI_F08=ON``
     (default), this may be either ``integer`` or ``type(MPI_Comm)``.
     With ``USE_MPI_F08=OFF``, only the ``integer`` communicator
     interface is available. In serial builds, ``comm`` is omitted from
     the interface.
   * ``mode``: one of ``VTKHDF_STATIC``, ``VTKHDF_FIXED_MESH``,
     ``VTKHDF_MOVING_MESH``, or ``VTKHDF_DYNAMIC_MESH``. The default is
     ``VTKHDF_STATIC``.
   * ``data_chunk_bytes``: target chunk size in bytes for chunked datasets
     storing mesh-sized data. The default is ``128*1024`` (128 KiB).
   * ``stat``: status code. ``stat == 0`` indicates success.
   * ``errmsg``: allocated error message when ``stat /= 0``.

   In the MPI build, ``stat`` and ``errmsg`` are collective outputs;
   all ranks return identical values.

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

* ``VTKHDF_STATIC`` (default): fixed mesh and static point, cell, and field
  data; no time stepping.
* ``VTKHDF_FIXED_MESH``: fixed mesh and static or temporal point, cell, and
  field data.
* ``VTKHDF_MOVING_MESH``: fixed mesh topology with temporal geometry;
  point/cell datasets must be temporal.
* ``VTKHDF_DYNAMIC_MESH``: temporal mesh topology and geometry;
  point/cell datasets must be temporal.

Field data may be static or temporal in all temporal modes.

The mode constants are exported by ``vtkhdf_ug_file_type``. The original
names ``UG_STATIC``, ``UG_FIXED_MESH``, ``UG_MOVING_MESH``, and
``UG_DYNAMIC_MESH`` remain available as compatibility aliases with the same
values.

Temporal datasets are registered before time stepping; writes through their
handles supply values for the current step. Static data remains independent
of time stepping.

Mesh Data
---------

These procedures write the portion of an unstructured mesh provided by the
calling MPI rank.

.. code-block:: fortran

   use vtkhdf_ug_file_type
   type(vtkhdf_ug_file) :: file

.. method-doc::

   call file%write_mesh(points, cnode, xcnode, types)

   Write the full mesh, geometry and topology.

.. method-doc::

   call file%write_mesh_topology(cnode, xcnode, types)

   Write only topology. This is mainly used with ``VTKHDF_MOVING_MESH`` or
   ``VTKHDF_DYNAMIC_MESH`` when geometry is written separately.

.. method-doc::

   call file%write_mesh_geometry(points)

   Write only geometry.

The array arguments are the same in all cases:

* ``points``: ``real32`` or ``real64`` array of shape ``(3,npoints)``
  containing point coordinates. Supply all three coordinate components,
  including for 1D or 2D meshes; unused components may be set to zero.

* ``cnode``, ``xcnode``: ``int32`` or ``int64`` rank-1 arrays describing
  the mesh topology in compressed sparse row (CSR) format:

  - ``cnode`` contains the concatenated cell-node connectivity lists.
  - ``xcnode`` contains the starting index of each cell's connectivity in
    ``cnode``. For cell ``i``, ``cnode(xcnode(i):xcnode(i+1)-1)`` is its
    node list. ``size(xcnode)`` must equal ``ncells + 1``, and its final
    element equals ``size(cnode)`` + 1.

  CSR indexing is 1-based, following Fortran conventions. Conversion to the
  0-based VTKHDF representation is handled internally.

* ``types``: an ``int8`` array of length ``ncells`` containing VTK cell type
  codes. Named constants such as ``VTK_TETRA`` and ``VTK_HEXAHEDRON`` are
  provided by ``vtkhdf_vtk_cell_types``.

Temporal Mesh Data
~~~~~~~~~~~~~~~~~~
The mesh must be written before any associated point or cell data is written.
In ``VTKHDF_STATIC`` and ``VTKHDF_FIXED_MESH`` modes, this means one initial call to
``write_mesh``. In ``VTKHDF_MOVING_MESH`` and ``VTKHDF_DYNAMIC_MESH`` modes the mesh
may be written, or not, every time step. This is subject to the following rules:

* The mesh must be written once before writing temporal point/cell
  datasets for the first step.

* Any mesh component (topology or geometry) written before the first call
  to ``start_time_step`` is regarded as being written for the first time
  step.

* It is an error to write any mesh component more than once for a single
  time step.

* If a mesh component is not written for a given time step, its most
  recently written value is used.

* Once point or cell data is written for a time step, no mesh component
  may be written until the next step.


Point and Cell Data
-------------------
The ``write_cell_data`` and ``write_point_data`` methods support two
interfaces: writes by name create static point or cell datasets, and writes
by handle supply temporal data. In both cases, write the associated mesh
before writing the data. See :doc:`common` for naming and handle conventions.

.. code-block:: fortran

   use vtkhdf_ug_file_type
   type(vtkhdf_ug_file) :: file
   type(vtkhdf_cell_data_handle) :: cell_var
   type(vtkhdf_point_data_handle) :: point_var

Static Data
~~~~~~~~~~~
These writes by name are supported in ``VTKHDF_STATIC`` and
``VTKHDF_FIXED_MESH`` modes.

.. method-doc::

   call file%write_cell_data(name, array [, attribute])
   call file%write_point_data(name, array [, attribute])

   Write the data ``array`` to a new cell or point dataset ``name``.
   Scalar data uses a rank-1 ``array``; multicomponent data uses a rank-2
   ``array`` whose first dimension indexes components. The last dimension
   indexes mesh entities and must have extent ``ncells`` (cell data) or
   ``npoints`` (point data).
   Supported types are ``real32``, ``real64``, ``int8``, ``int32`` and
   ``int64``.

   The optional ``attribute`` argument designates a VTK role, such as
   ``Tensors`` or ``GlobalIds``; see :doc:`common` for attribute conventions.

.. warning::

   ``int8`` data is written using VTK/ParaView-compatible 8-bit
   unsigned storage. Use only nonnegative values. Negative ``int8``
   values are not a supported portable interface.

.. note::

   A tensor value represented by a rank-2 Fortran array must be packed into
   a rank-1 component tuple. Supply tensor data for the mesh as a rank-2
   array of shape ``(ncomp,nentities)``, including for a single entity;
   separate array dimensions for tensor rows and columns are not supported.
   VTK and ParaView receive one multicomponent tuple per mesh entity.

Temporal Data
~~~~~~~~~~~~~
Register temporal point and cell datasets before the first call to
``start_time_step``. Write their values during an open time step.

.. method-doc::

   cell_var = file%register_temporal_cell_data(name, mold [, attribute])
   point_var = file%register_temporal_point_data(name, mold [, attribute])

   Register ``name`` as a time-dependent cell or point dataset. The type
   and kind of ``mold`` determine the dataset type. For scalar data,
   pass a scalar ``mold``; for multicomponent data, pass a rank-1 ``mold``
   whose size equals the number of components. The value of ``mold`` is
   never referenced. Supported types are ``real32``, ``real64``, ``int8``,
   ``int32`` and ``int64``.

   The optional ``attribute`` argument designates a VTK role at registration;
   see :doc:`common` for attribute conventions.

   The returned handle is opaque; user code should store it and pass it
   to later temporal writes.

   .. warning::

      ``int8`` data is written using VTK/ParaView-compatible 8-bit
      unsigned storage. Use only nonnegative values. Negative ``int8``
      values are not a supported portable interface.

.. method-doc::

   call file%start_time_step(time)

   Start a new time step with the scalar ``real32`` or ``real64`` time value
   ``time``. In MPI builds, only rank zero's value is used; the call is collective.

   After ``start_time_step`` is called for the first time step, every
   temporal point or cell dataset must be written once for that step.

.. method-doc::

   call file%write_cell_data(cell_var, array)
   call file%write_point_data(point_var, array)

   Write ``array`` to the temporal dataset identified by the handle,
   associating it with the current time step. Its type, kind, scalar or
   multicomponent interpretation, and component count must match the
   registered ``mold``. Its last dimension must have extent ``ncells``
   (cell data) or ``npoints`` (point data) for the current mesh.

   After the first time step, a temporal dataset need not be written
   at every time step. If omitted, its most recently written value is
   used, provided it is compatible with the current mesh; otherwise
   the dataset must be rewritten. Writing the same temporal dataset
   more than once within a single time step is an error.

.. method-doc::

   call file%finalize_time_step()

   Complete the current time step after all its temporal writes. Starting
   another step or closing the file also finalizes an open step.

Field Data
----------
See :doc:`field_data` for the shared static and temporal field data interface,
supported values, array interpretation, and MPI requirements.

Workflow by Mode
----------------

In every mode, static field data may be written by name at any point after
file creation and before closing the file, independently of mesh writes.

.. rubric:: ``VTKHDF_STATIC``

* Write the mesh once with ``write_mesh``.
* Then write static point and cell datasets by name:

  - ``call file%write_point_data(name, array [, attribute])``
  - ``call file%write_cell_data(name, array [, attribute])``

Time stepping and temporal registration are not allowed.

.. rubric:: ``VTKHDF_FIXED_MESH``

* Write the mesh once with ``write_mesh``.
* Write static cell and point data (may be done at any time after the mesh)

  - ``call file%write_cell_data(name, array [, attribute])``
  - ``call file%write_point_data(name, array [, attribute])``

* Register temporal datasets before the first time step:

  - ``cell_var = file%register_temporal_cell_data(name, mold [, attribute])``
  - ``point_var = file%register_temporal_point_data(name, mold [, attribute])``
  - ``field_var = file%register_temporal_field_data(name, mold)``

* Then, for each step:

  * ``call file%start_time_step(time)``
  * Write temporal datasets (or not, if reuse is possible):

    - ``call file%write_cell_data(cell_var, array)``
    - ``call file%write_point_data(point_var, array)``
    - ``call file%write_field_data(field_var, array [, as_vector])``

  * ``call file%finalize_time_step``

.. rubric:: ``VTKHDF_MOVING_MESH``

* Write topology once with ``write_mesh_topology``.

* Register temporal datasets before the first time step:

  - ``cell_var = file%register_temporal_cell_data(name, mold [, attribute])``
  - ``point_var = file%register_temporal_point_data(name, mold [, attribute])``
  - ``field_var = file%register_temporal_field_data(name, mold)``

* Then, for each step:

  * ``call file%start_time_step(time)``
  * ``call file%write_mesh_geometry(points)`` first step, or when geometry changes
  * Write temporal datasets (or not, if reuse is possible):

    - ``call file%write_cell_data(cell_var, array)``
    - ``call file%write_point_data(point_var, array)``
    - ``call file%write_field_data(field_var, array [, as_vector])``

  * ``call file%finalize_time_step``

The number of points must be constant across time steps.

Static point and cell datasets are not allowed in this mode.

.. rubric:: ``VTKHDF_DYNAMIC_MESH``

* Write initial mesh with ``write_mesh``, or wait until the first step.

* Register temporal datasets before the first time step:

  - ``cell_var = file%register_temporal_cell_data(name, mold [, attribute])``
  - ``point_var = file%register_temporal_point_data(name, mold [, attribute])``
  - ``field_var = file%register_temporal_field_data(name, mold)``

* Then, for each step:

  * ``call file%start_time_step(time)``

  * Write mesh for first step if not already written, or when mesh changes:

    - ``call file%write_mesh_topology(cnode, xcnode, types)`` to replace topology only;
    - ``call file%write_mesh_geometry(points)`` to replace geometry only;
    - ``call file%write_mesh(points, cnode, xcnode, types)`` to replace both.

  * Write temporal datasets (or not, if reuse is possible):

    - ``call file%write_cell_data(cell_var, array)``
    - ``call file%write_point_data(point_var, array)``
    - ``call file%write_field_data(field_var, array [, as_vector])``

  * ``call file%finalize_time_step``

Static point and cell datasets are not allowed in this mode.

This is the most general mode. Point counts, cell counts, connectivity, and
geometry may all vary from one step to the next.
