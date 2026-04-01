The vtkhdf_ug_file_type module
==============================

This module defines the ``vtkhdf_ug_file`` derived type for writing
VTKHDF UnstructuredGrid files.

In the MPI build, all type-bound procedures are collective over the
communicator passed to ``create``. Every rank must call the same method
in the same order and supply identical values for all non-distributed
arguments.

File Creation and Management
----------------------------

.. code-block:: fortran

   use vtkhdf_ug_file_type
   type(vtkhdf_ug_file) :: file

``call file%create(filename, [comm,] stat, errmsg [, mode] [, data_chunk_bytes])``
   Create a new VTKHDF UnstructuredGrid file.

   * ``filename``: path to the file to create. The recommended extension is
     ``.vtkhdf``.
   * ``comm``: MPI communicator. In MPI builds with ``USE_MPI_F08=ON``
     (default), this may be either ``integer`` or ``type(MPI_Comm)``.
     With ``USE_MPI_F08=OFF``, only the ``integer`` communicator
     interface is available. In serial builds, ``comm`` is omitted from
     the interface.
   * ``mode``: one of ``UG_STATIC``, ``UG_FIXED_MESH``,
     ``UG_MOVING_MESH``, or ``UG_DYNAMIC_MESH``. The default is
     ``UG_STATIC``.
   * ``data_chunk_bytes``: target chunk size in bytes for chunked datasets
     storing mesh-sized data. The default is ``128*1024`` (128 KiB).
   * ``stat``: status code. ``stat == 0`` indicates success.
   * ``errmsg``: allocated error message when ``stat /= 0``.

   In the MPI build, ``stat`` and ``errmsg`` are collective outputs;
   all ranks return identical values.

``call file%close()``
   Close the file and release internal resources. Users should always call
   this to "finalize" the object; automatic finalization cannot perform a
   proper collective cleanup.

``call file%flush()``
   Flush the file's HDF5 buffers and request the OS to flush the file buffers.

Writer Modes
~~~~~~~~~~~~
The ``mode`` argument selects one of four workflows:

* ``UG_STATIC``: single mesh with associated cell/point datasets; field
  datasets; no time stepping (static).
* ``UG_FIXED_MESH``: single fixed mesh; cell/point datasets may be static
  (independent of time step) or temporal (changing with time step); field
  datasets must be temporal.
* ``UG_MOVING_MESH``: fixed mesh topology with temporal geometry;
  all datasets must be temporal.
* ``UG_DYNAMIC_MESH``: temporal mesh topology and geometry;
  all datasets must be temporal.

In temporal modes, registration determines which datasets participate in time
stepping, and each write supplies the value for the current step.

Mesh Data
---------

These procedures write the portion of an unstructured mesh provided by the
calling MPI rank.

.. code-block:: fortran

   use vtkhdf_ug_file_type
   type(vtkhdf_ug_file) :: file

``call file%write_mesh(points, cnode, xcnode, types)``
   Write the full mesh, geometry and topology.

``call file%write_mesh_topology(cnode, xcnode, types)``
   Write only topology. This is mainly used with ``UG_MOVING_MESH`` or
   ``UG_DYNAMIC_MESH`` when geometry is written separately.

``call file%write_mesh_geometry(points)``
   Write only geometry.

The array arguments are the same in all cases:

* ``points``: ``real32`` or ``real64`` array of shape (3, `npoints`)
  containing node coordinates. Coordinates are always interpreted as 3-D;
  for 1-D or 2-D problems, the unused components must be set (e.g., to 0).

* ``cnode``, ``xcnode``: ``int32`` or ``int64`` rank-1 arrays describing
  the mesh topology in CSR format:

  - ``cnode`` contains the concatenated cell-node connectivity lists.
  - ``xcnode`` contains the starting index of each cell's connectivity in
    ``cnode``. For cell ``i``, ``cnode(xcnode(i):xcnode(i+1)-1)`` is its
    node list. ``size(xcnode)`` must equal `ncells` + 1, and its final
    element equals ``size(cnode)`` + 1.

  CSR indexing is 1-based, following Fortran conventions. Conversion to the
  0-based VTKHDF representation is handled internally.

* ``types``: an ``int8`` array of length `ncells` containing VTK cell type
  codes. Named constants such as ``VTK_TETRA`` and ``VTK_HEXAHEDRON`` are
  provided by ``vtkhdf_vtk_cell_types``.

Temporal Mesh Data
~~~~~~~~~~~~~~~~~~
The mesh must be written before any associated mesh-centered data is written.
In ``UG_STATIC`` and ``UG_FIXED_MESH`` modes, this means one initial call to
``write_mesh``. In ``UG_MOVING_MESH`` and ``UG_DYNAMIC_MESH`` modes the mesh
may be written, or not, every time step. This is subject to the following rules:

* The mesh must be written once prior to writing temporal mesh-centered
  datasets for the first step.

* Any mesh component (topology or geometry) written before the first call
  to ``start_time_step`` is regarded as being written for the first time
  step.

* It is an error to write any mesh component more than once for a single
  time step.

* If a mesh component is not written for a given time step, its most
  recently written value is used.

* Once mesh-centered data is written for a time step, no mesh component
  may be written until the next step.


Cell and Point Data
-------------------
The ``write_cell_data`` and ``write_point_data`` methods support two
interfaces: `by-name` for writing static mesh-centered data, and `by-handle`
for writing temporal mesh-centered data. In both cases the mesh associated
with the data must be written prior to writing the data.

Cell datasets and point datasets have separate namespaces. Input names are
trimmed, empty names are replaced by a default, the characters ``/``, ``.``,
and space are replaced by an underscore, and duplicate internal names are made
unique by appending a suffix such as ``_1``. This normalized internal name is
what VTK and ParaView see. For user reference, the original name is preserved
as the value of the dataset ``Name`` attribute.

.. code-block:: fortran

   use vtkhdf_ug_file_type
   type(vtkhdf_ug_file) :: file
   type(vtkhdf_cell_data_handle) :: cell_var
   type(vtkhdf_point_data_handle) :: point_var

Static Data
~~~~~~~~~~~
These by-name writes are allowed only where the selected ``mode`` permits
static mesh-centered data.

.. glossary::

   ``call file%write_cell_data(name, array [, attribute])``
   ``call file%write_point_data(name, array [, attribute])``
      Write the data ``array`` to a new cell or point dataset ``name``.
      Scalar data (rank-1 ``array``) and vector data (rank-2 ``array``)
      are supported. The last dimension of ``array`` indexes the mesh entity
      and must have extent `ncells` (cell data) or `npoints` (point data).
      Supported types are ``real32``, ``real64``, ``int8``, ``int32`` and
      ``int64``.

      If present, ``attribute`` writes the VTKHDF dataset ``Attribute``
      HDF5 attribute for arrays with a specific VTK role, for example
      ``Tensors`` or ``GlobalIds``. Blank values are ignored. For the
      full set of allowed values, see the VTKHDF specification section
      `Attribute Data <https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/vtkhdf_specifications.html#attribute-data>`_.

.. warning::
   ``int8`` data is written using VTK/ParaView-compatible 8-bit
   unsigned storage. Use only nonnegative values. Negative ``int8``
   values are not a supported portable interface.

.. note::
   VTK and ParaView only support scalar and vector-valued mesh-centered data.
   Other kinds such as tensor values must be packed into a vector.

Temporal Data
~~~~~~~~~~~~~
Temporal mode files support time-dependent datasets. Temporal datasets must
be registered before the first call to ``start_time_step``. After the first
time step is started, no further registrations are allowed. The first time
step must be started before temporal datasets are written.

.. glossary::

   ``cell_var = file%register_temporal_cell_data(name, mold [, attribute])``
   ``point_var = file%register_temporal_point_data(name, mold [, attribute])``
      Register ``name`` as a time-dependent cell or point dataset. The type
      and kind of ``mold`` determines the dataset type. For scalar data,
      pass a scalar ``mold``, and for vector data, pass a rank-1 ``mold``
      whose size equals the number of components. The value of ``mold`` is
      never referenced. Supported types are ``real32``, ``real64``, ``int8``,
      ``int32`` and ``int64``.

      If present, ``attribute`` writes the VTKHDF dataset ``Attribute``
      HDF5 attribute for the registered array, for example ``Tensors`` or
      ``GlobalIds``. For the full set of allowed values, see the VTKHDF
      specification section `Attribute Data <https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/vtkhdf_specifications.html#attribute-data>`_.

      The returned handle is opaque; user code should store it and pass it
      to later temporal writes.

      .. warning::
         ``int8`` data is written using VTK/ParaView-compatible 8-bit
         unsigned storage. Use only nonnegative values. Negative ``int8``
         values are not a supported portable interface.

   ``call file%start_time_step(time)``
      Start a new time step with time value ``time``.

      After ``start_time_step`` is called for the first time step, every
      temporal mesh-centered dataset must be written once for that step.

   ``call file%write_cell_data(cell_var, array)``
   ``call file%write_point_data(point_var, array)``
      Write ``array`` to the temporal dataset identified by the handle,
      associating it with the current time step. ``array`` must conform
      to the shape and type implied by the registered ``mold`` and must
      have extent `ncells` (cell data) or `npoints` (point data) in its
      last dimension.

      After the first time step, a temporal dataset need not be written
      at every time step. If omitted, its most recently written value is
      used, provided it is compatible with the current mesh; otherwise
      the dataset must be rewritten. Writing the same temporal dataset
      more than once within a single time step is an error.

   ``call file%finalize_time_step()``
      Finalize the current time step. Until this call, the file is in an
      in-progress state for that step. Best practice is to call
      ``finalize_time_step`` immediately after all temporal writes for
      the step are complete.

      As a safeguard, ``finalize_time_step`` is called automatically if
      ``start_time_step`` or ``close`` is called while the previous time
      step is still open, but users are discouraged from relying on this.

Field data
----------
Writing field data is similar to cell and point data, and follows the
same static/temporal pattern, but with some notable differences:

* Field data is not tied to mesh entities, and its writing is not ordered
  with respect to mesh writing. *Scalars*, rank-1 arrays, and rank-2
  arrays are supported for ``int32``, ``int64``, ``real32``, and
  ``real64`` types.

* In MPI builds, all field-data calls remain collective. However, only
  rank 0 (the root) contributes the written payload, but all ranks must
  call the same method in the same order and pass ``array`` arguments
  with matching type, kind, and shape. On non-root ranks ``array`` values
  are never referenced.

* For rank-1 ``array``, ``as_vector=.true.`` causes the data to interpreted
  as a rank-2 array with shape (`n`,1). This conforms to how Paraview
  naturally understands geometric vectors.

* The shape (rank and extents) of temporal field data may vary from step
  to step.

.. code-block:: fortran

   use vtkhdf_ug_file_type
   type(vtkhdf_ug_file) :: file
   type(vtkhdf_field_data_handle) :: field_var

``call file%write_field_data(name, array [, as_vector])``
   Write ``array`` to a new static field dataset ``name``.
   This is allowed only when the selected ``mode`` is ``UG_STATIC``.

``field_var = file%register_temporal_field_data(name, mold)``
   Register ``name`` as a temporal field dataset and return an opaque handle.
   The type and kind of the *scalar* ``mold`` determine the dataset type; its
   value is ignored. The returned handle is passed to later temporal writes.

``call file%write_field_data(field_var, array [, as_vector])``
   Write ``array`` to the temporal field dataset identified by the handle
   ``field_var``, and associate it with the current time step.

   If a temporal field dataset is not written for a time step its most recently
   written value, if any, is used. Before the first write, ommission yields
   an empty field data for that step. Writing the same temporal dataset more
   than once within a single time step is an error.

.. note::
   String-valued field data is not currently supported by VTKHDF/ParaView.

Workflow by Mode
----------------

.. rubric:: ``UG_STATIC``

* Write the mesh once with ``write_mesh``.
* Then write static datasets by name:

  - ``call file%write_point_data(name, array [, attribute])``
  - ``call file%write_cell_data(name, array [, attribute])``
  - ``call file%write_field_data(name, array [, as_vector])``

Time stepping and temporal registration are not allowed.

.. rubric:: ``UG_FIXED_MESH``

* Write the mesh once with ``write_mesh``.
* Write static cell and point data (may be done at any time after the mesh)

  - ``call file%write_cell_data(name, array [, attribute])``
  - ``call file%write_point_data(name, array [, attribute])``

  Static field data is not allowed in this mode.

* Register temporal datasets before the first time step:

  - ``cell_var = file%register_temporal_cell_data(name, mold [, attribute])``
  - ``point_var = file%register_temporal_point_data(name, mold [, attribute])``
  - ``field_var = file%register_temporal_field_data(name, mold)``

* Then, for each step:

  * ``call file%start_time_step(time)``
  * Write temporal datasets (or not, if allowed):

    - ``call file%write_cell_data(cell_var, array)``
    - ``call file%write_point_data(point_var, array)``
    - ``call file%write_field_data(field_var, array [, as_vector])``

  * ``call file%finalize_time_step``

.. rubric:: ``UG_MOVING_MESH``

* Write topology once with ``write_mesh_topology``.

* Register temporal datasets before the first time step:

  - ``cell_var = file%register_temporal_cell_data(name, mold [, attribute])``
  - ``point_var = file%register_temporal_point_data(name, mold [, attribute])``
  - ``field_var = file%register_temporal_field_data(name, mold)``

* Then, for each step:

  * ``call file%start_time_step(time)``
  * ``call file%write_mesh_geometry(points)`` first step, or when geometry changes
  * Write temporal datasets (or not, if allowed):

    - ``call file%write_cell_data(cell_var, array)``
    - ``call file%write_point_data(point_var, array)``
    - ``call file%write_field_data(field_var, array [, as_vector])``

  * ``call file%finalize_time_step``

The number of points must be constant across time steps.

Static datasets are not allowed in this mode.

.. rubric:: ``UG_DYNAMIC_MESH``

* Write initial mesh with ``write_mesh``, or wait until the first step.

* Register temporal datasets before the first time step:

  - ``cell_var = file%register_temporal_cell_data(name, mold [, attribute])``
  - ``point_var = file%register_temporal_point_data(name, mold [, attribute])``
  - ``field_var = file%register_temporal_field_data(name, mold)``

* Then, for each step:

  * ``call file%start_time_step(time)``

  * Write mesh for first step if not already written, or when mesh changes:

    - ``write_mesh_topology(cnode, xcnode, types)`` to replace topology only;
    - ``write_mesh_geometry(points)`` to replace geometry only;
    - ``write_mesh(points, cnode, xcnode, types)`` to replace both.

  * Write temporal datasets (or not, if allowed):

    - ``call file%write_cell_data(cell_var, array)``
    - ``call file%write_point_data(point_var, array)``
    - ``call file%write_field_data(field_var, array [, as_vector])``

  * ``call file%finalize_time_step``

Static datasets are not allowed in this mode.

This is the most general mode. Point counts, cell counts, connectivity, and
geometry may all vary from one step to the next.
