Common Conventions
==================

Writers for simple VTK datasets use the same conventions for data associations,
names, handles, and collective calls. Mesh-specific array shapes and mesh-writing
requirements are described in the corresponding reference chapters.
These conventions also apply to the typed blocks of a :doc:`composite file
<composite>` (MB or PDC); each block is a separate simple VTK dataset.

The caller is responsible for supplying all required mesh and data writes,
including complete coverage of datasets written in regions. Closing a file
releases resources and finalizes any open time step; it does not certify that
the output is complete.

Data Associations and Components
--------------------------------

* **Point data** supplies one tuple per mesh point.
* **Cell data** supplies one tuple per mesh cell.
* **Field data** supplies values associated with the VTK dataset as a whole,
  independently of mesh points and cells. See :doc:`field_data`.

A tuple contains one or more numeric components. Scalar data has one component
per tuple; multicomponent data may represent vectors, tensors, or other grouped
values. For multicomponent Fortran arrays, the first dimension indexes the
components. The remaining dimensions index tuples, using a mesh-dependent
layout for point and cell data.

The optional ``attribute`` argument to point and cell data creation or
registration designates a VTK role, such as ``Scalars``, ``Vectors``,
``Tensors``, or ``GlobalIds``. This is separate from the array's component
count. Blank values are ignored. See the VTK dataset chapters for supported
shapes and the VTKHDF specification for the full set of
`attribute roles <https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/vtkhdf_specifications.html#attribute-data>`_.
These designations are supported by VTK 9.7.0 and later and ParaView 6.2.0
and later; earlier versions ignore them.

Names and Handles
-----------------

Composite block names form a namespace within the composite file. Point, cell,
and field data have separate namespaces within each simple VTK dataset,
including each block of a composite. The following naming rules apply to both
blocks and data arrays.

Names are trimmed; blank names are replaced by a default. The characters
``/``, ``.``, and space are replaced by underscores. Duplicate names, including
names that become identical after these substitutions, are made unique by
appending a suffix such as ``_1``. The resulting name is what VTK and ParaView
see. The original name is written as the ``Name`` attribute on the HDF5
dataset or block group for user reference.

Each write by name or temporal registration creates a new dataset.
Repeating a name does not select an existing dataset. Where subsequent writes
are supported, they identify the dataset by an opaque handle:
``vtkhdf_cell_data_handle``, ``vtkhdf_point_data_handle``, or
``vtkhdf_field_data_handle``. These types are exported by the file-type
modules. Store the returned handle and use it only with the standalone writer
or composite block that created it.

Collective Calls
----------------

In MPI builds, writer methods are collective over the communicator passed
to ``create``. The in-memory ``add_node`` and ``add_block_ref`` methods of
:doc:`composite files <composite>` perform no communication, but still require
matching construction calls on all ranks. All ranks must call the same methods
in the same order.
Non-distributed arguments, including filenames, modes, names, and attribute
designations, must have identical values on all ranks. Distributed arrays and
placement arguments may differ as specified by the writer interface.
For ``start_time_step``, only rank zero's time value is used.

Field data has a distinct rule: all ranks supply arguments with matching type,
kind, and shape, but only rank zero contributes values.

In serial builds, the same methods are used with the communicator argument
omitted from ``create``. Always call ``close`` explicitly; automatic Fortran
finalization cannot perform collective cleanup.

Writer Modes and Time Stepping
------------------------------

Each file-type module exports the mode constants supported by its writer.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Mode
     - Behavior
   * - ``VTKHDF_STATIC``
     - Fixed mesh and static data; no time stepping. This is the default.
   * - ``VTKHDF_FIXED_MESH``
     - Fixed mesh and static or temporal data.
   * - ``VTKHDF_MOVING_MESH``
     - Fixed topology and temporal geometry; data may be static or temporal
       as supported by the writer.
   * - ``VTKHDF_DYNAMIC_MESH``
     - Temporal topology, geometry, and point/cell data.

ImageData supports ``VTKHDF_STATIC`` and ``VTKHDF_FIXED_MESH``.
RectilinearGrid and StructuredGrid additionally support ``VTKHDF_MOVING_MESH``.
UnstructuredGrid and HyperTreeGrid support all four modes.
In moving-mesh mode, point/cell data may be static or temporal.
Static field data is supported in every writer mode.

Static data is independent of time stepping and has the same value at every step.
Temporal data must be registered before the first ``start_time_step`` call.
Writes by handle then supply data for the current step. A step is completed
by ``finalize_time_step``; starting another step or closing the file also
finalizes an open step. Rules for required writes and reuse of omitted values
depend on the data association and VTK dataset type and are specified in their
reference chapters.
