Common Conventions
==================

Writers for simple VTK datasets use the same conventions for data associations,
names, handles, and collective calls. Mesh-specific array shapes and mesh-writing
requirements are described in the corresponding reference chapters.

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

Point, cell, and field data have separate namespaces within each VTK dataset.
Names are trimmed; blank names are replaced by a default. The characters
``/``, ``.``, and space are replaced by underscores. Duplicate names, including
names that become identical after these substitutions, are made unique by
appending a suffix such as ``_1``. The resulting name is what VTK and ParaView
see. The original name is written as the ``Name`` attribute on the HDF5
dataset for user reference.

Each write by name or temporal registration creates a new dataset.
Repeating a name does not select an existing array. Where subsequent writes
are supported, they identify the array by an opaque handle:
``vtkhdf_cell_data_handle``, ``vtkhdf_point_data_handle``, or
``vtkhdf_field_data_handle``. These types are exported by the file-type
modules. Store the returned handle and use it only with its owning file.

Collective Calls
----------------

In MPI builds, all writer methods are collective over the communicator passed
to ``create``. All ranks must call the same methods in the same order.
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
UnstructuredGrid supports all four modes.
In moving-mesh mode, the structured writers also allow static point/cell data;
UnstructuredGrid requires temporal point/cell data.
Static field data is supported in every writer mode.

Temporal data must be registered before the first ``start_time_step`` call.
Writes by handle then supply data for the current step. A step is completed
by ``finalize_time_step``; starting another step or closing the file also
finalizes an open step. Rules for required writes and reuse of omitted values
depend on the data association and VTK dataset type and are specified in their
reference chapters.
