Field Data
==========

Field data consists of values associated with a simple VTK dataset as a whole
rather than with its mesh points or cells. All writers for simple VTK datasets
provide the interface described here. Field data writes need not be ordered
with respect to mesh writes.

In the signatures below, ``file`` denotes a file writer object for a simple
VTK dataset; see :doc:`overview` for the available writer types and their
modules.

Values and Interpretation
-------------------------

The ``array`` argument may be a scalar, a rank-1 array, or a rank-2 array of
``integer(int32)``, ``integer(int64)``, ``real(real32)``, or ``real(real64)``.
The kind constants are available from ``iso_fortran_env``.

VTK and ParaView interpret these values as tuples with one or more components:

.. list-table::
   :header-rows: 1
   :widths: 45 25 30

   * - Fortran argument
     - Number of tuples
     - Components per tuple
   * - Scalar
     - 1
     - 1
   * - Rank-1 ``array(n)``
     - ``n``
     - 1
   * - Rank-1 ``array(n)``, ``as_vector=.true.``
     - 1
     - ``n``
   * - Rank-2 ``array(ncomp,ntuples)``
     - ``ntuples``
     - ``ncomp``

``as_vector`` is optional and defaults to ``.false.``. It may be ``.true.``
only for rank-1 arguments. It changes the grouping into components and tuples;
it does not require three components or assign a geometric vector role.
For a rank-2 argument, each Fortran column is one tuple.

In MPI builds, all field data calls are collective. All ranks must supply
matching argument types, kinds, and shapes. Only rank 0 contributes the
written values; values on other ranks are never referenced.

.. note::

   String-valued field data is not currently supported by fVTKHDF due to
   incompatibilities between parallel HDF5 writing and VTK's string reading.

Static Field Data
-----------------

``call file%write_field_data(name, array [, as_vector])``
   Write a complete static field dataset.

   * ``name``: character name, handled according to :doc:`common`.
   * ``array``: complete scalar or array value, interpreted as described above.
   * ``as_vector``: optional logical controlling rank-1 interpretation.

   Static field data is allowed only in ``VTKHDF_STATIC`` mode. Each call
   creates a new field, including when a name is repeated. No handle is
   returned, and writes are not incremental.

Temporal Field Data
-------------------

.. code-block:: fortran

   ! The handle type is exported by each file-type module.
   type(vtkhdf_field_data_handle) :: field_var

``field_var = file%register_temporal_field_data(name, mold)``
   Register a time-dependent field dataset and return an opaque handle.

   * ``name``: character name, handled according to :doc:`common`.
   * ``mold``: scalar determining the numeric type and kind. Its value is
     ignored. It does not constrain the rank, shape, or component count of
     later writes.

   Registration is allowed only in a temporal mode and must precede the first
   call to ``start_time_step``.

``call file%write_field_data(field_var, array [, as_vector])``
   Write the complete value for the current time step to the registered field.

   * ``field_var``: handle returned by registration.
   * ``array``: scalar or array with the registered type and kind.
   * ``as_vector``: optional logical controlling rank-1 interpretation.

   A time step must be open. Each field may be written at most once per step.
   The rank, extents, and number of components may change between steps.

If a field is omitted from a step, its most recently written value is reused.
Before the first write, omission gives an empty field for that step. These
rules are independent of the first-step requirements for point and cell data.

Example
-------

The following calls write three scalar tuples and one three-component tuple,
respectively, in a static file:

.. code-block:: fortran

   call file%write_field_data('values', [1.0,2.0,3.0])
   call file%write_field_data('direction', [1.0,0.0,0.0], as_vector=.true.)
