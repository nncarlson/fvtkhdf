.. toctree::
   :hidden:

.. role:: f(code)
   :language: fortran

===================
The fVTKHDF Library
===================

The vtkhdf_ug_file_type module
==============================

This module defines the `vtkhdf_ug_file` derived type for exporting mesh-based
solution data to a VTKHDF file readable by ParaView. It produces an
UnstructuredGrid (UG) dataset and supports both static and time-dependent cell
and point data on a single static mesh.

.. code-block:: fortran

   use vtkhdf_ug_file_type
   type(vtkhdf_ug_file) :: file

File procedures
---------------

:f:`call file%create(filename, comm, stat, errmsg [,is_temporal])`
   Initializes a new VTKHDF unstructured grid file.

   - `filename` (character(*), intent(in)): The name of the file to create.
   - `comm` (integer, intent(in)): MPI communicator.
   - `stat` (integer, intent(out)): Status code (0 for success).
   - `errmsg` (character(:), allocatable, intent(out)): Error message if `stat` is non-zero.
   - `is_temporal` (logical, intent(in), optional): Set to `.true.` if the file will contain time-dependent data.

``call file%close()``
    Close the file.

``call file%flush()``
    Flush the internal HDF5 buffers associated with the file and request
    the operating system to flush the system buffers for the open file.

Mesh data
-----------
Writes the unstructured mesh geometry and topology to the file. This procedure **must** be called before any data writing or registration procedures.

``call file%write_mesh(x, cnode, xcnode, types)``
   TEXT

Static mesh data
----------------

.. glossary::

   ``call file%write_cell_data(name, array)``
   ``call file%write_point_data(name, array)``
      Write the cell/point-based data `array` to a new cell/point dataset
      `name`. Scalar-valued data (rank-1 `array`) and vector-valued data
      (rank-2 `array`) are supported. For temporal files, this dataset is
      static and not associated with any time step.
   

Time-dependent mesh data
------------------------

.. glossary::

   ``call file%register_temporal_cell_data(name, mold)``
   ``call file%register_temporal_point_data(name, mold)``
      TEXT


   ``call file%write_time_step(time)``
      TEXT

   ``call file%write_temporal_cell_data(name, array)``
   ``call file%write_temporal_point_data(name, array)``
      TEXT

