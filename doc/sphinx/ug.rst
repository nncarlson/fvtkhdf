The vtkhdf_ug_file_type module
==============================

This module defines the ``vtkhdf_ug_file`` derived type for writing
VTKHDF UnstructuredGrid files. It supports:

* A static mesh
* Static point and cell datasets
* Optional time-dependent point and cell datasets

In the MPI build, all type-bound procedures are collective over the
communicator passed to ``create``. Every rank must call the same method
in the same order and supply identical values for all non-distributed
arguments.

.. code-block:: fortran

   use vtkhdf_ug_file_type
   type(vtkhdf_ug_file) :: file
   type(vtkhdf_cell_data_handle) :: cell_var
   type(vtkhdf_point_data_handle) :: point_var

File Creation and Management
----------------------------

``call file%create(filename, [comm,] stat, errmsg [,is_temporal])``
   Create a new VTKHDF "UnstructuredGrid" file.
   
   * ``filename``: path to the file to create. The recommended file extension
     is ``.vtkhdf``.
   * ``comm``: the MPI communicator: either ``integer`` or ``type(MPI_Comm)``.
     In serial builds ``comm`` is omitted from the interface.
   * ``is_temporal`` (optional): set ``.true.`` to enable time-dependent
     datasets. The default is ``.false.``
   * ``stat``: status code. ``stat == 0`` indicates success.
   * ``errmsg``: allocated error message when ``stat /= 0``.

   In the MPI build, ``stat`` and ``errmsg`` are collective outputs:
   all ranks return identical values.

``call file%close()``
    Close the file and release internal resources. Users should *always* call
    this to "finalize" the object; automatic finalization cannot perform a
    proper collective cleanup. 
    

``call file%flush()``
    Flush the file's HDF5 buffers and request the OS to flush the file buffers.

Mesh Data
---------
Writes the portion of the unstructured mesh provided by the calling MPI rank.
The mesh must be written before any mesh-centered data is written.

``call file%write_mesh(points, cnode, xcnode, types)``
  * ``points``: ``real32`` or ``real64`` array of shape (3, `npoints`)
    containing the node coordinates. Coordinates are always interpreted
    as 3D; for 1D or 2D geometries, the unused components must be set
    (e.g., to 0.0).
  * ``cnode``, ``xcnode``: ``integer32`` or ``integer64`` arrays
    describing the mesh topology in CSR format:

    - ``cnode`` contains the concatenated cell-node connectivity data.
    - ``xcnode`` contains the starting index of each cell's connectivity
      in ``cnode``. For cell ``i``, ``cnode(xcnode(i):xcnode(i+1)-1)`` gives
      its node list. ``size(xcnode)`` must equal `ncells`\ +1 and
      its final element equals ``size(cnode)+1``.
      
    CSR indexing is 1-based (Fortran style). Conversion to
    0-based indexing required by VTKHDF is handled internally.
  * ``types``: an ``int8`` array of size `ncells` containing VTK cell type
    codes. Named constants such as ``VTK_TETRA`` and ``VTK_HEXAHEDRON``
    are provided by ``vtkhdf_vtk_cell_types``.

Static mesh-centered data
-------------------------
Writes static cell or point data arrays. These procedures may be called
after ``write_mesh``. For temporal files, this data is not associated with
any time step.

.. glossary::

   ``call file%write_cell_data(name, array)``
   ``call file%write_point_data(name, array)``
      Write the data ``array`` to a new cell or point dataset ``name``.
      Scalar data (rank-1 ``array``) and vector data
      (rank-2 ``array``) are supported. The last dimension of ``array``
      indexes the mesh entity and must have extent `ncells` (cell data)
      or `npoints` (point data).

      Cell and point datasets each have their own namespace. Input names are
      trimmed, empty input is replaced by a default name, the characters
      ``/``, ``.``, and space are replaced by ``_``, and duplicate internal
      names are made unique by appending a suffix such as ``_1``. The
      sanitized internal dataset name is what the VTKHDF reader and ParaView
      use. The original user-facing name is still written to the dataset
      ``Name`` attribute for informational purposes.

.. note::
   VTK only supports scalar and vector-valued mesh-centered data. Other kinds
   such as tensor values must be packed into a vector.

Time-dependent mesh-centered data
---------------------------------
Temporal files support time-dependent datasets. Temporal datasets must be
registered before the first call to ``start_time_step``. After the first
time step is started, no further registrations are allowed. The first time
step must be started before temporal datasets are written.

.. glossary::

   ``cell_var = file%register_temporal_cell_data(name, mold)``
   ``point_var = file%register_temporal_point_data(name, mold)``
      Register ``name`` as a time-dependent cell or point dataset.
      The type and kind of ``mold`` determines the dataset type.
      For scalar data, pass a scalar ``mold``,
      and for vector data, pass a rank-1 ``mold`` whose size equals
      the number of components. The value of ``mold`` is never referenced.

      Naming follows the same normalization and disambiguation rules as
      ``write_cell_data`` and ``write_point_data``. The returned handle is
      opaque; user code should store it and pass it to later temporal writes.


   ``call file%start_time_step(time)``
      Start a new time step with time value ``time``.

      After ``start_time_step`` is called for the first time step, every
      registered temporal dataset must be written once for that step.


   ``call file%write_temporal_cell_data(cell_var, array)``
   ``call file%write_temporal_point_data(point_var, array)``
      Write the ``array`` to the temporal dataset identified by the handle,
      associating it
      with the current time step.
      ``array`` must conform to the shape and type implied by the
      registered ``mold`` and must have extent `ncells` (cell data) or
      `npoints` (point data) in its last dimension.

      After the first time step, a temporal dataset need not be written at
      every time step; if omitted,
      its most recently written value is used. Writing the same temporal
      dataset more than once within a single time step is an error.

   ``call file%finalize_time_step()``
      Finalize the current time step. Until this call, the file is in an
      in-progress state for that step.
      Best practice is to call ``finalize_time_step`` immediately after all
      temporal writes for the step are complete.

      ``finalize_time_step`` is called implicitly when ``start_time_step``
      begins a new step while one is still open, and when ``close`` is called.
