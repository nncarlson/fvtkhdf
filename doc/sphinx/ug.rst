The vtkhdf_ug_file_type module
==============================

This module defines a ``vtkhdf_ug_file`` derived type for writing
VTKHDF UnstructuredGrid files. It supports three temporal mesh levels:

* **VTKHDF_STATIC_MESH**: Static mesh (fixed geometry and topology) with temporal data
* **VTKHDF_DEFORMED_MESH**: Deformed mesh (varying geometry but fixed topology) with temporal data
* **VTKHDF_TEMPORAL_MESH**: Fully temporal mesh (geometry and topology vary over time) - *not yet implemented*

Additionally supports:

* Static point and cell datasets
* Optional time-dependent point and cell datasets

In MPI build, all type-bound procedures are collective over the
communicator passed to ``create``. Every rank must call the same method
in the same order and supply identical values for all non-distributed
arguments.

Several methods return a status code ``stat`` and an allocatable
character ``errmsg``.

* ``stat == 0`` indicates success.
* ``stat /= 0`` indicates failure and ``errmsg`` is allocated.

In MPI build, these return values are collective: all ranks return
identical values.

.. code-block:: fortran

   use vtkhdf_ug_file_type
   use vtkhdf_temporal_level
   type(vtkhdf_ug_file) :: file

File Creation and Management
----------------------------

``call file%create(filename, [comm,] stat, errmsg [,temporal_level])``
   Create a new VTKHDF "UnstructuredGrid" file.
   
   * ``filename``: path to the file to create. The recommended file extension
     is ``.vtkhdf``.
   * ``comm``: MPI communicator: either ``integer`` or ``type(MPI_Comm)``.
     In serial builds ``comm`` is omitted from the interface.
   * ``temporal_level`` (optional): specifies the temporal mesh level:
     
     - ``VTKHDF_STATIC_MESH``: static mesh with temporal data (default)
     - ``VTKHDF_DEFORMED_MESH``: deformed mesh with temporal geometry
     - ``VTKHDF_TEMPORAL_MESH``: fully temporal mesh (not yet implemented)

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
  Write the mesh geometry and topology. Use this for static meshes.
  
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

``call file%write_mesh_topology(mold, cnode, xcnode, types)``
  Write mesh topology only (for deformed and temporal meshes). Use this together with
  ``write_temporal_points`` for meshes where geometry changes over time but
  topology remains fixed. Requires ``temporal_level=VTKHDF_DEFORMED_MESH``
  or ``temporal_level=VTKHDF_TEMPORAL_MESH`` at file creation.  
  * ``mold``: ``real32`` or ``real64`` array of shape (3, `npoints`) used
    to determine data type and shape for the temporal Points dataset.
    The actual coordinate values are written later via ``write_temporal_points``.
  * ``cnode``, ``xcnode``, ``types``: same as for ``write_mesh``.

Static mesh-centered data
-------------------------
Writes static cell or point data arrays. These procedures may be called
after ``write_mesh``. For temporal files, this data is not associated with
any time step.

.. glossary::

   ``call file%write_cell_data(name, array)``
   ``call file%write_point_data(name, array)``
      Write the data ``array`` to a new cell or point dataset ``name``;
      ``name`` must not already be defined for the indicated entity type.
      Scalar data (rank-1 ``array``) and vector data
      (rank-2 ``array``) are supported. The last dimension of ``array``
      indexes the mesh entity and must have extent `ncells` (cell data)
      or `npoints` (point data).

.. note::
   VTK only supports scalar and vector-valued mesh-centered data. Other kinds
   such as tensor values must be packed into a vector.

Time-dependent mesh-centered data
---------------------------------
Temporal files support time-dependent datasets. Temporal datasets must be
registered before the first call to ``write_time_step``. After the first
time step is started, no further registrations are allowed. The first time
step must be started before temporal datasets are written.

.. glossary::

   ``call file%register_temporal_cell_data(name, mold)``
   ``call file%register_temporal_point_data(name, mold)``
      Register ``name`` as a time-dependent cell or point dataset.
      The type and kind of ``mold`` determines the dataset type.
      For scalar data, pass a scalar ``mold``,
      and for vector data, pass a rank-1 ``mold`` whose size equals
      the number of components. The value of ``mold`` is never referenced.
      ``name`` must not be defined as a static or temporal dataset of the
      same entity type.


   ``call file%write_time_step(time)``
      Start a new time step with time value ``time``.

   ``call file%write_temporal_points(points)``
      Write mesh geometry for the current time step (deformed and temporal meshes).
      This procedure is used together with ``write_mesh_topology`` for meshes
      where geometry changes over time but topology remains fixed. Must be called after
      ``write_time_step`` for each time step. Requires ``temporal_level=VTKHDF_DEFORMED_MESH``
      or ``temporal_level=VTKHDF_TEMPORAL_MESH`` at file creation.
      
      * ``points``: ``real32`` or ``real64`` array of shape (3, `npoints`)
        containing node coordinates for this time step.

   ``call file%write_temporal_cell_data(name, array)``
   ``call file%write_temporal_point_data(name, array)``
      Write the ``array`` to the temporal dataset ``name``, associating it
      with the current time step.
      ``array`` must conform to the shape and type implied by the
      registered ``mold`` and must have extent `ncells` (cell data) or
      `npoints` (point data) in its last dimension.
      
      A temporal dataset need not be written at every time step; if omitted,
      its most recently written value is used. Writing the same temporal
      dataset more than once within a single time step is an error.
