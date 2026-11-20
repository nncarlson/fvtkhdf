MultiBlockDataSet
=================

The ``vtkhdf_mb_file_type`` module defines ``vtkhdf_mb_file`` for writing a
flat collection of named simple VTK datasets to one VTKHDF file. Blocks may
have different dataset types and writer modes. Each has its own mesh and
point, cell, and field data; all temporal blocks share one timeline.
Hierarchical nesting and PartitionedDataSetCollection output are not supported.

The file owns its blocks. ``add_block`` returns a typed pointer through which
you write the block's mesh and data directly. File creation, time stepping,
flushing, and closing are performed on the owning ``vtkhdf_mb_file``.

Supported Blocks and Modes
--------------------------

All five block types and the shared ``VTKHDF_*`` mode constants are exported
by ``vtkhdf_mb_file_type``. The pointer type selects the VTK dataset type:

.. list-table::
   :header-rows: 1

   * - VTK dataset type
     - Block pointer type
     - Supported modes
   * - :doc:`ImageData <image>`
     - ``vtkhdf_image_block``
     - Static, fixed mesh
   * - :doc:`RectilinearGrid <rectilinear>`
     - ``vtkhdf_rectilinear_block``
     - Static, fixed mesh, moving mesh
   * - :doc:`StructuredGrid <structured_grid>`
     - ``vtkhdf_structured_block``
     - Static, fixed mesh, moving mesh
   * - :doc:`UnstructuredGrid <unstructured>`
     - ``vtkhdf_unstructured_block``
     - Static, fixed mesh, moving mesh, dynamic mesh
   * - :doc:`HyperTreeGrid <hypertree>`
     - ``vtkhdf_hypertree_block``
     - Static, fixed mesh, moving mesh, dynamic mesh

The modes correspond to ``VTKHDF_STATIC``, ``VTKHDF_FIXED_MESH``,
``VTKHDF_MOVING_MESH``, and ``VTKHDF_DYNAMIC_MESH``; see :doc:`common`.
A block has the same geometry, data-writing interfaces, and reuse rules as
the standalone writer for its dataset type. The corresponding chapters give
its mesh arguments, array shapes, region placement, and mode restrictions.

File Creation and Management
----------------------------

.. code-block:: fortran

   use vtkhdf_mb_file_type
   type(vtkhdf_mb_file) :: file

.. method-doc::

   call file%create(filename, [comm,] stat, errmsg [, data_chunk_bytes])

   Create a new VTKHDF MultiBlockDataSet file. Modes are selected per block;
   there is no file-level ``mode`` argument.

   * ``filename``: output file path, normally ending in ``.vtkhdf``.
   * ``comm``: MPI communicator, omitted in serial builds. With
     ``USE_MPI_F08=ON`` (default), either ``integer`` or ``type(MPI_Comm)`` is
     accepted. With ``USE_MPI_F08=OFF``, only ``integer`` is accepted.
   * ``stat``: output status code; zero indicates success.
   * ``errmsg``: allocated error message when ``stat /= 0``.
   * ``data_chunk_bytes``: optional positive ``int64`` target chunk size,
     defaulting to ``128*1024`` (128 KiB), used by every block.

.. method-doc::

   call file%flush()

   Flush the shared file's HDF5 buffers.

.. method-doc::

   call file%close()

   Finalize any open time step, close every block and the shared file, and
   release their resources. All borrowed block pointers and their data handles
   become invalid. Call this explicitly; automatic Fortran finalization cannot
   perform collective cleanup.

In MPI builds, file and block operations are collective over the communicator
passed to ``create``. All ranks must create the same blocks, with the same
types, names, and modes, and make calls in the same order. Distributed mesh
and point/cell data follow the corresponding simple writer's rules. Field
data is supplied on every rank, with matching type and shape; only rank zero's
values are written. See :doc:`common`.

Block Creation and Ownership
----------------------------

.. method-doc::

   call file%add_block(name, dataset [, mode])

   Create a block and return a pointer to it in ``dataset``.

   * ``name``: block name. Names are normalized and disambiguated within the
     file according to :doc:`common`; the original name is preserved in the
     block group's ``Name`` attribute.
   * ``dataset``: pointer of one of the five concrete block types above.
     The file allocates and initializes its target; the caller does not
     allocate it.
   * ``mode``: optional mode supported by that block type; defaults to
     ``VTKHDF_STATIC``.

Pointers to existing blocks remain valid when more blocks are added. The file
retains ownership of their targets: use the pointers for mesh and data writes,
and let the file manage their lifecycle. After ``file%close()``, discard or
nullify the pointers rather than accessing or deallocating their old targets.

Temporal blocks must be added before the first ``file%start_time_step`` call.
Static blocks may also be added later. Static and temporal blocks may coexist
in the writer, subject to the reader limitations below.

For example, create a RectilinearGrid block named ``fluid`` and an
UnstructuredGrid block named ``solid``:

.. code-block:: fortran

   type(vtkhdf_rectilinear_block), pointer :: fluid
   type(vtkhdf_unstructured_block), pointer :: solid

   call file%add_block('fluid', fluid, mode=VTKHDF_MOVING_MESH)
   call file%add_block('solid', solid, mode=VTKHDF_FIXED_MESH)

Mesh and Data Writes
--------------------

Call mesh and data methods on the block pointer, with the same arguments as
for its standalone writer. For example, the RectilinearGrid block uses
``fluid%write_mesh(point_bounds, x, y, z)``; the UnstructuredGrid block
uses ``solid%write_mesh(points, cnode, xcnode, types)``. RectilinearGrid and
StructuredGrid also provide their usual ``write_mesh_geometry`` methods.

The :doc:`structured` interface applies to point/cell data in ImageData,
RectilinearGrid, and StructuredGrid blocks. The :doc:`unstructured` interface
applies to UnstructuredGrid blocks. These four types use :doc:`field_data`,
including the ``data=`` payload keyword and ``as_vector`` interpretation.
HyperTreeGrid blocks use the cell and field methods described in
:doc:`hypertree`, with the ``array=`` payload keyword; they have no point data.

The ``vtkhdf_cell_data_handle``, ``vtkhdf_point_data_handle``, and
``vtkhdf_field_data_handle`` types are exported by ``vtkhdf_mb_file_type``.
Each handle belongs to the block that created it. Supply it directly to that
block's write method.

Each block has independent namespaces for point, cell, and field datasets.
Within each namespace, creation by name reserves a new name; repeating a name
does not select an existing dataset. Subsequent writes use a data handle where
supported by the simple writer. Field data belongs to an individual block;
the MB file has no composite-level field-data write interface.

When several blocks use the same point or cell variable name, keep its
quantity, numeric type, and component count consistent across those blocks so
readers can treat it as one variable across the collection.

Shared Time Stepping
--------------------

The file is temporal when at least one block uses a temporal mode. A collection
with only static blocks has no time stepping.

Register each block's temporal point, cell, and field datasets through that
block before the first file-level time step. Perform time stepping on the file, so every
temporal block receives the same steps and time values.

.. method-doc::

   call file%start_time_step(time)

   Start a step in every temporal block. ``time`` is a scalar ``real32`` or
   ``real64``; only rank zero's value is used in MPI builds. Starting another
   step first finalizes an open step. Static blocks are unaffected.

Within a step, update each block's geometry and data as required by its mode.
Omitted writes follow that block type's reuse rules. For example, a moving
RectilinearGrid block can reuse its coordinate vectors independently of its
point/cell arrays. Static field data may be written at any point after block
creation and before closing the file, including after time stepping.

.. method-doc::

   call file%finalize_time_step()

   Complete the open step in every temporal block. Closing the file also
   finalizes an open step.

Example
-------

This serial example writes the rectilinear ``fluid`` and unstructured ``solid``
blocks with a shared temperature time series. The fluid coordinates change
between steps; the solid mesh remains fixed.

.. code-block:: fortran

   program multiblock_example
     use,intrinsic :: iso_fortran_env, only: int8
     use vtkhdf_mb_file_type
     implicit none
     type(vtkhdf_mb_file) :: file
     type(vtkhdf_rectilinear_block), pointer :: fluid
     type(vtkhdf_unstructured_block), pointer :: solid
     type(vtkhdf_point_data_handle) :: fluid_temperature, solid_temperature
     integer :: bounds(3,2), stat, step
     real :: fluid_values(2,2,2), solid_values(4), solid_points(3,4)
     character(:), allocatable :: errmsg

     call file%create('multiblock.vtkhdf',stat,errmsg)
     if (stat /= 0) error stop errmsg
     call file%add_block('fluid',fluid,mode=VTKHDF_MOVING_MESH)
     call file%add_block('solid',solid,mode=VTKHDF_FIXED_MESH)
     bounds(:,1) = [0,0,0]
     bounds(:,2) = [1,1,1]
     call fluid%write_mesh(bounds,[0.,1.],[0.,1.],[0.,1.])
     solid_points = reshape([2.,0.,0., 3.,0.,0., 2.,1.,0., 2.,0.,1.],[3,4])
     call solid%write_mesh(solid_points,[1,2,3,4],[1,5],[10_int8])
     call fluid%write_field_data('material',data=1)
     call solid%write_field_data('material',data=2)
     fluid_temperature = fluid%register_temporal_point_data('temperature',0.)
     solid_temperature = solid%register_temporal_point_data('temperature',0.)

     do step = 0,2
       call file%start_time_step(real(step))
       if (step > 0) call fluid%write_mesh_geometry([0.,1.+0.1*step],[0.,1.],[0.,1.])
       fluid_values = real(step)
       call fluid%write_point_data(fluid_temperature,fluid_values,[0,0,0])
       solid_values = real(step)+10.
       call solid%write_point_data(solid_temperature,solid_values)
       call file%finalize_time_step()
     end do
     call file%close()
     nullify(fluid,solid)
   end program

:repo:`example/mb_example1.F90` demonstrates UnstructuredGrid blocks, distributed
mesh partitions, and temporal point/cell data in serial and MPI builds.

Reader Compatibility
--------------------

Each block retains the reader version requirements of its simple dataset type;
see :doc:`image`, :doc:`rectilinear`, :doc:`structured_grid`, and :doc:`hypertree`.

The VTKHDF reader in VTK currently has several issues that may limit the use
of composite files. These issues are being actively addressed. Use the latest
stable release of ParaView or VTK for best results.

The following limitations are known in ParaView 6.2.1 (VTK 9.7.2):

.. list-table::
   :header-rows: 1
   :widths: 25 45 30

   * - Feature
     - Known limitation
     - Workaround
   * - Mixed static and temporal blocks
     - Combining static and temporal blocks can lose time information, omit
       blocks or data, or terminate the reader.
     - Make every block temporal if any block is temporal. Use
       ``VTKHDF_FIXED_MESH`` for otherwise static blocks; these need no
       temporal-data registration or repeated mesh or data writes.
   * - Block field data
     - Static and temporal field arrays stored on composite blocks are
       omitted by the reader.
     - Use standalone output when field data is required. For hypertree files,
       also observe the :ref:`field-cache limitation <hypertree-reader-limitations>`.
   * - Temporal data read directly with VTK
     - The reader can lose time-step information and fail to read point/cell
       arrays in ImageData, RectilinearGrid, and StructuredGrid blocks.
     - Write the affected datasets as standalone files.
