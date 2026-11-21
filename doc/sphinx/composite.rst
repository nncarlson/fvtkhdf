Composite Datasets
==================

A composite dataset is a collection of named simple VTK datasets in one file.
Its blocks may have different dataset types and writer modes. Each block has
its own mesh and point, cell, and field data; all temporal blocks share one
timeline.

Choosing the VTK Representation
-------------------------------

Choose the file type according to the composite
representation you want VTK's ``vtkHDFReader`` to produce:

.. list-table::
   :header-rows: 1

   * - Output of VTK's ``vtkHDFReader``
     - Fortran file type
     - Module
   * - ``vtkMultiBlockDataSet`` (MB)
     - ``vtkhdf_mb_file``
     - ``vtkhdf_mb_file_type``
   * - ``vtkPartitionedDataSetCollection`` (PDC)
     - ``vtkhdf_pdc_file``
     - ``vtkhdf_pdc_file_type``

The file type specifies the composite dataset type recorded in the VTKHDF file.

The two file types have identical APIs. Switching composite file types changes only
the module and file type declarations; all subsequent block creation, mesh and
data writes, time stepping, flushing, and closing calls stay the same. The
choice does not change how mesh and data contributions are distributed among
MPI ranks.

For MultiBlockDataSet output, declare:

.. code-block:: fortran

   use vtkhdf_mb_file_type
   type(vtkhdf_mb_file) :: file

For PartitionedDataSetCollection output, declare instead:

.. code-block:: fortran

   use vtkhdf_pdc_file_type
   type(vtkhdf_pdc_file) :: file

In the interfaces below, ``file`` denotes either concrete file type. Both
modules export the same block types, data handle types, ``VTKHDF_*`` mode
constants, and ``vtkhdf_version``.

The file owns its blocks. ``add_block`` returns a typed pointer through which
you write the block's mesh and data directly. The owning composite file manages
file creation, time stepping, flushing, and closing.

By default, a flat assembly is created; see :ref:`composite-assembly-construction`
for optional hierarchical organization.

Supported Blocks and Modes
--------------------------

The pointer type passed to ``add_block`` selects the simple VTK dataset type,
with the same choices and modes for MB and PDC:

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

.. method-doc::

   call file%create(filename, [comm,] stat, errmsg [, data_chunk_bytes])

   Create a new VTKHDF file with the composite dataset type selected by the
   file type declaration. Modes are selected per block; there is no file-level
   ``mode`` argument.

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
   perform collective cleanup. It releases local heap storage only; explicit
   ``close`` is required to release MPI/HDF5 resources.

In MPI builds, file and block I/O operations are collective over the communicator
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

In flat output, VTK's zero-based block indices follow the order of ``add_block``
calls.

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
applies to UnstructuredGrid blocks. HyperTreeGrid blocks use the cell methods
described in :doc:`hypertree`; they have no point data. All block types use
:doc:`field_data`, including the ``data=`` payload keyword and ``as_vector``
interpretation.

The ``vtkhdf_cell_data_handle``, ``vtkhdf_point_data_handle``, and
``vtkhdf_field_data_handle`` types are exported by both composite modules.
Each handle belongs to the block that created it. Supply it directly to that
block's write method.

Each block has independent namespaces for point, cell, and field datasets.
Within each namespace, creation by name reserves a new name; repeating a name
does not select an existing dataset. Subsequent writes use a data handle where
supported by the simple writer. Field data belongs to an individual block;
the composite file has no field-data write interface of its own.

When several blocks use the same point or cell variable name, keep its
quantity, numeric type, and component count consistent across those blocks so
readers can treat it as one variable across the collection.

Shared Time Stepping
--------------------

The file is temporal when at least one block uses a temporal mode. A collection
with only static blocks has no time stepping.

Register each block's temporal point, cell, and field datasets through that
block before the first file-level time step. Perform time stepping on the file,
so every temporal block receives the same steps and time values.

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

This serial example selects PDC output and writes rectilinear ``fluid`` and
unstructured ``solid`` blocks with a shared temperature time series. The fluid
coordinates change between steps; the solid mesh remains fixed. For MB output,
substitute the module and file type declarations shown above; the rest of the
program is unchanged.

.. code-block:: fortran

   program composite_example
     use,intrinsic :: iso_fortran_env, only: int8
     use vtkhdf_pdc_file_type
     implicit none
     type(vtkhdf_pdc_file) :: file
     type(vtkhdf_rectilinear_block), pointer :: fluid
     type(vtkhdf_unstructured_block), pointer :: solid
     type(vtkhdf_point_data_handle) :: fluid_temperature, solid_temperature
     integer :: bounds(3,2), stat, step
     real :: fluid_values(2,2,2), solid_values(4), solid_points(3,4)
     character(:), allocatable :: errmsg

     call file%create('composite.vtkhdf',stat,errmsg)
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
:repo:`example/pdc_example1.F90` demonstrates UnstructuredGrid and ImageData
blocks with a shared timeline in serial and MPI builds. Both use the composite
API described here.

.. _composite-assembly-construction:

Assembly Construction
---------------------

An assembly organizes blocks into a hierarchy of nodes and block references.
Its root node is created automatically and has index 0. In the default flat
assembly, all block references are direct children of the root.

The optional assembly construction methods group blocks under additional nodes.
A block may appear under more than one node while retaining one mesh, one set
of data, and one writer object. Construction changes how the blocks are
organized; it does not affect their mesh or data-writing interfaces.

For MB output, VTK's reader represents ordinary assembly nodes as
``vtkMultiBlockDataSet`` objects, nested according to the hierarchy. For PDC
output, it represents the hierarchy as a ``vtkDataAssembly`` associated with
the ``vtkPartitionedDataSetCollection`` object created by the reader.

.. method-doc::

   node = file%add_node(name [, parent])

   Add an ordinary hierarchy node and return its stable positive integer index.
   ``parent`` is the index of a previously added ordinary node; omitted or zero
   means the root. Node indices are local to the owning file.

   ``name`` follows the normalization and disambiguation rules in :doc:`common`.
   A blank name defaults to ``Unnamed_Node``. Names are resolved when the
   assembly is committed; disambiguation does not change node indices.

.. method-doc::

   call file%add_block_ref(dataset [, parent])

   Add an explicit reference to an existing block under ``parent``. Omitted or
   zero ``parent`` means the root. ``dataset`` is a typed block pointer returned
   by this file's ``add_block``. Standalone writers and blocks belonging to
   another file are invalid.

   Adding the same block under the same parent again is a no-op and retains
   the first reference's position. References under different parents are
   allowed, including an explicit root reference alongside deeper placements.

An explicit reference suppresses that block's automatic root placement. Blocks
without explicit references remain at the root, so no stored block is omitted.
Explicit children appear in insertion order, including interleaved ordinary
nodes and block references. Implicit root references follow all explicit root
children in block-creation order.

Every reference uses its block's established, globally disambiguated name.
Ordinary node names need uniqueness only among siblings. Block-reference names
have priority within each parent; colliding ordinary node names are disambiguated
automatically. This includes collisions with blocks implicitly retained at the
root. Identical ordinary node names under different parents are allowed.

Block creation, mesh/data writes, and assembly construction may be interleaved.
Construction changes only the pending description: the automatic flat assembly
remains in the output until commitment. ``flush`` does not commit a custom
assembly.

.. method-doc::

   call file%set_assembly()

   Commit the assembly once while the file is open. This operation is collective
   in MPI builds. Calling it without custom construction is valid and preserves
   the default flat organization, with all blocks at the assembly root.

   Successful commitment freezes block creation, node creation, and explicit
   reference creation. A second ``set_assembly`` call is invalid. Existing block
   pointers and data handles remain valid for their usual mesh, data, and time
   operations; ``flush`` and ``close`` remain available.

Normal ``close`` implicitly commits pending construction before closing block
resources. With no custom construction, it preserves the automatic flat
organization. Automatic Fortran finalization does not commit an assembly;
explicit collective cleanup still requires ``close``.

Commitment may occur before or after time stepping starts, including during an
open step. The usual restriction on adding temporal blocks after the first step
still applies. The construction API provides no deletion or move operations.

In MPI builds, ``add_node`` and ``add_block_ref`` perform no communication, but
all ranks must make matching construction calls with the same names, parent
indices, and corresponding block objects. Rank-zero-only construction is not
supported.

Invalid parents, foreign block objects, and operations forbidden after
commitment are errors.

For example, reference a block in two groups and explicitly commit:

.. code-block:: fortran

   regions = file%add_node('Regions')
   call file%add_block_ref(fluid,parent=regions)
   call file%add_block_ref(solid,parent=regions)
   materials = file%add_node('Materials')
   call file%add_block_ref(fluid,parent=materials)
   call file%set_assembly()

:repo:`example/composite_assembly.F90` provides a complete serial/MPI example
with two placements of the fluid block.

VTK Reader Compatibility
------------------------

These compatibility notes concern VTK's ``vtkHDFReader``, including its use in
ParaView. Each block retains the VTK reader version requirements of its simple
dataset type;
see :doc:`image`, :doc:`rectilinear`, :doc:`structured_grid`, and :doc:`hypertree`.

The writers target VTKHDF version 2.8. Reader shortcomings do not change the
supported block types or modes. The following known reader limitations concern
both composite representations; they have been observed with VTK 9.7.2,
including the VTK build bundled with ParaView 6.2.1:

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
   * - Composite temporal state
     - The reader can lose time-step information and fail to read point/cell
       arrays in ImageData, RectilinearGrid, and StructuredGrid blocks.
     - Write the affected datasets as standalone files. Changing array
       selections after metadata discovery and before the first data request
       can avoid the state reset in the tested reader; this is a reader-specific
       workaround, not a change to the file format.
