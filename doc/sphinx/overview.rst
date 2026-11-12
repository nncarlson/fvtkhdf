========
Overview
========

fVTKHDF is a modern Fortran library for writing
`VTKHDF <https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html>`_
files readable by VTK and ParaView. VTKHDF is an HDF5-based file format designed
for scalable scientific data storage and visualization.

fVTKHDF is designed primarily for MPI-parallel high-performance computing (HPC)
applications. A serial build is also available, preserving the same high-level
API.

Dataset Types and Writers
=========================

The VTKHDF format distinguishes *simple* VTK dataset types, such as ImageData
and UnstructuredGrid, from *composite* types, such as MultiBlockDataSet and
PartitionedDataSetCollection. Composite types are containers for simple VTK
datasets. Each simple VTK dataset has its own mesh and associated point, cell,
and field data.

Currently supported VTK dataset types:

* **ImageData, RectilinearGrid, and StructuredGrid**: logically rectangular
  grids. ImageData uses uniform spacing; RectilinearGrid is defined by the
  tensor product of three coordinate vectors; StructuredGrid gives each point
  explicit physical coordinates. RectilinearGrid and StructuredGrid also
  support moving geometry with fixed logical grids.

* **UnstructuredGrid (UG)**: meshes with arbitrary cell connectivity,
  supporting moving geometry and changing topology.

* **HyperTreeGrid (HTG)**: tree-based adaptive mesh refinement on a
  rectilinear grid of root cells, supporting moving geometry and changing
  refinement.

* **MultiBlockDataSet (MB)**: a composite type containing a flat collection
  of UnstructuredGrid blocks.

fVTKHDF provides a file writer for each supported VTK dataset type. Each
writer is a Fortran derived type exported by its corresponding module.
Writers for simple VTK datasets provide methods for defining the mesh and
writing point, cell, and field data. Composite writers manage collections
of simple VTK datasets.

.. list-table::
   :header-rows: 1

   * - VTK dataset type
     - File writer type
     - Module
   * - ImageData
     - ``vtkhdf_image_file``
     - ``vtkhdf_image_file_type``
   * - RectilinearGrid
     - ``vtkhdf_rectilinear_file``
     - ``vtkhdf_rectilinear_file_type``
   * - StructuredGrid
     - ``vtkhdf_structured_file``
     - ``vtkhdf_structured_file_type``
   * - UnstructuredGrid
     - ``vtkhdf_ug_file``
     - ``vtkhdf_ug_file_type``
   * - HyperTreeGrid
     - ``vtkhdf_htg_file``
     - ``vtkhdf_htg_file_type``
   * - MultiBlockDataSet
     - ``vtkhdf_mb_file``
     - ``vtkhdf_mb_file_type``

The ``vtkhdf_vtk_cell_types`` module additionally provides VTK cell type
integer codes relevant to UnstructuredGrid.

Data and Writer Modes
=====================

For simple VTK datasets, point and cell data associate scalar or multicomponent
values with mesh entities; field data associates values with the VTK dataset
as a whole. Their writers share conventions for static writes by name and
temporal registration followed by writes through handles. Field data uses a
common interface and array interpretation across simple VTK dataset types.
Structured writers additionally support rectangular region placement and
incremental point/cell writes.

Writers for simple VTK datasets support subsets of the shared modes
``VTKHDF_STATIC``, ``VTKHDF_FIXED_MESH``, ``VTKHDF_MOVING_MESH``, and
``VTKHDF_DYNAMIC_MESH``. The last three modes produce temporal VTK datasets.
See :doc:`common` for mode definitions and shared conventions, and
:doc:`field_data` for the field data interface.

Reader Compatibility
====================

This library targets VTKHDF Version 2.8.
Unless otherwise indicated, reader support targets ParaView 5.13 and later.
Older releases may work, but compatibility is not guaranteed. Additional
version requirements and known reader limitations are described in the
individual VTK dataset chapters.

Parallel and Serial I/O
=======================

In MPI builds, parallel HDF5 uses MPI-IO to write mesh and data collectively
to a single shared file. The library manages HDF5 operations and MPI
coordination internally. Writer methods are collective over the communicator
supplied at file creation; all participating ranks must make the same sequence
of calls. See :doc:`common` for the collective-call requirements.

A serial build uses a standard HDF5 installation without MPI. The API is the
same, except that file creation omits the communicator argument.

Data Placement
--------------

UnstructuredGrid stores each MPI rank's mesh contribution as a separate
partition, with its own points and cell connectivity. HyperTreeGrid likewise
stores each rank's trees as a separate partition. Each structured writer
assembles caller-owned rectangular regions into one global structured VTK
dataset without recording rank partitions. The library performs no
repartitioning or data redistribution.
