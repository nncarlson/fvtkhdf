========
Overview
========

fVTKHDF is a modern Fortran library for writing
`VTKHDF <https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html>`_
files readable by ParaView. VTKHDF is an HDF5-based file format designed
for scalable scientific data storage and visualization.

.. note::

   In this manual, an unqualified *dataset* means an HDF5 dataset. A *VTK
   dataset* refers to a VTK data object, such as ImageData, UnstructuredGrid,
   or MultiBlockDataSet.

fVTKHDF is designed primarily for MPI-parallel high-performance computing (HPC)
applications. It uses parallel HDF5 (via MPI-IO) to enable all MPI ranks to
collectively write to a single shared file.

A serial build is also available, removing the MPI dependency while preserving
the same high-level API.

The VTKHDF format distinguishes *simple* VTK dataset types, such as ImageData
and UnstructuredGrid, from *composite* types, such as MultiBlockDataSet and
PartitionedDataSetCollection. Composite types are containers for simple VTK
datasets. Each simple VTK dataset has its own mesh and associated point, cell,
and field data.

fVTKHDF provides a file writer for each supported VTK dataset type. Each
writer is a Fortran derived type exported by its corresponding module.
Writers for simple VTK datasets provide methods for defining the mesh and
writing point, cell, and field data. Composite writers manage collections
of simple VTK datasets.

Writer modes use the shared names
``VTKHDF_STATIC``, ``VTKHDF_FIXED_MESH``, ``VTKHDF_MOVING_MESH``, and
``VTKHDF_DYNAMIC_MESH``. The last three modes produce temporal VTK datasets.
Each file-type module exports its supported subset: ImageData supports
the first two, while UnstructuredGrid supports all four.

Currently supported VTK dataset types:

* **UnstructuredGrid (UG)**: four writer modes
  (``VTKHDF_STATIC``, ``VTKHDF_FIXED_MESH``, ``VTKHDF_MOVING_MESH``,
  ``VTKHDF_DYNAMIC_MESH``) covering static, moving-mesh, and dynamic-mesh
  workflows for point, cell, and field data.

* **ImageData**: fixed regular geometry and rectangular point/cell block
  writes, with static data and temporal replacement or reuse. See
  :doc:`structured` for the data-writing interface and :doc:`image` for
  geometry and reader compatibility.

* **MultiBlockDataSet (MB)**: a composite type containing a flat collection
  of UnstructuredGrid blocks.

For simple VTK datasets, point and cell data associate scalar or multicomponent
values with mesh entities; field data associates values with the VTK dataset
as a whole. Their writers create static data by name and register temporal
data before time stepping, returning handles for subsequent writes. Names are
sanitized and disambiguated automatically. Field data uses a common interface and array
interpretation across simple VTK dataset types. Structured writers additionally
support rectangular block placement and incremental point/cell writes. See
:doc:`common` for shared conventions and :doc:`field_data` for the field data
interface.

The reference is organized into structured, unstructured, and composite VTK
dataset types. RectilinearGrid and StructuredGrid are planned additions to
the structured family; PartitionedDataSetCollection is a planned addition
to the composite family. These writers are not yet available.

This library targets VTKHDF Version 2.8.
See the individual VTK dataset chapters for reader requirements and verified
compatibility, including the ImageData reader limitations.

The available file writers are:

.. list-table::
   :header-rows: 1

   * - VTK dataset type
     - File writer type
     - Module
   * - UnstructuredGrid
     - ``vtkhdf_ug_file``
     - ``vtkhdf_ug_file_type``
   * - ImageData
     - ``vtkhdf_image_file``
     - ``vtkhdf_image_file_type``
   * - MultiBlockDataSet
     - ``vtkhdf_mb_file``
     - ``vtkhdf_mb_file_type``

The ``vtkhdf_vtk_cell_types`` module additionally provides VTK cell type
integer codes.

Parallel I/O and Data Model
===========================
The library uses Parallel HDF5 (via MPI-IO) for parallel output.
Although the VTKHDF specification does not require parallel writing,
fVTKHDF is designed for HPC use, allowing all MPI ranks to write
concurrently to a shared file.

Key Architectural Features
--------------------------
* **Single-File Output:** All simulation data, including geometry,
  connectivity, and temporal fields, is stored within a single `.vtkhdf` file.
* **Data placement:** UG and MB preserve mesh partitions. ImageData uses
  one global structured VTK dataset, assembled from caller-owned rectangular
  regions without recording rank partitions. The library performs no
  repartitioning or data redistribution.
* **Library Encapsulation:**  All parallel coordination and HDF5 logic is
  managed internally. The user interacts with a high-level API that requires
  only an MPI communicator.

Build Variants (Parallel vs. Serial)
====================================
The library can be configured for either parallel or serial environments:

* MPI Build: Requires MPI and a parallel-enabled HDF5 installation. The file
  creation method accepts an MPI communicator to coordinate collective I/O
  across the cluster.

* Serial Build: Requires only a standard HDF5 installation. The API remains
  identical, with the exception that the file creation method does not require
  a communicator argument.


Visualization Compatibility
===========================
For UG and MB, ParaView 5.13+ is the existing compatibility target. ImageData
reader limitations are described in :doc:`image`. By default, ParaView
treats internal partitions as independent mesh entities. Its filters therefore
process each partition independently.
This can lead to discontinuities in spatial filters (e.g., Contour, Slice) and
'Z-fighting' artifacts in the renderer if overlapping ghost cells are present.
Users can resolve these issues by including special VTK data arrays in the
VTKHDF file which allow ParaView to "stitch" the partitions into a single
logical mesh.
For more information on these arrays see the
:doc:`Parallel Stitching and Ghost Cells <stitching>` guide.
