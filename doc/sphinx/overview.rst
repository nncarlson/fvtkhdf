========
Overview
========

fVTKHDF is a modern Fortran library for writing
`VTKHDF <https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html>`_
files readable by ParaView. VTKHDF is an HDF5-based file format designed
for scalable scientific data storage and visualization.

fVTKHDF is designed primarily for MPI-parallel high-performance computing (HPC)
applications. It uses parallel HDF5 (via MPI-IO) to enable all MPI ranks to
collectively write to a single shared file.

A serial build is also available, removing the MPI dependency while preserving
the same high-level API.

Writer modes are defined in ``vtkhdf_file_modes`` using the shared names
``VTKHDF_STATIC``, ``VTKHDF_FIXED_MESH``, ``VTKHDF_MOVING_MESH``, and
``VTKHDF_DYNAMIC_MESH``. Each file-type module also exports its supported
subset: ImageData supports the first two, while UnstructuredGrid and the
current MultiBlock blocks support all four.

Currently supported VTKHDF dataset types:

* **UnstructuredGrid (UG)**: four writer modes
  (``VTKHDF_STATIC``, ``VTKHDF_FIXED_MESH``, ``VTKHDF_MOVING_MESH``,
  ``VTKHDF_DYNAMIC_MESH``) covering static, moving-mesh, and dynamic-mesh
  workflows for point, cell, and field data.

* **ImageData**: fixed regular geometry and rectangular point/cell block
  writes, with static data and temporal replacement or reuse. See
  :doc:`image` for reader compatibility and FieldData semantics.

* **MultiBlockDataSet (MB)**: a flat collection of UnstructuredGrid blocks
  with per-block UG writer modes and static or temporal point, cell, and
  field data.

Dataset writers follow the same interface conventions and data semantics as
UnstructuredGrid wherever their mesh representations permit. Field naming,
static versus temporal data, and FieldData behavior are shared conventions;
structured writers add rectangular block placement and incremental writes.

This library targets VTKHDF Version 2.8.
See the individual dataset chapters for reader requirements and verified
compatibility, including the ImageData reader limitations.

Public modules documented here:

* ``vtkhdf_ug_file_type`` — write UnstructuredGrid files
* ``vtkhdf_image_file_type`` — write ImageData files
* ``vtkhdf_mb_file_type`` — write MultiBlockDataSet files
* ``vtkhdf_vtk_cell_types`` — VTK cell type integer codes

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
  one global structured dataset, assembled from caller-owned rectangular
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
