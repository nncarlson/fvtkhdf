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

Currently supported VTKHDF dataset types:

* **UnstructuredGrid (UG)**  
  - Static mesh  
  - Static or time-dependent point, cell, and field data

* **MultiBlockDataSet (MB)**  
  - Static mesh  
  - Flat collection of UnstructuredGrid blocks  
  - Static or time-dependent point, cell, and field data per block

This library targets VTKHDF Version 2.5.
While basic files may open in ParaView 5.10, ParaView 5.13 or newer is
recommended for full feature support.

Public modules documented here:

* ``vtkhdf_ug_file_type`` — write UnstructuredGrid files
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
* **Preserved Partitioning:** The file structure records data according to
  the application's native MPI partitioning. Each MPI rank contributes
  exactly one VTKHDF partition. No repartitioning or data redistribution
  is performed by the library.
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
The generated files are compatible with ParaView 5.13+. By default, ParaView
treats internal partitions as independent mesh entities. Its filters therefore
process each partition independently.
This can lead to discontinuities in spatial filters (e.g., Contour, Slice) and
'Z-fighting' artifacts in the renderer if overlapping ghost cells are present.
Users can resolve these issues by including special VTK data arrays in the
VTKHDF file which allow ParaView to "stitch" the partitions into a single
logical mesh.
For more information on these arrays see the
:doc:`Parallel Stitching and Ghost Cells <stitching>` guide.
