===================
The fVTKHDF Library
===================

fVTKHDF is a modern Fortran library for writing VTKHDF files readable by ParaView.
It targets MPI-parallel HPC workflows and uses parallel HDF5 for collective I/O.
A serial build is also available.

.. toctree::
   :maxdepth: 1
   :caption: Reference Manual

   overview
   common
   field_data

.. toctree::
   :maxdepth: 1
   :caption: Structured VTK Dataset Types

   Common Interface <structured>
   ImageData <image>

.. toctree::
   :maxdepth: 1
   :caption: Unstructured VTK Dataset Types

   UnstructuredGrid <ug>

.. toctree::
   :maxdepth: 1
   :caption: Composite VTK Dataset Types

   MultiBlockDataSet <mb>

.. toctree::
   :maxdepth: 1
   :caption: Reference and Guides

   vtk_cell_types
   stitching
