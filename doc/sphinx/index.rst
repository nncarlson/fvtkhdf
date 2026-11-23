===================
The fVTKHDF Library
===================

fVTKHDF is a modern Fortran library for writing VTKHDF files readable by VTK
and ParaView.
It targets MPI-parallel HPC workflows and uses parallel HDF5 for collective I/O.
A serial build is also available.

* :doc:`installation`: prerequisites, serial and MPI builds, examples, and
  integration with your project.
* :doc:`reference`: writer interfaces, dataset types, and reader compatibility.

.. toctree::
   :hidden:
   :maxdepth: 1

   installation
   reference
