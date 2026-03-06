Examples
========

A minimal usage pattern for writing an UnstructuredGrid file is:

.. code-block:: fortran

   use vtkhdf_ug_file_type
   use vtkhdf_temporal_level
   type(vtkhdf_ug_file) :: file
   integer :: stat
   character(:), allocatable :: errmsg

   call file%create("mesh.vtkhdf", comm, stat, errmsg, temporal_level=VTKHDF_STATIC_MESH)
   call file%write_mesh(points, cnode, xcnode, types)
   call file%write_cell_data("pressure", pressure)
   call file%close()

This illustrates the required ordering of calls: create the file,
write the mesh, write any static or temporal datasets, and finally close.

For deformed meshes (varying geometry but fixed topology), use
``temporal_level=VTKHDF_DEFORMED_MESH`` and call
``write_mesh_topology`` followed by ``write_temporal_points`` at each
time step:

.. code-block:: fortran

   use vtkhdf_ug_file_type
   use vtkhdf_temporal_level
   type(vtkhdf_ug_file) :: file
   integer :: stat
   character(:), allocatable :: errmsg

   call file%create("deformed.vtkhdf", comm, stat, errmsg, temporal_level=VTKHDF_DEFORMED_MESH)
   call file%write_mesh_topology(points, cnode, xcnode, types)
   call file%register_temporal_cell_data("pressure", mold)
   call file%write_time_step(0.0)
   call file%write_temporal_points(points)
   call file%write_temporal_cell_data("pressure", pressure)
   call file%close()

More complete examples (serial and MPI-parallel, UnstructuredGrid and
MultiBlockDataSet) are provided in the ``examples`` directory of
the project repository.
