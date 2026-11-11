Parallel Stitching and Ghost Cells
==================================

UnstructuredGrid files store each MPI rank's mesh as a separate partition.
To produce a seamless visualization and support correct calculations across
partition boundaries, additional point and cell arrays identify shared
entities and ghost copies. Without this information, rendering artifacts
and incorrect results can occur. VTK and ParaView filters may use these
arrays when processing the global mesh. Write them using the standard point
and cell data procedures.

* ``"vtkGhostType"``: an ``int8`` array that may be written for cells,
  points, or both, using the same name in their respective data groups.

  For cells, value 0 identifies a real cell owned by the local rank and
  value 1 a ghost (duplicate) cell. Marking duplicate cells allows supporting
  VTK filters to exclude them from rendering or calculations where appropriate.
  While not required, the inclusion of ghost cells in mesh output provides
  the context needed for filters to interpolate values across partition
  boundaries.

  For points, value 0 identifies the owned copy and value 1 a duplicate copy
  of a point shared between partitions. Shared boundary points can be
  duplicated even when no ghost cells are included. Global point IDs identify
  which copies represent the same point; the ghost flags mark their duplicate
  status.

* ``"vtkGlobalPointIds"``: an ``int64`` point-centered array providing a
  unique global index for every node.

* ``"vtkGlobalCellIds"``: an ``int64`` cell-centered array providing a
  unique global index for every cell.

Copies of the same point or cell in different partitions must have the same
global ID.

.. note::

   With VTK 9.7.0 / ParaView 6.2.0 and later, pass ``attribute='GlobalIds'``
   when writing static global point or cell IDs or registering temporal global
   IDs. The VTKHDF reader uses this designation to recognize the array as
   global IDs; the name alone is insufficient. Both point and cell IDs use
   the same attribute value, and their array names may be arbitrary.
   ``vtkGhostType`` continues to be recognized by name and does not require
   an ``attribute`` argument.
