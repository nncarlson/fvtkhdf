Parallel Stitching and Ghost Cells
==================================

When visualizing partitioned unstructured grids, ParaView treats each
VTKHDF partition as a separate piece of geometry. To produce a seamless
visualization and ensure mathematical correctness across the global domain,
specifically named mesh-centered arrays may be written to the file.
When ParaView encounters these arrays, it interprets them as metadata
rather than ordinary point or cell data. These arrays are written using
the standard point and cell data procedures.

* ``"vtkGhostType"``: an ``int8`` array that may be written for cells,
  points, or both, using the same name in their respective data groups.

  For cells, value 0 identifies a real cell owned by the local rank and
  value 1 a ghost (duplicate) cell. This prevents duplicate cells from being
  rendered as overlapping geometry. While not required, the inclusion of
  ghost cells in mesh output provides the context needed for filters to
  interpolate values across partition boundaries.

  For points, value 0 identifies the owned copy and value 1 a duplicate copy
  of a point shared between partitions. Shared boundary points can be
  duplicated even when no ghost cells are included. Global point IDs identify
  which copies represent the same point; the ghost flags mark their duplicate
  status.

* ``"vtkGlobalPointIds"``: an ``int64`` point-centered array providing a
  unique global index for every node. This prevents visual "cracks" in
  filters such as Contour or Slice.

* ``"vtkGlobalCellIds"``: an ``int64`` cell-centered array providing a
  unique global index for every cell. This is required for accurate data
  integration (e.g., the Integrate Variables filter), ensuring that ghost
  or duplicate cells are not double-counted.

.. note::

   With VTK 9.7.0 / ParaView 6.2.0 and later, pass ``attribute='GlobalIds'``
   when writing static global point or cell IDs or registering temporal global
   IDs. The VTKHDF reader uses this designation to recognize the array as
   global IDs; the name alone is insufficient. Both point and cell IDs use
   the same attribute value, and their array names may be arbitrary.
   ``vtkGhostType`` continues to be recognized by name and does not require
   an ``attribute`` argument.
