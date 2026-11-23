# Examples

These programs demonstrate the ImageData, RectilinearGrid, StructuredGrid,
UnstructuredGrid, HyperTreeGrid, MultiBlockDataSet, and PartitionedDataSetCollection
writers. All support serial and MPI builds. Configure CMake with `-DBUILD_EXAMPLES=ON` to compile them.

From the `build/example` directory:

**Parallel Build**

```bash
mpiexec -np 3 ./unstructured_fixed_mesh
paraview example.vtkhdf
```

**Serial Build**

```bash
./unstructured_fixed_mesh
paraview example.vtkhdf
```

The UnstructuredGrid and MultiBlockDataSet examples write `example.vtkhdf`.
Each structured and HyperTreeGrid example writes a file named after its
executable, such as `structured_moving_mesh.vtkhdf`.

## Structured Dataset Types

`image_fixed_mesh.F90` writes a regular 2D grid with one disjoint cell slab
per MPI rank (or one slab in a serial build). It demonstrates `VTKHDF_FIXED_MESH`,
global logical bounds, static cell data, temporal point and field data, and
shared-point ownership using `write_bounds`. See the [ImageData manual](../doc/sphinx/image.rst) for multiple blocks, halo exclusion, empty contributions, and reader
limitations.

`rectilinear_fixed_mesh.F90` follows the ImageData example with nonuniform
coordinate vectors supplied by rank zero. The data and shared-point ownership
are identical; only the mesh geometry differs. It supports serial and MPI builds.

`structured_fixed_mesh.F90` follows the fixed-mesh examples, mapping the same
logical grid to an annular strip with uniform angular and radial increments.
Height increases with azimuth so successive turns do not overlap. Distributed
point coordinates and temporal point data use the same shared-point ownership;
static cell data and temporal field data match the other examples. It supports
serial and MPI builds.

`rectilinear_moving_mesh.F90` follows the fixed-mesh example, shrinking x and
expanding y at each successive step while keeping the logical grid fixed.
Static diffusivity and temporal point and field data are written as in the
fixed-mesh example. It supports serial and MPI builds.

`structured_moving_mesh.F90` follows the fixed-mesh example, rotating the strip
about z by one angular increment at each successive step while keeping each
point's height fixed. Static diffusivity and temporal point and field data are
written as in the fixed-mesh example. It supports serial and MPI builds.

## HyperTreeGrid

`hypertree_dynamic_mesh.F90` writes a 2D grid of quadtrees that are refined around a
circle moving across the grid. It demonstrates `VTKHDF_DYNAMIC_MESH`, the
breadth-first cell order of each tree, a new topology at every step, and
temporal cell and field data. It registers and writes a material interface
that approximates the circle by straight segments in boundary cells. The full
grid is retained on both sides of the interface, with no masked cells.
Root cells are dealt to the MPI ranks in turn.
See the [HyperTreeGrid manual](../doc/sphinx/hypertree.rst). It supports serial and
MPI builds.

To show the full grid with the interface overlaid, run the supplied ParaView
script from the directory containing the output:

```bash
pvpython --force-offscreen-rendering \
  /path/to/fvtkhdf/example/hypertree_dynamic_mesh_view.py hypertree_dynamic_mesh.vtkhdf
```

This command writes the state and exits; it does not start an interactive
viewer. The offscreen option prevents a temporary render window from opening.
In ParaView, use **File > Load State** to open the resulting
`hypertree_dynamic_mesh.pvsm`, then use the time controls to follow the circle.
The view colors the full grid by `level`, shows its cell edges, and overlays
the interface in red. Color by `distance` to see signed distance from the circle.

The script uses **Hyper Tree Grid To Unstructured Grid** for the full cells,
and **Hyper Tree Grid Geometry Filter** with **Fill Material** disabled for
the interface lines. Opening the data file directly uses the usual HTG
material surface, which clips boundary cells to one side of the interface.

## UnstructuredGrid and MultiBlockDataSet

`unstructured_static.F90`, `unstructured_fixed_mesh.F90`, `unstructured_moving_mesh.F90`, and
`unstructured_dynamic_mesh.F90` demonstrate the four UnstructuredGrid writer modes.
`mb_example1.F90` writes a collection of named UnstructuredGrid blocks with a
shared timeline, using typed block pointers for mesh and data writes. The MB
writer also accepts ImageData, RectilinearGrid, StructuredGrid, and HyperTreeGrid blocks;
the [composite dataset manual](../doc/sphinx/composite.rst) includes a heterogeneous example.
These examples write `example.vtkhdf`, replacing any file
with that name in the working directory.

## PartitionedDataSetCollection

`pdc_example1.F90` writes `pdc_example1.vtkhdf` with an UnstructuredGrid fluid
block and an ImageData solid block. Both use fixed-mesh mode and share a
three-step timeline. Each rank supplies one fluid quad and one solid cell.
The [composite dataset manual](../doc/sphinx/composite.rst) describes the shared composite API,
flat assembly, borrowed block pointers, and reader limitations.

The MB and PDC examples use the same composite API. Changing the imported
module and file type selects which composite representation VTK's
`vtkHDFReader` (also used by ParaView) should produce; the block and data-writing calls remain the same.

## Composite Assembly

`composite_assembly.F90` uses the shared assembly API to place fluid and solid
blocks under `Regions`, and reference the same fluid block under `Materials`.
It interleaves block and node creation, commits the hierarchy, and then writes
a time series through the original block pointers. It writes
`composite_assembly.vtkhdf` and supports serial and MPI builds. The example
selects PDC; substitute the MB module and file type for MB output.
