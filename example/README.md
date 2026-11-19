# Examples

These programs demonstrate the ImageData, RectilinearGrid, StructuredGrid,
UnstructuredGrid, HyperTreeGrid, and MultiBlockDataSet writers. All support serial and MPI
builds. Configure CMake with `-DBUILD_EXAMPLES=ON` to compile them.

From the `build/example` directory:

**Parallel Build**

```bash
mpiexec -np 3 ./ug_fixed_mesh
paraview example.vtkhdf
```

**Serial Build**

```bash
./ug_fixed_mesh
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

`htg_dynamic_mesh.F90` writes a 2D grid of quadtrees that are refined around a
circle moving across the grid. It demonstrates `VTKHDF_DYNAMIC_MESH`, the
breadth-first cell order of each tree, a new topology at every step, and
temporal cell and field data. Root cells are dealt to the MPI ranks in turn.
See the [HyperTreeGrid manual](../doc/sphinx/htg.rst). It supports serial and
MPI builds.

## UnstructuredGrid and MultiBlockDataSet

`ug_static.F90`, `ug_fixed_mesh.F90`, `ug_moving_mesh.F90`, and
`ug_dynamic_mesh.F90` demonstrate the four UnstructuredGrid writer modes.
`mb_example1.F90` writes a collection of named UnstructuredGrid blocks with a
shared timeline. These examples write `example.vtkhdf`, replacing any file
with that name in the working directory.
