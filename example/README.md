# Examples

Here are some example programs for writing UnstructuredGrid (UG) and
MultiBlockDataSet (MB) VTKHDF datasets. Use `-DBUILD_EXAMPLES=YES` when
running `cmake` to compile the examples.

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

All examples produce an `example.vtkhdf` output file.

## ImageData

`image_fixed_mesh.F90` writes a regular 2D grid with one disjoint cell slab
per MPI rank (or one slab in a serial build). It demonstrates `VTKHDF_FIXED_MESH`,
global logical bounds, static cell data, temporal point and field data, and
shared-point ownership using `write_bounds`. See the ImageData
manual for multiple blocks, halo exclusion, empty contributions, and reader
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
