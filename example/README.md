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
global logical bounds, explicit array placement, full temporal replacement,
omitted-write reuse, and independent rank-zero FieldData. See the ImageData
manual for multiple blocks, halo exclusion, empty contributions, and reader
limitations.
