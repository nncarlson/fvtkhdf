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
