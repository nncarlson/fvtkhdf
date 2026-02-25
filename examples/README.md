# Examples

There are parallel and serial example programs for writing UnstructuredGrid (UG)
and MultiBlockDataSet (MB) VTKHDF datasets. Use `-DBUILD_EXAMPLES=YES` when
running `cmake` to compile the examples.

From the `build/examples` directory:

**Parallel Build**
```bash
mpiexec -np 3 ./vtkhdf_ug_demo
paraview ug_demo.vtkhdf
```

```bash
mpiexec -np 3 ./vtkhdf_mb_demo
paraview mb_demo.vtkhdf
```

**Serial Build**
```bash
./vtkhdf_ug_demo
paraview ug_demo.vtkhdf
```

```bash
./vtkhdf_mb_demo
paraview mb_demo.vtkhdf
```
