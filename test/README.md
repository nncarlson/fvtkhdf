# Test Notes

## Optional VTK Python Reader Checks

The UG and MB writer tests can optionally be followed by Python-based readback
checks using VTK's `vtkHDFReader`. These checks are enabled only when the
configured Python interpreter can import the `vtk` package.

A repo-local virtual environment is a convenient way to keep this optional
dependency isolated:

```bash
python3 -m venv .venv
.venv/bin/pip install vtk
```

Configure CMake to use that interpreter:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DENABLE_MPI=OFF \
  -DPython3_EXECUTABLE=$PWD/.venv/bin/python
```

If `vtk` is not available, the build still succeeds and the reader checks are
skipped.
