# fpm Usage

This repository provides fpm support as a secondary build path. The regression
test suite requires CMake/CTest.

fVTKHDF generates some of its Fortran sources from Fypp templates. Because fpm
does not natively support Fypp preprocessing, run `generate_fpm_src.py` before
building with fpm. The script expands the `.F90.fypp` templates into `.F90`
files and assembles an fpm project in `fpm-src/`, including the remaining
sources, examples, and manifest. This step requires Python 3 and the `fypp`
executable on your `PATH`.

Use matching serial or MPI options for generation and the subsequent fpm build,
as shown below. Regenerate after changing the library sources or switching
between serial and MPI builds. Make source changes in the original files;
the script recreates the generated tree each time.

## Serial build

```bash
cd fpm
export PKG_CONFIG_PATH=/path/to/your/hdf5/lib/pkgconfig
./generate_fpm_src.py --profile release
cd fpm-src
fpm build --profile release
```

## MPI build

```bash
cd fpm
export PKG_CONFIG_PATH=/path/to/your/parallel-hdf5/lib/pkgconfig
./generate_fpm_src.py --features mpi-release
cd fpm-src
fpm build --features mpi-release
```

For MPI builds using **HDF5 older than 1.13.0**, add the preprocessor define
`HDF5_NULL_MEMSPACE_BUG` to work around collective-write corruption involving
null memory dataspaces:

```bash
fpm build --features mpi-release --flag "-DHDF5_NULL_MEMSPACE_BUG"
```

The same flag is required with `mpi-nof08-release`. Serial builds and builds
using HDF5 1.13.0 or newer do not need it. CMake builds select the workaround
automatically.

`PKG_CONFIG_PATH` must point to the `pkgconfig` directory for the HDF5 installation you want `fpm` to use. For a serial build, use the serial HDF5 installation. For an MPI build, use the parallel HDF5 installation that matches your MPI compiler wrappers.

Useful checks:

```bash
pkg-config --cflags --libs hdf5
which h5cc
which h5pcc
```

## GNU compiler workaround

All GNU compiler features in `fpm.toml` define `GNU_BUG20261103` to work around
incorrect handling of noncontiguous arrays passed through `select rank`.
With **GNU Fortran 16.1.0 or newer**, this define can be removed from the
`gfortran.preprocess.cpp.macros` lists. Leaving it enabled is safe.
CMake builds select the workaround automatically based on the compiler version.
