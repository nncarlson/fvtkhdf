# fpm Usage

This repository includes limited `fpm` support as a secondary build path.

## Serial build

```bash
export PKG_CONFIG_PATH=</path/to/hdf5/lib/pkgconfig>
./generate_fpm_src.py --profile release
cd fpm-src
fpm build --profile release
```

## MPI build

```bash
export PKG_CONFIG_PATH=</path/to/parallel-hdf5/lib/pkgconfig>
./generate_fpm_src.py --features mpi-release
cd fpm-src
fpm build --features mpi-release
```

`PKG_CONFIG_PATH` must point to the `pkgconfig` directory for the HDF5 installation you want `fpm` to use. For a serial build, use the serial HDF5 installation. For an MPI build, use the parallel HDF5 installation that matches your MPI compiler wrappers.

Useful checks:

```bash
pkg-config --cflags --libs hdf5
which h5cc
which h5pcc
```
