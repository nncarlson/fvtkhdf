[![CI](https://github.com/nncarlson/fvtkhdf/actions/workflows/ci.yml/badge.svg)](https://github.com/nncarlson/fvtkhdf/actions/workflows/ci.yml)
[![Docs](https://github.com/nncarlson/fvtkhdf/actions/workflows/docs.yml/badge.svg)](https://nncarlson.github.io/fvtkhdf/)

# fVTKHDF

A modern Fortran library for writing **VTKHDF** format files.

fVTKHDF provides a Fortran 2018 interface for writing HDF5-based VTKHDF files
readable by VTK and ParaView. It supports collective output to a single file
using parallel HDF5 and MPI, as well as serial builds without MPI.

fVTKHDF targets version 2.8 of the
[VTKHDF File Format Specification](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html)
for the dataset types it supports. VTK is used for reader validation; the
library itself depends on HDF5, with MPI required only for parallel builds.

## Dataset Support

* **ImageData, RectilinearGrid, and StructuredGrid:** Logically rectangular
  grids defined by uniform spacing, the tensor product of coordinate vectors,
  or explicit point coordinates, respectively.

  - Point and cell data are assembled from rectangular regions, with support
    for halos, multiple regions per rank, and empty collective contributions.
  - RectilinearGrid and StructuredGrid support moving geometry with fixed
    logical grids; ImageData geometry remains fixed.
  - Static point/cell data may be written incrementally. Temporal data supports
    complete replacement or reuse of the previous value.

* **UnstructuredGrid (UG):** Meshes with arbitrary cell connectivity, including
  moving geometry and changing topology. Each MPI rank supplies a separate
  mesh partition.

* **MultiBlockDataSet (MB):** A flat collection of UnstructuredGrid blocks.
  Each block selects its own writer mode, and all temporal blocks share a
  common timeline. Structured block types and hierarchical nesting are not
  supported.

All simple dataset writers support point, cell, and field data. They use the
shared mode names `VTKHDF_STATIC`, `VTKHDF_FIXED_MESH`, `VTKHDF_MOVING_MESH`, and
`VTKHDF_DYNAMIC_MESH`, with each writer supporting the applicable subset.
The original `UG_*` mode names remain available as compatibility aliases.

Field data uses a common interface for scalar and multicomponent numeric
values. Static field data is supported in every mode, including within MB
blocks. In temporal files, a static value applies to every step. Registered
temporal fields may change between steps or reuse their preceding values.

## Reader Compatibility

Unless otherwise indicated, reader support targets ParaView 5.13 and later.
Older releases may work, but compatibility is not guaranteed.

| Feature | Reader requirement |
| --- | --- |
| RectilinearGrid and StructuredGrid | ParaView 6.2+ (VTK 9.7+) |
| Dataset attribute designations (`Scalars`, `Vectors`, `Tensors`, `GlobalIds`, etc.) | ParaView 6.2+ (VTK 9.7+); older readers ignore these designations |

Readers earlier than VTK 9.7.2 (ParaView 6.2.1) require trailing singleton axes
for structured grids: use `(nx,ny,1)` for 2D and `(nx,1,1)` for 1D. For
**ImageData only**, these readers also require zero-based global logical point
lower bounds. A nonzero physical origin and nonzero local write bounds remain
supported.

VTK 9.7.1 also has a reader-cache defect affecting out-of-order time requests.
It is fixed in VTK 9.7.2 (ParaView 6.2.1); disable **Use Cache** with an affected
reader. See the [reference manual](https://nncarlson.github.io/fvtkhdf/) for
format-specific interfaces and reader limitations.

## Documentation

See the [Reference Manual](https://nncarlson.github.io/fvtkhdf/).

For local documentation work, install `sphinx`, `furo`, and `sphinx-autobuild`,
then run from the repository root:

```sh
sphinx-autobuild doc/sphinx build-doc/doc/sphinx/html
```

Open the local URL printed by the command. Changes to the documentation rebuild
and reload the preview automatically. No CMake configuration or compiler is
required. For a single build, use:

```sh
sphinx-build -b html doc/sphinx build-doc/doc/sphinx/html
```

CMake builds configured with `-DBUILD_HTML=ON` can also use
`cmake --build <build-dir> --target docs`.

The [developer guide](doc/developer.md) describes source organization,
ownership, runtime checks, and validation practices.

## Installation

### Prerequisites

* Fortran 2018 Compiler (GCC 13+, Intel oneAPI, etc.)
* CMake (3.28+)
* HDF5 C Library (1.10+) (Parallel HDF5 for MPI-enabled build)
* Python 3 + `fypp` (`pip install fypp`)
* MPI (Optional)

Limited `fpm` build instructions are available in [fpm/README.md](./fpm/README.md).

### Build and Install

By default, the library builds with MPI support enabled. Use
`-DENABLE_MPI=OFF` for a serial-only build. If your MPI installation does
not provide the `mpi_f08` module, use `-DUSE_MPI_F08=OFF` to disable the
`mpi_f08` communicator overloads while keeping the integer-communicator MPI
interface enabled.

```bash
git clone https://github.com/nncarlson/fvtkhdf
cd fvtkhdf
# Configure the build and set install location
cmake -B build --install-prefix /path/to/install
# Build and install
cmake --build build --parallel
cmake --install build
```

Run the tests for your chosen build configuration with:

```bash
ctest --test-dir build --output-on-failure
```

### Spack

fVTKHDF is available through Spack:

```bash
spack install fvtkhdf
```

The Spack package builds the MPI-enabled configuration and currently
targets parallel HDF5 1.14.x.

## Using fVTKHDF in your Project

Once installed, you can use **fVTKHDF** in your own CMake-based project by
adding the following to your `CMakeLists.txt`:

```CMake
find_package(fVTKHDF REQUIRED)

add_executable(my_simulation main.f90)
target_link_libraries(my_simulation PRIVATE fVTKHDF::fvtkhdf)
```

## Quick Start

Here is a minimal serial example that writes a **VTKHDF** UnstructuredGrid
dataset for an unstructured mesh consisting of a single tetrahedral cell
with a scalar field. More complete examples (serial and MPI-parallel,
UnstructuredGrid, ImageData, RectilinearGrid, StructuredGrid, and
MultiBlockDataSet) are provided in the [example directory](example/README.md).

```fortran
program quick_start

  use,intrinsic :: iso_fortran_env, only: int8
  use vtkhdf_ug_file_type
  use vtkhdf_vtk_cell_types, only: VTK_TETRA
  implicit none

  ! 1. Declare the file object (Unstructured Grid)
  type(vtkhdf_ug_file) :: myfile

  ! Data arrays
  real, allocatable :: points(:,:)
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  real, allocatable :: temperature(:)

  integer :: stat
  character(:), allocatable :: errmsg

  ! 2. Define a simple Tetrahedron (4 points, 1 cell)
  points = reshape([real :: 0,0,0,  1,0,0,  0,1,0,  0,0,1], [3,4])
  cnode = [1, 2, 3, 4] ! cell-node connectivity
  xcnode = [1, 5]      ! cnode start index of connectivity list
  types = [VTK_TETRA]
  temperature = [100.0, 200.0, 300.0, 400.0] ! Point data

  ! 3. Create and Write
  call myfile%create("simple.vtkhdf", stat, errmsg, mode=VTKHDF_STATIC)
  if (stat /= 0) error stop errmsg

  ! Write the mesh topology (Static Mesh)
  call myfile%write_mesh(points, cnode, xcnode, types)

  ! Write a data field
  call myfile%write_point_data("Temperature", temperature)

  call myfile%close()

end program quick_start
```

## Development Roadmap

Planned and prospective features include:

* Completing UnstructuredGrid (UG) feature coverage

  - Add support for [polyhedron cells](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/vtkhdf_specifications.html#polyhedron-support)

* Composite VTK datasets
  - Add structured block types to MultiBlockDataSet.
  - Replace the legacy MB with PartitionedDataSetCollection (PDC) ([issue #26](https://github.com/nncarlson/fvtkhdf/issues/26))
  - Implement a general `Assembly` tree for PDC to replace the existing flat tree.

### Additional Dataset Types

The [VTKHDF format specification](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html)
defines the format for additional VTK dataset types (e.g., HyperTreeGrid,
OverlappingAMR, etc.) that are not currently supported.

Contributions adding support for these formats are very welcome.
I am happy to collaborate on design and integration so that new
functionality fits naturally within the existing API and internal
infrastructure.

If you are interested in contributing support for one of these
dataset types, please open an issue to discuss design and implementation
details.

## License

**fVTKHDF** is distributed under the 2-clause BSD license.
See [LICENSE.md](./LICENSE.md) for details.
