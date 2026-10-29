[![CI](https://github.com/nncarlson/fvtkhdf/actions/workflows/ci.yml/badge.svg)](https://github.com/nncarlson/fvtkhdf/actions/workflows/ci.yml)
[![Docs](https://github.com/nncarlson/fvtkhdf/actions/workflows/docs.yml/badge.svg)](https://nncarlson.github.io/fvtkhdf/)

# fVTKHDF

A modern Fortran library for writing **VTKHDF** format files.

This library provides a high-level, object-oriented Fortran interface for
generating **VTKHDF** files, a relatively new and evolving HDF5-based VTK
file format used by ParaView. The library is designed for high-performance
computing (HPC) applications, offering robust MPI-parallel output via HDF5,
while also supporting serial workflows. By utilizing HDF5 as the underlying
storage mechanism, **fVTKHDF** provides a more scalable alternative to older
VTK ASCII or XML formats.

fVTKHDF targets version 2.8 of the
[VTKHDF File Format Specification](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html)
for the dataset types it supports.

**ParaView Compatibility**

**ParaView 6.2 or later is recommended. ParaView 5.13 is the oldest
compatibility target.** Older ParaView releases support the existing dataset
types but may ignore newer metadata.

| Output feature | ParaView support |
| --- | --- |
| UnstructuredGrid and MultiBlock, including temporal FieldData | 5.13+ expected |
| ImageData | 5.13+ expected, with the reader limitations below |
| Dataset attribute designations (`Scalars`, `Vectors`, etc.) | 6.2+ |
| RectilinearGrid and StructuredGrid (planned) | Will require 6.2+ |

**ImageData reader limitations:** ParaView 6.2.0 includes VTK 9.7.1 and still
has two geometry reader defects. For affected readers, use **zero-based global
point bounds** on every axis and trailing singleton axes: `(nx,ny,1)` for 2D
or `(nx,1,1)` for 1D. Positive local block indices and a nonzero physical origin
are supported. Both defects are fixed in the
[tested upstream development source](https://github.com/Kitware/VTK/commit/9f0588e2fe0d10573555834a66ac3dab95be4114)
and the fixes have been backported and queued for VTK 9.7.2. The corresponding
ParaView release has not yet been identified. fVTKHDF preserves the supplied
global bounds and physical origin. See the [ImageData manual](doc/sphinx/image.rst)
for details.

**Dataset Support**

**fVTKHDF** currently supports the following VTK data models:

* **UnstructuredGrid (UG):** For meshes with arbitrary cell types
  (tetrahedra, hexahedra, etc.).

  * Supports a **static** workflow and three **transient** workflows:
    fixed, moving, and dynamic mesh.

  * Supports point-centered, cell-centered, and field data.

* **ImageData:** Fixed regular-grid geometry with static and temporal
  point/cell arrays assembled from rectangular blocks, including halo
  exclusion, multiple blocks per rank, and empty collective contributions.

  * Uses explicit three-dimensional spatial array shapes and component-first
    multicomponent arrays.
  * Temporal fields support complete replacement or whole-field reuse.
  * FieldData retains the existing UG semantics.

* **MultiBlockDataSet (MB):** Supports a flat assembly/collection of
  UnstructuredGrid blocks.

  * Mixed-Mode Support: Each block can independently utilize any of the
    UG workflows within the same file.

  * All transient blocks share a common set of time steps.

  * Useful for a logical decomposition of a model into distinct components
    (e.g., piston, cylinder, valves), independent of any parallel domain
    decomposition.

  * Hierarchical assembly of MultiBlock datasets is not supported.


## Documentation
See the [Reference Manual](https://nncarlson.github.io/fvtkhdf/).


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
UnstructuredGrid, ImageData, and MultiBlockDataSet) are provided in the `example`
directory.

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
