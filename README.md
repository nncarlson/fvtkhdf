[![CI](https://github.com/nncarlson/fvtkhdf/actions/workflows/ci.yml/badge.svg)](https://github.com/nncarlson/fvtkhdf/actions/workflows/ci.yml)
[![Docs on Read the Docs](.github/docs-badge.svg)](https://fvtkhdf.readthedocs.io/)

# fVTKHDF

A modern Fortran library for writing **VTKHDF** format files.

fVTKHDF provides a Fortran 2018 interface for writing HDF5-based VTKHDF files
readable by VTK and ParaView. It supports collective output to a single file
using parallel HDF5 and MPI, as well as serial builds without MPI.

fVTKHDF targets version 2.8 of the
[VTKHDF File Format Specification](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html)
for the dataset types it supports. The library depends on HDF5, with MPI
required only for parallel builds.

## Documentation

* [Installation Guide](https://fvtkhdf.readthedocs.io/page/installation.html):
  prerequisites, building, testing, and project integration.
* [Reference Manual](https://fvtkhdf.readthedocs.io/page/reference.html):
  writer interfaces, dataset types, and reader compatibility.

## Dataset Support

* **ImageData, RectilinearGrid, and StructuredGrid:** Logically rectangular
  grids defined by uniform spacing, coordinate vectors, or explicit point
  coordinates, respectively. RectilinearGrid and StructuredGrid support moving
  geometry; ImageData geometry remains fixed.

* **UnstructuredGrid:** Meshes with arbitrary cell connectivity, including
  moving geometry and changing topology. Each MPI rank supplies a separate
  mesh partition.

* **HyperTreeGrid:** Tree-based adaptive mesh refinement on a rectilinear
  grid of root cells, including moving geometry, changing refinement, and
  masked cells. Each MPI rank supplies a separate partition of the trees.

* **MultiBlockDataSet and PartitionedDataSetCollection:** A flat collection of
  ImageData, RectilinearGrid, StructuredGrid, UnstructuredGrid, and HyperTreeGrid
  blocks in one file. Typed block pointers
  provide each simple writer's mesh and data methods. Each block selects its
  own mode, and all temporal blocks share a common timeline. Blocks form a flat
  collection by default; assembly construction organizes them hierarchically
  and allows multiple references to a block. The two file types share one composite API; choose
  the file type according to the composite representation VTK's `vtkHDFReader`
  (also used by ParaView) should produce. See the [composite dataset manual](doc/sphinx/composite.rst).

All simple dataset writers support static and time-dependent cell and field
data. All except HyperTreeGrid also support point data. Depending on the dataset
type, the mesh can remain fixed, move, or change topology between time steps.

## Reader Compatibility

Unless otherwise indicated, reader support targets ParaView 5.13 and later.
Older releases may work, but compatibility is not guaranteed. For best results,
use the latest stable release of ParaView.

| Feature | Reader requirement |
| --- | --- |
| RectilinearGrid and StructuredGrid | ParaView 6.2+ (VTK 9.7+) |
| HyperTreeGrid | ParaView 6.2+ (VTK 9.7+) |
| Dataset attribute designations (`Scalars`, `Vectors`, `Tensors`, `GlobalIds`, etc.) | ParaView 6.2+ (VTK 9.7+); older readers ignore these designations |

See the [reference manual](https://fvtkhdf.readthedocs.io/page/reference.html) for
format-specific reader limitations and workarounds.

## Installation

Select a Fortran 2018 compiler, CMake 3.28+, the HDF5 C library, and Python 3
with `fypp`. The default build enables MPI and requires parallel HDF5 built
with a compatible MPI installation.

```bash
git clone https://github.com/nncarlson/fvtkhdf
cd fvtkhdf
cmake -B build --install-prefix /path/to/install
cmake --build build --parallel
ctest --test-dir build --output-on-failure
cmake --install build
```

Use `-DENABLE_MPI=OFF` for a serial build. See the
[Installation Guide](https://fvtkhdf.readthedocs.io/page/installation.html)
for prerequisites, MPI options, examples, fpm, Spack, and
[project integration](https://fvtkhdf.readthedocs.io/page/installation.html#using-fvtkhdf-in-your-project).

## Quick Start

This minimal example requires a serial build configured with `-DENABLE_MPI=OFF`.
It writes a **VTKHDF** UnstructuredGrid dataset containing a single tetrahedral
cell with scalar point data. For MPI usage, see [unstructured_static.F90](example/unstructured_static.F90).
More complete examples for all supported dataset types are provided in the
[example directory](example/README.md).

```fortran
program quick_start

  use,intrinsic :: iso_fortran_env, only: int8
  use vtkhdf_unstructured_file_type
  use vtkhdf_vtk_cell_types, only: VTK_TETRA
  implicit none

  ! 1. Declare the file object (Unstructured Grid)
  type(vtkhdf_unstructured_file) :: myfile

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

* Completing UnstructuredGrid feature coverage

  - Add support for [polyhedron cells](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/vtkhdf_specifications.html#polyhedron-support)

### Additional Dataset Types

The [VTKHDF format specification](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html)
defines the format for additional VTK dataset types, such as OverlappingAMR,
that are not currently supported.

Contributions supporting additional dataset types are welcome. Please open an
issue to discuss design and integration before starting an implementation.

## Questions and Bug Reports

Please use [GitHub Issues](https://github.com/nncarlson/fvtkhdf/issues) for bug
reports, feature requests, and usage questions. For bug reports, include a small
reproducing example and your compiler, HDF5, and MPI versions, as applicable.

## License

**fVTKHDF** is distributed under the 2-clause BSD license.
See [LICENSE.md](./LICENSE.md) for details.
See [CONTRIBUTORS.md](./CONTRIBUTORS.md) for contributor credits.
