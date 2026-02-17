## fVTKHDF

fVTKHDF is a modern Fortran library for writing files in the VTKHDF format,
designed for high-performance scientific simulations. VTKHDF is a relatively
new VTK file format designed for high-performance and parallel data storage.
It utilizes the HDF5 standard as its underlying storage mechanism, providing
a more scalable alternative to older VTK ASCII or XML formats.

### Documentation

### Compiling

### Installation

### Testing

### License
fVTKHDF distributed under the 2-clause BSD license.
See [LICENSE.md](./LICENSE.md) for details.

# fVTKHDF

A modern Fortran library for writing **VTKHDF** format files. 

This library provides a high-level, object-oriented Fortran interface to
generating VTKHDF files, a standard widely supported by visualization tools
like ParaView (5.10+) and VisIt. It is designed for High-Performance Computing
(HPC) applications, offering robust MPI-parallel output via HDF5, while fully
supporting serial workflows.

## Features

* Targets version **2.5** of the [VTKHDF File Format Specification](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html).
* **Dataset Support:**
    * **UnstructuredGrid (UG):** Full support for arbitrary cell types and topologies.
    * **MultiBlockDataSet (MB):** Supports flat collections of UnstructuredGrid blocks.
* **Parallel & Serial:** * Built on top of **HDF5**, allowing for efficient, collective MPI-parallel
    I/O.
    * Can be built as a purely serial library (no MPI dependency required) for simpler workflows.
* **Time-Dependent Data:** Supports writing temporal data (transient fields on a static mesh) using the VTKHDF `Steps` group mechanism.
* **Type Polymorphism:**
    * **Coordinates/Data:** Supports both `real32` and `real64`.
    * **Connectivity/Ids:** Supports `int32` and `int64`.
    * **Mesh Size:** Supports both "small" (standard integer) and "large" (64-bit integer) mesh addressing.

## Installation

### Prerequisites
* Fortran Compiler (GCC 13+, Intel oneAPI, etc.)
* CMake (3.30+)
* HDF5 Library (1.10+)
* Python 3 + `fypp` (`pip install fypp`)
* MPI (Optional)

### Build Instructions

**Standard Parallel Build:**
```bash
git clone [https://github.com/nncarlson/fvtkhdf](https://github.com/nncarlson/fvtkhdf)
cd fvtkhdf
cmake -B build -DENABLE_MPI=ON
cmake --build build --parallel

**Standard Serial Build:**
```bash
cmake -B build -DENABLE_MPI=OFF
cmake --build build
```

## Quck Start

Here is a minimal serial example that writes a VTKHDF UstructuredGrid dataset
for an unstructured mesh consisting of a single tetrahedral cell with a scalar
field.

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
  integer, allocatable :: connectivity(:), cell_ptr(:)
  integer(int8), allocatable :: types(:)
  real, allocatable :: temperature(:)

  integer :: stat
  character(:), allocatable :: errmsg

  ! 2. Define a simple Tetrahedron (4 points, 1 cell)
  points = reshape([0,0,0,  1,0,0,  0,1,0,  0,0,1], shape=[3,4])
  connectivity = [1, 2, 3, 4] ! Point IDs
  cell_ptr = [1, 5]           ! starting index of each cell in connectivity
  types = [VTK_TETRA]
  temperature = [100.0, 200.0, 300.0, 400.0] ! Point data

  ! 3. Create and Write
  call myfile%create("simple.vtkhdf", stat, errmsg)
  if (stat /= 0) error stop errmsg
  
  ! Write the mesh topology (Static Mesh)
  call myfile%write_mesh(points, connectivity, offsets, types)
  
  ! Write a data field
  call myfile%write_point_data("Temperature", temperature)

  call myfile%close()

end program quick_start

```

 simple Unstructured Grid (a single tetrahedon) with a scalar field.
