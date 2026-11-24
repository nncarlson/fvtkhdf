==================
Installation Guide
==================

Prerequisites
=============

* A Fortran 2018 compiler and companion C compiler.
* CMake 3.28 or newer.
* The HDF5 C library, version 1.10 or newer, including development headers.
  The HDF5 Fortran library is not required. MPI builds require parallel HDF5.
* Python 3 and ``fypp`` (``python3 -m pip install fypp``), with the ``fypp``
  executable on your ``PATH``.
* MPI for parallel builds, with C and Fortran development interfaces.

Tested compiler families include GNU Fortran (GCC), Intel oneAPI, NAG,
LLVM Flang, and LFortran.

The ``h5dump`` tool, provided by HDF5, and the Python packages ``h5py``,
``numpy``, and ``vtk`` are **strongly recommended** for testing. They verify
generated file layouts, stored data, and VTK reader compatibility. Without
them, the Fortran writer and API-contract tests still run, but CTest reports
most file-layout and reader-validation checks as skipped, substantially
reducing coverage of output correctness and interoperability.

Before configuring, select the intended compiler and HDF5 environment, along
with a compatible MPI installation for parallel builds. Use the same Fortran
compiler for fVTKHDF and your application. For parallel builds, MPI's Fortran
interface must also have been built with that compiler.

Configure, build, test, and install
===================================

Get the source:

.. code-block:: console

   git clone https://github.com/nncarlson/fvtkhdf
   cd fvtkhdf

The default branch contains the latest development sources. To build a
release, check out its release tag before configuring.

Choose one of the following configurations. Both enable the supplied examples
and set the installation prefix; replace ``/path/to/install`` with your desired
location. Separate build directories keep serial and MPI configurations apart.
CMake defaults to a Release build; add ``-DCMAKE_BUILD_TYPE=Debug`` to enable
additional runtime diagnostics.

To select compilers explicitly, set the ``FC`` and ``CC`` environment variables
before the first CMake configuration. Use a new build directory when changing
compilers.

.. _serial-build:

Serial build
------------

Disable MPI explicitly. This configuration works with serial HDF5:

.. code-block:: console

   cmake -B build-serial -DENABLE_MPI=OFF -DBUILD_EXAMPLES=ON \
     --install-prefix /path/to/install
   cmake --build build-serial --parallel
   ctest --test-dir build-serial --output-on-failure
   cmake --install build-serial

MPI build (default)
-------------------

MPI support is enabled by default. Parallel HDF5 built with a compatible MPI
installation is required; configuration fails if only serial HDF5 is found.

.. code-block:: console

   cmake -B build-mpi -DENABLE_MPI=ON -DBUILD_EXAMPLES=ON \
     --install-prefix /path/to/install
   cmake --build build-mpi --parallel
   ctest --test-dir build-mpi --output-on-failure
   cmake --install build-mpi

CTest launches the parallel tests through the MPI launcher, normally with four
ranks. Run it in an environment that permits MPI jobs with that many ranks.

The ``USE_MPI_F08`` option is enabled by default. If your MPI installation does
not provide the ``mpi_f08`` module, add ``-DUSE_MPI_F08=OFF`` to the MPI configure
command. This disables the ``mpi_f08`` communicator overloads while retaining
the integer-communicator MPI interface. Serial builds disable these overloads
automatically.

Dependency discovery
--------------------

If CMake cannot find the intended HDF5 installation or selects the wrong one,
add ``-DHDF5_ROOT=/path/to/hdf5`` to the configure command. Similarly, use
``-DPython3_EXECUTABLE=/path/to/python3`` only if CMake cannot find the intended
Python interpreter or selects one without the required packages. These are
optional overrides; normal automatic discovery needs neither argument.

Building and running examples
=============================

The configurations above set ``BUILD_EXAMPLES=ON``. To enable examples in an
existing build, reconfigure it with that option and rebuild. Executables are
placed in the build's ``example/`` directory.

From the repository root, run the static UnstructuredGrid example using the
configuration you built:

.. code-block:: console

   # Serial configuration
   cd build-serial/example
   ./unstructured_static

For the MPI configuration, instead use:

.. code-block:: console

   cd build-mpi/example
   mpiexec -n 3 ./unstructured_static

Each run writes ``example.vtkhdf`` in its working directory, replacing an
existing file with that name. Open it in ParaView with
``paraview example.vtkhdf``. See the :repo:`example directory <example/README.md>`
for the full set of serial and MPI examples.

Using fVTKHDF in your project
=============================

The installation includes Fortran module files, the library, and CMake package
configuration files. In your application's ``CMakeLists.txt``, enable both C
and Fortran so the package can discover its HDF5 C dependency:

.. code-block:: cmake

   cmake_minimum_required(VERSION 3.28)
   project(my_simulation LANGUAGES C Fortran)

   find_package(fVTKHDF REQUIRED)
   add_executable(my_simulation main.f90)
   target_link_libraries(my_simulation PRIVATE fVTKHDF::fvtkhdf)

Configure and build your application with the installation prefix:

.. code-block:: console

   cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/install
   cmake --build build --parallel

Keep the compatible compiler, HDF5, and MPI environment selected while building
your application. The imported target carries the library's dependencies. To
require an MPI-enabled installation explicitly, use
``find_package(fVTKHDF REQUIRED COMPONENTS MPI)``.

Other installation routes
=========================

Fortran Package Manager (fpm)
-----------------------------

See the :repo:`fpm instructions <fpm/README.md>` for the secondary fpm build
path, including serial and MPI configurations and compiler workarounds.

Spack
-----

fVTKHDF is also available through Spack:

.. code-block:: console

   spack install fvtkhdf

Spack installs a packaged release, which may differ from the development
version described in this guide.

The Spack package builds the MPI-enabled configuration and currently targets
parallel HDF5 1.14.x.
