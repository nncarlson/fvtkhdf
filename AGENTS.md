# Repository Guidelines

## Project Structure & Module Organization
`src/` contains the library sources. Core modules live in `.F90` files, while repeated interfaces are generated from `.F90.fypp` templates through `cmake/Fypp.cmake`. `test/` holds CTest-driven regression programs plus `test/packaging/` for the installed-package consumer check. `examples/` contains serial and MPI demo programs, and `doc/sphinx/` builds the HTML reference manual. Treat `build-gcc/` and `build-nag/` as local build artifacts, not source.

## Build, Test, and Development Commands
Load the default MPI toolchain first with `ml gcc mpich`. Configure a default build with `cmake -B build -DCMAKE_BUILD_TYPE=Debug -DENABLE_MPI=YES`. Use `-DENABLE_MPI=NO` for serial work; enable docs or examples with `-DBUILD_HTML=ON` or `-DBUILD_EXAMPLES=ON`. Build with `cmake --build build --parallel`. Run serial tests with `ctest --test-dir build --output-on-failure`. MPI-parallel tests also use `ctest --test-dir build --output-on-failure`, but must be launched in escalated mode so the 4-rank test executables can run. Install locally with `cmake --install build --prefix /tmp/fvtkhdf-stage`. Validate the installed package with `cmake -S test/packaging -B build-consumer -DCMAKE_PREFIX_PATH=/tmp/fvtkhdf-stage -DENABLE_MPI=YES`.

## Coding Style & Naming Conventions
Follow the existing Fortran style: two-space indentation, `implicit none`, lowercase keywords, and compact `use,intrinsic :: ...` imports. Keep module and file names aligned, for example `vtkhdf_ctx_type` in `src/vtkhdf_ctx_type.F90`. Use descriptive derived-type and procedure names with `vtkhdf_` prefixes. Edit fypp templates in `src/*.F90.fypp` when changing generated overload families; do not hand-edit generated files in build directories. No formatter is configured, so match nearby code exactly.

## Testing Guidelines
Add or update CTest executables in `test/CMakeLists.txt` when behavior changes. MPI-enabled tests run under 4 ranks; serial coverage uses the `*-serial.F90` variants. Keep test names tied to the module or behavior under test, such as `vtkhdf_ug_test.F90` or `test_leak.F90`. At minimum, run `ctest --test-dir build --output-on-failure` for the MPI or serial mode you changed.

## Commit & Pull Request Guidelines
Recent commits use short, imperative subjects such as `Add examples` or `Redo MB block data structure`. Keep the first line focused on the user-visible change and avoid mixing unrelated refactors. Pull requests should state whether the change targets MPI, serial, packaging, or docs, list the exact CMake and `ctest` commands run, and link any relevant issue. Include screenshots only for documentation or rendered output changes.
