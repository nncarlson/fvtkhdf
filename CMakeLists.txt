cmake_minimum_required(VERSION 3.30)

project(fVTKHDF Fortran C)

option(ENABLE_MPI "Enable MPI support" ON)
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Default to Release build type on single-configurator generators
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
endif()

# Disable assertions unless build type is Debug
add_compile_definitions("$<$<NOT:$<CONFIG:Debug>>:NDEBUG>")

if(ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS C Fortran)
  set(HDF5_PREFER_PARALLEL True)
  add_compile_definitions("$<$<COMPILE_LANGUAGE:C>:USE_MPI>")
else()
  set(HDF5_PREFER_PARALLEL False)
endif()

find_package(HDF5 REQUIRED COMPONENTS C HL)

if (ENABLE_MPI AND NOT HDF5_IS_PARALLEL)
  message(FATAL_ERROR "Parallel HDF5 is required when MPI is enabled, but only serial HDF5 was found")
endif()

# Compiling requires the Fypp Fortran preprocessor
find_package(Python3 REQUIRED COMPONENTS Interpreter)
find_program(FYPP_EXECUTABLE fypp)
if (NOT FYPP_EXECUTABLE)
  message(FATAL_ERROR
    "fypp not found. Install with:\n"
    "  ${Python3_EXECUTABLE} -m pip install --user fypp\n"
    "or use your package manager/conda.\n")
endif()

# Required compiler flags

add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffree-line-length-none>")

# Force Fortran standard behavior for Intel (not the default!). This meta-option enables a
# collection of individual compiler setting, including the viral `-assume std_mod_proc_name`
# which changes the way module procedures are named internally. This matters to any consumer
# of the module, which must be compiled with the same option. So we revert back to the default
# non-standard naming scheme unless explicitly instructed not to by the `ENABLE_STD_MOD_PROC_NAME`
# cmake option.
option(ENABLE_STD_MOD_PROC_NAME "Build with -assume std_mod_proc_name when using Intel" OFF)
add_compile_options(
  "$<$<COMPILE_LANG_AND_ID:Fortran,Intel,IntelLLVM>:-standard-semantics>"
  "$<$<AND:$<NOT:$<BOOL:${ENABLE_STD_MOD_PROC_NAME}>>,$<COMPILE_LANG_AND_ID:Fortran,Intel,IntelLLVM>>:SHELL:-assume nostd_mod_proc_name>"
)

# Compiler bug workarounds
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel*)
  add_compile_definitions(INTEL_BUG20240327)
endif()

# Fypp preprocessing modules
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(Fypp)

add_subdirectory(src)

include(CTest)
add_subdirectory(test)
