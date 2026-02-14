cmake_minimum_required(VERSION 3.30)

project(fVTKHDF VERSION 0.1.0 LANGUAGES Fortran C)

# --- Options ---
option(ENABLE_MPI "Enable MPI support" ON)
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(ENABLE_STD_MOD_PROC_NAME "Build with -assume std_mod_proc_name when using Intel" OFF)

# Default to Release build if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
endif()

# Standardize install paths (lib64, bin, include)
include(GNUInstallDirs)

# --- Dependencies ---

# 1. MPI
if(ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS C Fortran)
  set(HDF5_PREFER_PARALLEL True)
  add_compile_definitions(USE_MPI)
else()
  set(HDF5_PREFER_PARALLEL False)
endif()

# 2. HDF5
find_package(HDF5 REQUIRED COMPONENTS C HL)

if(ENABLE_MPI AND NOT HDF5_IS_PARALLEL)
  message(FATAL_ERROR "Parallel HDF5 is required when MPI is enabled, but only serial HDF5 was found.")
endif()

# 3. Python & Fypp
find_package(Python3 REQUIRED COMPONENTS Interpreter)
find_program(FYPP_EXECUTABLE fypp)
if(NOT FYPP_EXECUTABLE)
  message(FATAL_ERROR "fypp not found.\nInstall with: ${Python3_EXECUTABLE} -m pip install fypp")
endif()

# --- Compiler Settings ---

# Disable assertions in Release
add_compile_definitions("$<$<NOT:$<CONFIG:Debug>>:NDEBUG>")

# Fortran Flags
add_compile_options(
  "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffree-line-length-none>"
  "$<$<COMPILE_LANG_AND_ID:Fortran,NAG>:-w=longlines>"
)

# Intel Standard Semantics
add_compile_options(
  "$<$<COMPILE_LANG_AND_ID:Fortran,Intel,IntelLLVM>:-standard-semantics>"
  "$<$<AND:$<NOT:$<BOOL:${ENABLE_STD_MOD_PROC_NAME}>>,$<COMPILE_LANG_AND_ID:Fortran,Intel,IntelLLVM>>:SHELL:-assume nostd_mod_proc_name>"
)

# Compiler Bug Workarounds
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  add_compile_definitions(INTEL_BUG20240327)
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  add_compile_definitions(GNU_PR123899)
endif()

# --- Subdirectories ---
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(Fypp)

add_subdirectory(src)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()
