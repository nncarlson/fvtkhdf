cmake_minimum_required(VERSION 3.20)
project(fvtkhdf_mb_compatibility LANGUAGES C Fortran)

find_package(fVTKHDF CONFIG REQUIRED)
if(fVTKHDF_MPI_SUPPORT)
  message(FATAL_ERROR "This initial compatibility harness requires a serial fVTKHDF installation")
endif()

add_executable(generate generate.F90)
target_link_libraries(generate PRIVATE fVTKHDF::fvtkhdf)
