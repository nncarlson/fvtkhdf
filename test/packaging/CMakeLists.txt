cmake_minimum_required(VERSION 3.20)
project(fvtkhdf_consumer LANGUAGES C Fortran)

option(ENABLE_MPI "Build MPI consumer" ON)

find_package(fVTKHDF CONFIG REQUIRED)

if(ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran C) # ensures mpi_f08.mod is available
  add_executable(consumer consumer_mpi.f90)
else()
  add_executable(consumer consumer_serial.f90)
endif()

target_link_libraries(consumer PRIVATE fVTKHDF::fvtkhdf)
