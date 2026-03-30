cmake_minimum_required(VERSION 3.20)
project(fvtkhdf_consumer LANGUAGES C Fortran)

option(ENABLE_MPI "Build MPI consumer" ON)

find_package(fVTKHDF CONFIG REQUIRED)

if(ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran C)
  add_executable(consumer consumer_mpi.F90)
  if(fVTKHDF_MPI_F08_SUPPORT)
    target_compile_definitions(consumer PRIVATE USE_MPI_F08)
  endif()
else()
  add_executable(consumer consumer_serial.f90)
endif()

target_link_libraries(consumer PRIVATE fVTKHDF::fvtkhdf)
