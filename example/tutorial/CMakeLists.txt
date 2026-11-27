cmake_minimum_required(VERSION 3.28)
project(fvtkhdf_tutorial LANGUAGES C Fortran)

find_package(fVTKHDF REQUIRED)

set(examples image_static image_temporal triangles_static)
if(fVTKHDF_MPI_SUPPORT)
  list(APPEND examples image_parallel image_parallel_points
                       rectilinear_parallel rectilinear_moving
                       structured_parallel structured_moving
                       triangles_parallel triangles_temporal)
endif()

foreach(example IN LISTS examples)
  add_executable(${example} ${example}.F90)
  set_target_properties(${example} PROPERTIES Fortran_PREPROCESS ON)
  target_link_libraries(${example} PRIVATE fVTKHDF::fvtkhdf)
  if(fVTKHDF_MPI_SUPPORT)
    target_compile_definitions(${example} PRIVATE USE_MPI)
  endif()
endforeach()
