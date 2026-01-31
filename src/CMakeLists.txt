set(SRC
  f90_assert.F90
  hdf5_c_binding.F90
  hl_hdf5.F90
  vtkhdf_file_type.F90
  vtkhdf_ug_type.F90
)

set(MPI_USERS hl_hdf5.F90)
set_property(SOURCE ${MPI_USERS} APPEND PROPERTY COMPILE_OPTIONS
  "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fallow-argument-mismatch;-w>"
  "$<$<COMPILE_LANG_AND_ID:Fortran,NAG>:-mismatch;-w>")

# Separate C target to prevent HDF5/MPI C build details from contaminating Fortran targets
add_library(cshim STATIC hdf5_c_binding_ext.c)
target_link_libraries(cshim PUBLIC MPI::MPI_C HDF5::HDF5)
if (BUILD_SHARED_LIBS)
  set_target_properties(cshim PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

add_library(fvtkhdf ${SRC})

set_target_properties(fvtkhdf PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)

target_include_directories(fvtkhdf PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)

target_link_libraries(fvtkhdf MPI::MPI_Fortran cshim)
