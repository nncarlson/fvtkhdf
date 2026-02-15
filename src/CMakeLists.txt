set(FYPP_SRC
  vtkhdf_h5.F90.fypp
  vtkhdf_ug_type.F90.fypp
  vtkhdf_mb_file_type.F90.fypp
  vtkhdf_ug_file_type.F90.fypp
)

add_fypp_sources(GEN_SRC SOURCES ${FYPP_SRC} FLAGS --line-numbering)
#add_fypp_sources(GEN_SRC SOURCES ${FYPP_SRC})

set(SRC
  vtkhdf_assert.F90
  vtkhdf_vtk_cell_types.F90
  vtkhdf_h5_c_binding.F90
  vtkhdf_version_param.F90
  vtkhdf_h5e_context_type.F90
)

if(ENABLE_MPI)
  list(APPEND SRC vtkhdf_ctx_type.F90)
else()
  list(APPEND SRC vtkhdf_ctx_type-serial.F90)
endif()

set(MPI_USERS vtkhdf_ctx_type.F90)
set_property(SOURCE ${MPI_USERS} APPEND PROPERTY COMPILE_OPTIONS
  "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fallow-argument-mismatch;-w>"
  "$<$<COMPILE_LANG_AND_ID:Fortran,NAG>:-mismatch;-w>")

# Separate C target to prevent HDF5/MPI C build details from contaminating Fortran targets
add_library(cshim STATIC vtkhdf_h5_c_shim.c)
if(USE_MPI)
  target_link_libraries(cshim PUBLIC MPI::MPI_C)
endif()
target_link_libraries(cshim PUBLIC HDF5::HDF5)
if (BUILD_SHARED_LIBS)
  set_target_properties(cshim PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

add_library(fvtkhdf ${SRC} ${GEN_SRC})

set_target_properties(fvtkhdf PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)

target_include_directories(fvtkhdf PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)

if(ENABLE_MPI)
  target_link_libraries(fvtkhdf MPI::MPI_Fortran cshim)
endif()
target_link_libraries(fvtkhdf cshim)
