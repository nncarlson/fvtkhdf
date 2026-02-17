include(CMakePackageConfigHelpers)

# --- Source Definitions ---
set(FYPP_SRC
  vtkhdf_h5.F90.fypp
  vtkhdf_ug_type.F90.fypp
  vtkhdf_ug_file_type.F90.fypp
  vtkhdf_mb_file_type.F90.fypp
)

add_fypp_sources(GEN_SRC SOURCES ${FYPP_SRC} FLAGS --line-numbering)

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

# Generate the library info file and add to source list
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/fvtkhdf_lib_info.F90.in"
  "${CMAKE_CURRENT_BINARY_DIR}/fvtkhdf_lib_info.F90"
  @ONLY
)
list(APPEND SRC "${CMAKE_CURRENT_BINARY_DIR}/fvtkhdf_lib_info.F90")

# Suppress mismatch warnings for legacy MPI interfaces (GNU/NAG)
set(MPI_USERS vtkhdf_ctx_type.F90)
set_property(SOURCE ${MPI_USERS} APPEND PROPERTY COMPILE_OPTIONS
  "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fallow-argument-mismatch;-w>"
  "$<$<COMPILE_LANG_AND_ID:Fortran,NAG>:-mismatch;-w>")

# --- Targets ---

# 1. C Shim Library (Internal) Note: Separate target to prevent
# HDF5/MPI C build details from contaminating Fortran targets
add_library(cshim OBJECT vtkhdf_h5_c_shim.c)

# Note: Object libraries cannot generally link libraries in the traditional sense.
# We just need the include paths here.
target_link_libraries(cshim PUBLIC HDF5::HDF5)
if(ENABLE_MPI)
  target_link_libraries(cshim PUBLIC MPI::MPI_C)
endif()
if(BUILD_SHARED_LIBS)
  set_target_properties(cshim PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

# 2. Main Fortran Library
add_library(fvtkhdf ${SRC} ${GEN_SRC})
target_sources(fvtkhdf PRIVATE $<TARGET_OBJECTS:cshim>)

set_target_properties(fvtkhdf PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

# Include Directories
# Note: We install headers to 'include/fvtkhdf' to keep the namespace clean
target_include_directories(fvtkhdf PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/fvtkhdf>
)

# Linking
if(ENABLE_MPI)
  target_link_libraries(fvtkhdf PUBLIC MPI::MPI_Fortran)
endif()
target_link_libraries(fvtkhdf PRIVATE cshim)

# --- Installation & Export ---

set_target_properties(fvtkhdf PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
    SKIP_BUILD_RPATH FALSE
    BUILD_WITH_INSTALL_RPATH FALSE
)

# 1. Install Targets (Library & Shim)
install(TARGETS fvtkhdf
    EXPORT fVTKHDFTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/fvtkhdf
)

# 2. Install Fortran Modules
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/fvtkhdf
)

# 3. Generate & Install Config File (fVTKHDFConfig.cmake)
set(CONFIG_FILE "${CMAKE_CURRENT_BINARY_DIR}/fVTKHDFConfig.cmake")
set(VERSION_FILE "${CMAKE_CURRENT_BINARY_DIR}/fVTKHDFConfigVersion.cmake")

file(WRITE "${CONFIG_FILE}.in" "
@PACKAGE_INIT@

include(CMakeFindDependencyMacro)

# Always require HDF5
find_dependency(HDF5 COMPONENTS C HL)

# Capture build state
set(fVTKHDF_MPI_SUPPORT @ENABLE_MPI@)

# If built with MPI, require MPI
if(fVTKHDF_MPI_SUPPORT)
    find_dependency(MPI COMPONENTS C Fortran)
endif()

# Safety check: Prevent linking Serial lib if MPI component requested
if(MPI IN_LIST fVTKHDF_FIND_COMPONENTS AND NOT fVTKHDF_MPI_SUPPORT)
    set(fVTKHDF_FOUND FALSE)
    set(fVTKHDF_NOT_FOUND_MESSAGE \"Requested MPI component, but fVTKHDF was built without MPI.\")
endif()

include(\"\${CMAKE_CURRENT_LIST_DIR}/fVTKHDFTargets.cmake\")
check_required_components(fVTKHDF)
")

configure_package_config_file("${CONFIG_FILE}.in" "${CONFIG_FILE}"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/fVTKHDF"
)

write_basic_package_version_file("${VERSION_FILE}"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(FILES "${CONFIG_FILE}" "${VERSION_FILE}"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/fVTKHDF"
)

# 4. Install Export Set (Target definitions)
install(EXPORT fVTKHDFTargets
    FILE fVTKHDFTargets.cmake
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/fVTKHDF"
)
