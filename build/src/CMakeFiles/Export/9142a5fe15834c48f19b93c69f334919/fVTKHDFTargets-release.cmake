#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "fvtkhdf" for configuration "Release"
set_property(TARGET fvtkhdf APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(fvtkhdf PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libfvtkhdf.so.0.1.0"
  IMPORTED_SONAME_RELEASE "libfvtkhdf.so.0"
  )

list(APPEND _cmake_import_check_targets fvtkhdf )
list(APPEND _cmake_import_check_files_for_fvtkhdf "${_IMPORT_PREFIX}/lib/libfvtkhdf.so.0.1.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
