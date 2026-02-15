

####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was fVTKHDFConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

include(CMakeFindDependencyMacro)

# Always require HDF5
find_dependency(HDF5 COMPONENTS C HL)

# Capture build state
set(fVTKHDF_MPI_SUPPORT OFF)

# If built with MPI, require MPI
if(fVTKHDF_MPI_SUPPORT)
    find_dependency(MPI COMPONENTS C Fortran)
endif()

# Safety check: Prevent linking Serial lib if MPI component requested
if(MPI IN_LIST fVTKHDF_FIND_COMPONENTS AND NOT fVTKHDF_MPI_SUPPORT)
    set(fVTKHDF_FOUND FALSE)
    set(fVTKHDF_NOT_FOUND_MESSAGE "Requested MPI component, but fVTKHDF was built without MPI.")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/fVTKHDFTargets.cmake")
check_required_components(fVTKHDF)
