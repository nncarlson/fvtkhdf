function(add_fypp_sources out_var)
  # Usage:
  #   add_fypp_sources(GEN_SRCS
  #     SOURCES a.F90.fypp b.F90.fypp
  #     INCLUDES fypp_macros.fypp f90_assert.fpp
  #     FLAGS --line-numbering
  #   )

  set(options)
  set(oneValueArgs)
  set(multiValueArgs SOURCES INCLUDES FLAGS)
  cmake_parse_arguments(FYPP "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT FYPP_SOURCES)
    message(FATAL_ERROR "add_fypp_sources: SOURCES is required")
  endif()

  set(gen_files "")
  foreach(infile IN LISTS FYPP_SOURCES)
    get_filename_component(in_abs "${infile}" ABSOLUTE)
    get_filename_component(fname "${infile}" NAME)

    # Output name: strip trailing ".fypp" if present
    set(outname "${fname}")
    string(REGEX REPLACE "\\.fypp$" "" outname "${outname}")

    set(outfile "${CMAKE_CURRENT_BINARY_DIR}/${outname}")
    #set_source_files_properties(${outfile} PROPERTIES GENERATED TRUE LANGUAGE Fortran Fortran_PREPROCESS ON)

    add_custom_command(
      OUTPUT  "${outfile}"
      COMMAND "${FYPP_EXECUTABLE}" ${FYPP_FLAGS} "${in_abs}" "${outfile}"
      DEPENDS "${in_abs}" ${FYPP_INCLUDES}
      VERBATIM
      COMMENT "Fypp: ${fname} -> ${outname}"
    )

    list(APPEND gen_files "${outfile}")
  endforeach()

  set(${out_var} ${gen_files} PARENT_SCOPE)
endfunction()
