add_custom_target(docs
  COMMAND ${SPHINX_EXECUTABLE} -b html
    "${CMAKE_CURRENT_SOURCE_DIR}"     # Source directory for .rst files
    "${CMAKE_CURRENT_BINARY_DIR}/html" # Output directory
  COMMENT "Building Sphinx documentation"
  VERBATIM
)
