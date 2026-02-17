configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
  "${CMAKE_CURRENT_BINARY_DIR}/conf.py"
  @ONLY
)

add_custom_target(docs
  COMMAND ${SPHINX_EXECUTABLE} -b singlehtml
    -c "${CMAKE_CURRENT_BINARY_DIR}"   # Use the generated conf.py
    "${CMAKE_CURRENT_SOURCE_DIR}"     # Source directory for .rst files
    "${CMAKE_CURRENT_BINARY_DIR}/html" # Output directory
  COMMENT "Building single-page Sphinx documentation"
)
