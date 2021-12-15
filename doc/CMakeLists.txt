find_package(Doxygen)
if(DOXYGEN_FOUND)
  set(DOXYGEN_EXCLUDE_PATTERNS "${PROJECT_SOURCE_DIR}/build/*")
  set(DOXYGEN_FULL_PATH_NAMES NO)
  set(DOXYGEN_USE_MDFILE_AS_MAINPAGE README.md)
  set(DOXYGEN_WARN_NO_PARAMDOC YES)
  doxygen_add_docs(
    doc ${CMAKE_SOURCE_DIR}
    COMMENT "Generating doxygen documentation for ${PROJECT_NAME}")
endif(DOXYGEN_FOUND)