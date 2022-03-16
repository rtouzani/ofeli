# ==============================================================================
#
#                                    O  F  E  L  I
#
#                           Object  Finite  Element  Library
#
# ==============================================================================
#
#   Copyright (C) 1998 - 2022 Rachid Touzani
#
#   This file is part of OFELI.
#
#   OFELI is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   OFELI is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with OFELI. If not, see <http://www.gnu.org/licenses/>.
#
# ==============================================================================

message (STATUS "Building documentation ...")
message (STATUS "")


find_package (Doxygen)
if (DOXYGEN_FOUND)

    add_subdirectory (css)
    add_subdirectory (demos)
    add_subdirectory (im)
    add_subdirectory (presentations)

    file (COPY citing.html
               copyright.html
               ofeli.cfg
               demos.html
               download.html
               extlib.html
               fformats.html
               footer.html
               header.html
               index.html
               installation.html
               news.html
               presentations.html
               doxygen.css
               head.css
               tabs.css
               logo.png
               ofeli.pdf
               userguide.pdf
          DESTINATION .
         )

    set (DOXYGEN_PROJECT_NAME "OFELI")
    set (DOXYGEN_PROJECT_BRIEF "The OFELI Library")
    set (DOXYGEN_PROJECT_LOGO logo.png)
    set (DOXYGEN_EXCLUDE_PATTERNS "${PROJECT_SOURCE_DIR}/build/*")
    set (DOXYGEN_FULL_PATH_NAMES NO)
    set (DOXYGEN_USE_MDFILE_AS_MAINPAGE README.md)
    set (DOXYGEN_WARN_NO_PARAMDOC YES)
    set (DOXYGEN_INPUT "${PROJECT_SOURCE_DIR}/include")
    set (DOXYGEN_GENERATE_HTML YES)
    set (DOXYGEN_HTML_OUTPUT html)
    set (DOXYGEN_HIDE_UNDOC_MEMBERS YES)
    set (DOXYGEN_HIDE_UNDOC_CLASSES YES)
    set (DOXYGEN_HIDE_FRIEND_COMPOUNDS YES)
    set (DOXYGEN_HIDE_SCOPE_NAMES YES)
    set (DOXYGEN_INLINE_INFO NO)
    set (DOXYGEN_IMAGE_PATH im)
    set (DOXYGEN_HTML_HEADER header.html)
    set (DOXYGEN_HTML_FOOTER footer.html)
    set (DOXYGEN_BINARY_TOC YES)
    set (DOXYGEN_TOC_EXPAND YES)

    set (DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    message ("Doxygen build started")

    doxygen_add_docs (
                      docs
                      ${PROJECT_SOURCE_DIR}
                      COMMENT "Generate man pages"
                     )

    add_custom_target (doc_doxygen ALL
                       COMMAND ${DOXYGEN_EXECUTABLE} 
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       COMMENT "Generating API documentation with Doxygen"
                       VERBATIM 
                      )

    install (DIRECTORY html/ DESTINATION ${INSTALL_DOCDIR}/html)
else (DOXYGEN_FOUND)
      message (FATAL_ERROR "Could not find doxygen to build documentation")
endif (DOXYGEN_FOUND)
