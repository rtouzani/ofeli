# ==============================================================================
#
#                                    O  F  E  L  I
#
#                           Object  Finite  Element  Library
#
# ==============================================================================
#
#   Copyright (C) 1998 - 2024 Rachid Touzani
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
    add_subdirectory (js)
    add_subdirectory (presentations)

    file (COPY doxygen.css
               extra.css
               head.css
               tabs.css
               ofeli.pdf
               ofeli.cfg
               citing.html
               copyright.html
               demos.html
               download.html
               extlib.html
               fformats.html
               footer.html
               header.html
               index.html
               news.html
               presentations.html
               rita.html
          DESTINATION .
         )

#   Documentation
    set (DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    message ("Doxygen build started")

    doxygen_add_docs (
                      docs
                      ${PROJECT_SOURCE_DIR}
                      COMMENT "Generate documentation"
                     )

    add_custom_target (doc_doxygen ALL
                       COMMAND ${DOXYGEN_EXECUTABLE} ofeli.cfg
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       COMMENT "Generating API documentation with Doxygen"
                       VERBATIM 
                      )

#    install (DIRECTORY html/ DESTINATION ${INSTALL_DOCDIR}/html)
else (DOXYGEN_FOUND)
     message (FATAL_ERROR "Could not find doxygen to build documentation")
endif (DOXYGEN_FOUND)#
