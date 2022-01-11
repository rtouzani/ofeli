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
#
#                         cmake for the ofeli project
#
# ==============================================================================


cmake_minimum_required (VERSION 3.16.0)
message (STATUS "This is CMake ${CMAKE_VERSION}")
message (STATUS "")

# General information

set (OFELI_PACKAGE_VENDOR
     "The ofeli authors <http://ofeli.org/>"
    )
set (OFELI_PACKAGE_DESCRIPTION
     "An object finite element library"
    )

# Set up ofeli Version

file (STRINGS "${CMAKE_SOURCE_DIR}/VERSION" _version LIMIT_COUNT 2)
set (OFELI_VERSION "${_version}")
string (REGEX REPLACE "^([0-9]+)\\..*" "\\1"
        OFELI_VERSION_MAJOR "${OFELI_VERSION}"
       )
string (REGEX REPLACE "^[0-9]+\\.([0-9]+).*" "\\1"
        OFELI_VERSION_MINOR "${OFELI_VERSION}"
       )
string (REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1"
        OFELI_VERSION_SUBMINOR "${OFELI_VERSION}"
       )
#string (FIND "${OFELI_VERSION}" ";." NB)
string (SUBSTRING "${OFELI_VERSION}" 6 -1 OFELI_VERSION_DATE)
set (OFELI_VERSION ${OFELI_VERSION_MAJOR}.${OFELI_VERSION_MINOR}.${OFELI_VERSION_SUBMINOR})
string(REPLACE "^" "" "${OFELI_VERSION_DATE}" OFELI_VERSION_DATE)

project (ofeli)

message (STATUS "Building release ${PROJECT_NAME}-${OFELI_VERSION}, ${OFELI_VERSION_DATE}")

# Configuration

set (CMAKE_BUILD_TYPE "Release")
set (CMAKE_EXPORT_COMPILE_COMMANDS ON)
set (CMAKE_CXX_STANDARD 11)
enable_testing ()

option (BUILD_TESTS "Build test programs" OFF)

add_definitions (-Wall -fPIE -Ofast -mtune=native -D_REENTRANT -DNDEBUG -fstrict-aliasing -fno-strict-overflow)

# Define installation directories
#
include (GNUInstallDirs)
set (INSTALL_BINDIR ${CMAKE_INSTALL_PREFIX}/bin)
set (INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
set (INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/include/ofeli)
set (INSTALL_OFELIDIR ${CMAKE_INSTALL_PREFIX}/share/ofeli)
set (INSTALL_DOCDIR ${INSTALL_OFELIDIR}/doc)
set (INSTALL_MATDIR ${INSTALL_OFELIDIR}/material)
set (INSTALL_DEMODIR ${INSTALL_OFELIDIR}/demos)
set (INSTALL_UTILDIR ${INSTALL_OFELIDIR}/util)

# Define path for material data
#
option (MATERIAL_ENABLE "Enable Material Path" ON)
if (MATERIAL_ENABLE)
  set (PATH_MATERIAL "${INSTALL_MATDIR}")
endif ()
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/include/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)
file (APPEND ${CMAKE_CURRENT_BINARY_DIR}/config.h "\n\n#define OFELI_VERSION \"${OFELI_VERSION}\"\n")
file (APPEND ${CMAKE_CURRENT_BINARY_DIR}/config.h "#define OFELI_DATE \"${OFELI_VERSION_DATE}\"")

include_directories (include ${CMAKE_CURRENT_BINARY_DIR})
add_subdirectory (src)
add_subdirectory (include)
add_subdirectory (demos)
add_subdirectory (util)

# Directories to install
install (TARGETS ofeli LIBRARY DESTINATION ${INSTALL_LIBDIR})
install (DIRECTORY include/ DESTINATION ${INSTALL_INCLUDEDIR})
install (DIRECTORY doc/ DESTINATION ${INSTALL_DOCDIR})
install (DIRECTORY material/ DESTINATION ${INSTALL_MATDIR})
install (DIRECTORY demos/ DESTINATION ${INSTALL_DEMODIR})
install (DIRECTORY util/ DESTINATION ${INSTALL_UTILDIR})
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/config.h DESTINATION ${INSTALL_INCLUDEDIR})

# Preparing Cpack
include (${CMAKE_SOURCE_DIR}/cmake/setup_cpack.cmake)

# uninstall target
if (NOT TARGET uninstall)
   configure_file (
                   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/uninstall.cmake.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/uninstall.cmake"
                   IMMEDIATE @ONLY
                  )

  add_custom_target (uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake/uninstall.cmake)
endif ()
