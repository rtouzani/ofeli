# ==============================================================================
#
#                                 r  i  t  a
#
#            An environment for Modelling and Numerical Simulation
#
# ==============================================================================
#
#   Copyright (C) 2021 - 2026 Rachid Touzani
#
#   This file is part of rita.
#
#   rita is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   rita is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
# ==============================================================================
#
#                          cmake for the rita project
#
# ==============================================================================

# General information

set (RITA_PACKAGE_VENDOR
     "The rita authors <http://ofeli.org/rita.html>"
    )
set (RITA_PACKAGE_DESCRIPTION 
     "A console application for modelling and numerical simulation"
    )

# Set up rita Version

project (rita)

file (STRINGS "${CMAKE_SOURCE_DIR}/rita/src/VERSION" _version LIMIT_COUNT 2)
set (RITA_VERSION "${_version}")
string (REGEX REPLACE "^([0-9]+)\\..*" "\\1"
        RITA_VERSION_MAJOR "${RITA_VERSION}"
       )
string (REGEX REPLACE "^[0-9]+\\.([0-9]+).*" "\\1"
        RITA_VERSION_MINOR "${RITA_VERSION}"
       )
string (REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1"
        RITA_VERSION_SUBMINOR "${RITA_VERSION}"
       )

string (SUBSTRING "${RITA_VERSION}" 6 -1 RITA_VERSION_DATE)
set (RITA_VERSION ${RITA_VERSION_MAJOR}.${RITA_VERSION_MINOR}.${RITA_VERSION_SUBMINOR})
string (REPLACE "^" "" "${RITA_VERSION_DATE}" RITA_VERSION_DATE)
set (INSTALL_TUTORIALDIR ${INSTALL_RITADIR}/tutorial)
set (INSTALL_APPDIR ${INSTALL_RITADIR}/applications)

message (STATUS "Building release ${PROJECT_NAME}-${RITA_VERSION}, ${RITA_VERSION_DATE}")

# Check gmsh
include(FindPackageHandleStandardArgs)
find_library(GMSH_LIBRARY gmsh)
find_package_handle_standard_args(GMSH REQUIRED_VARS GMSH_LIBRARY)

if (GMSH_FOUND)
   message (STATUS "gmsh package found")  
   add_definitions (-DUSE_GMSH)
else (GMSH_FOUND)
   message (FATAL_ERROR "Could not find package gmsh\n")
endif (GMSH_FOUND)


# Configuration
set (CMAKE_BUILD_TYPE "Release")
set (CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Compiler specific flags
add_definitions (-DNDEBUG)

enable_testing ()
option (BUILD_TESTS "Build test programs" OFF)

set (OFELI_LIBDIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set (OFELI_INCDIR ${CMAKE_SOURCE_DIR}/include)

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/src/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)
file (APPEND ${CMAKE_CURRENT_BINARY_DIR}/config.h "\n\n#define RITA_VERSION \"${RITA_VERSION}\"\n")
file (APPEND ${CMAKE_CURRENT_BINARY_DIR}/config.h "#define RITA_DATE \"${RITA_VERSION_DATE}\"")
configure_file (${CMAKE_CURRENT_BINARY_DIR}/config.h ${CMAKE_CURRENT_SOURCE_DIR}/src/config.h)

include_directories (${OFELI_INCDIR} ${GMSH_INCDIR})

add_subdirectory (src)
add_subdirectory (muparserx)
add_subdirectory (tutorial)

install (TARGETS rita RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install (DIRECTORY tutorial/ DESTINATION ${INSTALL_TUTORIALDIR})
