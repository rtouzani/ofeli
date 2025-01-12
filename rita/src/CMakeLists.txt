# ==============================================================================
#
#                                 r  i  t  a
#
#            An environment for Modelling and Numerical Simulation
#
# ==============================================================================
#
#   Copyright (C) 2021 - 2025  Rachid Touzani
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
#   GNU General Public License for more details.
#
# ==============================================================================

# Look for the ofeli library
#find_library (OFELI_LIB NAMES ofeli PATHS ${OFELI_LIBDIR})
#if (NOT OFELI_LIB)
#   message (FATAL_ERROR "ofeli library not found")
#endif ()

# Look for the gmsh library
find_library (GMSH_LIB NAMES gmsh PATHS ${GMSH_LIBDIR})
if (NOT GMSH_LIB)
  message (FATAL_ERROR "gmsh package not found")
endif ()

add_executable (rita)
target_link_libraries (rita ofeli)

target_sources (${PROJECT_NAME} PRIVATE
                ae.cpp
                approximation.cpp
                calc.cpp
                cmd.cpp
                configure.cpp
                data.cpp
                eigen.cpp
                funct.cpp
                help.cpp
                integration.cpp
                ls.cpp
                main.cpp
                mesh.cpp
                ode.cpp
                optim.cpp
                pde.cpp
                plot.cpp
                rita.cpp
                solve.cpp
                transient.cpp
               )

target_link_libraries (rita ${OFELI_LIB};${GMSH_LIB})
#install (TARGETS rita RUNTIME DESTINATION ${INSTALL_BINDIR})

#
# Configure config.h
configure_file (
                ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
                ${CMAKE_CURRENT_BINARY_DIR}/config.h
               )
