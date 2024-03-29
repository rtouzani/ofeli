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

include (CTest)

add_subdirectory (1D)
add_subdirectory (linear_pde)
add_subdirectory (adapt)
add_subdirectory (electromagnetics)
add_subdirectory (thermal)
add_subdirectory (solid)
add_subdirectory (fluid)
add_subdirectory (interface)
add_subdirectory (laplace)
add_subdirectory (porous)
add_subdirectory (solvers)

add_custom_target (demos
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   COMMENT "Running demos..."
                  )

message (STATUS "Setting up testsuite - Done")
