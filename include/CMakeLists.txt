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

#
# Configure config.h
configure_file (
                ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
                ${CMAKE_CURRENT_BINARY_DIR}/config.h
               )
               project(ofeli)

file (COPY Acoustics.h
           CL.h
           Electromagnetics.h
           Fluid.h
           Laplace.h
           LinearPDE.h
           OFELI.h
           OFELIException.h
           OFELI_Config.h
           Porous.h
           Solid.h
           Therm.h
      DESTINATION .
     )

add_subdirectory (equations)
add_subdirectory (io)
add_subdirectory (linear_algebra)
add_subdirectory (mesh)
add_subdirectory (post)
add_subdirectory (shape_functions)
add_subdirectory (solvers)
add_subdirectory (util)
