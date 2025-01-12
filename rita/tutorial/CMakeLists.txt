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

include (CTest)

set (CMAKE_RITA_EXEC ${CMAKE_BINARY_DIR}/rita/src/rita)

add_subdirectory (ae)
add_subdirectory (approximation)
add_subdirectory (calc)
add_subdirectory (eigen)
add_subdirectory (integration)
add_subdirectory (ls)
add_subdirectory (ode)
add_subdirectory (optim)
add_subdirectory (pde)
add_subdirectory (plot)


add_custom_target (tutorial
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   COMMENT "Running tutorial..."
                  )

message (STATUS "Setting up testsuite - Done")
