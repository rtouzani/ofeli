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

project (1D)

add_executable (elliptic elliptic.cpp)
target_link_libraries (elliptic ofeli)

add_executable (heat heat.cpp)
target_link_libraries (heat ofeli)

add_executable (transport transport.cpp)
target_link_libraries (transport ofeli)

add_test (1D-1 elliptic)
add_test (1D-2 heat 10 0.1)
add_test (1D-3 transport 10 0.1)

install (FILES
         README.md
         elliptic.cpp
         heat.cpp
         transport.cpp
         DESTINATION ${INSTALL_DEMODIR}/${PROJECT_NAME}
        )
