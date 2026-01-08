# ==============================================================================
#
#                                    O  F  E  L  I
#
#                           Object  Finite  Element  Library
#
# ==============================================================================
#
#   Copyright (C) 1998 - 2026 Rachid Touzani
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

project(linear_1d)

add_executable(demo1D1 demo1.cpp)
target_link_libraries(demo1D1 ofeli)

add_executable(demo1D2 demo2.cpp)
target_link_libraries(demo1D2 ofeli)

file (COPY
      README.md
      DESTINATION .
     )

install (FILES README.md
               demo1.cpp
               demo2.cpp
         DESTINATION ${INSTALL_DEMODIR}/linear_pde/1d/${PROJECT_NAME}
        )

add_test (demo1D1 demo1D1 10)
add_test (demo1D2 demo1D2 10 0.1)
