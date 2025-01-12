# ==============================================================================
#
#                                    O  F  E  L  I
#
#                           Object  Finite  Element  Library
#
# ==============================================================================
#
#   Copyright (C) 1998 - 2025 Rachid Touzani
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

project(linear_2d)

add_executable(demo2D1 demo1.cpp)
target_link_libraries(demo2D1 ofeli)

add_executable(demo2D2 demo2.cpp)
target_link_libraries(demo2D2 ofeli)

add_executable(demo2D3 demo3.cpp)
target_link_libraries(demo2D3 ofeli)

file (COPY
      demo1.dat
      demo2.dat
      demo3.dat
      mesh1.m
      mesh2.m
      mesh3.m
      README.md
      DESTINATION .
     )

install (FILES README.md
               demo1.cpp
               demo2.cpp
               demo3.cpp
               mesh1.m
               mesh2.m
               mesh3.m
               demo1.dat
               demo2.dat
               demo3.dat
         DESTINATION ${INSTALL_DEMODIR}/linear_pde/2d/${PROJECT_NAME}
        )

add_test (demo2D1 demo2D1 demo1.dat)
add_test (demo2D2 demo2D2 demo2.dat)
add_test (demo2D3 demo2D3 demo3.dat)
