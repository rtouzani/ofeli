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

project(laplace)

add_executable (laplace_demo1 laplace_demo1.cpp)
target_link_libraries (laplace_demo1 ofeli)

add_executable (laplace_demo2 laplace_demo2.cpp)
target_link_libraries (laplace_demo2 ofeli)

add_executable (laplace_demo3 laplace_demo3.cpp)
target_link_libraries (laplace_demo3 ofeli)

add_executable (laplace_demo4 laplace_demo4.cpp)
target_link_libraries (laplace_demo4 ofeli)

file (COPY
      annulus.geo
      annulus.m
      cube.geo
      cube.m
      proj.dat
      proj2.dat
      proj3.dat
      square.geo
      square.m
      DESTINATION .
     )

add_test (laplace1 laplace_demo1 proj2.dat)
add_test (laplace2 laplace_demo2 proj2.dat)
add_test (laplace3 laplace_demo3 proj3.dat)
add_test (laplace4 laplace_demo4 proj.dat)

install (FILES README.md
               laplace_demo1.cpp
               laplace_demo2.cpp
               laplace_demo3.cpp
               laplace_demo4.cpp
               annulus.geo
               annulus.m
               cube.geo
               cube.m
               proj.dat
               proj2.dat
               proj3.dat
               square.geo
               square.m
         DESTINATION ${INSTALL_DEMODIR}/${PROJECT_NAME}
        )
