# ==============================================================================
#
#                                    O  F  E  L  I
#
#                           Object  Finite  Element  Library
#
# ==============================================================================
#
#   Copyright (C) 1998 - 2023 Rachid Touzani
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

project (OPTIM)

add_executable (opt_demo1 opt_demo1.cpp)
target_link_libraries (opt_demo1 ofeli)

add_executable (opt_demo2 opt_demo2.cpp)
target_link_libraries (opt_demo2 ofeli)

add_executable (opt_demo3 opt_demo3.cpp)
target_link_libraries (opt_demo3 ofeli)

add_executable (opt_demo4 opt_demo4.cpp)
target_link_libraries (opt_demo4 ofeli)

add_executable (opt_demo5 opt_demo5.cpp)
target_link_libraries (opt_demo5 ofeli)

add_executable (opt_demo6 opt_demo6.cpp)
target_link_libraries (opt_demo6 ofeli)

file (COPY test.dat test.m DESTINATION .)

add_test (opt1 opt_demo1)
add_test (opt1 opt_demo2)
add_test (opt1 opt_demo3 test.dat)
add_test (opt1 opt_demo4)
add_test (opt1 opt_demo5)
add_test (opt1 opt_demo6)

install (FILES README.md
               opt_demo1.cpp
               opt_demo2.cpp
               opt_demo3.cpp
               opt_demo4.cpp
               opt_demo5.cpp
               opt_demo6.cpp
               test.dat
               test.m
         DESTINATION ${INSTALL_DEMODIR}/solvers/${PROJECT_NAME}
        )
