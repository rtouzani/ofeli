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

project (ODE)

add_executable (ode_demo1 ode_demo1.cpp)
target_link_libraries (ode_demo1 ofeli)

add_executable (ode_demo2 ode_demo2.cpp)
target_link_libraries (ode_demo2 ofeli)

add_executable (ode_demo3 ode_demo3.cpp)
target_link_libraries (ode_demo3 ofeli)

add_executable (ode_demo4 ode_demo4.cpp)
target_link_libraries (ode_demo4 ofeli)

add_executable (ode_demo5 ode_demo5.cpp)
target_link_libraries (ode_demo5 ofeli)

add_executable (ode_demo6 ode_demo6.cpp)
target_link_libraries (ode_demo6 ofeli)

add_executable (ode_demo7 ode_demo7.cpp)
target_link_libraries (ode_demo7 ofeli)

add_executable (ts_demo1 ts_demo1.cpp)
target_link_libraries (ts_demo1 ofeli)

add_executable (ts_demo2 ts_demo2.cpp)
target_link_libraries (ts_demo2 ofeli)

file (COPY mesh1.m mesh2.m DESTINATION .)

add_test (ode1 ode_demo1 0.1)
add_test (ode2 ode_demo2 0.1)
add_test (ode3 ode_demo3 0.1)
add_test (ode4 ode_demo4 0.1)
add_test (ode5 ode_demo5 0.1)
add_test (ode6 ode_demo6 0.1 1.0)
add_test (ode7 ode_demo7 0.1 1.0 output.dat)
add_test (ts1 ts_demo1 mesh1.m 0.1)
add_test (ts2 ts_demo2 mesh2.m 0.1)

install (FILES README.md
               ode_demo1.cpp
               ode_demo2.cpp
               ode_demo3.cpp
               ode_demo4.cpp
               ode_demo5.cpp
               ode_demo6.cpp
               ode_demo7.cpp
               ts_demo1.cpp
               ts_demo2.cpp
               mesh1.m
               mesh2.m
         DESTINATION ${INSTALL_DEMODIR}/solvers/${PROJECT_NAME}
        )
