# ==============================================================================
#
#                                    O  F  E  L  I
#
#                           Object  Finite  Element  Library
#
# ==============================================================================
#
#   Copyright (C) 1998 - 2022 Rachid Touzani
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

project (EIGEN)

add_executable (eigen_demo1 eigen_demo1.cpp)
target_link_libraries (eigen_demo1 ofeli)

add_executable (eigen_demo2 eigen_demo2.cpp)
target_link_libraries (eigen_demo2 ofeli)

file (COPY eigen_demo2.dat test.m DESTINATION .)

install (TARGETS eigen_demo1 eigen_demo2 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install (FILES README.md
               eigen_demo1.cpp
               eigen_demo2.cpp
               eigen_demo2.dat
               test.m
         DESTINATION ${INSTALL_DEMODIR}/solvers/${PROJECT_NAME}
        )
