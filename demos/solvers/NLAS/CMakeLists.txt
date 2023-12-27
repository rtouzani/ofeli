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

project (NLAS)

add_executable (nl_demo1 nl_demo1.cpp)
target_link_libraries (nl_demo1 ofeli)

add_executable (nl_demo2 nl_demo2.cpp)
target_link_libraries (nl_demo2 ofeli)

add_executable (nl_demo3 nl_demo3.cpp)
target_link_libraries (nl_demo3 ofeli)

add_executable (nl_demo4 nl_demo4.cpp)
target_link_libraries (nl_demo4 ofeli)

add_test (nl1 nl_demo1)
add_test (nl2 nl_demo2)
add_test (nl3 nl_demo3)
add_test (nl4 nl_demo4)

install (FILES README.md
               nl_demo1.cpp
               nl_demo2.cpp
               nl_demo3.cpp
               nl_demo4.cpp
         DESTINATION ${INSTALL_DEMODIR}/solvers/${PROJECT_NAME}
        )
