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

project (adapt)

add_executable (ad1 adapt_demo1.cpp)
target_link_libraries (ad1 ofeli)

add_executable (ad2 adapt_demo2.cpp)
target_link_libraries (ad2 ofeli)

file (COPY L.dom rect.dom DESTINATION .)

add_test (ad1 ad1 rect.dom)
add_test (ad2 ad2 rect.dom)

install (FILES
         README.md
         adapt_demo1.cpp
         adapt_demo2.cpp
         L.dom
         rect.dom
         DESTINATION ${INSTALL_DEMODIR}/${PROJECT_NAME}
        )
