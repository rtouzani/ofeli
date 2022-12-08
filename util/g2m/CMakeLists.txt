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

project(g2m)

add_executable (g2m main.cpp)
target_link_libraries (g2m ofeli)

file (COPY test.dom
           README.md
       DESTINATION .
      )

add_test (g2m g2m -d test.dom -o test.m)

install (FILES README.md test.dom
         DESTINATION ${INSTALL_UTILDIR}/${PROJECT_NAME}
        )
install (TARGETS g2m RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
