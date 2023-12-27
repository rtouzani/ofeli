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

project (cfield)

add_executable (cfield cfield.cpp)
target_link_libraries (cfield ofeli)

file (COPY cavity.m
           cavity.pos
           cavity.s
           spider_tecplot.dat
           README.md
       DESTINATION .
      )

add_test (cfield cfield -f gmsh -m cavity.m -i cavity.s -o cavity.pos)

install (TARGETS cfield RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install (FILES README.md  
              cavity.m
              cavity.pos
              cavity.s
              spider_tecplot.dat
         DESTINATION ${INSTALL_UTILDIR}/${PROJECT_NAME}
        )
