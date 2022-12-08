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

project (tins2)

add_executable (tins2 main.cpp)
target_link_libraries (tins2 ofeli)

file (COPY
      cavity-30.m
      cavity-60.dat
      cavity-60.m
      cavity.dat
      step.dat
      step.dom
      step.m
      DESTINATION .
     )

add_test (tins2 tins2 cavity.dat)

install (FILES
         README.md
         main.cpp
         step.dat
         cavity.dat
         cavity-30.m
         cavity-60.m
         cavity-60.dat
         step.dom
         step.m
         DESTINATION ${INSTALL_DEMODIR}/fluid/${PROJECT_NAME}
        )
