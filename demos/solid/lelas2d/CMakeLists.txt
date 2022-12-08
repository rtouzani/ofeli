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

project (lelas2d)

add_executable (lelas2d main.cpp)
target_link_libraries (lelas2d ofeli)

file (COPY
      beam-10x5.m
      beam-20x10.m
      beam.dat
      beam.m
      DESTINATION .
     )

add_test (lelas2d lelas2d beam.dat)

install (FILES README.md
               main.cpp
               beam-10x5.m
               beam-20x10.m
               beam.dat
               beam.m
         DESTINATION ${INSTALL_DEMODIR}/solid/${PROJECT_NAME}
        )
