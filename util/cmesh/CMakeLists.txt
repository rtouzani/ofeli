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

project (cmesh)

add_executable (cmesh cmesh.cpp)
target_link_libraries (cmesh ofeli)

install (TARGETS cmesh RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

file (COPY README.md
           A.ele
           A.node
           beam.m
           cube.ele
           cube.node
           cube.vol
           disk.bamg
           disk.m
           gear.ele
           gear.node
           l.bamg
           torus.neu
       DESTINATION .
      )

add_test (cmesh cmesh --input disk.bamg --from bamg --to ofeli)

install (FILES README.md
               A.ele
               A.node
               beam.m
               cube.ele
               cube.node
               cube.vol
               disk.bamg
               disk.m
               gear.ele
               gear.node
               l.bamg
               torus.neu
         DESTINATION ${INSTALL_UTILDIR}/${PROJECT_NAME}
        )
