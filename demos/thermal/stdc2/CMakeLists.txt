# ==============================================================================
#
#                                    O  F  E  L  I
#
#                           Object  Finite  Element  Library
#
#  ==============================================================================
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
#  ==============================================================================

project(stdc2)

add_executable(stdc2 main.cpp)
target_link_libraries(stdc2 ofeli)

file(COPY proj.dat proj.m proj.v DESTINATION .)

add_test (stdc2 stdc2 proj.dat)

install (FILES README.md
               main.cpp
               proj.dat
               proj.m
               proj.v
         DESTINATION ${INSTALL_DEMODIR}/thermal/${PROJECT_NAME}
        )

