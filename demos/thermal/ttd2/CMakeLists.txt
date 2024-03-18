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

project(ttd2)

add_executable(ttd2 main.cpp)
target_link_libraries(ttd2 ofeli)

target_sources(ttd2 PRIVATE error.cpp)

file(COPY proj-5x10.m proj-5x20.m proj-5x40.m proj-5x80.m proj.dat DESTINATION .)

add_test (ttd2 ttd2 proj.dat)

install (FILES README.md
               main.cpp
               error.cpp
               proj-5x10.m
               proj-5x20.m
               proj-5x40.m
               proj-5x80.m
               proj.dat
         DESTINATION ${INSTALL_DEMODIR}/thermal/${PROJECT_NAME}
        )

