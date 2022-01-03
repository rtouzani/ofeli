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

#
# Configure config.h
configure_file (
                ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
                ${CMAKE_CURRENT_BINARY_DIR}/config.h
               )

#
# Add a rule for how to install the header files:
#install (DIRECTORY ofeli
#         DESTINATION ${INSTALL_INCLUDE_DIR}
#         COMPONENT library
#         FILES_MATCHING PATTERN "*.h"
#        )

# Install all generated header files, too:
#
#install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ofeli
#         DESTINATION ${OFELI_INCLUDE_DIR}
#         COMPONENT library
#         FILES_MATCHING PATTERN "*.h"
#        )
