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

message (STATUS "Setting up CPack")
include (InstallRequiredSystemLibraries)

set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "rita")
set (CPACK_PACKAGE_VENDOR "Rachid Touzani")
set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set (CPACK_PACKAGE_VERSION_MAJOR "${RITA_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${RITA_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${RITA_VERSION_SUBMINOR}")
set (CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
if (WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full UNIX paths properly.
  # Make sure there is at least one set of four backlashes.
  set (CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/doc/images\\\\rita.icns")
  set (CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\rita")
  set (CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} A PDE Solver")
  set (CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.ofeli.org/rita")
  set (CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.rachidtouzani.com")
  set (CPACK_NSIS_CONTACT "rachid.touzani@uca.fr")
  set (CPACK_NSIS_MODIFY_PATH ON)
else ()
  set (CPACK_STRIP_FILES "bin/rita")
  set (CPACK_SOURCE_STRIP_FILES "")
endif ()
set (CPACK_PACKAGE_EXECUTABLES "rita")
include (CPack)
message (STATUS "Setting up CPack - Done")
