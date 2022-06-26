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
#   Provide a 'make uninstall' tool to uninstall the package
#
# ==============================================================================

if (NOT EXISTS "/Users/touzani/ofeli/build/install_manifest.txt")
  message (FATAL_ERROR "Cannot find install manifest: /Users/touzani/ofeli/build/install_manifest.txt")
endif()

file (READ "/Users/touzani/ofeli/build/install_manifest.txt" files)
string (REGEX REPLACE "\n" ";" files "${files}")
foreach (file ${files})
  message (STATUS "Uninstalling $ENV{DESTDIR}${file}")
  if (IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
     exec_program ("/opt/local/bin/cmake" ARGS "-E remove \"$ENV{DESTDIR}${file}\""
                   OUTPUT_VARIABLE rm_out
                   RETURN_VALUE rm_retval
                  )
    if (NOT "${rm_retval}" STREQUAL 0)
       message (FATAL_ERROR "Problem when removing $ENV{DESTDIR}${file}")
    endif ()
    else (IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
      message (STATUS "File $ENV{DESTDIR}${file} does not exist.")
    endif ()
endforeach ()
