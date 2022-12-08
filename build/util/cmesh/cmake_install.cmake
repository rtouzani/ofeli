# Install script for directory: /Users/touzani/ofeli/util/cmesh

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Library/Developer/CommandLineTools/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/touzani/ofeli/build/util/cmesh/cmesh")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cmesh" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cmesh")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cmesh")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/ofeli/util/cmesh/README.md;/usr/local/share/ofeli/util/cmesh/A.ele;/usr/local/share/ofeli/util/cmesh/A.node;/usr/local/share/ofeli/util/cmesh/beam.m;/usr/local/share/ofeli/util/cmesh/cube.ele;/usr/local/share/ofeli/util/cmesh/cube.node;/usr/local/share/ofeli/util/cmesh/cube.vol;/usr/local/share/ofeli/util/cmesh/disk.bamg;/usr/local/share/ofeli/util/cmesh/disk.m;/usr/local/share/ofeli/util/cmesh/gear.ele;/usr/local/share/ofeli/util/cmesh/gear.node;/usr/local/share/ofeli/util/cmesh/l.bamg;/usr/local/share/ofeli/util/cmesh/torus.neu")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/usr/local/share/ofeli/util/cmesh" TYPE FILE FILES
    "/Users/touzani/ofeli/util/cmesh/README.md"
    "/Users/touzani/ofeli/util/cmesh/A.ele"
    "/Users/touzani/ofeli/util/cmesh/A.node"
    "/Users/touzani/ofeli/util/cmesh/beam.m"
    "/Users/touzani/ofeli/util/cmesh/cube.ele"
    "/Users/touzani/ofeli/util/cmesh/cube.node"
    "/Users/touzani/ofeli/util/cmesh/cube.vol"
    "/Users/touzani/ofeli/util/cmesh/disk.bamg"
    "/Users/touzani/ofeli/util/cmesh/disk.m"
    "/Users/touzani/ofeli/util/cmesh/gear.ele"
    "/Users/touzani/ofeli/util/cmesh/gear.node"
    "/Users/touzani/ofeli/util/cmesh/l.bamg"
    "/Users/touzani/ofeli/util/cmesh/torus.neu"
    )
endif()

