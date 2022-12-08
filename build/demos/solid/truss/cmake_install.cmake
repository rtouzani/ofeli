# Install script for directory: /Users/touzani/ofeli/demos/solid/truss

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
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/ofeli/demos/solid/truss/README.md;/usr/local/share/ofeli/demos/solid/truss/main.cpp;/usr/local/share/ofeli/demos/solid/truss/bridge.dat;/usr/local/share/ofeli/demos/solid/truss/HoweBridge.m;/usr/local/share/ofeli/demos/solid/truss/PrattBridge.m;/usr/local/share/ofeli/demos/solid/truss/tower1.dat;/usr/local/share/ofeli/demos/solid/truss/tower1.m;/usr/local/share/ofeli/demos/solid/truss/tower2.dat;/usr/local/share/ofeli/demos/solid/truss/tower2.m;/usr/local/share/ofeli/demos/solid/truss/truss.dat;/usr/local/share/ofeli/demos/solid/truss/truss.m")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/usr/local/share/ofeli/demos/solid/truss" TYPE FILE FILES
    "/Users/touzani/ofeli/demos/solid/truss/README.md"
    "/Users/touzani/ofeli/demos/solid/truss/main.cpp"
    "/Users/touzani/ofeli/demos/solid/truss/bridge.dat"
    "/Users/touzani/ofeli/demos/solid/truss/HoweBridge.m"
    "/Users/touzani/ofeli/demos/solid/truss/PrattBridge.m"
    "/Users/touzani/ofeli/demos/solid/truss/tower1.dat"
    "/Users/touzani/ofeli/demos/solid/truss/tower1.m"
    "/Users/touzani/ofeli/demos/solid/truss/tower2.dat"
    "/Users/touzani/ofeli/demos/solid/truss/tower2.m"
    "/Users/touzani/ofeli/demos/solid/truss/truss.dat"
    "/Users/touzani/ofeli/demos/solid/truss/truss.m"
    )
endif()

