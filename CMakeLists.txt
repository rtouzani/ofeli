cmake_minimum_required(VERSION 3.16.0)
project(
  Ofeli
  VERSION 0.0.0
  DESCRIPTION "An Object Finite Element Library"
  LANGUAGES CXX)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_INSTALL_PREFIX "./")
set(CMAKE_INSTALL_LIBDIR "/home/martin/Documents/github/ofeli")
set(CMAKE_INSTALL_INCLUDEDIR "/home/martin/Documents/github/ofeli/includers")

option(BUILD_TESTS "Build test programs" OFF)

add_definitions(-Wall)

include_directories(include)

add_subdirectory(src)
add_subdirectory(demos)
add_subdirectory(util)

install(TARGETS ofeli LIBRARY DESTINATION CMAKE_INSTALL_LIBDIR)
install(TARGETS ofeli PUBLIC_HEADER DESTINATION CMAKE_INSTALL_INCLUDEDIR)