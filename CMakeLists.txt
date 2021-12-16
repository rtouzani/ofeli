cmake_minimum_required(VERSION 3.16.0)
project(
  Ofeli
  VERSION 0.0.0
  DESCRIPTION "An Object Finite Element Library"
  LANGUAGES CXX)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 11)

option(BUILD_TESTS "Build test programs" OFF)

add_definitions(-Wall)

include(GNUInstallDirs)

set(CMAKE_INSTALL_BINDIR ${CMAKE_INSTALL_PREFIX}/bin)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/include)
set(CMAKE_INSTALL_MATDIR ${CMAKE_INSTALL_PREFIX}/material)

include_directories(include)

add_subdirectory(src)
add_subdirectory(demos)
add_subdirectory(util)
add_subdirectory(doc)

file(COPY test_ofeli.sh DESTINATION .)

install(TARGETS ofeli LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY material/ DESTINATION ${CMAKE_INSTALL_MATDIR})