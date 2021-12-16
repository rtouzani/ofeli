project(conv)

add_executable(cfield cfield.cpp)
target_link_libraries(cfield ofeli)

add_executable(cmesh cmesh.cpp)
target_link_libraries(cmesh ofeli)

install(TARGETS cfield cmesh RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})