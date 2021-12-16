project(vmesh)

add_executable(vmesh main.cpp)
target_link_libraries(vmesh ofeli)

install(TARGETS vmesh RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})