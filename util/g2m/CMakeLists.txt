project(g2m)

add_executable(g2m main.cpp)
target_link_libraries(g2m ofeli)

install(TARGETS g2m RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})