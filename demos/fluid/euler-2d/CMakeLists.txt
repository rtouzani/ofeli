project(euler-2d)

add_executable(euler2d main.cpp)
target_link_libraries(euler2d ofeli)

file(COPY proj.dat rectangle.m DESTINATION .)

install(TARGETS euler2d RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})