project(lh2d)

add_executable(lh2d main.cpp)
target_link_libraries(lh2d ofeli)

file(COPY proj.dat rectangle.m rectangle.v DESTINATION .)

install(TARGETS lh2d RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})