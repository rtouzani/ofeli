project(lelas3d)

add_executable(lelas3d main.cpp)
target_link_libraries(lelas3d ofeli)

file(COPY beam.dat beam.m DESTINATION .)

install(TARGETS lelas3d RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})