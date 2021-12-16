project(lelas2d)

add_executable(lelas2d main.cpp)
target_link_libraries(lelas2d ofeli)

file(COPY beam-10x5.m beam-20x10.m beam.dat beam.m DESTINATION .)

install(TARGETS lelas2d RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})