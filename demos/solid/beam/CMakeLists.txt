project(beam)

add_executable(beam main.cpp)
target_link_libraries(beam ofeli)

file(COPY beam.dat DESTINATION .)

install(TARGETS beam RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})