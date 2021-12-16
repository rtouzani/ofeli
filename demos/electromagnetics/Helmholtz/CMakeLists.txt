project(Helmholtz)

add_executable(helmholtz main.cpp)
target_link_libraries(helmholtz ofeli)

file(COPY proj.dat proj.m DESTINATION .)

install(TARGETS helmholtz RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})