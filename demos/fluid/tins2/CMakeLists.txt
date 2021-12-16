project(tins2)

add_executable(tins2 main.cpp)
target_link_libraries(tins2 ofeli)

file(COPY cavity-30.m cavity-60.dat cavity-60.m cavity.dat step.dat step.dom step.m DESTINATION .)

install(TARGETS tins2 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})