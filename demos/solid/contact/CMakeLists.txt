project(contact)

add_executable(contact main.cpp)
target_link_libraries(contact ofeli)

file(COPY beam.dat beam.m MyMaterial.md DESTINATION .)

install(TARGETS contact RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})