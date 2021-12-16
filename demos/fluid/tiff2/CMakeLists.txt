project(tiff2)

add_executable(tiff2 main.cpp)
target_link_libraries(tiff2 ofeli)

file(COPY cavity.dat cavity.m DESTINATION .)

install(TARGETS tiff2 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})