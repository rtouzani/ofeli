project(laplace)

add_executable(laplace_demo1 laplace_demo1.cpp)
target_link_libraries(laplace_demo1 ofeli)

add_executable(laplace_demo2 laplace_demo2.cpp)
target_link_libraries(laplace_demo2 ofeli)

add_executable(laplace_demo3 laplace_demo3.cpp)
target_link_libraries(laplace_demo3 ofeli)

add_executable(laplace_demo4 laplace_demo4.cpp)
target_link_libraries(laplace_demo4 ofeli)

file(COPY annulus.geo annulus.m cube.geo cube.m proj.dat proj2.dat proj3.dat square.geo square.m DESTINATION .)

install(TARGETS laplace_demo1 laplace_demo2 laplace_demo3 laplace_demo4 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})