project(ODE)

add_executable(ode_demo1 ode_demo1.cpp)
target_link_libraries(ode_demo1 ofeli)

add_executable(ode_demo2 ode_demo2.cpp)
target_link_libraries(ode_demo2 ofeli)

add_executable(ode_demo3 ode_demo3.cpp)
target_link_libraries(ode_demo3 ofeli)

add_executable(ode_demo4 ode_demo4.cpp)
target_link_libraries(ode_demo4 ofeli)

add_executable(ode_demo5 ode_demo5.cpp)
target_link_libraries(ode_demo5 ofeli)

add_executable(ode_demo6 ode_demo6.cpp)
target_link_libraries(ode_demo6 ofeli)

add_executable(ode_demo7 ode_demo7.cpp)
target_link_libraries(ode_demo7 ofeli)

add_executable(ts_demo1 ts_demo1.cpp)
target_link_libraries(ts_demo1 ofeli)

add_executable(ts_demo2 ts_demo2.cpp)
target_link_libraries(ts_demo2 ofeli)

file(COPY mesh1.m mesh2.m DESTINATION .)

install(TARGETS ode_demo1 ode_demo2 ode_demo3 ode_demo4 ode_demo5 ode_demo6 ode_demo7 ts_demo1 ts_demo2 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})