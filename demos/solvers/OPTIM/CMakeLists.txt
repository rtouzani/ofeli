project(OPTIM)

add_executable(opt_demo1 opt_demo1.cpp)
target_link_libraries(opt_demo1 ofeli)

add_executable(opt_demo2 opt_demo2.cpp)
target_link_libraries(opt_demo2 ofeli)

add_executable(opt_demo3 opt_demo3.cpp)
target_link_libraries(opt_demo3 ofeli)

add_executable(opt_demo4 opt_demo4.cpp)
target_link_libraries(opt_demo4 ofeli)

add_executable(opt_demo5 opt_demo5.cpp)
target_link_libraries(opt_demo5 ofeli)

add_executable(opt_demo6 opt_demo6.cpp)
target_link_libraries(opt_demo6 ofeli)

file(COPY test.dat test.m DESTINATION .)

install(TARGETS opt_demo1 opt_demo2 opt_demo3 opt_demo4 opt_demo5 opt_demo6 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})
install(FILES Opt3.h Opt4.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})