project(NLAS)

add_executable(nl_demo1 nl_demo1.cpp)
target_link_libraries(nl_demo1 ofeli)

add_executable(nl_demo2 nl_demo2.cpp)
target_link_libraries(nl_demo2 ofeli)

add_executable(nl_demo3 nl_demo3.cpp)
target_link_libraries(nl_demo3 ofeli)

add_executable(nl_demo4 nl_demo4.cpp)
target_link_libraries(nl_demo4 ofeli)

install(TARGETS nl_demo1 nl_demo2 nl_demo3 nl_demo4 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})