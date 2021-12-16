project(LS)

add_executable(ls_demo1 ls_demo1.cpp)
target_link_libraries(ls_demo1 ofeli)

add_executable(ls_demo2 ls_demo2.cpp)
target_link_libraries(ls_demo2 ofeli)

file(COPY mesh.m DESTINATION .)

install(TARGETS ls_demo1 ls_demo2 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})