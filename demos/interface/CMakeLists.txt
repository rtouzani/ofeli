project(interface)

add_executable(fmm_demo fmm_demo.cpp)
target_link_libraries(fmm_demo ofeli)

install(TARGETS fmm_demo RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})