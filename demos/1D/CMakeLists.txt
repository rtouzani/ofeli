project(1D)

add_executable(elliptic elliptic.cpp)
target_link_libraries(elliptic ofeli)

add_executable(heat heat.cpp)
target_link_libraries(heat ofeli)

add_executable(transport transport.cpp)
target_link_libraries(transport ofeli)