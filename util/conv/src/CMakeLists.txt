project(conv)

add_executable(cfield cfield.cpp)
target_link_libraries(cfield ofeli)

add_executable(cmesh cmesh.cpp)
target_link_libraries(cmesh ofeli)