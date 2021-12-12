project(ttd2)

add_executable(ttd2 main.cpp)
target_link_libraries(ttd2 ofeli)

target_sources(ttd2 PRIVATE error.cpp)
