project(EIGEN)

add_executable(eigen_demo1 eigen_demo1.cpp)
target_link_libraries(eigen_demo1 ofeli)

add_executable(eigen_demo2 eigen_demo2.cpp)
target_link_libraries(eigen_demo2 ofeli)