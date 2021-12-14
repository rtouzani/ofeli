project(ofeli)

add_library(${PROJECT_NAME})
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "../include/")

add_subdirectory(equations)
add_subdirectory(io)
add_subdirectory(linear_algebra)
add_subdirectory(mesh)
add_subdirectory(post)
add_subdirectory(shape_functions)
add_subdirectory(solvers)
add_subdirectory(util)

