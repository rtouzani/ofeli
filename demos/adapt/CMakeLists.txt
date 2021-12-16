project(adapt)

add_executable(ad1 adapt_demo1.cpp)
target_link_libraries(ad1 ofeli)

add_executable(ad2 adapt_demo2.cpp)
target_link_libraries(ad2 ofeli)

file(COPY L.dom rect.dom DESTINATION .)

install(TARGETS ad1 ad2 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})