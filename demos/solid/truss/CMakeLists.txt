project(truss)

add_executable(truss main.cpp)
target_link_libraries(truss ofeli)

file(COPY bridge.dat HoweBridge.m PrattBridge.m tower1.dat tower1.m tower2-1.m tower2.dat tower2.m truss.dat truss.m DESTINATION .)

install(TARGETS truss RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})