project(std3)

add_executable(std3 main.cpp)
target_link_libraries(std3 ofeli)

file(COPY beam.dat beam.m DESTINATION .)

add_test (std3 std3 beam.dat)

install (FILES README.md
               main.cpp
               beam.dat
               beam.m
         DESTINATION ${INSTALL_DEMODIR}/thermal/${PROJECT_NAME}
        )
