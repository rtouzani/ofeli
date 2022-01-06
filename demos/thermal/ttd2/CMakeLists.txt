project(ttd2)

add_executable(ttd2 main.cpp)
target_link_libraries(ttd2 ofeli)

target_sources(ttd2 PRIVATE error.cpp)

file(COPY proj-5x10.m proj-5x20.m proj-5x40.m proj-5x80.m proj.dat DESTINATION .)

add_test (ttd2 ttd2 proj.dat)

install (FILES README.md
               main.cpp
               error.cpp
               proj-5x10.m
               proj-5x20.m
               proj-5x40.m
               proj-5x80.m
               proj.dat
         DESTINATION ${INSTALL_DEMODIR}/thermal/${PROJECT_NAME}
        )

