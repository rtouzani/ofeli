project(stdc2)

add_executable(stdc2 main.cpp)
target_link_libraries(stdc2 ofeli)

file(COPY proj.dat proj.m proj.v DESTINATION .)

add_test (stdc2 stdc2 proj.dat)

install (FILES README.md
               main.cpp
               proj.dat
               proj.m
               proj.v
         DESTINATION ${INSTALL_DEMODIR}/thermal/${PROJECT_NAME}
        )

