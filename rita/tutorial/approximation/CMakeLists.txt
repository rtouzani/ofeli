# ==============================================================================
#
#                                 r  i  t  a
#
#            An environment for Modelling and Numerical Simulation
#
# ==============================================================================
#
#   Copyright (C) 2021 - 2024  Rachid Touzani
#
#   This file is part of rita.
#
#   rita is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   rita is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
# ==============================================================================

project (approximation)

file (COPY example1.rita example2.rita example3.rita example4.rita example5.rita ex1.dat ex2.dat ex4.dat DESTINATION .)

add_test (approx1 ${CMAKE_RITA_EXEC} example1.rita)
add_test (approx2 ${CMAKE_RITA_EXEC} example2.rita)
add_test (approx3 ${CMAKE_RITA_EXEC} example3.rita)
add_test (approx4 ${CMAKE_RITA_EXEC} example4.rita)
add_test (approx5 ${CMAKE_RITA_EXEC} example5.rita)

install (FILES
         README.md
         example1.rita
         example2.rita
         example3.rita
         example4.rita
         example5.rita
         ex1.dat
         ex2.dat
         ex4.dat
         DESTINATION ${INSTALL_TUTORIALDIR}/${PROJECT_NAME}
        )
