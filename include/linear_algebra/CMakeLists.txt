# ==============================================================================
#
#                                    O  F  E  L  I
#
#                           Object  Finite  Element  Library
#
# ==============================================================================
#
#   Copyright (C) 1998 - 2025 Rachid Touzani
#
#   This file is part of OFELI.
#
#   OFELI is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   OFELI is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with OFELI. If not, see <http://www.gnu.org/licenses/>.
#
# ==============================================================================

file (COPY Assembly.h
           Assembly_impl.h
           BMatrix.h
           BMatrix_impl.h
           DMatrix.h
           DMatrix_impl.h
           DSMatrix.h
           DSMatrix_impl.h
           GraphOfMatrix.h
           LocalMatrix.h
           LocalMatrix_impl.h
           LocalVect.h
           LocalVect_impl.h
           Matrix.h
           Matrix_impl.h
           Point.h
           Point2D.h
           SkMatrix.h
           SkMatrix_impl.h
           SkSMatrix.h
           SkSMatrix_impl.h
           SpMatrix.h
           SpMatrix_impl.h
           SpaceTime.h
           TrMatrix.h
           TrMatrix_impl.h
           Vect.h
           Vect_impl.h
      DESTINATION .
     )

add_subdirectory (petsc)
