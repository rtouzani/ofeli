/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

   This file is part of OFELI.

   OFELI is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   OFELI is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with OFELI. If not, see <http://www.gnu.org/licenses/>.

  ==============================================================================

              Definition of class Simplex to solve Linear Programming
                     problems by the simplex method

  ==============================================================================*/

#ifndef __SIMPLEX_H
#define __SIMPLEX_H

#include "OFELI_Config.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/DMatrix.h"

namespace OFELI {

/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Simplex.h
 *  \brief Definition file for class Simplex.
 */

/*! \class Simplex
 *  \ingroup VectMat
 * \brief To run the Simplex method.
 *
 * This class uses the Simplex method to solve linear programming problems.
 * Namely, we consider solving the linear program
 *     Min C1*x1 + ... + Cn*xn
 * under the constraints:
 *     ai1*x1 + ai2*x2 + ... + ain*xn <= bi       i=1,...,m 
 *
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */

class Simplex {

 public:

/** Constructor using Linear Program data
 *  @param [in] A Matrix defining coefficients of constraints
 *  @param [in] b Vector of constraint values
 *  @param [in] c Vector containing coefficients of linear cost functional
 *  @param [out] x Vector of obtained solution: After executing the function run
 *  this vector contains solution 
 */
    Simplex(const DMatrix<real_t>& A, Vect<real_t>& b, const Vect<real_t>& c, Vect<real_t>& x);

/** \brief Run the Simplex algorithm
 *  @return Number of performed iterations to run the Simplex method
 */
    int run();

/// \brief Return obtained cost
    real_t getCost() const { return _max; }

 private:
    size_t _nr, _nc, _nb_it;
    DMatrix<real_t> _A;
    Vect<real_t>    *_b, _c, *_x;
    real_t _max;
    bool _isUnbounded;

/// if the table has further negative constraints,then it is not optimal
    bool checkOptimality();
    void doPivoting(size_t pr, size_t pc);
    size_t getPivotRow(size_t pivotColumn);
    size_t getPivotColumn();
    bool solve();
};
 
/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
