/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani

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

                          Definition of class 'LPSolver'

  ==============================================================================*/

#ifndef __LP_SOLVER_H
#define __LP_SOLVER_H

#include <stdlib.h>
#include <math.h>

#include <iostream>
using std::ostream;
using std::endl;
using std::cerr;

#include <iomanip>
using std::setw;

#include <string>
using std::string;

#include "OFELI_Config.h"
#include "solvers/simplex.h"

namespace OFELI {

/*! \file LPSolver.h
 *  \brief Definition file for class LPSolver.
 */

//-----------------------------------------------------------------------------
// Class LPSolver
//-----------------------------------------------------------------------------

/*! \class LPSolver
 * \ingroup Solver
 * \brief To solve a linear programming problem
 * \details The Linear Program reads:
 *
 *     Minimise: d(1)*x(1) + ... + d(n)*x(n) + e
 *     Subject to the constraints:
 *            A(i,1)*x(1) + ... + A(i,n)*x(n) <= a(i)  i=1,...,n_le
 *            B(i,1)*x(1) + ... + B(i,n)*x(n) >= b(i)  i=1,...,n_ge
 *            C(i,1)*x(1) + ... + C(i,n)*x(n) =  c(i)  i=1,...,n_eq
 *            x(i) >= 0, 1<=i<=n
 *
 * Solution is held by the Simplex method
 * Reference: "Numerical Recipes By W.H. Press, B. P. Flannery,
 *             S.A. Teukolsky and W.T. Vetterling, Cambridge
 *             University Press, 1986"
 *
 * C-implementation copied from J-P Moreau, Paris
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
class LPSolver
{

 public:

/*! \enum Setting
 * Selects setting option: Objective or Constraints
 */
enum Setting {
   OBJECTIVE     =  0,    /*!< Objective function coefficients               */
   LE_CONSTRAINT =  1,    /*!< 'Less or Equal' constraint coefficients       */
   GE_CONSTRAINT =  2,    /*!< 'Greater or Equal' constraint coefficients    */
   EQ_CONSTRAINT =  3,    /*!< 'Equality' constraint coefficients            */
};

//----------------------------   BASIC OPERATIONS   ----------------------------

/// \brief Default constructor
    LPSolver();

/** \brief Constructor using Linear Program data
 *  @param [in] nv Number of optimization variables
 *  @param [in] nb_le Number of '<=' inequality constraints
 *  @param [in] nb_ge Number of '>=' inequality constraints
 *  @param [in] nb_eq Number of '=' equality constraints
 */
    LPSolver(int nv,
             int nb_le,
             int nb_ge,
             int nb_eq);

/// \brief Destructor
    ~LPSolver();

//-------------------------------   MODIFIERS  ---------------------------------

/** \brief Set optimization parameters
 *  @param [in] nv Number of optimization variables
 *  @param [in] nb_le Number of '<=' inequality constraints
 *  @param [in] nb_ge Number of '>=' inequality constraints
 *  @param [in] nb_eq Number of '=' equality constraints
 */
    void setSize(int nv,
                 int nb_le,
                 int nb_ge,
                 int nb_eq);

/** \brief vector of optimization variables
 *  @param [in] x Vector of optimization variables. Its size must be at least equal
 *              to number of optimization variables
 */
    void set(Vect<real_t>& x);

/** \brief Set optimization data 
 *  \details This function enables providing all optimization data.
 *  It has to be used for the objectice
 *  function and once for each constraint. 
 *  @param [in] opt Option for data, to choose among enumerated values:
 *              <ul>
 *                <li> OBJECTIVE To set objective function to minimize
 *                <li> LE_CONSTRAINT To set a '<=' inequality constraint
 *                <li> GE_CONSTRAINT To set a '>=' inequality constraint
 *                <li> EQ_CONSTRAINT To set an equality constraint
 *              </ul>
 *  @param [in] a Vector coefficients if the chosen function. If \c opt==OBJECTIVE,
 *  vector components are the coefficients multiplying the variables in the objective function.
 *  if \c xx_CONSTRAINT, vector components are the coefficients multiplying the variables
 *  in the corresponding constraint.
 *  @param [in] b Constant value in the objective function or in a constraint. Its default value
 *  is 0.0 
 */
    void set(Setting             opt,
             const Vect<real_t>& a,
             real_t              b=0.0);

/** \brief Run the linear program solver
 *  \details This function runs the linear programming solver using the Simplex
 *  algorithm
 *  @return 0 if process is complete, >0 otherwise 
 */
    int run();

/// \brief Return objective
/// \details Once execution is complete, this function returns optimal value of objective 
    real_t getObjective() const { return _theSimplex.getObjective(); }

/// Output class information
    friend ostream & operator<<(ostream&        s,
                                const LPSolver& os);

 private:

   Vect<real_t> _A, *_x;
   simplex _theSimplex;
   int _it, _nv, _nb_le, _nb_ge, _nb_eq, _nc, _i_le, _i_ge, _i_eq;
   real_t _ret;
};

/// \fn ostream & operator<<(ostream& s, const LPSolver &os)
/// \brief Output solver information
/// \ingroup Solver
    ostream & operator<<(ostream&        s,
                         const LPSolver& os);

} /* namespace OFELI */

#endif
