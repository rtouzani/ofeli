#ifndef __SIMPLEX_H
#define __SIMPLEX_H

/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

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

                          Definition of class 'simplex'

  ==============================================================================*/

#include "OFELI_Config.h"
#include "linear_algebra/Vect.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/*! \file simplex.h
 *  \brief Definition file for class simplex.
 */

namespace OFELI {

//-----------------------------------------------------------------------------
// Class simplex
//-----------------------------------------------------------------------------

/*! \class simplex
 * \ingroup Solver
 * \brief To solve a linear programming problem using the simplex method
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
class simplex {

 public:

/** \brief Default constructor
 *  \details Useful to set principal constants
 */
    simplex();

  /*---------------------------------------------------------------------------------------- 
 Simplex method for linear programming. Input parameters a, m, n, mp, np, m1, m2, and m3, 
 and output parameters a, icase, izrov, and iposv are described above (see reference). 
 Parameters: MMAX is the maximum number of constraints expected; NMAX is the maximum number 
 of variables expected; EPS is the absolute precision, which should be adjusted to the 
 scale of your variables.                                    
 -----------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------
 Initial left-hand variables. m1 type constraints are represented by having their slackv ariable 
 initially left-hand, with no artificial variable. m2 type constraints have their slack 
 variable initially left-hand, with a minus sign, and their artificial variable handled implicitly 
 during their first exchange. m3 type constraints have their artificial variable initially 
 left-hand.     
------------------------------------------------------------------------------------------------*/ 
   simplex(Vect<real_t>& A,
           int           nv,
           int           nb_le,
           int           nb_ge,
           int           nb_eq,
           Vect<real_t>& x);

/// \brief Destructor
    ~simplex();

/**
  */
    void set(Vect<real_t>& A,
             int           nv,
             int           nb_le,
             int           nb_ge,
             int           nb_eq,
             Vect<real_t>& x);

/**
  */
    int run();

/** \brief Return objective value
 *  \details Return objective value corresponding to obtained optimal value
 *  \return Objective value for reached optimal solution
 */
    real_t getObjective() const;

 private:
   int _nv, _nb_le, _nb_ge, _nb_eq, _nb, _nl2, _ret;
   real_t **_A;
   Vect<real_t> *_x;
   vector<int> _l1, _l2, _l3, _zerov, _posv;
   real_t _eps;

// Determine the maximum of elements whose index is contained in the supplied list 
// _l1, either with or without taking the absolute value, as flagged by iabf. 
   void simp1(int mm, int iabf, int &kp, real_t &bmax);

// Locate a pivot element, taking degeneracy into account. 
   void simp2(int& ip, int kp, real_t& q1);

// Matrix operations to exchange a left-hand and right-hand variable (see text). 
   void simp3(int ii, int ip, int kp);

   void setSolution();
};

} /* namespace OFELI */

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif
