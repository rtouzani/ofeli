/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2018 Rachid Touzani

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

           Template function for iteration (relaxed) Gauss-Seidel Method

  ==============================================================================*/

#ifndef __GS_H
#define __GS_H

#include <iostream>
using std::ostream;
using std::endl;

#include <iomanip>
using std::setw;


#include "OFELI_Config.h"
#include "linear_algebra/SpMatrix.h"
#include "linear_algebra/Vect.h"
#include "util/util.h"
#include "io/output.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file GS.h
 *  \brief Function to solve a linear system of equations
 *  using the Gauss-Seidel method.
 *
 */

/** \fn int GS(const SpMatrix<T_>& A, const Vect<T_> &b, Vect<T_> &x, real_t omega, int max_it, real_t toler, int verbose)
 *  \ingroup Solver
 *  \brief Gauss-Seidel solver function.
 *  @param [in] A Problem matrix (Instance of class SpMatrix).
 *  @param [in] b Right-hand side vector (class Vect)
 *  @param [in,out] x Vect instance containing initial solution guess in input and solution
 *  of the linear system in output (If iterations have succeeded).
 *  @param [in] omega Relaxation parameter.
 *  @param [in] max_it Maximum number of iterations.
 *  @param [in] toler Tolerance for convergence (measured in relative weighted 2-Norm).
 *  @param [in] verbose Information output parameter
 *    - 0: No output
 *    - 1: Output iteration information
 *    - 2 and greater: Output iteration information and solution at each iteration.
 *  @return Number of performed iterations
 *
 * \tparam <T_> Data type (real_t, float, complex<real_t>, ...)
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
template< class T_>
int GS(const SpMatrix<T_>& A,
       const Vect<T_>&     b,
       Vect<T_>&           x,
       real_t              omega,
       int                 max_it,
       real_t              toler,
       int                 verbose)
{
   if (verbose>0)
      cout << "Running Gauss-Seidel method ..." << endl;
   size_t i, j, l;
   int it;
   size_t size = A.getNbRows();
   real_t nrm = x.getNorm2();
   if (nrm == 0)
      nrm = 1;
   Vect<T_> y(x);

   for (i=1; i<=size; i++) {
      if (Abs(A(i,i)) < OFELI_EPSMCH) {
         cerr << "Error in function GS : " << i << "-th Diagonal";
         cerr << " entry is zero." << endl;
         exit(1);
      }
   }

   for (it=1; it<=max_it; it++) {
      l = 1;
      T_ s;
      for (i=1; i<=size; ++i) {
         s = 0;
         for (j=1; j<=size; ++j)
            s += A(i,j)*y(j);
         s -= A(i,i)*y(i);
         y(i) = (b(i)-s)/A(i,i);
      }

      real_t err = 0;
      for (i=0; i<size; i++) {
         y[i] = omega*y[i] + (1-omega)*x[i];
         err += (y[i]-x[i])*(y[i]-x[i]);
      }
      err = sqrt(err/size)/nrm;
      x = y;

      if (verbose > 1)
         cout << "Iteration: " << setw(4) << it << "  ... Error: " << err << endl;

      nrm = x.getNorm2();
      if (verbose > 2)
         cout << x;
      if (err < toler) {
         toler = err;
         if (verbose)
            cout << "Convergence of the GS method after " << it << " iterations." << endl;
         return it;
      }
   }
   if (verbose)
      cout << "No Convergence after " << it << " iterations." << endl;
   return -it;
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
