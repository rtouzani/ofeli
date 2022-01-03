/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani

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
                   Template function for SSOR iteration Method
  ==============================================================================*/

#ifndef __SSOR_H
#define __SSOR_H

#include <iostream>
using std::ostream;
using std::endl;

#include <iomanip>
using std::setw;


#include "OFELI_Config.h"
#include "linear_algebra/Vect_impl.h"
#include "util/util.h"
#include "io/output.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file SSOR.h
 *  \brief Function to solve a linear system of equations
 *  using the Symmetric Successive Over Relaxation method.
 *
 */

/** \fn int SSOR(const M_ &A, const Vect<T_> &b, Vect<T_> &x, int max_it, real_t toler)
 *  \ingroup Solver
 *  \brief SSOR solver function.
 *  @param [in] A Problem matrix (Instance of abstract class \b M_).
 *  @param [in] b Right-hand side vector (class Vect)
 *  @param [in,out] x Vect instance containing initial solution guess in input and solution
 *  of the linear system in output (If iterations have succeeded).
 *  @param [in] max_it Maximum number of iterations.
 *  @param [in] toler Tolerance for convergence (measured in relative weighted 2-Norm).
 *  @return Number of performed iterations,
 *
 * \b Template \b Arguments:
 *
 *  \arg \a T_ data type (double, float, ...)\n
 *  \arg \a M_ %Matrix storage class\n
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template< class T_, class M_>
int SSOR(const M_&       A,
         const Vect<T_>& b,
         Vect<T_>&       x,
         int             max_it,
         real_t          toler)
{
   size_t i, l, k;
   int it;
   size_t size=x.Size();
   real_t nrm=x.Norm2();
   if (nrm == 0)
      nrm = 1;
   Vect<T_> y(x);

   for (i=1; i<=size; i++) {
      if (Abs(A(i,i)) < OFELI_EPSMCH) {
         cerr << "Error in function SSOR: " << i << "-th Diagonal entry is null." << endl;
         exit(1);
      }
   }

   for (it=1; it<=max_it; it++) {
      l = A.Length() - 1;
      for (i=size; i>=1; i--) {
         T_ s = 0;
         for (k=0; k<A.RowPtr(i+1)-A.RowPtr(i); k++) {
            if (A.ColInd(l+1) != i)
               s += A[l] * y(A.ColInd(l+1));
            l--;
         }
         y(i) = (b(i)-s)/A(i,i);
      }

      real_t err=0;
      for (i=1; i<=size; i++)
         err += (x(i)-y(i))*(x(i)-y(i));
      err = sqrt(err/size)/nrm;
      x = y;

      if (Verbosity > 1)
         cout << "Iteration: " << setw(4) << it << "  ... Error: " << err << endl;

      nrm = x.Norm2();
      if (Verbosity > 2)
         cout << x;
      if (err < toler) {
         if (Verbosity)
            cout << "Convergence of the SSOR method after " << it << " iterations." << endl;
         return it;
      }
   }
   if (Verbosity)
      cout << "No Convergence of the SSOR method after " << it << " iterations." << endl;
   return -it;
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
