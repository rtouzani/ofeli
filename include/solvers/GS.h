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
#include "linear_algebra/SpMatrix_impl.h"
#include "linear_algebra/Vect_impl.h"
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

/** \fn int GS(const SpMatrix<T_>& A, const Vect<T_> &b, Vect<T_> &x, real_t omega, int max_it, real_t toler)
 *  \ingroup Solver
 *  \details This function uses the relaxed Gauss-Seidel algorithm to solve a
 *  linear system with a sparse matrix.\n
 *  The global variable Verbosity enables choosing output message level
 *  <ul>
 *    <li> Verbosity < 2 : No output message
 *    <li> Verbosity > 1 : Notify executing the function GS
 *    <li> Verbosity > 2 : Notify convergence with number of performed iterations or divergence
 *    <li> Verbosity > 3 : Output each iteration number and residual
 *    <li> Verbosity > 6 : Print final solution if convergence
 *    <li> Verbosity > 10 : Print obtained solution at each iteration
 *  </ul>
 *  \brief Gauss-Seidel solver function.
 *  @param [in] A Problem matrix (Instance of class SpMatrix).
 *  @param [in] b Right-hand side vector (class Vect)
 *  @param [in,out] x Vect instance containing initial solution guess in input and solution
 *  of the linear system in output (If iterations have succeeded).
 *  @param [in] omega Relaxation parameter.
 *  @param [in] max_it Maximum number of iterations.
 *  @param [in] toler Tolerance for convergence (measured in relative weighted 2-Norm).
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
       real_t              toler)
{
   if (Verbosity>2)
      cout << "Running Gauss-Seidel method ..." << endl;
   for (i=1; i<=size; i++)
      if (Abs(A(i,i)) < OFELI_EPSMCH)
            throw OFELIException("In GS(A,b,x,omega,max_it,toler): null diagonal term: " + itos(i));
   size_t i, j, l;
   int it;
   size_t size = A.getNbRows();
   real_t nrm = x.getNorm2();
   if (nrm == 0)
      nrm = 1;
   Vect<T_> y(x);

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

      if (Verbosity>3)
         cout << "Iteration: " << setw(4) << it << "  ... Error: " << err << endl;

      nrm = x.getNorm2();
      if (Verbosity>10)
         cout << "Solution at iteration " << it << ": \n" << x;
      if (err < toler) {
         toler = err;
         if (Verbosity>2)
            cout << "Convergence of the GS method after " << it << " iterations." << endl;
         if (Verbosity>6)
            cout << "Solution:\n" << x;
         return it;
      }
   }
   if (Verbosity>2)
      cout << "No Convergence after " << it << " iterations." << endl;
   return -it;
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
int GS(const Matrix<T_>* A,
       const Vect<T_>&   b,
       Vect<T_>&         x,
       real_t            omega,
       int               max_it,
       real_t            toler)
{
   SpMatrix<T_> &AA = MAT(SpMatrix<T_>,A);
   return GS(AA,b,x,omega,max_it,toler);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
