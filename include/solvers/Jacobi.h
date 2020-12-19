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

             Template function for iteration (relaxed) Jacobi Method

  ==============================================================================*/

#ifndef __JACOBI_H
#define __JACOBI_H

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

/*! \file Jacobi.h
 *  \brief Function to solve a linear system of equations using the Jacobi method.
 */

/** \fn int Jacobi(const SpMatrix<T_>& A, const Vect<T_> &b, Vect<T_> &x, real_t omega, int max_it, real_t toler)
 *  \ingroup Solver
 *  \brief Jacobi solver function.
 *  @param [in] A Problem matrix (Instance of class SpMatrix).
 *  @param [in] b Right-hand side vector (class Vect)
 *  @param [in,out] x Vect instance containing initial solution guess in input and solution
 *  of the linear system in output (If iterations have succeeded).
 *  @param [in] omega Relaxation parameter.
 *  @param [in] max_it Maximum number of iterations.
 *  @param [in,out] toler Tolerance for convergence (measured in relative weighted 2-Norm).
 *  @return Number of performed iterations,
 *
 * \tparam <T_> Data type (real_t, float, complex<real_t>, ...)
 * \tparam <M_> %Matrix storage class
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
template<class T_>
int Jacobi(const SpMatrix<T_>& A,
           const Vect<T_>&     b,
           Vect<T_>&           x,
           real_t              omega,
           int                 max_it,
           real_t              toler)
{
   if (Verbosity)
      cout << "Running Jacobi method ..." << endl;
   size_t size = x.Size();
   for (size_t i=1; i<=size; i++)
      if (std::abs(A(i,i)) < OFELI_EPSMCH)
         throw OFELIException("In Jacobi(A,b,x,omega,max_it,toler): null diagonal term: " + std::to_string(i));
   size_t j;
   int it;
   real_t nrm = x.Norm2();
   if (nrm == 0)
      nrm = 1;
   Vect<T_> y(size);

   A.Mult(x,y);
   y -= b;

   for (it=1; it<=max_it; it++) {
      for (size_t i=1; i<=size; i++) {
         T_ s = 0;
         for (size_t j=1; j<i; j++)
            s += A(i,j)*x(j);
         for (j=i+1; j<=size; j++)
            s += A(i,j)*x(j);
         y(i) = omega*((b(i)-s)/A(i,i) - x(i));
      }
      real_t err = y.Norm2()/nrm;

      if (Verbosity > 1)
         cout << "Iteration: " << setw(4) << it << "  ... Error: " << err << endl;

      x += y;
      nrm = x.Norm2();
      if (Verbosity>10)
         cout << x;
      if (err < toler) {
         if (Verbosity>2)
            cout << "Convergence of the Jacobi method after " << it << " iterations." << endl;
         if (Verbosity>6)
         cout << "Solution at iteration " << it << ": \n" << x;
         return it;
      }
   }
   if (Verbosity>2)
      cout << "No Convergence of the Jacobi method after " << it << " iterations." << endl;
   return -it;
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
int Jacobi(const Matrix<T_>* A,
           const Vect<T_>&   b,
           Vect<T_>&         x,
           real_t            omega,
           int               max_it,
           real_t            toler)
{
   SpMatrix<T_> &AA = MAT(SpMatrix<T_>,A);
   return Jacobi(AA,b,x,omega,max_it,toler);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
