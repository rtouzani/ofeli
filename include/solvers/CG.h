/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzani

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

                 Template function for Conjugate Gradient Method
                   CG is inspired from the SIAM Templates book

  ==============================================================================*/

#ifndef __CG_H
#define __CG_H

#include <iostream>
using std::ostream;
using std::cout;
using std::endl;

#include <iomanip>
using std::setw;
using std::ios;
using std::setprecision;


#include "OFELI_Config.h"
#include "linear_algebra/SpMatrix.h"
#include "linear_algebra/Vect.h"
#include "solvers/Prec.h"
#include "util/util.h"
#include "io/output.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file CG.h
 *  \brief Functions to solve a symmetric positive definite linear system of equations
 *  using the Conjugate Gradient method.
 *
 *  Preconditioning is possible using a preconditioning class.
 *
 */

template<class T_> class SpMatrix;
template<class T_> class Prec;

/** \fn int CG(const SpMatrix<T_> &A, const Prec<T_> &P, const Vect<T_> &b, Vect<T_> &x, int max_it, real_t toler, int verbose)
 *  \ingroup Solver
 *  \brief Conjugate gradient solver function.
 *  @param [in] A Problem matrix (Instance of class SpMatrix).
 *  @param [in] P Preconditioner (Instance of class Prec).
 *  @param [in] b Right-hand side vector (class Vect)
 *  @param [in,out] x Vect instance containing initial solution guess in input and solution
 *  of the linear system in output (If iterations have succeeded).
 *  @param [in] max_it Maximum number of iterations.
 *  @param [in] toler Tolerance for convergence (measured in relative weighted 2-Norm).
 *  @param [in] verbose Information output parameter
 *    - 0: No output
 *    - 1: Output iteration information,
 *    - 2 and greater: Output iteration information and solution at each iteration.
 *  @return Number of performed iterations,
 *
 * \tparam <T_> Data type (double, float, complex<double>, ...)
 */
template<class T_>
int CG(const SpMatrix<T_>& A,
       const Prec<T_>&     P,
       const Vect<T_>&     b,
             Vect<T_>&     x,
             int           max_it,
             real_t        toler,
             int           verbose)
{
   if (verbose>0)
      cout << "Running preconditioned CG method ..." << endl;
   size_t size=x.size();
   real_t res, nrm=b.getNorm2();
   T_ rho=0, rho_1=T_(1), beta=0;
   if (nrm==0)
      nrm = 1;
   Vect<T_> r(size), p(size), q(size), z(size);
   r = b - A*x;
   if ((res=r.getNorm2()/nrm) <= toler) {
      if (verbose>1)
         cout << "Convergence after 0 iterations." << endl;
      return 0;
   }

   int it=0;
   for (it=1; it<=max_it; it++) {
      P.solve(r,z);
      rho = (r,z);
      if (it==1)
         p = z;
      else {
         beta = rho/rho_1;
         p = z + beta*p;
      }
      q = A*p;
      T_ alpha = rho/(p,q);
      x += alpha * p;
      r -= alpha * q;
      res = r.getNorm2()/nrm;

      if (verbose>1) {
         cout << "Iteration: " << setw(4) << it;
         cout << "  ... Residual: " << res << endl;
      }
      if (verbose>2)
         cout << x;
      if (res<=toler) {
         toler = res;
         if (verbose)
            cout << "Convergence after " << it << " iterations." << endl;
         return it;
      }
      rho_1 = rho;
   }
   if (verbose)
      cout << "No Convergence after " << it << " iterations." << endl;
   return -it;
}


/** \fn int CG(const SpMatrix<T_> &A, int prec, const Vect<T_> &b, Vect<T_> &x, int max_it, real_t toler, int verbose)
 *  \ingroup Solver
 *  \brief Conjugate gradient solver function.
 *  @param [in] A Problem matrix (Instance of abstract class SpMatrix).
 *  @param [in] prec Enum variable selecting a preconditioner, among the values <tt>IDENT_PREC</tt>,
 *  <tt>DIAG_PREC</tt>, <tt>ILU_PREC</tt> or <tt>SSOR_PREC</tt>
 *  @param [in] b Right-hand side vector (class Vect)
 *  @param [in,out] x Vect instance containing initial solution guess in input and solution
 *  of the linear system in output (If iterations have succeeded).
 *  @param [in] max_it Maximum number of iterations.
 *  @param [in] toler Tolerance for convergence (measured in relative weighted 2-Norm).
 *  @param [in] verbose Information output parameter
 *    - 0: No output
 *    - 1: Output iteration information,
 *    - 2 and greater: Output iteration information and solution at each iteration.
 *  @return Number of performed iterations,
 *
 * \tparam <T_> Data type (double, float, complex<double>, ...)
 */
template<class T_>
int CG(const SpMatrix<T_>& A,
             int           prec,
       const Vect<T_>&     b,
             Vect<T_>&     x,
             int           max_it,
             real_t        toler,
             int           verbose)
{
   return CG(A,Prec<T_>(A,prec),b,x,max_it,toler,verbose);
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
