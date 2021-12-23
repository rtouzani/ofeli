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
#include "solvers/Prec_impl.h"
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
template<class T_> class Vect;

/** \fn int CG(const SpMatrix<T_> &A, const Prec<T_> &P, const Vect<T_> &b, Vect<T_> &x, int max_it, real_t toler)
 *  \ingroup Solver
 *  \brief Conjugate gradient solver function.
 *  \details This function uses the preconditioned Conjugate Gradient algorithm to solve a
 *  linear system with a sparse matrix.\n
 *  The global variable Verbosity enables choosing output message level
 *  <ul>
 *    <li> Verbosity < 2 : No output message
 *    <li> Verbosity > 1 : Notify executing the function CG
 *    <li> Verbosity > 2 : Notify convergence with number of performed iterations or divergence
 *    <li> Verbosity > 3 : Output each iteration number and residual
 *    <li> Verbosity > 6 : Print final solution if convergence
 *    <li> Verbosity > 10 : Print obtained solution at each iteration
 *  </ul>
 *  @param [in] A Problem matrix (Instance of class SpMatrix).
 *  @param [in] P Preconditioner (Instance of class Prec).
 *  @param [in] b Right-hand side vector (class Vect)
 *  @param [in,out] x Vect instance containing initial solution guess in input and solution
 *  of the linear system in output (If iterations have succeeded).
 *  @param [in] max_it Maximum number of iterations.
 *  @param [in] toler Tolerance for convergence (measured in relative weighted 2-Norm).
 *  @return Number of performed iterations,
 *
 * \tparam <T_> Data type (double, float, complex<double>, ...)
 *
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
template<class T_>
int CG(const SpMatrix<T_>& A,
       const Prec<T_>&     P,
       const Vect<T_>&     b,
       Vect<T_>&           x,
       int                 max_it,
       real_t              toler)
{
   if (Verbosity>2)
      cout << "Running preconditioned CG method ..." << endl;
   size_t size=x.size();
   real_t res, nrm=b.getNorm2();
   T_ rho=0, rho_1=T_(1), beta=0;
   if (nrm==0)
      nrm = 1;
   Vect<T_> r(size), p(size), q(size), z(size);
   r = b - A*x;
   if ((res=r.getNorm2()/nrm) <= toler) {
      if (Verbosity>2)
         cout << "Convergence after 0 iterations." << endl;
      if (Verbosity>6)
         cout << "Solution:\n" << x;
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

      if (Verbosity>3)
         cout << "Iteration: " << setw(4) << it << ", ... Residual: " << res << endl;
      if (Verbosity>10)
         cout << "Solution at iteration " << it << ": \n" << x;
      if (res<=toler) {
         toler = res;
         if (Verbosity>2)
            cout << "Convergence after " << it << " iterations." << endl;
         if (Verbosity>6)
            cout << "Solution:\n" << x;
         return it;
      }
      rho_1 = rho;
   }
   if (Verbosity>2)
      cout << "No Convergence after " << it << " iterations." << endl;
   return -it;
}


/** \fn int CG(const SpMatrix<T_> &A, int prec, const Vect<T_> &b, Vect<T_> &x, int max_it, real_t toler)
 *  \ingroup Solver
 *  \brief Conjugate gradient solver function.
 *  \details This function uses the preconditioned Conjugate Gradient algorithm to solve a
 *  linear system with a sparse matrix.\n
 *  The global variable Verbosity enables choosing output message level
 *  <ul>
 *    <li> Verbosity < 2 : No output message
 *    <li> Verbosity > 1 : Notify executing the function CG
 *    <li> Verbosity > 2 : Notify convergence with number of performed iterations or divergence
 *    <li> Verbosity > 3 : Output each iteration number and residual
 *    <li> Verbosity > 6 : Print final solution if convergence
 *    <li> Verbosity > 10 : Print obtained solution at each iteration
 *  </ul>
 *  @param [in] A Problem matrix (Instance of abstract class SpMatrix).
 *  @param [in] prec Enum variable selecting a preconditioner, among the values <tt>IDENT_PREC</tt>,
 *  <tt>DIAG_PREC</tt>, <tt>ILU_PREC</tt> or <tt>SSOR_PREC</tt>
 *  @param [in] b Right-hand side vector (class Vect)
 *  @param [in,out] x Vect instance containing initial solution guess in input and solution
 *  of the linear system in output (If iterations have succeeded).
 *  @param [in] max_it Maximum number of iterations.
 *  @param [in] toler Tolerance for convergence (measured in relative weighted 2-Norm).
 *  @return Number of performed iterations,
 *
 * \tparam <T_> Data type (double, float, complex<double>, ...)
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
template<class T_>
int CG(const SpMatrix<T_>& A,
       int                 prec,
       const Vect<T_>&     b,
       Vect<T_>&           x,
       int                 max_it,
       real_t              toler)
{
   Prec<T_> p(A,prec);
   int nb_it = CG(A,p,b,x,max_it,toler);
   return nb_it;
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
int CG(const Matrix<T_>* A,
       int               prec,
       const Vect<T_>&   b,
       Vect<T_>&         x,
       int               max_it,
       real_t            toler)
{
   SpMatrix<T_> &AA = MAT(SpMatrix<T_>,A);
   return CG(AA,prec,b,x,max_it,toler);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
