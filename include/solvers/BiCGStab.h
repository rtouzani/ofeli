/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2024 Rachid Touzani

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

          Template function for BiConjugate Gradient Stabilized Method
               BiCGStab is inspired from the SIAM Templates book.

  ==============================================================================*/

#ifndef __BICG_STAB_H
#define __BICG_STAB_H

#include <iostream>
using std::ostream;
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

/*! \file BiCGStab.h
 *  \brief Solves an unsymmetric linear system of equations using the BiConjugate Gradient Stabilized method.
 *
 *  Preconditioning is possible using a preconditioning class.
 *
 */

/** \fn int BiCGStab(const SpMatrix<T_> &A, const Prec<T_> &P, const Vect<T_> &b, Vect<T_> &x, int max_it, real_t toler)
 *  \ingroup Solver
 *  \brief Biconjugate gradient stabilized solver function.
 *  \details This function uses the preconditioned Conjugate Gradient Stabilized algorithm to solve a
 *  linear system with a sparse matrix.\n
 *  The global variable Verbosity enables choosing output message level
 *  <ul>
 *    <li> Verbosity < 2 : No output message
 *    <li> Verbosity > 1 : Notify executing the function BiCGStab
 *    <li> Verbosity > 2 : Notify convergence with number of performed iterations or divergence
 *    <li> Verbosity > 3 : Output each iteration number and residual
 *    <li> Verbosity > 6 : Print final solution if convergence
 *    <li> Verbosity > 10 : Print obtained solution at each iteration
 *  </ul>
 *  @param [in] A Problem matrix (Instance of class SpMatrix).
 *  @param [in] P Preconditioner (Instance of class Prec).
 *  @param [in] b Right-hand side vector (class Vect)
 *  @param [in,out] x Vect instance containing initial solution guess on input and solution
 *  of the linear system on output (If iterations have succeeded).
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
int BiCGStab(const SpMatrix<T_>& A,
             const Prec<T_>&     P,
             const Vect<T_>&     b,
             Vect<T_>&           x,
             int                 max_it,
             real_t              toler)
{
   if (Verbosity>2)
      cout << "Running preconditioned BiCGStab method ..." << endl;
   int it;
   size_t size=b.size();
   T_ rho_1=0, rho_2=0, alpha=0, beta=0, omega=0;
   Vect<T_> p(size), pp(size), s(size), ss(size), t(size), v(size), r(size), rr(size);
   real_t res, nrm=b.getNorm2();
   r = b - A*x;
   rr = r;

   if (nrm==0.0)
      nrm = 1;
   if ((res=r.getNorm2()/nrm)<=toler) {
      if (Verbosity>2)
         cout << "Convergence of the BiCGStab method after 0 iterations." << endl;
      if (Verbosity>6)
         cout << "Solution\n" << x;
      return 0;
   }

   for (it=1; it<=max_it; it++) {
      rho_1 = (rr,r);
      if (rho_1==T_(0)) {
         toler = r.getNorm2()/nrm;
         return 2;
      }
      if (it==1)
         p = r;
      else {
         beta = (rho_1/rho_2) * (alpha/omega);
         p = r + beta*(p - omega*v);
      }
      P.solve(p,pp);
      v = A*pp;
      alpha = rho_1/(rr,v);
      s = r - alpha*v;
      if ((res=s.getNorm2()/nrm)<toler) {
         x += alpha*pp;
         if (Verbosity>3)
            cout << "Iteration: " << setw(4) << it << "  ... Residual: " << res << endl;
         toler = res;
         return it;
      }
      P.solve(s,ss);
      t = A*ss;
      omega = (t,s)/(t,t);
      x += alpha*pp + omega*ss;
      r = s - omega*t;
      rho_2 = rho_1;
      if (Verbosity>3)
         cout << "Iteration: " << setw(4) << it << "  ... Residual: " << res << endl;
      if (Verbosity>10)
         cout << "Solution at iteration " << it << ": \n" << x;
      if ((res=r.getNorm2()/nrm)<toler) {
         toler = res;
         if (Verbosity>2)
            cout << "Convergence of the BiCGStab method after " << it << " iterations." << endl;
         return it;
      }
      if (omega==T_(0)) {
         toler = r.getNorm2()/nrm;
         return 3;
      }
   }
   if (Verbosity>2)
      cout << "No Convergence after " << it << " iterations." << endl;
   return -it;
}


/** \fn int BiCGStab(const SpMatrix<T_> &A, int prec, const Vect<T_> &b, Vect<T_> &x, int max_it, real_t toler)
 *  \ingroup Solver
 *  \brief Biconjugate gradient stabilized solver function.
 *  \details This function uses the preconditioned Conjugate Gradient Stabilized algorithm to solve a
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
int BiCGStab(const SpMatrix<T_>& A,
             int                 prec,
             const Vect<T_>&     b,
             Vect<T_>&           x,
             int                 max_it,
             real_t              toler)
{
   return BiCGStab(A,Prec<T_>(A,prec),b,x,max_it,toler);
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
int BiCGStab(const Matrix<T_>* A,
             int               prec,
             const Vect<T_>&   b,
             Vect<T_>&         x,
             int               max_it,
             real_t            toler)
{
   SpMatrix<T_> &AA = MAT(SpMatrix<T_>,A);
   return BiCGStab(AA,prec,b,x,max_it,toler);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */

#endif
