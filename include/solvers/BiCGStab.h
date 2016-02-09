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
#include "linear_algebra/SpMatrix.h"
#include "linear_algebra/Vect.h"
#include "solvers/Prec.h"
#include "util/util.h"
#include "io/output.h"

namespace OFELI {

/*! \file BiCGStab.h
 *  \brief Solves an unsymmetric linear system of equations using the BiConjugate Gradient Stabilized method.
 *
 *  Preconditioning is possible using a preconditioning class.
 *
 */

/** \fn int BiCGStab(const SpMatrix<T_> &A, const Prec<T_> &P, const Vect<T_> &b, Vect<T_> &x, int max_it, real_t toler, int verbose)
 *  \ingroup Solver
 *  \brief Biconjugate gradient stabilized solver function.
 *  @param [in] A Problem matrix (Instance of class SpMatrix).
 *  @param [in] P Preconditioner (Instance of class Prec).
 *  @param [in] b Right-hand side vector (class Vect)
 *  @param [in,out] x Vect instance containing initial solution guess on input and solution
 *  of the linear system on output (If iterations have succeeded).
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
int BiCGStab(const SpMatrix<T_>& A,
             const Prec<T_>&     P,
             const Vect<T_>&     b,
                   Vect<T_>&     x,
                   int           max_it,
                   real_t        toler,
                   int           verbose)
{
   if (verbose>0)
      cout << "Running preconditioned BiCGStab method ..." << endl;
   int it;
   size_t size=x.size();
   T_ rho_1=0, rho_2=0, alpha=0, beta=0, omega=0;
   Vect<T_> p(size), pp(size), s(size), ss(size), t(size), v(size), r(size), rr(size);
   real_t res, nrm=b.getNorm2();
   r = b - A*x;
   rr = r;

   if (nrm==0.0)
      nrm = 1;
   if ((res=r.getNorm2()/nrm)<=toler) {
      if (verbose>1)
         cout << "Convergence of the BiCGStab method after 0 iterations." << endl;
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
         if (verbose>1)
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
      if (verbose>1)
         cout << "Iteration: " << setw(4) << it << "  ... Residual: " << res << endl;
      if (verbose>2)
         cout << x;
      if ((res=r.getNorm2()/nrm)<toler) {
         toler = res;
         if (verbose)
            cout << "Convergence of the BiCGStab method after " << it << " iterations." << endl;
         return it;
      }
      if (omega==T_(0)) {
         toler = r.getNorm2()/nrm;
         return 3;
      }
   }
   if (verbose)
      cout << "Convergence after " << it << " iterations." << endl;
   return -it;
}


/** \fn int BiCGStab(const SpMatrix<T_> &A, int prec, const Vect<T_> &b, Vect<T_> &x, int max_it, real_t toler, int verbose)
 *  \ingroup Solver
 *  \brief Biconjugate gradient stabilized solver function.
 *  @param [in] A Problem matrix (Instance of class SpMatrix).
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
int BiCGStab(const SpMatrix<T_>& A,
                   int           prec,
             const Vect<T_>&     b,
                   Vect<T_>&     x,
                   int           max_it,
                   real_t        toler,
                   int           verbose)
{
   return BiCGStab(A,Prec<T_>(A,prec),b,x,max_it,toler,verbose);
}

} /* namespace OFELI */

#endif
