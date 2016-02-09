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

                Template function for BiConjugate Gradient Method
                  BiCG is inspired from the SIAM Templates book

  ==============================================================================*/

#ifndef __BICG_H
#define __BICG_H

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

/*! \defgroup Solver Solver
 *  \brief Gathers Solver functions
 */

/*! \file BiCG.h
 *  \brief Solves an unsymmetric linear system of equations using the BiConjugate Gradient method.
 *
 */

/** \fn int BiCG(const SpMatrix<T_> &A, const Prec<T_> &P, const Vect<T_> &b, Vect<T_> &x, int max_it, real_t toler, int verbose)
 *  \ingroup Solver
 *  \brief Biconjugate gradient solver function.
 *  @param [in] A Problem matrix (Instance of class SpMatrix).
 *  @param [in] P Preconditioner (Instance of class Prec).
 *  @param [in] b Right-hand side vector (class Vect)
 *  @param [in,out] x Vect instance containing initial solution guess in input and solution
 *  of the linear system in output (If iterations have succeeded).
 *  @param [in] max_it Maximum number of iterations.
 *  @param toler [in] Tolerance for convergence (measured in relative weighted 2-Norm).
 *  @param verbose [in] Information output parameter
 *    - 0: No output
 *    - 1: Output iteration information,
 *    - 2 and greater: Output iteration information and solution at each iteration.
 *  @return Number of performed iterations,
 *
 * \tparam <T_> Data type (double, float, complex<double>, ...)
 */
template<class T_>
int BiCG(const SpMatrix<T_>& A,
         const Prec<T_>&     P,
         const Vect<T_>&     b,
               Vect<T_>&     x,
               int           max_it,
               real_t&       toler,
               int           verbose)
{
   if (verbose>0)
      cout << "Running preconditioned BiCG method ..." << endl;
   int it;
   size_t size=x.size();
   real_t res;
   T_ alpha=0, beta=0, rho_1=0, rho_2=0;
   Vect<T_> r(size), z(size), zz(size), p(size), pp(size), q(size), qq(size);

   real_t nrm=b.getNorm2();
   A.Mult(x,r);
   r = b - r;
   Vect<T_> rr(r);

   if (nrm == 0.0)
      nrm = 1;
   if ((res = r.getNorm2() / nrm) <= toler) {
      toler = res;
      if (verbose > 1)
         cout << "Convergence of the BiCG method after  0 iterations." << endl;
      return 0;
   }

   for (it=1; it<=max_it; it++) {
      P.solve(r,z);
      P.TransSolve(rr,zz);
      rho_1 = Dot(z,rr);
      if (rho_1 == T_(0)) {
         toler = r.getNorm2() / nrm;
         if (verbose > 1)
            cout << "Convergence of the BiCG method after " << it << " iterations." << endl;
         return it;
      }
      if (it == 1) {
         p = z;
         pp = zz;
      } else {
         beta = rho_1 / rho_2;
         p = z + beta*p;
         pp = zz + beta*pp;
      }
      A.Mult(p,q);
      A.TMult(pp,qq);
      alpha = rho_1/Dot(pp,q);
      Axpy( alpha,p,x);
      Axpy(-alpha,q,r);
      Axpy(-alpha,qq,rr);
      rho_2 = rho_1;

      if (verbose > 1)
         cout << "Iteration: " << setw(4) << it << "  ... Residual: " << res << endl;
      if (verbose > 2)
         cout << x;
      if ((res = r.getNorm2() / b.getNorm2()) < toler) {
         toler = res;
         if (verbose > 1)
            cout << "Convergence of the BiCG method after " << it << " iterations." << endl;
         return it;
      }
   }
   return it;
}


/** \fn int BiCG(const SpMatrix<T_> &A, int prec, const Vect<T_> &b, Vect<T_> &x, int max_it, real_t toler, int verbose)
 *  \ingroup Solver
 *  \brief Biconjugate gradient solver function.
 *  @param [in] A Problem matrix (Instance of class SpMatrix).
 *  @param [in] prec Enum variable selecting a preconditioner, among the values <tt>IDENT_PREC</tt>,
 *  <tt>DIAG_PREC</tt>, <tt>ILU_PREC</tt> or <tt>SSOR_PREC</tt>
 *  @param [in] b Right-hand side vector (class Vect)
 *  @param [in,out] x Vect instance containing initial solution guess in input and solution
 *  of the linear system in output (If iterations have succeeded).
 *  @param [in] max_it Maximum number of iterations.
 *  @param toler [in] Tolerance for convergence (measured in relative weighted 2-Norm).
 *  @param verbose [in] Information output parameter
 *    - 0: No output
 *    - 1: Output iteration information,
 *    - 2 and greater: Output iteration information and solution at each iteration.
 *  @return Number of performed iterations,
 *
 * \tparam <T_> Data type (double, float, complex<double>, ...)
 */
template<class T_>
int BiCG(const SpMatrix<T_>& A,
               int           prec,
         const Vect<T_>&     b,
               Vect<T_>&     x,
               int           max_it,
               real_t        toler,
               int           verbose)
{
   return BiCG(A,Prec<T_>(A,prec),b,x,max_it,toler,verbose);
}

} /* namespace OFELI */

#endif
