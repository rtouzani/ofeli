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

             Template function for Conjugate Gradient Squared Method
                 CGS is inspired from the SIAM Templates book

  ==============================================================================*/

#ifndef __CGS_H
#define __CGS_H

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
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file CGS.h
 *  \brief Solves an unsymmetric linear system of equations using the Conjugate Gradient Squared method.
 *
 *  Preconditioning is possible using a preconditioning class.
 */

/** \fn int CGS(const SpMatrix<T_> &A, const P_ &P, const Vect<T_> &b, Vect<T_> &x, int max_it, real_t toler, int verbose)
 * \ingroup Solver
 * \brief Conjugate Gradient Squared solver function.
 *  @param [in] A Problem matrix (Instance of class SpMatrix).
 *  @param [in] P Preconditioner (Instance of class Prec).
 *  @param [in] b Right-hand side vector (class Vect)
 *  @param [in,out] x Vect instance containing initial solution guess in input and solution
 *  of the linear system in output (If iterations have succeeded).
 *  @param [in] max_it Maximum number of iterations.
 *  @param [in] toler Tolerance for convergence (measured in relative weighted 2-Norm).
 *  @param [in] verbose Information output parameter
 *  <ul>
 *    <li><tt>0</tt>: No output
 *    <li><tt>1</tt>: Output iteration information,
 *    <li><tt>2</tt> and greater: Output iteration information and solution at each iteration.
 *  @return Number of performed iterations
 *
 * \tparam <T_> Data type (real_t, float, complex<real_t>, ...)
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
template<class T_>
int CGS(const SpMatrix<T_>& A,
        const Prec<T_>&     P,
        const Vect<T_>&     b,
        Vect<T_>&           x,
        int                 max_it,
        real_t              toler,
        int                 verbose)
{
   if (verbose>0)
      cout << "Running preconditioned CGS method ..." << endl;
   int it;
   size_t size=x.size();
   T_ rho_1=0, rho_2=0, alpha=0, beta=0;
   real_t res;
   Vect<T_> r(size), p(size), y(size), q(size), z(size), v(size), u(size), w(size);

   real_t nrm = b.getNorm2();
   Vect<T_> s = b - A*x;
   r = s;

   if (nrm==0)
      nrm = 1;
   if ((res=r.getNorm2()/nrm)<=toler) {
      toler = res;
      if (verbose)
         cout << "Convergence after 0 iterations." << endl;
      return 0;
   }

   for (it=1; it<=max_it; it++) {
      rho_1 = (s,r);
      if (rho_1 == T_(0)) {
         toler = r.getNorm2() / nrm;
         if (verbose)
            cout << "Convergence of the CGS method after 2 iterations." << endl;
         return 2;
      }
      if (it==1) {
         u = r;
         p = u;
      } else {
         beta = rho_1 / rho_2;
         u = r + beta * q;
         p = u + beta * (q + beta * p);
      }
      y = p;
      P.solve(y);
      v = A*y;
      alpha = rho_1/(s,v);
      w = u - alpha*v;
      w += u;
      P.solve(w);
      Axpy(alpha,w,x);
      A.Mult(w,z);
      Axpy(-alpha,z,r);
      rho_2 = rho_1;
      res = r.getNorm2()/nrm;
      if (verbose > 1) {
         cout << "Iteration: " << setw(4) << it;
         cout << "  ... Residual: " << res << endl;
      }
      if (verbose>2)
         cout << x;
      if (res<toler) {
         toler = res;
         if (verbose)
            cout << "Convergence of the CGS method after " << it << " iterations." << endl;
         return it;
      }
   }
   toler = res;
   if (verbose)
      cout << "No Convergence of the CGS method after " << it << " iterations." << endl;
   return -max_it;
}


/** \fn int CGS(const SpMatrix<T_> &A, int prec, const Vect<T_> &b, Vect<T_> &x, int max_it, real_t toler, int verbose)
 * \ingroup Solver
 * \brief Conjugate Gradient Squared solver function.
 *  @param [in] A Problem matrix (Instance of class SpMatrix).
 *  @param [in] prec Enum variable selecting a preconditioner, among the values <tt>IDENT_PREC</tt>,
 *  <tt>DIAG_PREC</tt>, <tt>ILU_PREC</tt> or <tt>SSOR_PREC</tt>
 *  @param [in] b Right-hand side vector (class Vect)
 *  @param [in,out] x Vect instance containing initial solution guess in input and solution
 *  of the linear system in output (If iterations have succeeded).
 *  @param [in] max_it Maximum number of iterations.
 *  @param [in] toler Tolerance for convergence (measured in relative weighted 2-Norm).
 *  @param [in] verbose Information output parameter
 *  <ul>
 *    <li><tt>0</tt>: No output
 *    <li><tt>1</tt>: Output iteration information,
 *    <li><tt>2</tt> and greater: Output iteration information and solution at each iteration.
 *  @return Number of performed iterations
 *
 * \tparam <T_> Data type (real_t, float, complex<real_t>, ...)
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
template<class T_>
int CGS(const SpMatrix<T_>& A,
        int                 prec,
        const Vect<T_>&     b,
        Vect<T_>&           x,
        int                 max_it,
        real_t              toler,
        int                 verbose)
{
   return CGS(A,Prec<T_>(A,prec),b,x,max_it,toler,verbose);
}
 
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
int CGS(const Matrix<T_>* A,
        int               prec,
        const Vect<T_>&   b,
        Vect<T_>&         x,
        int               max_it,
        real_t            toler,
        int               verbose)
{
   SpMatrix<T_> &AA = MAT(SpMatrix<T_>,A);
   return CGS(AA,prec,b,x,max_it,toler,verbose);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
