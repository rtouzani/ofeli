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

                 Template function for Generalized Minimum Residual
              GMRes is a transcription of the algorithm described in the
                                SIAM Templates book

  ==============================================================================*/

#ifndef __GMRES_H
#define __GMRES_H

#include <iostream>
using std::ostream;
using std::endl;
using std::cout;

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

/*! \file GMRes.h
 *  \brief Function to solve a linear system of equations
 *  using the Generalized Minimum Residual method.
 *
 *  Preconditioning is possible using a preconditioning class.
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template<class T_>
void rotate(T_& dx,
            T_& dy,
            T_& cs,
            T_& sn)
{
   T_ temp=cs*dx+sn*dy;
   dy = cs*dy - sn*dx;
   dx = temp;
}

template< class T_ >
void update(Vect<T_>& x,
            int       k,
            Vect<T_>& H,
            Vect<T_>& s,
            Vect<T_>  v[])
{
   Vect<T_> y(s);
   for (int i=k; i>=0; i--) {
      y[i] /= H(i+1,i+1);
      for (int j=i-1; j>=0; j--)
         y[j] -= H(j+1,i+1) * y[i];
   }
   for (int j=0; j<=k; j++)
      for (size_t l=0; l<x.size(); l++)
         x[l] += v[j][l] * y[j];
}


template<class T_>
void set_rotation(T_& dx,
                  T_& dy,
                  T_& cs,
                  T_& sn)
{
    if (dy == 0.0) {
       cs = 1.0;
       sn = 0.0;
    } else if (Abs(dy) > Abs(dx)) {
       T_ temp = dx/dy;
       sn = 1.0/sqrt(1.0 + temp*temp);
       cs = temp * sn;
    } else {
       T_ temp = dy/dx;
       cs = 1.0/sqrt(1.0 + temp*temp);
       sn = temp * cs;
   }
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \fn int GMRes(const SpMatrix<T_>& A, const Prec<T_>& P, const Vect<T_>& b, Vect<T_>& x, size_t m, int max_it, real_t toler)
 *  \ingroup Solver
 *  \brief GMRes solver function.
 *  \details This function uses the preconditioned GMRES algorithm to solve a
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
 *  @param  [in,out] x Vect instance containing initial solution guess in input and solution
 *  of the linear system in output (If iterations have succeeded).
 *  @param [in] m Number of subspaces to generate for iterations.
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
int GMRes(const SpMatrix<T_>& A,
          const Prec<T_>&     P,
          const Vect<T_>&     b,
          Vect<T_>&           x,
          size_t              m,
          int                 max_it,
          real_t              toler)
{
   if (Verbosity>2)
      cout << "Running preconditioned GMRes method ..." << endl;
   real_t res;
   int it=1;
   size_t n=b.size();
   Vect<T_> s(m+1), cs(m+1), sn(m+1), w(n), H(m+1,m+1);

   real_t normb=b.getNorm2();
   Vect<T_> r(b);
   r = b - A*x;
   P.solve(r);
   real_t beta=r.getNorm2();

   if (normb == 0.0)
      normb = 1;
   if ((res=r.getNorm2()/normb) <= toler) {
      toler = res;
      max_it = 0;
      if (Verbosity>2)
         cout << "Convergence of the GMRes method after " << it << " iterations." << endl;
      if (Verbosity>6)
         cout << "Solution:\n" << x;
      return it;
   }
   Vect<T_> *v = new Vect<T_> [m+1];
   for (size_t i=0; i<=m; i++) {
      v[i].resize(n);
      v[i] = T_(0);
   }

   while (it<=max_it) {
      Scale(T_(1.0/beta),r,v[0]);
      s = T_(0.);
      s[0] = beta;
      for (size_t i=0; i<m && it<=max_it; i++, it++) {
         w = A*v[i];
         P.solve(w);
         for (size_t k=0; k<=i; k++) {
            H(k+1,i+1) = (w,v[k]);
            for (size_t j=0; j<n; j++)
               w[j] -= H(k+1,i+1) * v[k][j];
         }
         H(i+2,i+1) = w.getNorm2();
         Scale((1.0/H(i+2,i+1)),w,v[i+1]);
         for (size_t k=0; k<i; k++)
            rotate(H(k+1,i+1), H(k+2,i+1), cs[k], sn[k]);
         set_rotation(H(i+1,i+1), H(i+2,i+1), cs[i], sn[i]);
         rotate(H(i+1,i+1), H(i+2,i+1), cs[i], sn[i]);
         rotate(s[i], s[i+1], cs[i], sn[i]);
         if (Verbosity>3)
            cout << "Iteration: " << setw(4) << it << ",  ... Residual: " << res << endl;
         if (Verbosity>10)
            cout << "Solution at iteration " << it << ": \n" << x;
         if ((res = Abs(s[i+1])/normb) < toler) {
            update(x,int(i),H,s,v);
            toler = res;
            delete [] v;
            if (Verbosity>2)
               cout << "Convergence of the GMRes method after " << it << " iterations." << endl;
            if (Verbosity>6)
               cout << "Solution:\n" << x;
            return it;
         }
      }

      update(x,int(m)-1,H,s,v);
      r = b - A*x;
      P.solve(r);
      beta = r.getNorm2();
      if ((res = beta/normb)< toler) {
         toler = res;
         delete [] v;
         if (Verbosity>2)
            cout << "Convergence of the GMRes method after " << it << " iterations." << endl;
         if (Verbosity>6)
            cout << "Solution:\n" << x;
         return it;
      }
   }
   toler = res;
   delete [] v;
   if (Verbosity>2)
      cout << "No Convergence after " << it << " iterations." << endl;
   return -it;
}


/** \fn int GMRes(const SpMatrix<T_>& A, int prec, const Vect<T_>& b, Vect<T_>& x, size_t m, int max_it, real_t toler)
 *  \ingroup Solver
 *  \brief GMRes solver function.
 *  \details This function uses the preconditioned GMRES algorithm to solve a
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
 *  @param  [in,out] x Vect instance containing initial solution guess in input and solution
 *  of the linear system in output (If iterations have succeeded).
 *  @param [in] m Number of subspaces to generate for iterations.
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
int GMRes(const SpMatrix<T_>& A,
          int                 prec,
          const Vect<T_>&     b,
          Vect<T_>&           x,
          size_t              m,
          int                 max_it,
          real_t              toler)
{
   return GMRes(A,Prec<T_>(A,prec),b,x,m,max_it,toler);
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
int GMRes(const Matrix<T_>* A,
          int               prec,
          const Vect<T_>&   b,
          Vect<T_>&         x,
          size_t            m,
          int               max_it,
          real_t            toler)
{
   SpMatrix<T_> &AA = MAT(SpMatrix<T_>,A);
   return GMRes(AA,prec,b,x,m,max_it,toler);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
