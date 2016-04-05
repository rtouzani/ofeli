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

/** \fn int GMRes(const SpMatrix<T_>& A, const Prec<T_>& P, const Vect<T_>& b, Vect<T_>& x, size_t m, int max_it, real_t toler, int verbose)
 *  \ingroup Solver
 *  \brief GMRes solver function.
 *  @param [in] A Problem matrix (Instance of class SpMatrix).
 *  @param [in] P Preconditioner (Instance of class Prec).
 *  @param [in] b Right-hand side vector (class Vect)
 *  @param  [in,out] x Vect instance containing initial solution guess in input and solution
 *  of the linear system in output (If iterations have succeeded).
 *  @param [in] m Number of subspaces to generate for iterations.
 *  @param [in] max_it Maximum number of iterations.
 *  @param [in] toler Tolerance for convergence (measured in relative weighted 2-Norm).
 *  @param [in] verbose Information output parameter (0: No output, 1: Output iteration information,
 *  2 and greater: Output iteration information and solution at each iteration.
 *  @return Number of performed iterations,
 *
 * \tparam <T_> Data type (double, float, complex<double>, ...)
 */
template<class T_>
int GMRes(const SpMatrix<T_>& A,
          const Prec<T_>&     P,
          const Vect<T_>&     b,
                Vect<T_>&     x,
                size_t        m,
                int           max_it,
                real_t        toler,
                int           verbose)
{
   if (verbose>0)
      cout << "Running preconditioned GMRes method ..." << endl;
   real_t res;
   int it=1;
   size_t n=b.size();
   Vect<T_> s(m+1), cs(m+1), sn(m+1), w(n);
   Vect<T_> H(m+1,m+1);

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
      if (verbose)
         cout << "Convergence of the GMRes method after " << it << " iterations." << endl;
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
         if (verbose>1)
            cout << "Iteration: " << setw(4) << it << ",  ... Residual: " << res << endl;
         if ((res = Abs(s[i+1])/normb) < toler) {
            update(x,int(i),H,s,v);
            toler = res;
            delete [] v;
            if (verbose)
               cout << "Convergence of the GMRes method after " << it << " iterations." << endl;
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
         if (verbose)
            cout << "Convergence of the GMRes method after " << it << " iterations." << endl;
         return it;
      }
   }
   toler = res;
   delete [] v;
   if (verbose)
      cout << "No Convergence of the GMRes method after " << it << " iterations." << endl;
   return -it;
}


/** \fn int GMRes(const SpMatrix<T_>& A, int prec, const Vect<T_>& b, Vect<T_>& x, size_t m, int max_it, real_t toler, int verbose)
 *  \ingroup Solver
 *  \brief GMRes solver function.
 *  @param [in] A Problem matrix (Instance of class SpMatrix).
 *  @param [in] prec Enum variable selecting a preconditioner, among the values <tt>IDENT_PREC</tt>,
 *  <tt>DIAG_PREC</tt>, <tt>ILU_PREC</tt> or <tt>SSOR_PREC</tt>
 *  @param [in] b Right-hand side vector (class Vect)
 *  @param  [in,out] x Vect instance containing initial solution guess in input and solution
 *  of the linear system in output (If iterations have succeeded).
 *  @param [in] m Number of subspaces to generate for iterations.
 *  @param [in] max_it Maximum number of iterations.
 *  @param [in] toler Tolerance for convergence (measured in relative weighted 2-Norm).
 *  @param [in] verbose Information output parameter (0: No output, 1: Output iteration information,
 *  2 and greater: Output iteration information and solution at each iteration.
 *  @return Number of performed iterations,
 *
 * \tparam <T_> Data type (double, float, complex<double>, ...)
 */
template<class T_>
int GMRes(const SpMatrix<T_>& A,
                int           prec,
          const Vect<T_>&     b,
                Vect<T_>&     x,
                size_t        m,
                int           max_it,
                real_t        toler,
                int           verbose)
{
   return GMRes(A,Prec<T_>(A,prec),b,x,m,max_it,toler,verbose);
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
