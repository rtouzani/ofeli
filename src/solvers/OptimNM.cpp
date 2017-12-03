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

                       Implementation of class 'OptimTN'

  ==============================================================================*/

#include "OFELI_Config.h"
#include "linear_algebra/Vect.h"
#include "solvers/OptSolver.h"
#include "util/util.h"

#include <iostream>
using std::cout;

namespace OFELI {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/** \brief OptimNM minimizes a function using the Nelder-Mead algorithm.
 * 
 * \details 
 * This function seeks the minimum value of a user-specified objective function.
 *
 * Simplex function minimization procedure due to Nelder+Mead(1965),
 * as implemented by O'Neill (1971, Appl. Statist. 20, 338-45), with
 * subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
 * 25, 97) and Hill(1978, 27, 380-2)
 *
 * This function does not include a termination test using the
 * fitting of a quadratic surface.
 * 
 * This code is distributed under the GNU LGPL license. 
 * 
 * 
 * Author:
 * 
 * Original FORTRAN77 version by R ONeill.
 * C version by John Burkardt.
 * Adaptation to the OFELI library by R. Touzani
 * 
 * Reference:
 * 
 * John Nelder, Roger Mead,
 * A simplex method for function minimization,
 * Computer Journal,
 * Volume 7, 1965, pages 308-313.
 * 
 * R ONeill,
 * Algorithm AS 47:
 * Function Minimization Using a Simplex Procedure,
 * Applied Statistics,
 * Volume 20, Number 3, 1971, pages 338-345.
 *
 * @param [in] Objective Name of the function to minimize. This function has as an argument
 * the vector of optimization variables and returns the objective
 * @param [in,out] x Vector containing the initial guess for the optimization variables
 * on input and on output the estimated values that mimimize
 * @param [out] ynewlo The minimum value of the objective function
 * @param [in] reqmin The terminating limit for the variance function values
 * @param [in] step determines the size and shape of the initial simplex. The
 * relative magnitudes of its elements should reflect the units of the variables 
 * @param [in] conv The convergence check is carried out every <tt>conv</tt> iterations
 * @param [int] max_eval The maximum number of function evaluations
 * @param [out] nb_eval The number of function evaluations used
 * @param [out] nb_restart The number of restarts
 *
 * @return Error indicator
 * <ul>
 *   <li>0: No errors detected
 *   <li>1: <tt>reqmin</tt> or <tt>conv</tt> has an illegal value, or the number
 *          of optimization variables is illegal
 *   <li>2: iteration terminated because <tt>max_eval</tt> was exceeded without 
 *          convergence
 * </ul>
 *
 */
int OptimNM(OptSolver&    opt,
            Vect<real_t>& x,
            real_t&       ynewlo,
            real_t        reqmin,
            Vect<real_t>& step,
            int           conv,
            int           max_eval,
            int&          nb_eval,
            int&          nb_restart)
{
   int ret=0;
   real_t ccoeff=0.5, ecoeff=2.0, eps=0.001, rcoeff=1.0, xx, y2star;

// Check input parameters
   if (reqmin <= 0.0)
      return 1;
   size_t n=x.size();
   if (n==0)
      return 1;
   if (conv<1)
      return 1;

   Vect<real_t> p(n*(n+1)), pstar(n), p2star(n), pbar(n), y(n+1), xmin(n);
   nb_eval = nb_restart = 0;
   int jcount=conv; 
   real_t del=1.0, rq=reqmin*n;

// Initial or restarted loop
   for (;;) {
      for (size_t i=0; i<n; i++)
         p[i+n*n] = x[i];
      y[n] = opt.Objective(x);
      nb_eval++;
      for (size_t j=0; j<n; j++) {
         xx = x[j];
         x[j] += step[j]*del;
         for (size_t i=0; i<n; i++)
            p[i+j*n] = x[i];
         y[j] = opt.Objective(x);
         nb_eval++;
         x[j] = xx;
      }

/*    The simplex construction is complete
      Find highest and lowest y values. ynewlo = y(ihi) indicates
      the vertex of the simplex to be replaced.
*/
      real_t ylo=y[0];
      int ilo=0;
      for (size_t i=1; i<=n; i++) {
         if (y[i]<ylo) {
            ylo = y[i];
            ilo = i;
         }
      }

//    Inner loop
      for (;;) {
         if (max_eval<=nb_eval)
            break;
         ynewlo = y[0];
         int ihi=0;
         for (size_t i=1; i<=n; i++) {
            if (ynewlo<y[i]) {
               ynewlo = y[i];
            ihi = i;
         }
      }

//    Calculate pbar, the centroid of the simplex vertices
//    excepting the vertex with y value ynewlo
      for (size_t i=0; i<n; i++) {
         real_t z=0.0;
         for (size_t  j=0; j<=n; j++)
            z += p[i+j*n];
         z -= p[i+ihi*n];  
         pbar[i] = z/n;
      }

//    Reflection through the centroid
      for (size_t i=0; i<n; i++)
         pstar[i] = pbar[i] + rcoeff*(pbar[i]-p[i+ihi*n]);
      real_t ystar = opt.Objective(pstar);
      nb_eval++;

//    Successful reflection, so extension
      if (ystar<ylo) {
         for (size_t i=0; i<n; i++)
            p2star[i] = pbar[i] + ecoeff*(pstar[i]-pbar[i]);
         y2star = opt.Objective(p2star);
         nb_eval++;

//       Check extension
         if (ystar<y2star) {
            for (size_t i=0; i<n; i++)
               p[i+ihi*n] = pstar[i];
            y[ihi] = ystar;
         }

//       Retain extension or contraction
         else {
            for (size_t i=0; i<n; i++)
               p[i+ihi*n] = p2star[i];
            y[ihi] = y2star;
         }
      }

//    No extension
      else {
         size_t l=0;
         for (size_t i=0; i<=n; i++) {
            if (ystar<y[i])
               l++;
         }

         if (1<l) {
            for (size_t i=0; i<n; i++)
               p[i+ihi*n] = pstar[i];
            y[ihi] = ystar;
         }

//       Contraction on the y(ihi) side of the centroid
         else if (l==0) {
            for (size_t i=0; i<n; i++)
               p2star[i] = pbar[i] + ccoeff*(p[i+ihi*n]-pbar[i]);
            y2star = opt.Objective(p2star);
            nb_eval++;

//          Contract the whole simplex
            if (y[ihi]<y2star) {
               for (size_t j=0; j<=n; j++) {
                  for (size_t i=0; i<n; i++) {
                     p[i+j*n] = (p[i+j*n] + p[i+ilo*n])*0.5;
                     xmin[i] = p[i+j*n];
                  }
                  y[j] = opt.Objective(xmin);
                  nb_eval++;
               }
               ylo = y[0];
               ilo = 0;

               for (size_t i=1; i<=n; i++) {
                  if (y[i]<ylo) {
                     ylo = y[i];
                     ilo = i;
                  }
               }
               continue;
            }

//          Retain contraction
            else {
               for (size_t i=0; i<n; i++)
                  p[i+ihi*n] = p2star[i];
               y[ihi] = y2star;
            }
         }

//       Contraction on the reflection side of the centroid
         else if (l==1) {
            for (size_t i=0; i<n; i++)
               p2star[i] = pbar[i] + ccoeff*(pstar[i]-pbar[i]);
            y2star = opt.Objective(p2star);
            nb_eval++;

//          Retain reflection ?
            if (y2star<=ystar) {
               for (size_t i=0; i<n; i++)
                  p[i+ihi*n] = p2star[i];
               y[ihi] = y2star;
            }
            else {
               for (size_t i=0; i<n; i++)
                  p[i+ihi*n] = pstar[i];
               y[ihi] = ystar;
            }
         }
      }

//    Check if ylo improved
      if (y[ihi]<ylo) {
         ylo = y[ihi];
         ilo = ihi;
      }
      jcount--;
      if (0<jcount)
         continue;

//    Check to see if minimum reached
      if (nb_eval<=max_eval) {
         jcount = conv;
         real_t z=0.0;
         for (size_t i=0; i<=n; i++)
            z += y[i];
         xx = z/(n+1);
         z = 0.0;
         for (size_t i=0; i<=n; i++)
            z += (y[i]-xx)*(y[i]-xx);
         if (z<=rq)
            break;
      }
   }

// Factorial tests to check that ynewlo is a local minimum
   for (size_t i=0; i<n; i++)
      xmin[i] = p[i+ilo*n];
      ynewlo = y[ilo];

      if (max_eval<nb_eval) {
         ret = 2;
         break;
      }
      ret = 0;

      for (size_t i=0; i<n; i++) {
         del = step[i] * eps;
         xmin[i] += del;
         real_t z = opt.Objective(xmin);
         nb_eval++;
         if (z<ynewlo) {
            ret = 2;
            break;
         }
         xmin[i] -= del + del;
         z = opt.Objective(xmin);
         nb_eval++;
         if (z<ynewlo) {
            ret = 2;
            break;
         }
         xmin[i] += del;
      }
      if (ret==0)
         break;

//    Restart the procedure
      x = xmin;
      del = eps;
      nb_restart++;
   }
   x = xmin;
   return ret;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */
