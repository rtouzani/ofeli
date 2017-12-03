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

                 Template function for iteration Richardson Method

  ==============================================================================*/

#ifndef __RICHARDSON_H
#define __RICHARDSON_H

#include <iostream>
using std::ostream;
using std::endl;

#include <iomanip>
using std::setw;

#include "OFELI_Config.h"
#include "linear_algebra/Vect.h"
#include "util/util.h"
#include "io/output.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Richardson.h
 *  \brief Function to solve a linear system of equations using the Richardson method.
 *
 */

/** \fn int Richardson(const M_ &A, const Vect<T_> &b, Vect<T_> &x, real_t omega, int max_it, real_t toler, int verbose)
 *  \ingroup Solver
 *  \brief Richardson solver function.
 *  @param [in] A Problem matrix problem (Instance of abstract class \b M_).
 *  @param [in] b Right-hand side vector (class Vect)
 *  @param x Vect instance containing initial solution guess in input and solution of the linear system in output (If iterations have succeeded).
 *  @param [in] omega Relaxation parameter.
 *  @param [in] max_it Maximum number of iterations.
 *  @param [in] toler Tolerance for convergence (measured in relative weighted 2-Norm).
 *  @param [in] verbose Information output parameter (<tt>0</tt>: No output,
 *  <tt>1</tt>: Output iteration information,
 *  <tt>2</tt> and greater: Output iteration information and solution at each iteration.
 *  @return nb_it Number of performed iterations,
 *
 * \tparam <T_> Data type (real_t, float, complex<real_t>, ...)
 * \tparam <M_> %Matrix storage class
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
template< class T_, class M_>
int Richardson(const M_&       A,
               const Vect<T_>& b,
               Vect<T_>&       x,
               real_t          omega,
               int             max_it,
               real_t          toler,
               int             verbose)
{
   if (verbose>0)
      cout << "Running Richardson method ..." << endl;
   size_t size = x.Size();
   real_t nrm = x.Norm2();
   if (nrm == 0)
      nrm = 1;
   Vect<T_> y(size);

   for (size_t it=1; it<=max_it; it++) {
      A.Mult(x,y);
      for (size_t i=1; i<=size; i++)
         y(i) = omega*(b(i)-y(i));
      real_t err = y.Norm2()/nrm;

      if (verbose > 1)
         cout << "Iteration: " << setw(4) << it << "  ... Error: " << err << endl;

      x += y;
      nrm = x.Norm2();
      if (verbose > 2)
         cout << x;
      if (err < toler) {
         if (verbose)
            cout << "Convergence of the Richardson method after " << it << " iterations." << endl;
         return it;
      }
   }
   if (verbose)
      cout << "No Convergence after " << it << " iterations." << endl;
   return -it;
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
