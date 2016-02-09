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

                             Definition of class Iter
                            Class to drive iterations

  ==============================================================================*/


#ifndef __ITER_H
#define __ITER_H

#include "OFELI_Config.h"

namespace OFELI {

/*! \class Iter
 *  \ingroup Solver
 * \brief Class to drive an iterative process
 *
 * \details This template class enables monitoring any iterative process.
 * It simply sets default values for tolerance, maximal number of iterations
 * and enables checking convergence using two successive iterates.
 */

template<class T_> class Iter
{

 public:

/** \brief Default Constructor
 *  \details This constructor set default values: the maximal number
 *  of iterations is set to <tt>100</tt> and the tolerance to <tt>1.e-8</tt>
 */
    Iter();

/** \brief Constructor with iteration parameters
 *  @param [in] max_it Maximum number of iterations
 *  @param [in] toler Tolerance value for convergence
 *  @param [in] verbose Verbosity parameter [default: 0]
 *  <tt>0</tt>: No message output, <tt> > 0</tt>: message output with increasing display.
 */
    Iter(int    max_it,
         real_t toler,
         int    verbose=0);

/// \brief Destructor
    ~Iter() { }

/// \brief Set maximal number of iterations
    void setMaxIter(int max_it) { _max_it = max_it; }

/// \brief Set tolerance value for convergence
    void setTolerance(real_t toler) { _toler = toler; }

/// \brief Set verbosity parameter
    void setVerbose(int v) { _verbose = v; }

/** \brief Check convergence
 *  @param [in,out] u Solution vector at previous iteration
 *  @param [in] v Solution vector at current iteration
 *  @param [in] opt Vector norm for convergence checking
 *  <tt>1</tt>: 1-norm, <tt>2</tt>: 2-norm, 0: Max. norm [default: <tt>2</tt>]
 *  @return true if convergence criterion is satisfied, false if not
 *  \details After checking, this function copied <tt>v</tt> into <tt>u</tt>.
 */
    bool check(      Vect<T_>& u,
               const Vect<T_>& v,
                     int       opt=2);

 private:

   bool _conv;
   real_t _toler;
   int _it, _verbose, _max_it;
};


template<class T_>
Iter<T_>::Iter()
{
   _it = 1;
   _max_it = 100;
   _toler = 1.e-8;
   _conv = false;
   _verbose = 0;
}


template<class T_>
Iter<T_>::Iter(int    max_it,
               real_t toler,
               int    verbose)
{
   _it = 1;
   _max_it = max_it;
   _toler = toler;
   _conv = false;
   _verbose = verbose;
}


template<class T_>
bool Iter<T_>::check(      Vect<T_>& u,
                     const Vect<T_>& v,
                           int       opt)
{
   real_t a, b=0;
   size_t n=u.size();
   _it++;
   if (opt==1) {
      a = u.getNorm1();
      if (fabs(a)<_toler)
         a = 1;
      for (size_t i=0; i<n; i++)
         b += Abs(u[i]-v[i]);
   }
   else if (opt==2) {
      a = u.getNorm2();
      if (fabs(a)<_toler)
         a = 1;
      for (size_t i=0; i<n; i++)
         b += Abs2(u[i]-v[i]);
      b = sqrt(b);
   }
   else {
      a = u.getNormMax();
      for (size_t i=0; i<n; i++)
         b = std::max(b,Abs(u[i]-v[i]));
   }
   b /= a;
   if (_verbose > 1)
      cout << "Iteration " << _it-1 << ", Discrepancy: " << b << endl;
   u = v;

   if (b>_toler && _it>_max_it) {
      if (_verbose > 0)
         cout << "No convergence after " << _it-1 << " iterations, Discrepancy: "
              << b << endl;
      return true;
   }
   if (b<_toler && _it<=_max_it) {
      if (_verbose > 0)
         cout << "Convergence after " << _it-1 << " iterations, Discrepancy: " 
              << b << endl;
      return true;
   }
   else
      return false;
}

} /* namespace OFELI */

#endif
