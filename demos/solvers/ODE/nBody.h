/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2026 Rachid Touzani

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

       User defined class to define the Ordinary Differential System
                        for the n-Body problem

  ==============================================================================*/

#ifndef __NBODY_H
#define __NBODY_H

#include "OFELI.h"
#include "solvers/MyODE.h"

using namespace OFELI;

class nBody : public MyODE
{

 public:

/// \brief Constructor using problem parameters
    nBody(size_t dim, real_t G, const Vect<real_t>& m) : _dim(dim), _nb(m.size()), _G(G), _m(&m) { }

/// \brief Destructor
    ~nBody() { }

/** \brief Virtual member function to define function defining ODE in case of a single ODE
 *  @param [in] t Time at which the ode function is evaluated
 *  @param [in] y Unknown vector
 *  @param [in] i Component of function 
 *  @return Value of function
 */
    real_t Function (real_t t, const Vect<real_t> &y, int i)
    {
       if (i<=_nb*_dim)
          return y(i+_nb*_dim);
       real_t z = 0., f = 0.;
       size_t k = (i%_dim==0) ? _dim : i%_dim;
       i = ceil(float(i)/float(_dim)) - _nb;
       for (size_t j=1; j<=_nb; ++j) {
          if (i!=j) {
             real_t N = 0;
             for (size_t l=1; l<=_dim; ++l) {
                z = y(j,l,1) - y(i,l,1);
                N += z*z;
             }
             f += (*_m)(j)*(y(j,k,1)-y(i,k,1))/sqrt(N*N*N);
          }
       }
       return f*_G;
    }

 private:
   size_t _dim, _nb;
   real_t _G;
   const Vect<real_t> *_m;
};

#endif