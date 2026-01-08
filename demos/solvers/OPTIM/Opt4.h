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

    User defined class to define the objective for the Brachistochrone problem

  ==============================================================================*/

#ifndef __OPT4_H
#define __OPT4_H

#include "OFELI.h"
#include "solvers/MyOpt.h"

using namespace OFELI;

class Opt4 : public MyOpt
{

 public:

/// Constructor using problem data
    Opt4(double g, double xm, double xM, Vect<double>& y)
    {
       N = y.size() - 1;
       h = (xM-xm)/N;
       for (size_t i=1; i<N; ++i)
          y[i] = y[i-1] + (y[N]-y[0])/N;
       cc = 2./sqrt(2*g);
    }

/*  Function to define the objective
 *  x (in) Vector of optimization variables
 *  Return value: Objective
 */
    double Objective(Vect<double>& y)
    {
       double f = 0.0;
       for (size_t i=0; i<N; ++i) {
          f += cc*sqrt(h*h + (y[i+1]-y[i])*(y[i+1]-y[i]))/(sqrt(y[0]-y[i])+sqrt(y[0]-y[i+1]));
       }
       return f; 
    }

/* Function to define the gradient
 * x (in) Vector of optimization variables
 * g (out) Vector of gradient
 */
    void Gradient(Vect<double>& y,
                  Vect<double>& G)
    {
       double a=0., b=0., u=0., v=0.;
       for (size_t j=0; j<=N; ++j) {
          G[j] = 0.;
          for (size_t i=0; i<N; ++i) {
             a = sqrt(y[0]-y[i]) + sqrt(y[0]-y[i+1]);
             b = sqrt(h*h + (y[i+1]-y[i])*(y[i+1]-y[i]));
             u = a*(y[i+1]-y[i])*(Kronecker(i+1,j)-Kronecker(i,j))/b;
             if (i==0)
                v = 0.5*b*(Kronecker(0,j)-Kronecker(i+1,j))/sqrt(y[0]-y[i+1]);
             else
                v = 0.5*b*((Kronecker(0,j)-Kronecker(i,j))/sqrt(y[0]-y[i]) +
                           (Kronecker(0,j)-Kronecker(i+1,j))/sqrt(y[0]-y[i+1]));
             G[j] += cc*(u-v)/(a*a);
          }
       }
    }

    size_t N;
    double cc, h;
};

#endif
