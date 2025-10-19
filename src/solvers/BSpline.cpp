/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2025 Rachid Touzani

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

                        Implementation of function BSpline

  ==============================================================================*/

#include "solvers/BSpline.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

void BSpline(size_t                n,
             size_t                t,
             Vect<Point<real_t> >& control, 
             Vect<Point<real_t> >& output,
             size_t                num_output)
{
   Point<real_t> calcxyz;
   Vect<size_t> u(n+t+1);
   compute_intervals(u,n,t);
// How much parameter goes up each time
   real_t inc = real_t(n-t+2)/(num_output-1);
   real_t interval = 0;

   for (size_t i=0; i<num_output-1; i++) {
      compute_point(u,n,t,interval,control,calcxyz);
      output[i] = calcxyz;
      interval += inc;
   }
// Put in the last point
   output[num_output-1] = control[n];
}


real_t blend(size_t              k,
             size_t              t,
             const Vect<size_t>& u,
             real_t              v)
// calculate the blending value
{
   real_t value;
   if (t==1) {
      if (u[k]<=v && v<u[k+1])
         value = 1.0;
      else
         value = 0.0;
   }
   else {
      if (u[k+t-1]==u[k] && u[k+t]==u[k+1])
         value = 0.0;
      else
         if (u[k+t-1]==u[k])
            value = (u[k+t]-v)/(u[k+t]-u[k+1])*blend(k+1,t-1,u,v);
         else
            if (u[k+t]==u[k+1])
               value = (v-u[k])/(u[k+t-1]-u[k])*blend(k,t-1,u,v);
            else
               value = (v-u[k])/(u[k+t-1]-u[k])*blend(k,t-1,u,v) +
                       (u[k+t]-v)/(u[k+t]-u[k+1])*blend(k+1,t-1,u,v);
   }
   return value;
}


void compute_intervals(Vect<size_t>& u,
                       size_t        n,
                       size_t        t)
{
   for (size_t j=0; j<=n+t; j++) {
      if (j<t)
         u[j] = 0;
      else
         if (t<=j && j<=n)
            u[j] = j-t+1;
         else
//          if n-t=-2 then we're screwed, everything goes to 0
            if (j>n)
               u[j] = n-t+2; 
   }
}


void compute_point(const Vect<size_t>&         u,
                   size_t                      n,
                   size_t                      t,
                   real_t                      v,
                   const Vect<Point<real_t> >& control,
                   Point<real_t>&              output)
{
   real_t temp;
// initialize the variables that will hold our outputted point
   output = 0.0;
   for (size_t k=0; k<=n; k++) {
      temp = blend(k,t,u,v);
//    same blend is used for each dimension coordinate
      output += control[k]*temp;
   }
}

} /* namespace OFELI */