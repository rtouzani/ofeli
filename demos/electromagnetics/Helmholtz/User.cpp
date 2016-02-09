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
                        User Defined functions to prescribe data
  ==============================================================================*/


#include "User.h"


complex<double> User::SurfaceForce(Point<double> x, int code, double time, size_t dof)
{
   time = 0;
   dof = 1;
   double z = -2*OFELI_PI*cos(OFELI_PI*x.y);
   if (code==1)
      return complex_t(0,z);
   else if (code==2)
      return complex_t(z,0);
   else
      return complex_t(0,0);
}


complex<double> User::ExactSolution(Point<double> x)
{
   double z = cos(OFELI_PI*x.y);
   return complex<double>(cos(2*OFELI_PI*x.x)*z,sin(2*OFELI_PI*x.x)*z);
}
