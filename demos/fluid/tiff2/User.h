/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2017 Rachid Touzani

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
                            Definition of class 'User'
  ==============================================================================*/

#ifndef __USER_H
#define __USER_H

#include "OFELI.h"
using namespace OFELI;

class User : public UserData<double>
{
 public:

   User(Mesh &mesh) : UserData<double>(mesh) { }

   double BoundaryCondition(Point<double> x, int code, double time=0., size_t dof=1)
   {
      time = 0;
      dof = 1;
      double ret = 0.0;
      if (code==1)
         ret = 0.0;
      if (code==3)
         ret = x.y*(1-x.y);
      if (code==2)
         ret = x.y;
      return ret;
   }
};

#endif
