/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzani

   This program is free software; you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free
   Software Foundation; Version 2 of the License.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
   details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the :

   Free Software Foundation
   Inc., 59 Temple Place - Suite 330
   Boston, MA  02111-1307, USA

  ==============================================================================

                            Prototypes for class 'User'

  ==============================================================================*/

#ifndef __USER_H
#define __USER_H

#include "OFELI.h"
using namespace OFELI;

class User : public UserData<double>  {

public :

   User(Mesh &mesh) : UserData<double> (mesh) { ; }

   double BoundaryCondition(Point<double> x, int code, double time=0., size_t dof=1)
   {
      time = 0;
      dof = 1;
      double ret = 0.0;
      if (code == 2)
         ret = 1.0;
      return ret;
   }

   double SurfaceForce(Point<double> x, int code, double time, size_t dof)
   {
      time = 0;
      dof = 1;
      if (code)
         return 1.0;
      else
         return 0.0;
   }
};

#endif
