/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2019 Rachid Touzani

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

  ===================================== =========================================
                            Definition of class 'User'
  ==============================================================================*/

#ifndef __USER_H
#define __USER_H

#include "io/UserData.h"
#include "linear_algebra/Point.h"
#include "mesh/Mesh.h"

using namespace OFELI;
using std::complex;

class User : virtual public UserData<complex_t> {

public :

   User(Mesh &mesh) : UserData<complex_t>(mesh) { ; }
   complex<double> SurfaceForce(Point<double> x, int code, double time, size_t dof);
   complex<double> ExactSolution(Point<double> x);
};

#endif
