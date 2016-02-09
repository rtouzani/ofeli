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

  ==============================================================================*/

#ifndef UPDATEMF_H__
#define UPDATEMF_H__

#include "OFELI.h"
#include "Electromagnetics.h"
#define MU0  1.2566370614359172953850573533118e-6

using namespace OFELI;

void UpdateMF(Mesh&            ms0,
              Mesh&            ms1,
              double           omega,
              complex_t&       volt,
              Vect<complex_t>& b0,
              Vect<complex_t>& b1);

#endif
