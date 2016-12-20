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

  ==============================================================================*/

#include "Electromagnetics.h"
#include "mesh/Mesh.h"
#include "linear_algebra/SkMatrix.h"
#include "linear_algebra/Vect.h"
#include <complex>

#ifndef LARGE
#define LARGE 1.e20
#endif

using namespace OFELI;

void ScaledMF(Mesh&                ms,
              SkMatrix<complex_t>& A,
              Vect<complex_t>&     b,
              complex_t&           current,
              double               omega,
              int                  flag)
//------------------------------------------------------------------------------
//              Calculate scaled magnetic field
//------------------------------------------------------------------------------
{
   flag = 0;
   MeshElements(ms) {
      EC2D1T3 eq(theElement);
      eq.Magnetic(omega,1.);
      eq.Electric();
      A.Assembly(TheElement,eq.eMat.get());
   }

   MeshNodes(ms) {
      int m = TheNode.getDOF(1);
      if (TheNode.getCode(1)==1) {
         A.set(m,m,A(m,m)*LARGE);
         b[m-1] = A(m,m)*current;
      }
      if (TheNode.getCode(1)==2) {
         A.set(m,m,A(m,m)*LARGE);
         b[m-1] = 0;
      }
   }
   A.solve(b);
}
