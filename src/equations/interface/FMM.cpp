/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2018 Rachid Touzani

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

                         Implementation of class 'FMM'

  ==============================================================================*/

#include "equations/interface/FMM.h"
#include "mesh/Grid.h"
#include <algorithm>

using std::max;

namespace OFELI {

FMM::FMM(const Grid&   g,
         Vect<real_t>& phi,
         bool          HA)
    : _phi(&phi), _high_accuracy(HA)
{
   _inf = std::numeric_limits<real_t>::max();
   _nx = g.getNx(); _ny = g.getNy(); _nz = g.getNz();
   Point<real_t> pMin(g.getXMin()), pMax(g.getXMax());
   if (_nz==0) {
      _gd = new Grid(pMin[0],pMax[0],pMin[1],pMax[1],_nx,_ny);
      _AlivePt.setSize(_nx+1,_ny+1);
   }
   else {
      _gd = new Grid(pMin[0],pMax[0],pMin[1],pMax[1],pMin[2],pMax[2],_nx,_ny,_nz);
      _AlivePt.setSize(_nx+1,_ny+1,_nz+1);
   }
   _AlivePt = _inf;
}


FMM::~FMM()
{
   if (_gd)
      delete _gd;
}


int FMM::MaxQuadratic(real_t  a,
                      real_t  b,
                      real_t  c,
                      real_t& x)
{
   if (a == 0) {
      if (b == 0)
         x = _inf;
      else
         x = -c/b;
      return 1;
   }
   
   real_t delta = b*b - 4.0*a*c;
   if (c != 0) {
      if (delta<0)
         return 0;
      else {
         if (delta==0)
            x = -b/(2.0*a);
         else
            x = max((-b-sqrt(delta))/(2.0*a),(-b+sqrt(delta))/(2.0*a));
      }
   }
   return 1;
}

} /* namespace OFELI */
