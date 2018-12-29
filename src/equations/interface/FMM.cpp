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

  ==============================================================================

                         Implementation of class 'FMM'

  ==============================================================================*/

#include "equations/interface/FMM.h"
#include "mesh/Grid.h"
#include <algorithm>

namespace OFELI {

FMM::FMM(const Grid&   g,
         Vect<real_t>& phi,
         bool          HA)
    : _u(&phi), _high_accuracy(HA)
{
   _inf = HUGE_VAL;
   _grid = &g;
   _dim = _grid->getDim();
   _nx = _grid->getNx(), _ny = _grid->getNy(), _nz = _grid->getNz();
   _hx = _grid->getHx(), _hy = _grid->getHy(), _hz = _grid->getHz();
   if (_dim==2) {
      _AlivePt.setSize(_nx+1,_ny+1);
      _Narrow.set((_nx+1)*(_ny+1));
   }
   else {
      _AlivePt.setSize(_nx+1,_ny+1,_nz+1);
      _Narrow.set((_nx+1)*(_ny+1)*(_nz+1));
   }
   _AlivePt = _inf;
}


FMM::~FMM()
{
}


int FMM::MaxQuad(real_t  a,
                 real_t  b,
                 real_t  c,
                 real_t& x)
{
   real_t delta = b*b - a*c;
   if (delta<0)
      return 1;
   real_t z = sqrt(delta);
   x = fmax((-b-z)/a,(-b+z)/a);
   return 0;
}

} /* namespace OFELI */
