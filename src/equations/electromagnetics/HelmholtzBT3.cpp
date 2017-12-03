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

                         Implementation of class HelmholtzBT3
                       for Helmholtz Equation in a Bounded Domain
                         using 3-node triangular finite element

  ==============================================================================*/


#include "equations/electromagnetics/HelmholtzBT3.h"

namespace OFELI {

HelmholtzBT3::HelmholtzBT3(Element* el)
{
   _nb_dof = 1;
   Init(el);
   _tr = new Triang3(el);
   _ln = NULL;
   _lx[0] = _ly[0] = _ly[1] = _lx[2] = 0.;
   _lx[1] = _ly[2] = 1.;
   ElementNodeCoordinates();
}


HelmholtzBT3::HelmholtzBT3(Side* sd)
{
   _nb_dof = 1;
   Init(sd);
   _ln = new Line2(sd);
   _tr = NULL;
   SideNodeCoordinates();
}


HelmholtzBT3::~HelmholtzBT3()
{
   if (_tr) { delete _tr; _tr = NULL; }
   if (_ln) { delete _ln; _ln = NULL; }
}


void HelmholtzBT3::LHS(real_t wave_nb)
{
   eMat = complex_t(0);
   eRHS = complex_t(0);
   real_t c = wave_nb*wave_nb*_tr->getArea()*OFELI_SIXTH*0.5;
   for (size_t i=1; i<=3; i++) {
      real_t aa = _tr->getArea()*_tr->DSh(i).x;
      real_t bb = _tr->getArea()*_tr->DSh(i).y;
      for (size_t j=1; j<=3; j++) {
         eMat(i,j) += aa*_tr->DSh(j).x + bb*_tr->DSh(j).y;
         eMat(i,j) -= c;
      }
      eMat(i,i) -= c;
   }
}


void HelmholtzBT3::BoundaryRHS(UserData<complex_t>& ud)
{
   sMat = complex_t(0);
   sRHS = complex_t(0);
   complex<real_t> f = ud.SurfaceForce(_x[0],_theSide->getCode(1));
   sRHS(1) += f*0.5*_ln->getLength();
   f = ud.SurfaceForce(_x[1],_theSide->getCode(1));
   sRHS(2) += f*0.5*_ln->getLength();
}

} /* namespace OFELI */
