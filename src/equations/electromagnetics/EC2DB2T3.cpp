/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani

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

                                 Class EC2DB2T3

  ==============================================================================*/


#include "equations/electromagnetics/EC2DB2T3.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"

namespace OFELI {

EC2DB2T3::EC2DB2T3(Element *el)
{
   _label = el->n();
   setMaterial();
   set(el);
}


EC2DB2T3::EC2DB2T3(Side *sd)
{
   _ns = sd->n();
   set(sd);
}


void EC2DB2T3::set(const Element *el)
{
   _theElement = el, _theSide = nullptr;
   setMaterial();
   Triang3 tr(_theElement);
   _center = tr.getCenter();
   _area = tr.getArea();
   _det = 2*_area;
   ElementNodeCoordinates();
   _dSh = tr.DSh();
   if (_omega_set)
      _omega = _omega_fct(_el_geo.center,0.);
   if (_Mu_set)
      _Mu = _Mu_fct(_el_geo.center,0.);
   if (_sigma_set)
      _sigma = _sigma_fct(_el_geo.center,0.);
   eMat = 0;
   eRHS = 0;
}


void EC2DB2T3::set(const Side *sd)
{
   _theElement = nullptr, _theSide = sd;
   Line2 ln(sd);
   SideNodeCoordinates();
   _dSh = ln.DSh();
   _center = ln.getCenter();
   _length = ln.getLength();
   sMat = 0;
   sRHS = 0;
}


void EC2DB2T3::RHS(real_t coef)
{
   eRHS(1) = eRHS(3) = eRHS(5) = coef*OFELI_THIRD*_area*_sigma;
   eRHS(2) = eRHS(4) = eRHS(6) = 0;
}


void EC2DB2T3::EMatr()
{
   real_t c = 0.25*_Mu*_omega*OFELI_THIRD*_area*_sigma;
   for (size_t i=1; i<=3; i++) {
      eMat(2*i-1,2*i  ) -= c;
      eMat(2*i  ,2*i-1) += c;
      for (size_t j=1; j<=3; j++) {
         real_t e = OFELI_THIRD*_area*(_dSh[i-1]*_dSh[j-1]);
         eMat(2*i-1,2*j-1) += e; eMat(2*i  ,2*j  ) += e;
         eMat(2*i-1,2*j  ) -= c; eMat(2*i  ,2*j-1) += c;
      }
   }
}


complex_t EC2DB2T3::Constant(const LocalVect<real_t,6>& u,
                             complex_t                  I)
{
   real_t ur = OFELI_THIRD*_area*(u(1)+u(3)+u(5));
   real_t ui = OFELI_THIRD*_area*(u(2)+u(4)+u(6));
   real_t c1 = (I.real() - _omega*ui*_sigma)/_area/_sigma;
   real_t c2 = (I.imag() + _omega*ur*_sigma)/_area/_sigma;
   return complex_t(c1,c2);
}


real_t EC2DB2T3::MagneticPressure(const LocalVect<real_t,6>& u)
{
   real_t c = 0.5*_Mu*_area;
   real_t dxr=0, dxi=0, dyr=0, dyi=0;
   for (size_t i=0; i<3; i++) {
      dxr += _dSh[i].x*u(2*i  );
      dxi += _dSh[i].x*u(2*i+1);
      dyr += _dSh[i].y*u(2*i  );
      dyi += _dSh[i].y*u(2*i+1);
   }
   return (c*(dxr*dxr+dxi*dxi+dyr*dyr+dyi*dyi));
}

} /* namespace OFELI */
