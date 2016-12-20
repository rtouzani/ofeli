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
   _nb_dof = 2;
   Init(el);
   set(el);
}


EC2DB2T3::EC2DB2T3(Side *sd)
{
   _ns = sd->n();
   _nb_dof = 2;
   set(sd);
}


void EC2DB2T3::set(const Element *el)
{
   _nb_dof = 1;
   Init(el);
   setMaterial();
   Triang3 tr(_theElement);
   _center = tr.getCenter();
   _area = tr.getArea();
   _det = 2*_area;
   ElementNodeCoordinates();
   _N[0] = tr.DSh(1);
   _N[1] = tr.DSh(2);
   _N[2] = tr.DSh(3);
   eMat = 0;
   eRHS = 0;
}


void EC2DB2T3::set(const Side *sd)
{
   Init(sd);
   Line2 ln(sd);
   SideNodeCoordinates();
   _N[0] = ln.DSh(1);
   _N[1] = ln.DSh(2);
   _center = ln.getCenter();
   _length = ln.getLength();
   sMat = 0;
   sRHS = 0;
}


void EC2DB2T3::RHS(real_t coef)
{
   eRHS(1) = eRHS(3) = eRHS(5) = coef*OFELI_THIRD*_area/_rho;
   eRHS(2) = eRHS(4) = eRHS(6) = 0;
}


void EC2DB2T3::EMatr(real_t omega)
{
   real_t c = 0.25*_mu*omega*OFELI_THIRD*_area/_rho;
   for (size_t i=1; i<=3; i++) {
      eMat(2*i-1,2*i  ) -= c;
      eMat(2*i  ,2*i-1) += c;
      for (size_t j=1; j<=3; j++) {
         real_t e = OFELI_THIRD*_area*(_N[i-1]*_N[j-1]);
         eMat(2*i-1,2*j-1) += e; eMat(2*i  ,2*j  ) += e;
         eMat(2*i-1,2*j  ) -= c; eMat(2*i  ,2*j-1) += c;
      }
   }
}


complex_t EC2DB2T3::Constant(      real_t               omega,
                             const LocalVect<real_t,6>& u,
                                   complex_t            I)
{
   real_t ur = OFELI_THIRD*_area*(u(1)+u(3)+u(5));
   real_t ui = OFELI_THIRD*_area*(u(2)+u(4)+u(6));
   real_t c1 = _rho*(I.real() - omega*ui/_rho)/_area;
   real_t c2 = _rho*(I.imag() + omega*ur/_rho)/_area;
   return complex_t (c1,c2);
}


real_t EC2DB2T3::MagneticPressure(const LocalVect<real_t,6>& u)
{
   real_t c = 0.5*_mu*_area;
   real_t dxr=0, dxi=0, dyr=0, dyi=0;
   for (size_t i=0; i<3; i++) {
      dxr += _N[i].x*u(2*i  );
      dxi += _N[i].x*u(2*i+1);
      dyr += _N[i].y*u(2*i  );
      dyi += _N[i].y*u(2*i+1);
   }
   return (c*(dxr*dxr+dxi*dxi+dyr*dyr+dyi*dyi));
}

} /* namespace OFELI */
