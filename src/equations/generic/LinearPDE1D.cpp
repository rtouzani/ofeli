/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2026 Rachid Touzani

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

                       Implementation of class LinearPDE1D

  ==============================================================================*/


#include "equations/generic/LinearPDE1D.h"
#include "shape_functions/Line2.h"
#include "equations/Equation_impl.h"
#include "linear_algebra/Vect_impl.h"
#include "util/constants.h"

namespace OFELI {

LinearPDE1D::LinearPDE1D()
{
   _equation_name = "Linear PDE";
   _finite_element = "1-D, 2-Node Lines (P1)";
   _lump = true;
   _stab = false;
}


LinearPDE1D::LinearPDE1D(Mesh& ms)
            : Equation<2,2,1,1>(ms)
{
   _equation_name = "Linear PDE";
   _finite_element = "1-D, 2-Node Lines (P1)";
   _theMesh = &ms;
   setMatrixType(TRIDIAGONAL);
   setSolver(DIRECT_SOLVER);
   _lump = true;
   _stab = false;
}


LinearPDE1D::LinearPDE1D(Mesh&         ms,
                         Vect<real_t>& u)
            : Equation<2,2,1,1>(ms,u)
{
   _equation_name = "Linear PDE";
   _finite_element = "1-D, 2-Node Lines (P1)";
   setMatrixType(TRIDIAGONAL);
   setSolver(DIRECT_SOLVER);
   _lump = true;
   _stab = false;
}


void LinearPDE1D::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   Line2 ln(_theElement);
   _el_geo.center = ln.getCenter();
   _el_geo.length = ln.getLength();
   _dSh = ln.DSh();
   ElementNodeCoordinates();
   ElementNodeVector(*_u,_eu);
   eA0.clear();
   eA1.clear();
   eA2.clear();
   eRHS.clear();
}


void LinearPDE1D::Mat_00(real_t coef)
{
   real_t c = OFELI_SIXTH*_el_geo.length*coef;
   c *= getPDECoef(PDECoefType::C00,SpaceTime(_el_geo.center.x,0.,0.,_TimeInt.time));
   eA0(1,1) += 2*c;
   eA0(2,2) += 2*c;
   eA0(1,2) +=   c;
   eA0(2,1) +=   c;
}


void LinearPDE1D::Mat_10(real_t coef)
{
   if (_lump) {
      real_t c = 0.5*_el_geo.length*coef;
      c *= getPDECoef(PDECoefType::C10,SpaceTime(_el_geo.center.x,0.,0.,_TimeInt.time));
      eA1(1,1) += c;
      eA1(2,2) += c;
   }
   else {
      real_t c = OFELI_SIXTH*_el_geo.length*coef;
      c *= getPDECoef(PDECoefType::C10,SpaceTime(_el_geo.center.x,0.,0.,_TimeInt.time));
      eA1(1,1) += 2*c;
      eA1(2,2) += 2*c;
      eA1(1,2) +=   c;
      eA1(2,1) +=   c;
   }
}


void LinearPDE1D::Mat_20(real_t coef)
{
   if (_lump) {
      real_t c = 0.5*_el_geo.length*coef;
      c *= getPDECoef(PDECoefType::C20,SpaceTime(_el_geo.center.x,0.,0.,_TimeInt.time));
      eA2(1,1) += c;
      eA2(2,2) += c;
   }
   else {
      real_t c = OFELI_SIXTH*_el_geo.length*coef;
      c *= getPDECoef(PDECoefType::C20,SpaceTime(_el_geo.center.x,0.,0.,_TimeInt.time));
      eA2(1,1) += 2*c;
      eA2(2,2) += 2*c;
      eA2(1,2) +=   c;
      eA2(2,1) +=   c;
   }
}


void LinearPDE1D::Mat_02(real_t coef)
{
   real_t c = coef*_el_geo.length;
   c *= getPDECoef(PDECoefType::C02,SpaceTime(_el_geo.center.x,0.,0.,_TimeInt.time));
   eA0(1,1) += c*_dSh[0].x*_dSh[0].x;
   eA0(1,2) += c*_dSh[0].x*_dSh[1].x;
   eA0(2,1) += c*_dSh[1].x*_dSh[0].x;
   eA0(2,2) += c*_dSh[1].x*_dSh[1].x;
}


void LinearPDE1D::Mat_01(real_t coef)
{
   real_t c01 = getPDECoef(PDECoefType::C01,SpaceTime(_el_geo.center.x,0.,0.,_TimeInt.time));
   real_t c = 0.25*_el_geo.length*c01*coef;
   eA0(1,1) += c*_dSh[0].x;
   eA0(1,2) += c*_dSh[1].x;
   eA0(2,1) += c*_dSh[0].x;
   eA0(2,2) += c*_dSh[1].x;
   if (_stab) {
      c = coef*_el_geo.length*_el_geo.length/c01;
      real_t dd1=c01*_dSh[0].x, dd2=c01*_dSh[1].x;
      eA0(1,1) += c*dd1*dd1;
      eA0(1,2) += c*dd1*dd2;
      eA0(2,1) += c*dd1*dd2;
      eA0(2,2) += c*dd2*dd2;
   }
}


void LinearPDE1D::BodyRHS(const Vect<real_t>& f)
{
   eRHS(1) += 0.5*_el_geo.length*f((*_theElement)(1)->n());
   eRHS(2) += 0.5*_el_geo.length*f((*_theElement)(2)->n());
}


real_t LinearPDE1D::Flux() const
{
   return _eu(1)*_dSh[0].x + _eu(2)*_dSh[1].x;
}

} /* namespace OFELI */
