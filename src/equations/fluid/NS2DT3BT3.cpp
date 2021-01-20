/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

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

                        Implementation of class NS2DT3BT3
      for 2-D Navier-Stokes equations using P1-Bubble/P1 (Mini) finite element

  ==============================================================================*/


#include "equations/fluid/NS2DT3BT3.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

NS2DT3BT3::NS2DT3BT3(Mesh& ms)
          : Equation<real_t,3,9,2,6>(ms)
{
   _equation_name = "Navier-Stokes";
   _finite_element = "2-D, P1+Bubble/P1";
   _Re = 1.;
}


NS2DT3BT3::NS2DT3BT3(Mesh&         ms,
                     Vect<real_t>& u)
          : Equation<real_t,3,9,2,6>(ms,u)
{
   _equation_name = "Navier-Stokes";
   _finite_element = "2-D, P1+Bubble/P1";
   _Re = 1.;
}


NS2DT3BT3::~NS2DT3BT3() { }


void NS2DT3BT3::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   setMaterial();
   Triang3 tr(_theElement);
   _el_geo.area = tr.getArea();
   _el_geo.center = tr.getCenter();
   ElementNodeCoordinates();
   _dSh = tr.DSh();
   if (_rho_set)
      _rho = _rho_fct(_el_geo.center,_TimeInt.time);
   if (_mu_set)
      _mu = _mu_fct(_el_geo.center,_TimeInt.time);
   if (_beta_set)
      _beta = _beta_fct(_el_geo.center,_TimeInt.time);
   eA0 = 0, eA1 = 0;
   eRHS = 0;
}


void NS2DT3BT3::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Line2 ln(_theSide);
   _el_geo.length = ln.getLength();
   SideNodeCoordinates();
   _dSh = ln.DSh();
   sRHS = 0;
}


void NS2DT3BT3::build()
{
   MESH_EL {
      set(the_element);
      Misc();
      Viscous();
      PressureGradient();
      if (AbsEqua<real_t>::_bc!=nullptr)
         Equation<real_t,3,9,2,6>::updateBC(The_element,*AbsEqua<real_t>::_bc);
      AbsEqua<real_t>::_A->Assembly(The_element,eA0.get());
      AbsEqua<real_t>::_b->Assembly(The_element,eRHS.get());
   }
}


void NS2DT3BT3::Misc()
{
   Point<real_t> x32=_x[2]-_x[1], x13=_x[0]-_x[2], x21=_x[1]-_x[0];
   real_t n13=(x13,x13), n21=(x21,x21), m32=(x13,x21);
   real_t a = 0.25*_mu/_el_geo.area;

   _aa(1,1) = a*x32.NNorm();
   _aa(2,1) = -a*(n13 + m32);
   _aa(1,2) = _aa(1,2);
   _aa(3,1) = -a*(n13 + m32);
   _aa(1,3) = _aa(3,1);
   _aa(2,2) = a*n13;
   _aa(3,2) = a*m32;
   _aa(2,3) = _aa(3,2);
   _aa(3,3) = a*n21;
   _aa(4,4) = a*(n13 + n21 + m32)/90.;

   _bb(1,1) = -OFELI_SIXTH*x32.y; 
   _bb(2,1) = _bb(1,1); 
   _bb(3,1) = _bb(1,1); 
   _bb(1,4) = 0.1*OFELI_TWELVETH*x32.y;
   _bb(1,2) = -OFELI_SIXTH*x13.y; 
   _bb(2,2) = _bb(1,2);
   _bb(3,2) = -OFELI_SIXTH*x13.y; 
   _bb(2,4) = 0.1*OFELI_TWELVETH*x13.y;
   _bb(1,3) = -OFELI_SIXTH*x21.y; 
   _bb(2,3) = -OFELI_SIXTH*x21.y; 
   _bb(3,3) = -OFELI_SIXTH*x21.y; 
   _bb(3,4) = 0.1*OFELI_TWELVETH*x21.y;

   _cc(1,1) = OFELI_SIXTH*x32.x;
   _cc(2,1) = _cc(1,1);
   _cc(3,1) = _cc(1,1);
   _cc(1,2) = OFELI_SIXTH*x13.x;
   _cc(2,2) = _cc(1,2);
   _cc(3,2) = _cc(1,2);
   _cc(1,3) = OFELI_SIXTH*x21.x;
   _cc(2,3) = _cc(1,3);
   _cc(3,3) = _cc(1,3);
   _cc(1,4) = -0.1*OFELI_TWELVETH*x32.x;
   _cc(2,4) = -0.1*OFELI_TWELVETH*x13.x;
   _cc(3,4) = -0.1*OFELI_TWELVETH*x21.x;
}


void NS2DT3BT3::LMass(real_t coef)
{
   real_t c = coef*_rho*OFELI_THIRD*_el_geo.area;
   for (size_t i=1; i<=3; i++) {
      eA1(3*i-2,3*i-2) += c;
      eA1(3*i-1,3*i-1) += c;
      eA1(3*i  ,3*i  ) += c;
   }
}


void NS2DT3BT3::Viscous(real_t coef)
{
   for (size_t i=0; i<3; i++) {
      eA0(3*i+1,1) = _aa(1,i+1)*coef;
      eA0(3*i+2,2) = _aa(1,i+1)*coef;
      eA0(3*i+1,4) = _aa(2,i+1)*coef;
      eA0(3*i+2,5) = _aa(2,i+1)*coef;
      eA0(3*i+1,7) = _aa(3,i+1)*coef;
      eA0(3*i+2,8) = _aa(3,i+1)*coef;
   }
}


void NS2DT3BT3::PressureGradient(real_t coef)
{
   for (size_t i=0; i<3; i++) {
      eA0(3*i+1,3) = -_bb(1,i+1);
      eA0(3*i+1,6) = -_bb(2,i+1);
      eA0(3*i+1,9) = -_bb(3,i+1);
      eA0(3*i+2,3) = -_cc(1,i+1);
      eA0(3*i+2,6) = -_cc(2,i+1);
      eA0(3*i+2,9) = -_cc(3,i+1);
      eA0(3*i+3,1) =  _bb(i+1,1);
      eA0(3*i+3,2) =  _cc(i+1,1);
      eA0(3*i+3,3) = (_bb(i+1,4)*_bb(1,4) + _cc(i+1,4)*_cc(1,4))/_aa(4,4);
      eA0(3*i+3,4) =  _bb(i+1,2);
      eA0(3*i+3,5) =  _cc(i+1,2);
      eA0(3*i+3,6) = (_bb(i+1,4)*_bb(2,4) + _cc(i+1,4)*_cc(2,4))/_aa(4,4);
      eA0(3*i+3,7) =  _bb(i+1,3);
      eA0(3*i+3,8) =  _cc(i+1,3);
      eA0(3*i+3,9) = (_bb(i+1,4)*_bb(3,4) + _cc(i+1,4)*_cc(3,4))/_aa(4,4);
   }
   for (size_t i=0; i<3; i++) {
      eA0(3*i+1,3) = -_bb(1,i+1);
      eA0(3*i+1,6) = -_bb(2,i+1);
      eA0(3*i+1,9) = -_bb(3,i+1);
      eA0(3*i+2,3) = -_cc(1,i+1);
      eA0(3*i+2,6) = -_cc(2,i+1);
      eA0(3*i+2,9) = -_cc(3,i+1);
      eA0(3*i+3,1) =  _bb(i+1,1);
      eA0(3*i+3,2) =  _cc(i+1,1);
      eA0(3*i+3,3) = (_bb(i+1,4)*_bb(1,4) + _cc(i+1,4)*_cc(1,4))/_aa(4,4);
      eA0(3*i+3,4) =  _bb(i+1,2);
      eA0(3*i+3,5) =  _cc(i+1,2);
      eA0(3*i+3,6) = (_bb(i+1,4)*_bb(2,4) + _cc(i+1,4)*_cc(2,4))/_aa(4,4);
      eA0(3*i+3,7) =  _bb(i+1,3);
      eA0(3*i+3,8) =  _cc(i+1,3);
      eA0(3*i+3,9) = (_bb(i+1,4)*_bb(3,4) + _cc(i+1,4)*_cc(3,4))/_aa(4,4);
   }
}


void NS2DT3BT3::RHS_Convection(real_t coef)
{
   Point<real_t> du = _eu(1)*_dSh[0] + _eu(3)*_dSh[1] + _eu(5)*_dSh[2];
   Point<real_t> dv = _eu(2)*_dSh[0] + _eu(4)*_dSh[1] + _eu(6)*_dSh[2];
   real_t c = coef*OFELI_THIRD*_el_geo.area;
   for (size_t i=1; i<=3; i++) {
      eRHS(3*i-2) -= c*(_eu(2*i-1)*du.x + _eu(2*i)*du.y);
      eRHS(3*i-1) -= c*(_eu(2*i-1)*dv.x + _eu(2*i)*dv.y);
   }
}


void NS2DT3BT3::BodyRHS(Vect<real_t>& f)
{
   for (size_t i=1; i<=3; i++) {
      eRHS(3*i-2) += OFELI_THIRD*_el_geo.area*f((*_theElement)(i)->n(),1);
      eRHS(3*i-1) += OFELI_THIRD*_el_geo.area*f((*_theElement)(i)->n(),2);
   }
}


void NS2DT3BT3::BoundaryRHS(Vect<real_t>& f)
{
   for (size_t i=1; i<=2; i++) {
      sRHS(3*i-2) += 0.5*_el_geo.length*f((*_theSide)(i)->n(),1);
      sRHS(3*i-1) += 0.5*_el_geo.length*f((*_theSide)(i)->n(),2);
   }
}

} /* namespace OFELI */
