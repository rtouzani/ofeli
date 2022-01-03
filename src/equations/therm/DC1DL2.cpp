/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani

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

	           Class DC1DL2 : Diffusion-Convection Element
                  using 2-Node Line Finite element in two dimensions

  ==============================================================================*/


#include "equations/therm/DC1DL2.h"
#include "shape_functions/Line2.h"
#include "equations/Equation_impl.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

DC1DL2::DC1DL2()
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "1-D, 2-Node Lines (P1)";
}


DC1DL2::DC1DL2(Mesh& ms)
       : Equation<2,2,1,1>(ms)
{
   _theMesh = &ms;
   setMatrixType(TRIDIAGONAL);
   setSolver(DIRECT_SOLVER);
}


DC1DL2::DC1DL2(Mesh&         ms,
               Vect<real_t>& u)
       : Equation<2,2,1,1>(ms,u)
{
   setMatrixType(TRIDIAGONAL);
   setSolver(DIRECT_SOLVER);
}


DC1DL2::~DC1DL2() { }


void DC1DL2::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   setMaterial();
   Line2 ln(_theElement);
   _el_geo.center = ln.getCenter();
   _dSh = ln.DSh();
   ElementNodeCoordinates();
   ElementNodeVector(*_u,_eu);
   if (_rho_set)
      _rho = _rho_fct(_el_geo.center,_TimeInt.time);
   if (_Cp_set)
      _cp = _Cp_fct(_el_geo.center,_TimeInt.time);
   if (_kappa_set)
      _diff = _kappa_fct(_el_geo.center,_TimeInt.time);
   eA0 = 0, eA1 = 0;
   eRHS = 0;
}


void DC1DL2::setInput(EqDataType    opt,
                      Vect<real_t>& u)
{
   Equa::setInput(opt,u);
   if (opt==VELOCITY_FIELD)
      _vel = &u;
}


void DC1DL2::LCapacity(real_t coef)
{
   real_t c = 0.5*_el_geo.length*_rho*_cp*coef;
   eA1(1,1) += c;
   eA1(2,2) += c;
}


void DC1DL2::Capacity(real_t coef)
{
   real_t c = OFELI_SIXTH*_el_geo.length*_rho*_cp*coef;
   eA1(1,1) += 2*c;
   eA1(2,2) += 2*c;
   eA1(1,2) +=   c;
   eA1(2,1) +=   c;
}


void DC1DL2::Diffusion(real_t coef)
{
   eA0(1,1) += coef*_diff*_el_geo.length*_dSh[0].x*_dSh[0].x;
   eA0(1,2) += coef*_diff*_el_geo.length*_dSh[0].x*_dSh[1].x;
   eA0(2,1) += coef*_diff*_el_geo.length*_dSh[1].x*_dSh[0].x;
   eA0(2,2) += coef*_diff*_el_geo.length*_dSh[1].x*_dSh[1].x;
}


void DC1DL2::Convection(const real_t& v,
                        real_t        coef)
{
   eA0(1,1) += 0.5*_el_geo.length*coef*v*_dSh[0].x;
   eA0(1,2) += 0.5*_el_geo.length*coef*v*_dSh[1].x;
   eA0(2,1) += 0.5*_el_geo.length*coef*v*_dSh[0].x;
   eA0(2,2) += 0.5*_el_geo.length*coef*v*_dSh[1].x;
}


void DC1DL2::Convection(real_t coef)
{
   real_t c = 0.25*_el_geo.length*coef*((*_vel)((*_theElement)(1)->n()) + (*_vel)((*_theElement)(2)->n()));
   eA0(1,1) += c*_dSh[0].x;
   eA0(1,2) += c*_dSh[1].x;
   eA0(2,1) += c*_dSh[0].x;
   eA0(2,2) += c*_dSh[1].x;
}


void DC1DL2::Convection(const Vect<real_t>& v,
                        real_t              coef)
{
   real_t c = 0.25*_el_geo.length*coef*(v((*_theElement)(1)->n()) + v((*_theElement)(2)->n()));
   eA0(1,1) += c*_dSh[0].x;
   eA0(1,2) += c*_dSh[1].x;
   eA0(2,1) += c*_dSh[0].x;
   eA0(2,2) += c*_dSh[1].x;
}


void DC1DL2::BodyRHS(const Vect<real_t>& f)
{
   eRHS(1) += 0.5*_el_geo.length*f((*_theElement)(1)->n());
   eRHS(2) += 0.5*_el_geo.length*f((*_theElement)(2)->n());
}


real_t DC1DL2::Flux() const
{
   return _diff*(_eu(1)*_dSh[0].x + _eu(2)*_dSh[1].x);
}

} /* namespace OFELI */
