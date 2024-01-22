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

                     Implementation of class Pres1DL2

  ==============================================================================*/


#include "equations/acoustics/Pres1DL2.h"
#include "shape_functions/Line2.h"
#include "equations/Equation_impl.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

Pres1DL2::Pres1DL2()
{
   _equation_name = "Acoustics";
   _finite_element = "1-D, 2-Node Lines (P1)";
}


Pres1DL2::Pres1DL2(Mesh& ms)
         : Equation<2,2,1,1>(ms)
{
   _theMesh = &ms;
   setMatrixType(TRIDIAGONAL);
   setSolver(DIRECT_SOLVER);
}


Pres1DL2::Pres1DL2(Mesh&         ms,
                  Vect<real_t>& u)
         : Equation<2,2,1,1>(ms,u)
{
   setMatrixType(TRIDIAGONAL);
   setSolver(DIRECT_SOLVER);
}


Pres1DL2::~Pres1DL2() { }


void Pres1DL2::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   Line2 ln(_theElement);
   _el_geo.center = ln.getCenter();
   _dSh = ln.DSh();
   ElementNodeCoordinates();
   ElementNodeVector(*_u,_eu);
   eA0 = 0, eA1 = 0, eA2 = 0;
   eRHS = 0;
}


void Pres1DL2::setInput(EType         opt,
                        Vect<real_t>& u)
{
   Equa::setInput(opt,u);
}


void Pres1DL2::LMass(real_t coef)
{
   real_t c = 0.5*_el_geo.length*coef/(_speed*_speed);
   eA1(1,1) += c;
   eA1(2,2) += c;
}


void Pres1DL2::Mass(real_t coef)
{
   real_t c = OFELI_SIXTH*_el_geo.length/(_speed*_speed)*coef;
   eA2(1,1) += 2*c;
   eA2(2,2) += 2*c;
   eA2(1,2) +=   c;
   eA2(2,1) +=   c;
}


void Pres1DL2::Diffusion(real_t coef)
{
   eA0(1,1) += coef*_el_geo.length*_dSh[0].x*_dSh[0].x;
   eA0(1,2) += coef*_el_geo.length*_dSh[0].x*_dSh[1].x;
   eA0(2,1) += coef*_el_geo.length*_dSh[1].x*_dSh[0].x;
   eA0(2,2) += coef*_el_geo.length*_dSh[1].x*_dSh[1].x;
}


void Pres1DL2::BodyRHS(const Vect<real_t>& f)
{
   eRHS(1) += 0.5*_el_geo.length*f((*_theElement)(1)->n());
   eRHS(2) += 0.5*_el_geo.length*f((*_theElement)(2)->n());
}


real_t Pres1DL2::Flux() const
{
   return _eu(1)*_dSh[0].x + _eu(2)*_dSh[1].x;
}

} /* namespace OFELI */
