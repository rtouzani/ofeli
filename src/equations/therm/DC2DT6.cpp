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

                       Class DC2DT6: Diffusion-Convection Element
                using 6-Node Triangular Finite element in two dimensions

  ==============================================================================*/


#include "equations/therm/DC2DT6.h"
#include "equations/Equation_impl.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

DC2DT6::DC2DT6()
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 6-Node Triangles (P2)";
}


DC2DT6::DC2DT6(Mesh& ms) 
       : Equation<6,6,3,3>(ms)
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 6-Node Triangles (P2)";
   _stab = false;
   setMatrixType(SPARSE);
   setSolver(BICG_STAB_SOLVER,DILU_PREC);
}


DC2DT6::DC2DT6(Mesh&         ms,
               Vect<real_t>& u)
       : Equation<6,6,3,3>(ms,u)
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 6-Node Triangles (P2)";
   _stab = false;
   setMatrixType(SPARSE);
   setSolver(BICG_STAB_SOLVER,DILU_PREC);
}


void DC2DT6::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   setMaterial();
   Triang6S tr(el);
   _el_geo.area = tr.getArea();
   _el_geo.center = tr.getCenter();
   ElementNodeCoordinates();
   tr.atMidEdges(_dSh,_wg);
   ElementNodeVector(*_u,_eu);
   if (_rho_set)
      _rho = _rho_fct(_el_geo.center,_TimeInt.time);
   if (_Cp_set)
      _cp = _Cp_fct(_el_geo.center,_TimeInt.time);
   if (_kappa_set)
      _diff = _kappa_fct(_el_geo.center,_TimeInt.time);
   eA0 = 0, eA1 = 0;
   eRHS = 0;
   _a3 = OFELI_THIRD*_el_geo.area;
}


void DC2DT6::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   _label = sd->n();
   SideNodeCoordinates();
   sRHS = 0;
}


void DC2DT6::LCapacity(real_t coef)
{
   real_t c = 0.5*_a3*_rho*_cp*coef;
   eA1(1,1) += c; eA1(2,2) += c;
   eA1(3,3) += c; eA1(4,4) += c;
   eA1(5,5) += c; eA1(6,6) += c;
}


void DC2DT6::Capacity(real_t coef)
{
   real_t c = OFELI_SIXTH*_el_geo.area*_rho*_cp*coef;
   real_t d = 0.5*c;
   eA1(1,1) += c;
   eA1(2,2) += c;
   eA1(3,3) += c;
   eA1(4,4) += c;
   eA1(5,5) += c;
   eA1(6,6) += c;
   eA1(1,2) += d; eA1(1,3) += d; eA1(1,4) += d; eA1(1,5) += d; eA1(1,6) += d;
   eA1(2,1) += d; eA1(2,3) += d; eA1(2,4) += d; eA1(2,5) += d; eA1(2,6) += d;
   eA1(3,1) += d; eA1(3,2) += d; eA1(3,4) += d; eA1(3,5) += d; eA1(3,6) += d;
   eA1(4,1) += d; eA1(4,2) += d; eA1(4,3) += d; eA1(4,5) += d; eA1(4,6) += d;
   eA1(5,1) += d; eA1(5,2) += d; eA1(5,3) += d; eA1(5,4) += d; eA1(5,6) += d;
   eA1(6,1) += d; eA1(6,2) += d; eA1(6,3) += d; eA1(6,4) += d; eA1(6,5) += d;
}


void DC2DT6::Diffusion(real_t coef)
{
   real_t a = coef*_diff;
   for (size_t k=0; k<3; k++) {
      for (size_t i=1; i<=6; i++) {
         for (size_t j=1; j<=6; j++)
            eA0(i,j) += a*(_dSh[3*(i-1)+k],_dSh[3*(j-1)+k]);
      }
   }
}


void DC2DT6::Convection(Point<real_t>& v,
                        real_t         coef)
{
   for (size_t j=1; j<=6; j++) {
      eA0(4,j) += _a3*coef*(v,_dSh[3*(j-1)  ]);
      eA0(5,j) += _a3*coef*(v,_dSh[3*(j-1)+1]);
      eA0(6,j) += _a3*coef*(v,_dSh[3*(j-1)+2]);
   }
}


void DC2DT6::BodyRHS(const Vect<real_t>& f)
{
   eRHS(4) += _a3*f((*_theElement)(4)->n());
   eRHS(5) += _a3*f((*_theElement)(5)->n());
   eRHS(6) += _a3*f((*_theElement)(6)->n());
}


void DC2DT6::BoundaryRHS(const Vect<real_t>& f)
{
   Line3 ln(_theSide);
   real_t c = OFELI_THIRD*ln.getDet();
   ln.setLocal(-1.0);
   sRHS(1) += c*f((*_theSide)(1)->n());
   ln.setLocal(0.0);
   sRHS(2) += c*f((*_theSide)(2)->n());
   ln.setLocal(1.0);
   sRHS(3) += 4*c*f((*_theSide)(3)->n());
}

} /* namespace OFELI */
