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

                       Class DC3DAT3 : Diffusion-Convection Element
                         using 3-Node Triangular Finite element
                             in Axisymmetric geometries

  ==============================================================================*/


#include "equations/Equation_impl.h"
#include "equations/therm/DC3DAT3.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/LocalMatrix_impl.h"

namespace OFELI {


DC3DAT3::DC3DAT3()
{
}


DC3DAT3::DC3DAT3(Mesh& ms) 
       : Equation<3,3,2,2>(ms)
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 3-Node Axisymmetric Triangles (P1)";
   _stab = false;
   setMatrixType(SPARSE);
   setSolver(GMRES_SOLVER,DILU_PREC);
}


DC3DAT3::DC3DAT3(Mesh&         ms,
                 Vect<real_t>& u)
        : Equation<3,3,2,2>(ms,u)
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 3-Node Axisymmetric Triangles (P1)";
   _stab = false;
   setMatrixType(SPARSE);
   setSolver(GMRES_SOLVER,DILU_PREC);
}


DC3DAT3::~DC3DAT3() { }


void DC3DAT3::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   setMaterial();
   Triang3 tr(_theElement);
   _el_geo.area = tr.getArea();
   _el_geo.center = tr.getCenter();
   _dSh = tr.DSh();
   ElementNodeCoordinates();
   ElementNodeVector(*_u,_eu);
   _r[0] = _x[0].x, _r[1] = _x[1].x, _r[2] = _x[2].x;
   if (_rho_set)
      _rho = _rho_fct(_el_geo.center,_TimeInt.time);
   if (_Cp_set)
      _cp = _Cp_fct(_el_geo.center,_TimeInt.time);
   if (_kappa_set)
      _diff = _kappa_fct(_el_geo.center,_TimeInt.time);
   eA0 = 0, eA1 = 0;
   eRHS = 0;
}


void DC3DAT3::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Line2 ln(sd);
   SideNodeCoordinates();
   _el_geo.length = ln.getLength();
   _r[0] = _x[0].x, _r[1] = _x[1].x;
   sMat = 0;
   sRHS = 0;
}


void DC3DAT3::LCapacity(real_t coef)
{
   real_t c = OFELI_THIRD*_el_geo.area*_rho*_cp*coef;
   eA1(1,1) += c*_r[0];
   eA1(2,2) += c*_r[1];
   eA1(3,3) += c*_r[2];
}


void DC3DAT3::Capacity(real_t coef)
{
   real_t c = 0.5*OFELI_TWELVETH*_el_geo.area*_rho*_cp*coef;
   eA1(1,1) += c*(2*_r[0] + _r[1] + _r[2]);
   eA1(2,2) += c*(_r[0] + 2*_r[1] + _r[2]);
   eA1(3,3) += c*(_r[0] + _r[1] + 2*_r[2]);
   eA1(1,2) += c*(_r[0]+_r[1]); eA1(2,1) += c*(_r[0]+_r[1]);
   eA1(2,3) += c*(_r[1]+_r[2]); eA1(3,2) += c*(_r[1]+_r[2]);
   eA1(1,3) += c*(_r[0]+_r[2]); eA1(3,1) += c*(_r[0]+_r[2]);
}


void DC3DAT3::Diffusion(real_t coef)
{
   real_t c = coef*OFELI_THIRD*_el_geo.area*_diff*(_r[0]+_r[1]+_r[2]);
   for (size_t i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eA0(i,j) += c*(_dSh[i-1],_dSh[j-1]);
}


void DC3DAT3::Diffusion(const LocalMatrix<real_t,2,2>& diff,
                        real_t                         coef)
{
   real_t c = coef*OFELI_THIRD*_el_geo.area*(_r[0]+_r[1]+_r[2]);
   for (size_t i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eA0(i,j) += c*(diff(1,1)*_dSh[i-1].x*_dSh[j-1].x + diff(2,2)*_dSh[i-1].y*_dSh[j-1].y
                      + diff(1,2)*_dSh[i-1].y*_dSh[j-1].x + diff(2,1)*_dSh[i-1].x*_dSh[j-1].y);
}


void DC3DAT3::BodyRHS(const Vect<real_t>& f)
{
   real_t c = 0.125*OFELI_THIRD*_el_geo.area;
   real_t f1 = f((*_theElement)(1)->n()) + f((*_theElement)(3)->n());
   real_t f2 = f((*_theElement)(1)->n()) + f((*_theElement)(2)->n());
   real_t f3 = f((*_theElement)(2)->n()) + f((*_theElement)(3)->n());
   eRHS(1) += c*(f1*(_r[0]+_r[2]) + f2*(_r[0]+_r[1]));
   eRHS(2) += c*(f2*(_r[0]+_r[1]) + f3*(_r[1]+_r[2]));
   eRHS(3) += c*(f1*(_r[0]+_r[2]) + f3*(_r[1]+_r[2]));
}


void DC3DAT3::BoundaryRHS(real_t flux)
{
   sRHS(1) += 0.5*flux*_el_geo.length*_r[0];
   sRHS(2) += 0.5*flux*_el_geo.length*_r[1];
}


void DC3DAT3::BoundaryRHS(const Vect<real_t>& f)
{
   if (_theSide->getCode(1)>0) {
      real_t c = 0.5*_el_geo.length*_r[1];
      if (f.getDOFType()==NODE_DOF) {
         sRHS(1) += c*f((*_theSide)(1)->n());
         sRHS(2) += c*f((*_theSide)(2)->n());
      }
      else if (f.getDOFType()==SIDE_DOF) {
         real_t ff = c*f(_theSide->n());
         sRHS(1) += ff;
         sRHS(2) += ff;
      }
   }
}


Point<real_t> &DC3DAT3::Grad(const Vect<real_t>& u)
{
   _grad = u[0]*_dSh[0] + u[1]*_dSh[1] + u[2]*_dSh[2];
   return _grad;
}

} /* namespace OFELI */
