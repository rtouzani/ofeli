/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

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

                         Implementation of class ELAS3DT4
  for 3-D Linear Elasticity Equations using 4-node tetrahedral P1 finite element

  ==============================================================================*/


#include "equations/solid/Elas3DT4.h"
#include "shape_functions/Tetra4.h"
#include "shape_functions/Triang3.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {


Elas3DT4::Elas3DT4(Mesh& ms)
         : Equation<real_t,4,12,3,9>(ms)
{
   _equation_name = "Linearized elasticity";
   _finite_element = "3-D, 4-Node tetrahedrals (P1)";
   setMatrixType(SPARSE|SYMMETRIC);
   setSolver(CG_SOLVER,DIAG_PREC);
   _terms = DEVIATORIC|DILATATION|BODY_FORCE|TRACTION;
}


Elas3DT4::Elas3DT4(Mesh&         ms,
                   Vect<real_t>& u)
         : Equation<real_t,4,12,3,9>(ms,u)
{
   _equation_name = "Linearized elasticity";
   _finite_element = "3-D, 4-Node tetrahedrals (P1)";
   setMatrixType(SPARSE|SYMMETRIC);
   setSolver(CG_SOLVER,DIAG_PREC);
   _terms = DEVIATORIC|DILATATION|BODY_FORCE|TRACTION;
}


Elas3DT4::~Elas3DT4() { }


void Elas3DT4::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   setMaterial();
   Tetra4 tetra(el);
   _el_geo.det = tetra.getDet();
   _el_geo.center = tetra.getCenter();
   ElementNodeCoordinates();
   if (AbsEqua<real_t>::_u!=nullptr)
      ElementNodeVector(*_u,_eu);
   if (AbsEqua<real_t>::_bf!=nullptr)
      ElementNodeVector(*_bf,_ebf);
   _dSh = tetra.DSh();
   _lambda = _nu*_E/((1+_nu)*(1-2*_nu));
   _G = 0.5*_E/(1+_nu);
   eA0 = 0, eA1 = 0, eA2 = 0;
   eRHS = 0;
}


void Elas3DT4::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Triang3 tr(sd);
   _el_geo.area = tr.getArea();
   SideNodeCoordinates();
   if (AbsEqua<real_t>::_u!=nullptr)
      SideNodeVector(*_u,_su);
   if (AbsEqua<real_t>::_sf!=nullptr)
      SideVector(*_sf,_ssf);
   sA0 = 0;
   sRHS = 0;
}


void Elas3DT4::Media(real_t E,
                     real_t nu,
                     real_t rho)
{
   _E = E;
   _nu = nu;
   _rho = rho;
}


void Elas3DT4::LMass(real_t coef)
{
   real_t c = 0.5*coef*OFELI_TWELVETH*_el_geo.det*_rho;
   for (size_t i=1; i<=4; i++) {
      eA2(3*i-2,3*i-2) += c;
      eA2(3*i-1,3*i-1) += c;
      eA2(3*i  ,3*i  ) += c;
   }
}


void Elas3DT4::Deviator(real_t coef)
{
   real_t c=_G*_el_geo.det*coef;
   for (size_t j=1; j<=4; j++) {
      Point<real_t> db=c*_dSh[j-1];
      for (size_t i=1; i<=4; i++) {
         eA0(3*i-2,3*j-2) += 2*_dSh[i-1].x*db.x + _dSh[i-1].z*db.z + _dSh[i-1].y*db.y;
         eA0(3*i-2,3*j-1) += _dSh[i-1].y*db.x;
         eA0(3*i-2,3*j  ) += _dSh[i-1].z*db.x;
         eA0(3*i-1,3*j-2) += _dSh[i-1].x*db.y;
         eA0(3*i-1,3*j-1) += 2*_dSh[i-1].y*db.y + _dSh[i-1].z*db.z + _dSh[i-1].x*db.x;
         eA0(3*i-1,3*j  ) += _dSh[i-1].z*db.y;
         eA0(3*i  ,3*j-2) += _dSh[i-1].x*db.z;
         eA0(3*i  ,3*j-1) += _dSh[i-1].y*db.z;
         eA0(3*i  ,3*j  ) += 2*_dSh[i-1].z*db.z + _dSh[i-1].y*db.y + _dSh[i-1].x*db.x;
      }
   }
}


void Elas3DT4::Dilatation(real_t coef)
{
   real_t c = _lambda*_el_geo.det*coef;
   for (size_t j=1; j<=4; j++) {
      for (size_t i=1; i<=4; i++) {
         eA0(3*i-2,3*j-2) += c*_dSh[i-1].x*_dSh[j-1].x;
         eA0(3*i-2,3*j-1) += c*_dSh[i-1].x*_dSh[j-1].y;
         eA0(3*i-2,3*j  ) += c*_dSh[i-1].x*_dSh[j-1].z;
         eA0(3*i-1,3*j-2) += c*_dSh[i-1].y*_dSh[j-1].x;
         eA0(3*i-1,3*j-1) += c*_dSh[i-1].y*_dSh[j-1].y;
         eA0(3*i-1,3*j  ) += c*_dSh[i-1].y*_dSh[j-1].z;
         eA0(3*i  ,3*j-2) += c*_dSh[i-1].z*_dSh[j-1].x;
         eA0(3*i  ,3*j-1) += c*_dSh[i-1].z*_dSh[j-1].y;
         eA0(3*i  ,3*j  ) += c*_dSh[i-1].z*_dSh[j-1].z;
      }
   }
}


void Elas3DT4::BodyRHS(const Vect<real_t>& f)
{
   real_t c = 0.5*OFELI_TWELVETH*_el_geo.det;
   for (size_t i=1; i<=4; ++i) {
      eRHS(3*i-2) += c*f((*_theElement)(i)->n(),1);
      eRHS(3*i-1) += c*f((*_theElement)(i)->n(),2);
      eRHS(3*i  ) += c*f((*_theElement)(i)->n(),3);
   }
}


void Elas3DT4::BodyRHS()
{
   real_t c = _el_geo.det/24.0;
   for (size_t i=1; i<=4; ++i) {
      eRHS(3*i-2) += c*_ebf(3*i-2);
      eRHS(3*i-1) += c*_ebf(3*i-1);
      eRHS(3*i  ) += c*_ebf(3*i  );
   }
}


void Elas3DT4::BoundaryRHS(const Vect<real_t>& f)
{
   for (size_t i=1; i<=3; ++i) {
      if (_theSide->getCode(1) != CONTACT)
         sRHS(3*i-2) += _el_geo.area*f(_theSide->n(),1);
      if (_theSide->getCode(2) != CONTACT)
         sRHS(3*i-1) += _el_geo.area*f(_theSide->n(),2);
      if (_theSide->getCode(3) != CONTACT)
         sRHS(3*i  ) += _el_geo.area*f(_theSide->n(),3);
   }
}


void Elas3DT4::BoundaryRHS()
{
   for (size_t i=1; i<=3; ++i) {
      if (_theSide->getCode(1) != CONTACT)
         sRHS(3*i-2) += _el_geo.area*_ssf[0];
      if (_theSide->getCode(2) != CONTACT)
         sRHS(3*i-1) += _el_geo.area*_ssf[1];
      if (_theSide->getCode(3) != CONTACT)
         sRHS(3*i  ) += _el_geo.area*_ssf[2];
   }
}

} /* namespace OFELI */
