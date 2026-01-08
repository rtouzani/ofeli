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

                       Implementation of class LinearPDE3D

  ==============================================================================*/


#include "equations/generic/LinearPDE3D.h"
#include "shape_functions/Tetra4.h"
#include "shape_functions/Triang3.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/LocalMatrix_impl.h"
#include "linear_algebra/LocalVect_impl.h"
#include "equations/Equation_impl.h"


namespace OFELI {

LinearPDE3D::LinearPDE3D()
{
   _equation_name = "Linear PDE";
   _finite_element = "3-D, 4-Node Tetrahedra (P1)";
   _lump = true;
   _stab = false;
}


LinearPDE3D::LinearPDE3D(Mesh& ms) 
            : Equation<4,4,3,3>(ms)
{
   _equation_name = "Linear PDE";
   _finite_element = "3-D, 4-Node Tetrahedra (P1)";
   setMatrixType(SPARSE);
   setSolver(CG_SOLVER,DILU_PREC);
   _lump = true;
   _stab = false;
}


LinearPDE3D::LinearPDE3D(Mesh&         ms,
                         Vect<real_t>& u)
            : Equation<4,4,3,3>(ms,u)
{
   _equation_name = "Linear PDE";
   _finite_element = "3-D, 4-Node Tetrahedra (P1)";
   setMatrixType(SPARSE);
   setSolver(CG_SOLVER,DILU_PREC);
   _lump = true;
   _stab = false;
}


void LinearPDE3D::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   Tetra4 t(_theElement);
   _el_geo.volume = t.getVolume();
   _el_geo.det = t.getDet();
   _el_geo.center = t.getCenter();
   _dSh = t.DSh();
   ElementNodeCoordinates();
   ElementNodeVector(*_u,_eu);
   eA0.clear();
   eA1.clear();
   eA2.clear();
   eRHS.clear();
}


void LinearPDE3D::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Triang3 t(sd);
   _el_geo.area = t.getArea();
   SideNodeCoordinates();
   sA0 = 0;
   sRHS = 0;
}


void LinearPDE3D::Mat_00(real_t coef)
{
   real_t c = 0.1*_el_geo.volume*coef;
   real_t d = 0.5*c;
   eA0(1,1) += c; eA0(2,2) += c; eA0(3,3) += c; eA0(4,4) += c;
   eA0(1,2) += d; eA0(2,1) += d; eA0(1,3) += d; eA0(1,4) += d;
   eA0(3,1) += d; eA0(2,3) += d; eA0(3,2) += d; eA0(2,4) += d;
   eA0(4,1) += d; eA0(4,2) += d; eA0(4,3) += d; eA0(3,4) += d;
}


void LinearPDE3D::Mat_10(real_t coef)
{
   if (_lump) {
      real_t c = coef*0.25*_el_geo.volume;
      for (size_t i=1; i<=4; i++)
         eA1(i,i) += c;
   }
   else {
      real_t c = 0.1*_el_geo.volume*coef;
      real_t d = 0.5*c;
      eA1(1,1) += c; eA1(2,2) += c; eA1(3,3) += c; eA1(4,4) += c;
      eA1(1,2) += d; eA1(2,1) += d; eA1(1,3) += d; eA1(1,4) += d;
      eA1(3,1) += d; eA1(2,3) += d; eA1(3,2) += d; eA1(2,4) += d;
      eA1(4,1) += d; eA1(4,2) += d; eA1(4,3) += d; eA1(3,4) += d;
   }
}


void LinearPDE3D::Mat_20(real_t coef)
{
   if (_lump) {
      real_t c = coef*0.25*_el_geo.volume;
      for (size_t i=1; i<=4; i++)
         eA2(i,i) += c;
   }
   else {
      real_t c = 0.1*_el_geo.volume*coef;
      real_t d = 0.5*c;
      eA2(1,1) += c; eA2(2,2) += c; eA2(3,3) += c; eA2(4,4) += c;
      eA2(1,2) += d; eA2(2,1) += d; eA2(1,3) += d; eA2(1,4) += d;
      eA2(3,1) += d; eA2(2,3) += d; eA2(3,2) += d; eA2(2,4) += d;
      eA2(4,1) += d; eA2(4,2) += d; eA2(4,3) += d; eA2(3,4) += d;
   }
}


void LinearPDE3D::Mat_01(real_t coef)
{
   LocalVect<real_t,4> dd;
   if (_type01==2) {
      for (size_t i=1; i<=4; i++)
         dd(i) = (_d01,_dSh[i-1]);
      for (size_t i=1; i<=4; i++)
         for (size_t j=1; j<=4; j++)
            eA0(i,j) += 0.25*coef*_el_geo.volume*dd(j);
      if (_stab) {
         real_t c=coef*_el_geo.volume*_el_geo.size/_d01.Norm(), d;
         for (size_t i=1; i<=4; i++) {
            d = c*dd(i);
            for (size_t j=1; j<=4; j++)
               eA0(i,j) += d*dd(j);
         }
      }
   }
}


void LinearPDE3D::Mat_02(real_t coef)
{
   real_t c = coef*_el_geo.volume;
   if (_type02==1) {
      for (size_t i=1; i<=4; i++)
         for (size_t j=1; j<=4; j++)
            eA0(i,j) += c*_c02*(_dSh[i-1],_dSh[j-1]);
   }
   else if (_type02==2) {
      for (size_t i=1; i<=4; i++)
         for (size_t j=1; j<=4; j++)
            eA0(i,j) += c*(_d02.x*_dSh[i-1].x*_dSh[j-1].x + _d02.y*_dSh[i-1].y*_dSh[j-1].y + _d02.z*_dSh[i-1].z*_dSh[j-1].z);
   }
   else if (_type02==3) {
      c *= _f02(_el_geo.center.x,_el_geo.center.y,_el_geo.center.z,_TimeInt.time);
      for (size_t i=1; i<=4; i++)
         for (size_t j=1; j<=4; j++)
            eA0(i,j) += c*(_dSh[i-1],_dSh[j-1]);
   }
   else if (_type02==4) {
      for (size_t i=1; i<=4; i++)
         for (size_t j=1; j<=4; j++)
          eA0(i,j) += c*(_e02(1,1)*_dSh[i-1].x*_dSh[j-1].x + _e02(1,2)*_dSh[i-1].y*_dSh[j-1].x +
                         _e02(1,3)*_dSh[i-1].z*_dSh[j-1].x + _e02(2,1)*_dSh[i-1].x*_dSh[j-1].y +
                         _e02(2,2)*_dSh[i-1].y*_dSh[j-1].y + _e02(2,3)*_dSh[i-1].z*_dSh[j-1].y +
                         _e02(3,1)*_dSh[i-1].x*_dSh[j-1].z + _e02(3,2)*_dSh[i-1].y*_dSh[j-1].z +
                         _e02(3,3)*_dSh[i-1].z*_dSh[j-1].z);
   }
}


void LinearPDE3D::BodyRHS(const Vect<real_t>& f)
{
   for (size_t i=1; i<=4; ++i)
      eRHS(i) += f((*_theElement)(i)->n())*0.25*_el_geo.volume;
}


void LinearPDE3D::BoundaryRHS(real_t flux)
{
   real_t c = flux*_el_geo.area*OFELI_THIRD;
   sRHS(1) += c;
   sRHS(2) += c;
   sRHS(3) += c;
}


void LinearPDE3D::BoundaryRHS(const Vect<real_t>& f)
{
   if (_theSide->getCode(1)>0) {
      real_t c = _el_geo.area*OFELI_THIRD;
      if (f.getDOFType()==NODE_DOF) {
         for (size_t i=1; i<=3; ++i)
            sRHS(i) += c*f((*_theSide)(i)->n());
      }
      else if (f.getDOFType()==SIDE_DOF) {
         real_t ff = c*f(_theSide->n());
         for (size_t i=1; i<=3; ++i)
            sRHS(i) += ff;
      }
   }
}


Point<real_t> LinearPDE3D::Flux() const
{
   Point<real_t> f;
   f = _eu(1)*_dSh[0] + _eu(2)*_dSh[1] +
       _eu(3)*_dSh[2] + _eu(4)*_dSh[3];
   return f;
}


void LinearPDE3D::Grad(Vect<Point<real_t> >& g)
{
   MESH_EL {
      set(the_element);
      g(element_label) = _eu(1)*_dSh[0] + _eu(2)*_dSh[1] + _eu(3)*_dSh[2] + _eu(4)*_dSh[3];
   }
}

} /* namespace OFELI */
