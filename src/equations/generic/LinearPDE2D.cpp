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

                     Implementation of class LinearPDE2D

  ==============================================================================*/


#include "equations/generic/LinearPDE2D.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/LocalMatrix_impl.h"
#include "linear_algebra/LocalVect_impl.h"
#include "equations/Equation_impl.h"

namespace OFELI {


LinearPDE2D::LinearPDE2D()
{
   _equation_name = "Linear PDE";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _lump = true;
   _stab = false;
}


LinearPDE2D::LinearPDE2D(Mesh& ms)
            : Equation<3,3,2,2>(ms)
{
   _equation_name = "Linear PDE";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   setMatrixType(SPARSE);
   setSolver(CG_SOLVER,DILU_PREC);
   _lump = true;
   _stab = false;
}


LinearPDE2D::LinearPDE2D(Mesh&         ms,
                         Vect<real_t>& u)
            : Equation<3,3,2,2>(ms,u)
{
   _equation_name = "Linear PDE";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   setMatrixType(SPARSE);
   setSolver(CG_SOLVER,DILU_PREC);
   _lump = true;
   _stab = false;
}


void LinearPDE2D::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   Triang3 tr(_theElement);
   _el_geo.area = tr.getArea();
   _el_geo.center = tr.getCenter();
   _dSh = tr.DSh();
   _el_geo.size = 2*tr.getCircumRadius();
   ElementNodeCoordinates();
   ElementNodeVector(*_u,_eu);
   eA0.clear();
   eA1.clear();
   eA2.clear();
   eRHS.clear();
}


void LinearPDE2D::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Line2 ln(sd);
   SideNodeCoordinates();
   _el_geo.center = ln.getCenter();
   _el_geo.length = ln.getLength();
   sMat = 0; sRHS = 0;
}


void LinearPDE2D::setInput(EqDataType    opt,
                           Vect<real_t>& u)
{
   Equa::setInput(opt,u);
}


void LinearPDE2D::Mat_00(real_t coef)
{
   real_t c = OFELI_SIXTH*_el_geo.area*coef;
   c *= getPDECoef(C00,SpaceTime(_el_geo.center.x,0.,0.,_TimeInt.time));
   eA0(1,1) = 2*c; eA0(2,2) = 2*c; eA0(3,3) = 2*c;
   eA0(1,2) =   c; eA0(1,3) =   c; eA0(2,3) =   c;
   eA0(2,1) =   c; eA0(3,1) =   c; eA0(3,2) =   c;
}


void LinearPDE2D::Mat_10(real_t coef)
{
   if (_lump) {
      real_t c = OFELI_THIRD*_el_geo.area*coef;
      c *= getPDECoef(C10,SpaceTime(_el_geo.center.x,0.,0.,_TimeInt.time));
      eA1(1,1) += c;
      eA1(2,2) += c;
      eA1(3,3) += c;
   }
   else {
      real_t c = OFELI_SIXTH*_el_geo.area*coef;
      c *= getPDECoef(C10,SpaceTime(_el_geo.center.x,0.,0.,_TimeInt.time));
      real_t d=0.5*c;
      eA1(1,1) += c; eA1(2,2) += c; eA1(3,3) += c;
      eA1(1,2) += d; eA1(2,1) += d; eA1(1,3) += d;
      eA1(3,1) += d; eA1(2,3) += d; eA1(3,2) += d;
   }
}


void LinearPDE2D::Mat_20(real_t coef)
{
   if (_lump) {
      real_t c = OFELI_THIRD*_el_geo.area*coef;
      c *= getPDECoef(C20,SpaceTime(_el_geo.center.x,0.,0.,_TimeInt.time));
      eA2(1,1) += c;
      eA2(2,2) += c;
      eA2(3,3) += c;
   }
   else {
      real_t c = OFELI_SIXTH*_el_geo.area*coef;
      c *= getPDECoef(C20,SpaceTime(_el_geo.center.x,0.,0.,_TimeInt.time));
      real_t d=0.5*c;
      eA2(1,1) += c; eA2(2,2) += c; eA2(3,3) += c;
      eA2(1,2) += d; eA2(2,1) += d; eA2(1,3) += d;
      eA2(3,1) += d; eA2(2,3) += d; eA2(3,2) += d;
   }
}


void LinearPDE2D::Mat_01(real_t coef)
{
   LocalVect<real_t,3> dd;
   if (_type01==2) {
      for (size_t i=1; i<=3; i++)
         dd(i) = (_d01,_dSh[i-1]);
      for (size_t i=1; i<=3; i++)
         for (size_t j=1; j<=3; j++)
            eA0(i,j) += OFELI_THIRD*coef*_el_geo.area*dd(j);
      if (_stab) {
         real_t c=coef*_el_geo.area*_el_geo.size/_d01.Norm(), d;
         for (size_t i=1; i<=3; i++) {
            d = c*dd(i);
            for (size_t j=1; j<=3; j++)
               eA0(i,j) += d*dd(j);
         }
      }
   }
}


void LinearPDE2D::Mat_02(real_t coef)
{
   real_t c = coef*_el_geo.area;
   if (_type02==1) {
      for (size_t i=1; i<=3; i++)
         for (size_t j=1; j<=3; j++)
            eA0(i,j) += c*_c02*(_dSh[i-1],_dSh[j-1]);
   }
   else if (_type02==2) {
      for (size_t i=1; i<=3; i++)
         for (size_t j=1; j<=3; j++)
            eA0(i,j) += c*(_d02.x*_dSh[i-1].x*_dSh[j-1].x + _d02.y*_dSh[i-1].y*_dSh[j-1].y);
   }
   else if (_type02==3) {
      c *= _f02(_el_geo.center.x,_el_geo.center.y,0.,_TimeInt.time);
      for (size_t i=1; i<=3; i++)
         for (size_t j=1; j<=3; j++)
            eA0(i,j) += c*(_dSh[i-1],_dSh[j-1]);
   }
   else if (_type02==4) {
      for (size_t i=1; i<=3; i++)
         for (size_t j=1; j<=3; j++)
            eA0(i,j) += c*(_e02(1,1)*_dSh[i-1].x*_dSh[j-1].x + _e02(1,2)*_dSh[i-1].y*_dSh[j-1].x +
                           _e02(2,1)*_dSh[i-1].x*_dSh[j-1].y + _e02(2,2)*_dSh[i-1].y*_dSh[j-1].y);
   }
}


void LinearPDE2D::BodyRHS(real_t bf)
{
   real_t c=OFELI_THIRD*_el_geo.area*bf;
   eRHS(1) += c;
   eRHS(2) += c;
   eRHS(3) += c;
}


void LinearPDE2D::BodyRHS(const Vect<real_t>& f)
{
   real_t c = OFELI_THIRD*_el_geo.area;
   for (size_t i=1; i<=3; i++)
      eRHS(i) += c*f((*_theElement)(i)->n());
}


void LinearPDE2D::BoundaryRHS(real_t flux)
{
   real_t c = 0.5*_el_geo.length*flux;
   for (size_t i=1; i<=2; i++)
       sRHS(i) += c;
}


void LinearPDE2D::BoundaryRHS(const Vect<real_t>& f)
{
   if (_theSide->getCode(1)>0) {
      real_t c = 0.5*_el_geo.length;
      if (f.getDOFType()==NODE_DOF) {
         sRHS(1) += c*f((*_theSide)(1)->n());
         sRHS(2) += c*f((*_theSide)(2)->n());
      }
      else if (f.getDOFType()==SIDE_DOF) {
         sRHS(1) += c*f(_theSide->n());
         sRHS(2) += c*f(_theSide->n());
      }
   }
}


Point<real_t> &LinearPDE2D::Flux() const
{
   _f = (_eu(1)*_dSh[0] + _eu(2)*_dSh[1] + _eu(3)*_dSh[2]);
   return _f;
}


real_t LinearPDE2D::Energy(Vect<real_t>& u)
{
   _u = &u;
   real_t E = 0.;
   MESH_EL {
      set(the_element);
      Point<real_t> du = _eu[0]*_dSh[0] + _eu[1]*_dSh[1] + _eu[2]*_dSh[2];
      E += 0.5*_el_geo.area*(du,du);
      if (_bf!=nullptr)
         E -= OFELI_THIRD*_el_geo.area*(_eu(1)*(*_bf)(The_element(1)->n()) +
                                        _eu(2)*(*_bf)(The_element(2)->n()) +
                                        _eu(3)*(*_bf)(The_element(3)->n())); 
   }
   return E;
}


void LinearPDE2D::EnergyGrad(Vect<real_t>& u,
                             Vect<real_t>& g)
{
   real_t f = 0.;
   g.clear();
   _u = &u;
   MESH_EL {
      set(the_element);
      Mat_02();
      for (size_t i=1; i<=3; ++i) {
         if (_bf!=nullptr)
            f = OFELI_THIRD*_el_geo.area*(*_bf)(The_element(i)->n());
         g(The_element(i)->n()) += eA0(i,1)*_eu(1) + eA0(i,2)*_eu(2) + eA0(i,3)*_eu(3) - f;
      }
   }
}


void LinearPDE2D::Grad(Vect<Point<real_t> >& g)
{
   MESH_EL {
      set(the_element);
      g(element_label) = _eu(1)*_dSh[0] + _eu(2)*_dSh[1] + _eu(3)*_dSh[2];
   }
}

} /* namespace OFELI */