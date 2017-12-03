/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2018 Rachid Touzani

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


#include "equations/therm/DC3DAT3.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"

namespace OFELI {


DC3DAT3::DC3DAT3()
{
}


DC3DAT3::DC3DAT3(const Element* el)
{
   set(el);
}


DC3DAT3::DC3DAT3(const Side* sd)
{
   set(sd);
}


DC3DAT3::DC3DAT3(const Element*      el,
                 const Vect<real_t>& u,
                       real_t        time)
{
   set(el);
   _time = time;
   ElementVector(u);
}


DC3DAT3::DC3DAT3(const Element*      el,
                 const Vect<real_t>& u,
                       real_t        time,
                       real_t        deltat,
                       int           scheme)
{
   set(el);
   _time = time;
   ElementVector(u);
   _time_step = deltat;
   setTimeIntegration(scheme);
}


DC3DAT3::DC3DAT3(const Side*         sd,
                 const Vect<real_t>& u,
                       real_t        time)
{
   set(sd);
   _time = time;
   SideVector(u);
}


DC3DAT3::DC3DAT3(const Side*         sd,
                 const Vect<real_t>& u,
                       real_t        time,
                       real_t        deltat,
                       int           scheme)
{
   set(sd);
   _time = time;
   SideVector(u);
   _time_step = deltat;
   setTimeIntegration(scheme);
}


DC3DAT3::~DC3DAT3() { }


void DC3DAT3::set(const Element* el)
{
   _nb_dof = 1;
   Init(el);
   setMaterial();
   Triang3 tr(_theElement);
   _area = tr.getArea();
   _center = tr.getCenter();
   _dSh(1) = tr.DSh(1);
   _dSh(2) = tr.DSh(2);
   _dSh(3) = tr.DSh(3);
   _h = 2*tr.getCircumRadius();
   ElementNodeCoordinates();
   eMat = 0;
   eRHS = 0;
}


void DC3DAT3::set(const Side* sd)
{
   _nb_dof = 1;
   Init(sd);
   Line2 ln(sd);
   SideNodeCoordinates();
   _center = ln.getCenter();
   _length = ln.getLength();
   _r[0] = _x[0].x; _r[1] = _x[1].x;
   sMat = 0;
   sRHS = 0;
}


void DC3DAT3::build()
{
   Equa_Therm<real_t,3,3,2,2>::build();
}


void DC3DAT3::LCapacityToLHS(real_t coef)
{
   real_t c = OFELI_THIRD*_area*_rhocp*coef;
   eMat(1,1) += c*_r[0];
   eMat(2,2) += c*_r[1];
   eMat(3,3) += c*_r[2];
}


void DC3DAT3::CapacityToLHS(real_t coef)
{
   real_t c = 0.5*OFELI_TWELVETH*_area*_rhocp*coef;
   eMat(1,1) += c*(2*_r[0] + _r[1] + _r[2]);
   eMat(2,2) += c*(_r[0] + 2*_r[1] + _r[2]);
   eMat(3,3) += c*(_r[0] + _r[1] + 2*_r[2]);
   eMat(1,2) += c*(_r[0]+_r[1]); eMat(2,1) += c*(_r[0]+_r[1]);
   eMat(2,3) += c*(_r[1]+_r[2]); eMat(3,2) += c*(_r[1]+_r[2]);
   eMat(1,3) += c*(_r[0]+_r[2]); eMat(3,1) += c*(_r[0]+_r[2]);
}


void DC3DAT3::LCapacityToRHS(real_t coef)
{
   real_t c = OFELI_THIRD*_area*_rhocp*coef;
   eRHS(1) += c*_r[0]*ePrev(1);
   eRHS(2) += c*_r[1]*ePrev(2);
   eRHS(3) += c*_r[2]*ePrev(3);
}


void DC3DAT3::CapacityToRHS(real_t coef)
{
   real_t m11, m22, m33, m12, m23, m13;
   real_t c = 0.5*OFELI_TWELVETH*_area*_rhocp*coef;
   m11 = c*(2*_r[0] + _r[1] + _r[2]);
   m22 = c*(_r[0] + 2*_r[1] + _r[2]);
   m33 = c*(_r[0] + _r[1] + 2*_r[2]);
   m12 = c*(_r[0]+_r[1]);
   m23 = c*(_r[1]+_r[2]);
   m13 = c*(_r[0]+_r[2]);
   eRHS(1) += m11*ePrev(1) + m12*ePrev(2) + m13*ePrev(3);
   eRHS(2) += m12*ePrev(1) + m22*ePrev(2) + m23*ePrev(3);
   eRHS(3) += m13*ePrev(1) + m23*ePrev(2) + m33*ePrev(3);
}


void DC3DAT3::Diffusion(real_t coef)
{
   real_t c = coef*OFELI_THIRD*_area*_diff*(_r[0]+_r[1]+_r[2]);
   for (size_t i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eMat(i,j) += c*(_dSh(i)*_dSh(j));
}


void DC3DAT3::DiffusionToRHS(real_t coef)
{
   real_t c = coef*OFELI_THIRD*_area*_diff*(_r[0]+_r[1]+_r[2]);
   for (size_t i=1; i<=3; i++) {
      Point<real_t> z = c*_dSh(i);
      eRHS(i) += z*(ePrev(1)*_dSh(1) + ePrev(2)*_dSh(2) + ePrev(3)*_dSh(3));
   }
}


void DC3DAT3::Diffusion(const LocalMatrix<real_t,2,2>& diff,
                              real_t                   coef)
{
   real_t c = coef*OFELI_THIRD*_area*(_r[0]+_r[1]+_r[2]);
   for (size_t i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eMat(i,j) += c*(diff(1,1)*_dSh(i).x*_dSh(j).x + diff(2,2)*_dSh(i).y*_dSh(j).y
                       + diff(1,2)*_dSh(i).y*_dSh(j).x + diff(2,1)*_dSh(i).x*_dSh(j).y);
}


void DC3DAT3::BodyRHS(UserData<real_t>& ud)
{
   real_t c = OFELI_SIXTH*_area;
   Point<real_t> x = 0.5*(_x[0]+_x[1]);
   Point<real_t> y = 0.5*(_x[1]+_x[2]);
   Point<real_t> z = 0.5*(_x[0]+_x[2]);
   eRHS(1) += c*(ud.BodyForce(z,_time)*(_r[0]+_r[2]) + ud.BodyForce(x,_time)*(_r[0]+_r[1]));
   eRHS(2) += c*(ud.BodyForce(x,_time)*(_r[0]+_r[1]) + ud.BodyForce(y,_time)*(_r[2]+_r[1]));
   eRHS(3) += c*(ud.BodyForce(z,_time)*(_r[0]+_r[2]) + ud.BodyForce(y,_time)*(_r[1]+_r[2]));
}


void DC3DAT3::BodyRHS(const Vect<real_t>& bf,
                            int           opt)
{
   real_t c = 0.125*OFELI_THIRD*_area;
   if (opt==LOCAL_ARRAY) {
      eRHS(1) += c*((bf[0]+bf[2])*(_r[0]+_r[2]) + (bf[0]+bf[1])*(_r[0]+_r[1]));
      eRHS(2) += c*((bf[0]+bf[1])*(_r[0]+_r[1]) + (bf[1]+bf[2])*(_r[1]+_r[2]));
      eRHS(3) += c*((bf[0]+bf[2])*(_r[0]+_r[2]) + (bf[1]+bf[2])*(_r[1]+_r[2]));
   }
   else {
      real_t f1 = bf(_theElement->getNodeLabel(1)) + bf(_theElement->getNodeLabel(3));
      real_t f2 = bf(_theElement->getNodeLabel(1)) + bf(_theElement->getNodeLabel(2));
      real_t f3 = bf(_theElement->getNodeLabel(2)) + bf(_theElement->getNodeLabel(3));
      eRHS(1) += c*(f1*(_r[0]+_r[2]) + f2*(_r[0]+_r[1]));
      eRHS(2) += c*(f2*(_r[0]+_r[1]) + f3*(_r[1]+_r[2]));
      eRHS(3) += c*(f1*(_r[0]+_r[2]) + f3*(_r[1]+_r[2]));
   }
}


void DC3DAT3::BoundaryRHS(real_t flux)
{
   sRHS(1) += 0.5*flux*_length*_r[0];
   sRHS(2) += 0.5*flux*_length*_r[1];
}


void DC3DAT3::BoundaryRHS(const Vect<real_t>& sf,
                                int           opt)
{
  if (opt==LOCAL_ARRAY)
      for (size_t i=1; i<=2; i++)
         sRHS(i) += sf(i)*0.5*_length*_r[0];
   else
      for (size_t i=1; i<=2; i++)
         sRHS(i) += sf((*_theSide)(i)->n())*0.5*_length*_r[1];
}


Point<real_t> &DC3DAT3::Grad(const Vect<real_t>& u)
{
   _grad = u[0]*_dSh(1) + u[1]*_dSh(2) + u[2]*_dSh(3);
   return _grad;
}

} /* namespace OFELI */
