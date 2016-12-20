/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2017 Rachid Touzani

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

                       Class DC3DT4 : Diffusion-Convection Element
             using 4-Node tetrahedral Finite element in three dimensions

  ==============================================================================*/


#include "equations/therm/DC3DT4.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Tetra4.h"

namespace OFELI {

DC3DT4::DC3DT4(const Element* el)
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "3-D, 4-Node Tetrahedra (P1)";
   set(el);
}


DC3DT4::DC3DT4(const Side* sd)
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "3-D, 4-Node Tetrahedra (P1)";
   set(sd);
}


DC3DT4::DC3DT4(const Element*      el,
               const Vect<real_t>& u,
                     real_t        time)
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "3-D, 4-Node Tetrahedra (P1)";
   set(el);
   _time = time;
   ElementVector(u);
}


DC3DT4::DC3DT4(const Element*      el,
               const Vect<real_t>& u,
                     real_t        time,
                     real_t        deltat,
                     int           scheme)
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "3-D, 4-Node Tetrahedra (P1)";
   set(el);
   _time = time;
   ElementVector(u);
   _time_step = deltat;
   setTimeIntegration(scheme);
}


DC3DT4::DC3DT4(const Side*         sd,
               const Vect<real_t>& u,
                     real_t        time)
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "3-D, 4-Node Tetrahedra (P1)";
   set(sd);
   _time = time;
   SideVector(u);
}


DC3DT4::DC3DT4(const Side*         sd,
               const Vect<real_t>& u,
                     real_t        time,
                     real_t        deltat,
                     int           scheme)
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "3-D, 4-Node Tetrahedra (P1)";
   set(sd);
   _time = time;
   SideVector(u);
   _time_step = deltat;
   setTimeIntegration(scheme);
}


DC3DT4::~DC3DT4() { }


void DC3DT4::set(const Element* el)
{
   _nb_dof = 1;
   Init(el);
   setMaterial();
   Tetra4 tetra(_theElement);
   _center = tetra.getCenter();
   _volume = tetra.getVolume();
   _det = tetra.getDet();
   _dSh(1) = tetra.DSh(1);
   _dSh(2) = tetra.DSh(2);
   _dSh(3) = tetra.DSh(3);
   _dSh(4) = tetra.DSh(4);
   ElementNodeCoordinates();
   eMat = 0;
   eRHS = 0;
}


void DC3DT4::set(const Side* sd)
{
   _nb_dof = 1;
   Init(sd);
   Triang3 triang(sd);
   _area = triang.getArea();
   _center = triang.getCenter();
   SideNodeCoordinates();
   _center = triang.getCenter();
   sMat = 0;
   sRHS = 0;
}


void DC3DT4::build()
{
   Equa_Therm<real_t,4,4,3,3>::build();
}


void DC3DT4::LCapacityToLHS(real_t coef)
{
   real_t c = coef*0.25*_volume*_rhocp;
   for (size_t i=1; i<=4; i++)
      eMat(i,i) += c;
}


void DC3DT4::LCapacityToRHS(real_t coef)
{
   real_t c = coef*0.25*_volume*_rhocp;
   eRHS(1) += c*ePrev(1);
   eRHS(2) += c*ePrev(2);
   eRHS(3) += c*ePrev(3);
   eRHS(4) += c*ePrev(4);
}


void DC3DT4::CapacityToLHS(real_t coef)
{
   real_t c = 0.1*_volume*_rhocp*coef;
   real_t d = 0.5*c;
   eMat(1,1) += c; eMat(2,2) += c; eMat(3,3) += c; eMat(4,4) += c;
   eMat(1,2) += d; eMat(2,1) += d; eMat(1,3) += d; eMat(1,4) += d;
   eMat(3,1) += d; eMat(2,3) += d; eMat(3,2) += d; eMat(2,4) += d;
   eMat(4,1) += d; eMat(4,2) += d; eMat(4,3) += d; eMat(3,4) += d;
}


void DC3DT4::CapacityToRHS(real_t coef)
{
   real_t c = 0.1*_volume*_rhocp*coef;
   real_t d = 0.5*c;
   eRHS(1) += c*ePrev(1) + d*(ePrev(2) + ePrev(3)+ ePrev(4));
   eRHS(2) += c*ePrev(2) + d*(ePrev(1) + ePrev(3)+ ePrev(4));
   eRHS(3) += c*ePrev(3) + d*(ePrev(1) + ePrev(2)+ ePrev(4));
   eRHS(4) += c*ePrev(4) + d*(ePrev(1) + ePrev(2)+ ePrev(3));
}


void DC3DT4::Diffusion(real_t coef)
{
   real_t c = coef*_diff*_volume;
   for (size_t i=1; i<=4; i++)
      for (size_t j=1; j<=4; j++)
         eMat(i,j) += c*(_dSh(i)*_dSh(j));
}


void DC3DT4::Diffusion(const DMatrix<real_t>& diff,
                             real_t           coef)
{
   real_t c = coef*_volume;
   for (size_t i=1; i<=4; i++)
       for (size_t j=1; j<=4; j++)
          eMat(i,j) += c*(diff(1,1)*_dSh(i).x*_dSh(j).x + diff(1,2)*_dSh(i).y*_dSh(j).x
                        + diff(1,3)*_dSh(i).z*_dSh(j).x + diff(2,1)*_dSh(i).x*_dSh(j).y
                        + diff(2,2)*_dSh(i).y*_dSh(j).y + diff(2,3)*_dSh(i).z*_dSh(j).y
                        + diff(3,1)*_dSh(i).x*_dSh(j).z + diff(3,2)*_dSh(i).y*_dSh(j).z
                        + diff(3,3)*_dSh(i).z*_dSh(j).z);
}


void DC3DT4::DiffusionToRHS(real_t coef)
{
   Point<real_t> u = _dSh(1)*ePrev(1) + _dSh(2)*ePrev(2) + _dSh(3)*ePrev(3) + _dSh(4)*ePrev(4);
   eRHS(1) -= coef*_diff*OFELI_SIXTH*_det*(_dSh(1)*u);
   eRHS(2) -= coef*_diff*OFELI_SIXTH*_det*(_dSh(2)*u);
   eRHS(3) -= coef*_diff*OFELI_SIXTH*_det*(_dSh(3)*u);
   eRHS(4) -= coef*_diff*OFELI_SIXTH*_det*(_dSh(4)*u);
}


void DC3DT4::Convection(const Point<real_t>& v,
                              real_t         coef)
{
   for (size_t i=1; i<=4; i++)
      for (size_t j=1; j<=4; j++)
         eMat(i,j) += coef*0.25*_volume*(v*_dSh(i));
}


void DC3DT4::Convection(const Vect<Point<real_t> >& v,
                              real_t                coef)
{
   size_t i;
   LocalMatrix<real_t,4,3> ve;
   for (i=1; i<=4; i++) {
      size_t n = _theElement->getNodeLabel(i);
      ve(i,1) = v(n).x;
      ve(i,2) = v(n).y;
      ve(i,3) = v(n).z;
   }
   Point<real_t> w;
   w.x = 0.25*(ve(1,1)+ve(2,1)+ve(3,1)+ve(4,1));
   w.y = 0.25*(ve(1,2)+ve(2,2)+ve(3,2)+ve(4,2));
   w.z = 0.25*(ve(1,3)+ve(2,3)+ve(3,3)+ve(4,3));
   for (i=1; i<=4; i++)
      for (size_t j=1; j<=4; j++)
         eMat(i,j) += coef*0.25*_volume*(w*_dSh(i));
}


void DC3DT4::Convection(real_t coef)
{
   LocalMatrix<real_t,4,3> ve;
   size_t i;
   for (i=1; i<=4; i++) {
      size_t n = _theElement->getNodeLabel(i);
      ve(i,1) = (*_vel)(3*n-2);
      ve(i,2) = (*_vel)(3*n-1);
      ve(i,3) = (*_vel)(3*n  );
   }
   Point<real_t> v;
   v.x = 0.25*(ve(1,1) + ve(2,1) + ve(3,1) + ve(4,1));
   v.y = 0.25*(ve(1,2) + ve(2,2) + ve(3,2) + ve(4,2));
   v.z = 0.25*(ve(1,3) + ve(2,3) + ve(3,3) + ve(4,3));
   real_t c = 0.25*_volume*coef;
   for (i=1; i<=4; i++)
      for (size_t j=1; j<=4; j++)
         eMat(i,j) += c*(v*_dSh(j));
}


void DC3DT4::RHS_Convection(const Point<real_t>& v,
                                  real_t         coef)
{
   Point<real_t> u = _dSh(1)*ePrev(1) + _dSh(2)*ePrev(2) + _dSh(3)*ePrev(3) + _dSh(4)*ePrev(4);
   eRHS(1) += coef*0.25*_volume*(v*u);
   eRHS(2) += coef*0.25*_volume*(v*u);
   eRHS(3) += coef*0.25*_volume*(v*u);
   eRHS(4) += coef*0.25*_volume*(v*u);
}


void DC3DT4::BodyRHS(UserData<real_t>& ud,
                     real_t            coef)
{
   real_t c = coef*0.25*_volume*ud.BodyForce(_center,_time);
   eRHS(1) += c;
   eRHS(2) += c;
   eRHS(3) += c;
   eRHS(4) += c;
}


void DC3DT4::BodyRHS(const Vect<real_t>& bf,
                           int           opt)
{
   if (opt==LOCAL_ARRAY) {
      eRHS(1) += bf(1)*0.25*_volume;
      eRHS(2) += bf(2)*0.25*_volume;
      eRHS(3) += bf(3)*0.25*_volume;
      eRHS(4) += bf(4)*0.25*_volume;
   }
   else {
      eRHS(1) += bf(_theElement->getNodeLabel(1))*0.25*_volume;
      eRHS(2) += bf(_theElement->getNodeLabel(2))*0.25*_volume;
      eRHS(3) += bf(_theElement->getNodeLabel(3))*0.25*_volume;
      eRHS(4) += bf(_theElement->getNodeLabel(4))*0.25*_volume;
   }
}


void DC3DT4::BoundaryRHS(UserData<real_t>& ud,
                         real_t            coef)
{
   if (_theSide->getCode(1)>0) {
      real_t c = OFELI_THIRD*coef*_area*ud.SurfaceForce(_center,_theSide->getCode(1),_time);
      sRHS(1) += c;
      sRHS(2) += c;
      sRHS(3) += c;
   }
}


void DC3DT4::BoundaryRHS(real_t flux)
{
   real_t c = flux*_area*OFELI_THIRD;
   sRHS(1) += c;
   sRHS(2) += c;
   sRHS(3) += c;
}


void DC3DT4::BoundaryRHS(const Vect<real_t>& sf,
                               int           opt)
{
   if (opt==LOCAL_ARRAY) {
      sRHS(1) += sf(1)*_area*OFELI_THIRD;
      sRHS(2) += sf(2)*_area*OFELI_THIRD;
      sRHS(3) += sf(3)*_area*OFELI_THIRD;
   }
   else {
      sRHS(1) += sf((*_theSide)(1)->n()) * _area*OFELI_THIRD;
      sRHS(2) += sf((*_theSide)(2)->n()) * _area*OFELI_THIRD;
      sRHS(3) += sf((*_theSide)(3)->n()) * _area*OFELI_THIRD;
   }
}


Point<real_t> DC3DT4::Flux() const
{
   Point<real_t> f;
   f = _diff*(ePrev(1)*_dSh(1) + ePrev(2)*_dSh(2) +
              ePrev(3)*_dSh(3) + ePrev(4)*_dSh(4));
   return f;
}


Point<real_t> DC3DT4::Grad(const Vect<real_t>& u) const
{
   Point<real_t> grad;
   grad = u[0]*_dSh(1) + u[1]*_dSh(2) + u[2]*_dSh(3) + u[3]*_dSh(4);
   return grad;
}


void DC3DT4::Periodic(real_t coef)
{
   for (size_t i=1; i<=3; i++) {
      real_t c = OFELI_THIRD*_area*coef;
      if (_theSide->getCode(1) == PERIODIC_A)
         eMat(i,i) += c;
      else if (_theSide->getCode(1) == PERIODIC_B)
         eMat(i,i) -= c;
   }
}

} /* namespace OFELI */
