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

                       Class DC2DT6: Diffusion-Convection Element
                using 6-Node Triangular Finite element in two dimensions

  ==============================================================================*/


#include "equations/therm/DC2DT6.h"

namespace OFELI {

DC2DT6::DC2DT6()
{
   _tr = NULL;
   _ln = NULL;
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 6-Node Triangles (P2)";
}


DC2DT6::DC2DT6(Mesh&              ms,
               SkSMatrix<real_t>& a,
               Vect<real_t>&      b)
{
   _tr = NULL;
   _ln = NULL;
   _theMesh = &ms;
   _A = &a;
   if (b.size()==0) { }
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 6-Node Triangles (P2)";
}


DC2DT6::DC2DT6(const Element* el)
{
   set(el);
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 6-Node Triangles (P2)";
   eMat = 0;
   eRHS = 0;
}


DC2DT6::DC2DT6(const Side* sd)
{
   set(sd);
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 6-Node Triangles (P2)";
   sMat = 0;
   sRHS = 0;
}


DC2DT6::DC2DT6(const Element*      el,
               const Vect<real_t>& u,
                     real_t        time)
{
   _time = time;
   set(el);
   ElementVector(u);
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 6-Node Triangles (P2)";
}


DC2DT6::DC2DT6(const Element*      el,
               const Vect<real_t>& u,
                     real_t        time,
                     real_t        deltat,
                     int           scheme)
{
   _time = time;
   set(el);
   ElementVector(u);
   _time_step = deltat;
   setTimeIntegration(scheme);
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 6-Node Triangles (P2)";
}


DC2DT6::DC2DT6(const Side*         sd,
               const Vect<real_t>& u,
                     real_t        time)
{
   _time = time;
   set(sd);
   SideVector(u);
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 6-Node Triangles (P2)";
}


DC2DT6::DC2DT6(const Side*         sd,
               const Vect<real_t>& u,
                     real_t        time,
                     real_t        deltat,
                     int           scheme)
{
   set(sd);
   SideVector(u);
   _time = time;
   _label = sd->n();
   _time_step = deltat;
   setTimeIntegration(scheme);
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 6-Node Triangles (P2)";
}


DC2DT6::~DC2DT6()
{
   if (_tr)
      delete _tr;
   if (_ln)
      delete _ln;
}


void DC2DT6::set(const Element* el)
{
   _nb_dof = 1;
   Init(el);
   setMaterial();
   _ln = NULL;
   _tr = new Triang6S(el);
   _area = _tr->getArea();
   ElementNodeCoordinates();
   eMat = 0;
   eRHS = 0;
   _s[0] = Point<real_t>(0.5,0);
   _s[1] = Point<real_t>(0.5,0.5);
   _s[2] = Point<real_t>(0,0.5);
   _a3 = OFELI_THIRD*_area;
}


void DC2DT6::set(const Side* sd)
{
   _label = sd->n();
   _nb_dof = 1;
   Init(sd);
   _tr = NULL;
   _ln = new Line3(sd);
   _center = _ln->getCenter();
   SideNodeCoordinates();
   eMat = 0;
   eRHS = 0;
}


void DC2DT6::build()
{
   Equa_Therm<real_t,6,6,3,3>::build();
}

/*
void DC2DT6::LCapacityToLHS(real_t coef)
{
   real_t c = 0.5*a3_*rhocp_*coef;
   eMat(1,1) += c; eMat(2,2) += c;
   eMat(3,3) += c; eMat(4,4) += c;
   eMat(5,5) += c; eMat(6,6) += c;
}


void DC2DT6::LCapacityToRHS(real_t coef)
{
   real_t c = 0/5*a3_*rhocp_*coef;
   eRHS(1) += c*ePrev(1);
   eRHS(2) += c*ePrev(2);
   eRHS(3) += c*ePrev(3);
   eRHS(4) += c*ePrev(4);
   eRHS(5) += c*ePrev(5);
   eRHS(6) += c*ePrev(6);
}


void DC2DT6::CapacityToLHS(real_t coef)
{
   real_t c = OFELI_SIXTH*_tr->Area()*rhocp_*coef;
   real_t d = 0.5*c;
   //eMat(1,1) += c; eMat(2,2) += c; eMat(3,3) += c;
    eMat(1,1) += c; eMat(2,2) += c; eMat(3,3) += c;eMat(4,4)+=c;eMat(5,5)+=c;eMat(6,6)+=c;

   //eMat(1,2) += d; eMat(2,1) += d; eMat(1,3) += d;
   //eMat(3,1) += d; eMat(2,3) += d; eMat(3,2) += d;

    eMat(1,2) += d; eMat(1,3) += d;eMat(1,4) += d; eMat(1,5) += d; eMat(1,6) += d;
    eMat(2,1) += d; eMat(2,3) += d;eMat(2,4) += d; eMat(2,5) += d; eMat(2,6) += d;
    eMat(3,1) += d; eMat(3,2) += d;eMat(3,4) += d; eMat(3,5) += d; eMat(3,6) += d;
    eMat(4,1) += d; eMat(4,2) += d;eMat(4,3) += d; eMat(4,5) += d; eMat(4,6) += d;
    eMat(5,1) += d; eMat(5,2) += d;eMat(5,3) += d; eMat(5,4) += d; eMat(5,6) += d;
    eMat(6,1) += d; eMat(6,2) += d;eMat(6,3) += d; eMat(6,4) += d; eMat(6,5) += d;
}


void DC2DT6::CapacityToRHS(real_t coef)
{
   real_t c = OFELI_SIXTH*_tr->Area()*rhocp_*coef;
   real_t d = 0.5*c;
   eRHS(1) += c*ePrev(1) + d*(ePrev(2) + ePrev(3));
   eRHS(2) += c*ePrev(2) + d*(ePrev(1) + ePrev(3));
   eRHS(3) += c*ePrev(3) + d*(ePrev(1) + ePrev(2));
}*/


void DC2DT6::Diffusion(real_t coef)
{
   real_t a = _a3*coef*_diff;
   for (size_t i=1; i<=6; i++)
      for (size_t j=1; j<=6; j++) {
         real_t c = 0;
         for (size_t k=0; k<3; k++)
            c +=  _tr->DSh(i,_s[k])*_tr->DSh(j,_s[k]);
         eMat(i,j) += a*c;
      }
}


void DC2DT6::Convection(Point<real_t>& v,
                        real_t         coef)
{
   for (size_t j=1; j<=6; j++) {
      eMat(4,j) += _a3*coef*(v*_tr->DSh(j,_s[0]));
      eMat(5,j) += _a3*coef*(v*_tr->DSh(j,_s[1]));
      eMat(6,j) += _a3*coef*(v*_tr->DSh(j,_s[2]));
   }
}


void DC2DT6::Convection(const Vect<real_t>& v,
                              real_t        coef)
{
   LocalMatrix<real_t,6,2> ve;
   for (size_t i=1; i<=6; i++) {
      size_t n = _theElement->getNodeLabel(i);
      ve(i,1) = v(2*n-1);
      ve(i,2) = v(2*n  );
   }
   for (size_t j=1; j<=6; j++) {
      eMat(4,j) += _a3*coef*(ve(4,1)*(_tr->DSh(j,_s[0])).x+ve(4,2)*(_tr->DSh(j,_s[0])).y);
      eMat(5,j) += _a3*coef*(ve(5,1)*(_tr->DSh(j,_s[1])).x+ve(5,2)*(_tr->DSh(j,_s[1])).y);
      eMat(6,j) += _a3*coef*(ve(6,1)*(_tr->DSh(j,_s[2])).x+ve(6,2)*(_tr->DSh(j,_s[2])).y);
   }
}


void DC2DT6::BodyRHS(UserData<real_t>& ud)
{
   eRHS(4) += _a3*ud.BodyForce(_x[3],_time);
   eRHS(5) += _a3*ud.BodyForce(_x[4],_time);
   eRHS(6) += _a3*ud.BodyForce(_x[5],_time);
}


void DC2DT6::BodyRHS(const Vect<real_t>& bf,
                           int           opt)
{
   if (opt==LOCAL_ARRAY) {
      eRHS(4) += _a3*bf(4);
      eRHS(5) += _a3*bf(5);
      eRHS(6) += _a3*bf(6);
   }
   else {
      eRHS(4) += _a3*bf(_theElement->getNodeLabel(4));
      eRHS(5) += _a3*bf(_theElement->getNodeLabel(5));
      eRHS(6) += _a3*bf(_theElement->getNodeLabel(6));
   }
}


void DC2DT6::BoundaryRHS(UserData<real_t>& ud)
{
   _ln->setLocal(-1.0);
   sRHS(1) += OFELI_THIRD*_ln->getDet()*ud.SurfaceForce((*_theSide)(1)->getCoord(),_theSide->getCode(1),_time);
   _ln->setLocal(0.0);
   sRHS(2) += OFELI_THIRD*_ln->getDet()*ud.SurfaceForce((*_theSide)(2)->getCoord(),_theSide->getCode(1),_time);
   _ln->setLocal(1.0);
   sRHS(3) += 4*OFELI_THIRD*_ln->getDet()*ud.SurfaceForce((*_theSide)(3)->getCoord(),_theSide->getCode(1),_time);
}


void DC2DT6::BoundaryRHS(const Vect<real_t>& sf,
                               int           opt)
{
   real_t c = OFELI_THIRD*_ln->getDet();
   if (opt==LOCAL_ARRAY) {
      _ln->setLocal(-1.0);
      sRHS(1) += c*sf(1);
      _ln->setLocal(0.0);
      sRHS(2) += c*sf(2);
      _ln->setLocal(1.0);
      sRHS(3) += 4*c*sf(3);
   }
   else {
      _ln->setLocal(-1.0);
      sRHS(1) += c*sf((*_theSide)(1)->n());
      _ln->setLocal(0.0);
      sRHS(2) += c*sf((*_theSide)(2)->n());
      _ln->setLocal(1.0);
      sRHS(3) += 4*c*sf((*_theSide)(3)->n());
   }
}

} /* namespace OFELI */
