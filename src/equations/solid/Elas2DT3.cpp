/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2019 Rachid Touzani

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

                         Implementation of class Elas2DT3
             for 2-D Linear Elasticity Equations with plane deformations
                      using 3-node triangular finite element

  ==============================================================================*/


#include "equations/solid/Elas2DT3.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

Elas2DT3::Elas2DT3()
{
}


Elas2DT3::Elas2DT3(Mesh& ms) 
         : Equation<real_t,3,6,2,4>(ms)
{
   _equation_name = "Linearized Elasticity";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _terms = DEVIATORIC|DILATATION|BODY_FORCE;
}


Elas2DT3::Elas2DT3(Mesh&         ms,
                   Vect<real_t>& u) 
         : Equation<real_t,3,6,2,4>(ms,u)
{
   _equation_name = "Linearized Elasticity";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _terms = DEVIATORIC|DILATATION|BODY_FORCE;
}


Elas2DT3::~Elas2DT3() { }


void Elas2DT3::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   setMaterial();
   PlaneStrain();
   Triang3 tr(el);
   _el_geo.area = tr.getArea();
   _dSh = tr.DSh();
   ElementNodeCoordinates();
   ElementNodeVector(*_u,_eu);
   eA0 = 0, eA1 = 0, eA2 = 0;
   eRHS = 0;
}


void Elas2DT3::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Line2 ln(sd);
   _el_geo.length = ln.getLength();
   SideNodeCoordinates();
   SideNodeVector(*_u,_su);
   sA0 = 0;
   sRHS = 0;
}


void Elas2DT3::Media(real_t E,
                     real_t nu,
                     real_t rho)
{
   _E = E;
   _nu = nu;
   _rho = rho;
   PlaneStrain(E,nu);
}


void Elas2DT3::PlaneStrain(real_t E,
                           real_t nu)
{
   _E = E;
   _nu = nu;
   _lambda = _nu*_E/((1+_nu)*(1-2*_nu));
   _G = 0.5*_E/(1+_nu);
   _E1 = _lambda*(1-_nu)/_nu;
   _E2 = _lambda;
   _E3 = _E1;
   _E6 = _G;
}


void Elas2DT3::PlaneStrain()
{
   _lambda = _nu*_E/((1+_nu)*(1-2*_nu));
   _G = 0.5*_E/(1+_nu);
   _E1 = _lambda*(1-_nu)/_nu;
   _E2 = _lambda;
   _E3 = _E1;
   _E6 = _G;
}


void Elas2DT3::PlaneStress(real_t E,
                           real_t nu)
{
   PlaneStrain(E,nu);
   _lambda = 2*_lambda*_G/(_lambda+2*_G);
   _E1 = _E/(1-_nu*_nu);
   _E2 = _E1*_nu;
   _E3 = _E1;
   _E6 = _G;
}


void Elas2DT3::PlaneStress()
{
   PlaneStrain();
   _lambda = 2*_lambda*_G/(_lambda+2*_G);
   _E1 = _E/(1-_nu*_nu);
   _E2 = _E1*_nu;
   _E3 = _E1;
   _E6 = _G;
}


void Elas2DT3::LMass(real_t coef)
{
   real_t c=_rho*coef*_el_geo.area*OFELI_THIRD;
   for (size_t i=1; i<=3; i++) {
      eA2(2*i-1,2*i-1) += c;
      eA2(2*i  ,2*i  ) += c;
   }
}


void Elas2DT3::Mass(real_t coef)
{
   real_t c=0.5*OFELI_SIXTH*_el_geo.area*_rho*coef;
   for (size_t i=1; i<=3; i++) {
      for (size_t j=1; j<=3; j++) {
         eA2(2*i-1,2*j-1) += c;
         eA2(2*i  ,2*j  ) += c;
      }
      eA2(2*i-1,2*i-1) += c;
      eA2(2*i  ,2*i  ) += c;
   }
}


void Elas2DT3::Deviator(real_t coef)
{
   real_t c=_G*_el_geo.area*coef;
   for (size_t i=1; i<=3; i++) {
      Point<real_t> a=c*_dSh[i-1];
      for (size_t j=1; j<=3; j++) {
         eA0(2*i-1,2*j-1) += 2*a.x*_dSh[j-1].x + a.y*_dSh[j-1].y;
         eA0(2*i-1,2*j  ) += a.y*_dSh[j-1].x;
         eA0(2*i  ,2*j-1) += a.x*_dSh[j-1].y;
         eA0(2*i  ,2*j  ) += 2*a.y*_dSh[j-1].y + a.x*_dSh[j-1].x;
      }
   }
}


void Elas2DT3::Dilatation(real_t coef)
{
   real_t c=_lambda*_el_geo.area*coef;
   for (size_t i=1; i<=3; i++) {
      Point<real_t> a=c*_dSh[i-1];
      for (size_t j=1; j<=3; j++) {
         eA0(2*i-1,2*j-1) += a.x*_dSh[j-1].x;
         eA0(2*i-1,2*j  ) += a.x*_dSh[j-1].y;
         eA0(2*i  ,2*j-1) += a.y*_dSh[j-1].x;
         eA0(2*i  ,2*j  ) += a.y*_dSh[j-1].y;
      }
   }
}


void Elas2DT3::BodyRHS(const Vect<real_t>& f)
{
   for (size_t k=1; k<=3; k++) {
      eRHS(2*k-1) += OFELI_THIRD*_el_geo.area*f(k,1);
      eRHS(2*k  ) += OFELI_THIRD*_el_geo.area*f(k,2);
   }
}


void Elas2DT3::BoundaryRHS(const Vect<real_t>& f)
{
   real_t c = 0.5*_el_geo.length;
   real_t fx = c*f(_theSide->n(),1);
   real_t fy = c*f(_theSide->n(),2);
   if (_theSide->getCode(1) != CONTACT) {
      sRHS(1) += fx;
      sRHS(3) += fx;
   }
   if (_theSide->getCode(2) != CONTACT) {
      sRHS(2) += fy;
      sRHS(4) += fy;
   }
}


void Elas2DT3::Periodic(real_t coef)
{
   for (size_t i=1; i<=2; i++) {
      real_t c = 0.5*_el_geo.length*coef;
      if (_theSide->getCode(1) == PERIODIC_A)
         sA0(2*i-1,2*i-1) += c;
      else if (_theSide->getCode(1) == PERIODIC_B)
         sA0(2*i-1,2*i-1) -= c;
      if (_theSide->getCode(2) == PERIODIC_A)
         sA0(2*i  ,2*i  ) += c;
      else if (_theSide->getCode(2) == PERIODIC_B)
         sA0(2*i  ,2*i  ) -= c;
   }
}


int Elas2DT3::Contact(real_t coef)
{
   int ret = 0;
   real_t c=0.5*coef*_el_geo.length;
   if (_cd==nullptr)
      throw OFELIException("In Elas2DT3::Contact(coef): No contact distance provided.");
   real_t dx=(*_cd)(_theSide->n(),1), dy=(*_cd)(_theSide->n(),2);
   if (_theSide->getCode(1) == CONTACT_BC) {
      if (dx > (*_u)((*_theSide)(1)->n(),1)) {
         ret = 1;
         sA0(1,1) += c;
         sRHS(1) += c*dx;
      }
      if (dx > (*_u)((*_theSide)(2)->n(),1)) {
         ret = 1;
         sA0(3,3) += c;
         sRHS(3) += c*dx;
      }
   }

   if (_theSide->getCode(2) == CONTACT_BC) {
      if (dy > (*_u)((*_theSide)(1)->n(),2)) {
         ret = 1;
         sA0(2,2) += c;
         sRHS(2) += c*dy;
      }
      if (dy > (*_u)((*_theSide)(2)->n(),2)) {
         ret = 1;
         sA0(4,4) += c;
         sRHS(4) += c*dy;
     }
   }
   return ret;
}


void Elas2DT3::Reaction(Vect<real_t>& r)
{
   r.setSize(_nb_sides,2);
   size_t n1, n2;
   mesh_elements(*_theMesh) {
      set(the_element);
      Deviator();
      Dilatation();
      for (size_t s=1; s<=3; s++) {
         Side *sd = _theElement->getPtrSide(s);
         size_t t = (s+1)%3;
         if (sd) {
            n1 = (*sd)(1)->n(), n2 = (*sd)(2)->n();
            real_t u1 = (*_u)((*_theSide)(1)->n(),1), u2 = (*_u)((*_theSide)(2)->n(),1);
            real_t v1 = (*_u)((*_theSide)(1)->n(),2), v2 = (*_u)((*_theSide)(2)->n(),2);
            if (sd->getCode(1)==CONTACT) {
               r((*sd)(1)->n(),1) += eA0(2*s-1,2*s-1)*u1 + eA0(2*s-1,2*t-1)*u2;
               r((*sd)(2)->n(),1) += eA0(2*t-1,2*s-1)*u1 + eA0(2*t-1,2*t-1)*u2;
            }
            if (sd->getCode(2)==CONTACT) {
               r((*sd)(1)->n(),2) += eA0(2*s  ,2*s  )*v1 + eA0(2*s  ,2*t  )*v2;
               r((*sd)(2)->n(),2) += eA0(2*t  ,2*s  )*v1 + eA0(2*t  ,2*t  )*v2;
            }
         }
      }
   }
}


void Elas2DT3::ContactPressure(const Vect<real_t>& f,
                               real_t              penal,
                               Point<real_t>&      p)
{
   p = 0;
   real_t ff = f(_theSide->n(),1);
   if (_theSide->getCode(1)==CONTACT) {
      if (ff>=_eu(1) && ff>=_eu(3))
         p.x = penal*(f(_theSide->n(),1)-0.5*(_eu(1)+_eu(3)));
   }
   ff = f(_theSide->n(),2);
   if (_theSide->getCode(2)==CONTACT) {
     if (ff>=_su(2) && ff>=_su(4))
         p.y = penal*(f(_theSide->n(),2)-0.5*(_su(2)+_su(4)));
   }
}


void Elas2DT3::Strain(Vect<real_t>& eps)
{
   eps.setSize(_nb_el,3);
   mesh_elements(*_theMesh) {
     size_t ne = The_element.n();
      size_t n1=The_element(1)->n(), n2=The_element(2)->n(), n3=The_element(3)->n();
      Triang3 tr(the_element);
      _dSh = tr.DSh();
      eps(ne,1) = (*_u)(n1,1)*_dSh[0].x + (*_u)(n2,1)*_dSh[1].x + (*_u)(n3,1)*_dSh[2].x;
      eps(ne,2) = (*_u)(n1,2)*_dSh[0].y + (*_u)(n2,2)*_dSh[1].y + (*_u)(n3,2)*_dSh[2].y;
      eps(ne,3) = (*_u)(n1,1)*_dSh[0].y + (*_u)(n2,1)*_dSh[1].y + (*_u)(n3,1)*_dSh[2].y +
                  (*_u)(n1,2)*_dSh[0].x + (*_u)(n2,2)*_dSh[1].x + (*_u)(n3,2)*_dSh[2].x;
   }
}


void Elas2DT3::Stress(Vect<real_t>& s,
                      Vect<real_t>& vm)
{
   s.setSize(_nb_el,2);
   vm.setSize(_nb_el);
   mesh_elements(*_theMesh) {
      size_t ne = The_element.n();
      size_t n1=The_element(1)->n(), n2=The_element(2)->n(), n3=The_element(3)->n();
      Triang3 tr(the_element);
      _dSh = tr.DSh();
      real_t e1 = (*_u)(n1,1)*_dSh[0].x + (*_u)(n2,1)*_dSh[1].x + (*_u)(n3,1)*_dSh[2].x;
      real_t e2 = (*_u)(n1,2)*_dSh[0].y + (*_u)(n2,2)*_dSh[1].y + (*_u)(n3,2)*_dSh[2].y;
      real_t e3 = (*_u)(n1,1)*_dSh[0].y + (*_u)(n2,1)*_dSh[1].y + (*_u)(n3,1)*_dSh[2].y +
                  (*_u)(n1,2)*_dSh[0].x + (*_u)(n2,2)*_dSh[1].x + (*_u)(n3,2)*_dSh[2].x;
      real_t sx=_E1*e1+_E2*e2, sy=_E2*e2+_E3*e3;
      real_t txy=_E6*e3;
      real_t delta=sqrt((sx+sy)*(sx+sy)+4*(txy*txy-sx*sy));
      real_t s1=0.5*(sx+sy-delta), s2=0.5*(sx+sy-delta);
      //      real_t s3=0.5*(sx+sy-delta);
      s(ne,1) = s1;
      s(ne,2) = s2;
      vm(ne) = sqrt(0.5*((s1-s2)*(s1-s2) + s2*s2 + s1*s2));
   }
}

} /* namespace OFELI */
