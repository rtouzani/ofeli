/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2024 Rachid Touzani

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
         : Equation<3,6,2,4>(ms)
{
   _equation_name = "Linearized Elasticity";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _terms = int(PDE_Terms::DEVIATORIC)|int(PDE_Terms::DILATATION)|int(PDE_Terms::BODY_RHS);
   bool contact = false;
   MESH_SD {
      if (The_side.getGlobalCode() == CONTACT_BC) {
         contact = true;
         break;
      }
   }
   if (contact)
      _terms = int(PDE_Terms::DEVIATORIC)|int(PDE_Terms::DILATATION)|int(PDE_Terms::BODY_RHS)|
               int(PDE_Terms::CONTACT);
}


Elas2DT3::Elas2DT3(Mesh&         ms,
                   Vect<real_t>& u) 
         : Equation<3,6,2,4>(ms,u)
{
   _equation_name = "Linearized Elasticity";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _terms = int(PDE_Terms::DEVIATORIC)|int(PDE_Terms::DILATATION)|int(PDE_Terms::BODY_RHS);
   bool contact = false;
   MESH_SD {
      if (The_side.getGlobalCode() == CONTACT_BC) {
         contact = true;
         break;
      }
   }
   if (contact)
      _terms = int(PDE_Terms::DEVIATORIC)|int(PDE_Terms::DILATATION)|int(PDE_Terms::BODY_RHS)|
               int(PDE_Terms::CONTACT);
}


Elas2DT3::~Elas2DT3() { }


void Elas2DT3::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   setMaterial();
   PlaneStrain();
   Triang3 tr(el);
   _el_geo.area = tr.getArea();
   _el_geo.center = tr.getCenter();
   _dSh = tr.DSh();
   ElementNodeCoordinates();
   ElementNodeVector(*_u,_eu);
   if (_terms&int(PDE_Terms::SOURCE) && Equa::_bf!=nullptr)
      ElementNodeVector(*_bf,_ebf);
   if (_rho_set)
      _rho = _rho_fct(_el_geo.center,_TimeInt.time);
   if (_young_set)
      _young = _young_fct(_el_geo.center,_TimeInt.time);
   if (_poisson_set)
      _poisson = _poisson_fct(_el_geo.center,_TimeInt.time);
   eA0 = 0, eA1 = 0, eA2 = 0;
   eMat = 0;
   eRHS = 0;
}


void Elas2DT3::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Line2 ln(sd);
   _el_geo.length = ln.getLength();
   SideNodeCoordinates();
   SideNodeVector(*_u,_su);
   if (Equa::_sf!=nullptr)
      SideVector(*_sf,_ssf);
   sA0 = 0;
   sMat = 0;
   sRHS = 0;
}


void Elas2DT3::Media(real_t E,
                     real_t nu,
                     real_t rho)
{
   _young = E;
   _poisson = nu;
   _rho = rho;
   PlaneStrain(E,nu);
}


void Elas2DT3::PlaneStrain(real_t E,
                           real_t nu)
{
   _young = E;
   _poisson = nu;
   _lambda = _poisson*_young/((1+_poisson)*(1-2*_poisson));
   _G = 0.5*_young/(1+_poisson);
   _E1 = _lambda*(1-_poisson)/_poisson;
   _E2 = _lambda;
   _E3 = _E1;
   _E6 = _G;
}


void Elas2DT3::PlaneStrain()
{
   _lambda = _poisson*_young/((1+_poisson)*(1-2*_poisson));
   _G = 0.5*_young/(1+_poisson);
   _E1 = _lambda*(1-_poisson)/_poisson;
   _E2 = _lambda;
   _E3 = _E1;
   _E6 = _G;
}


void Elas2DT3::PlaneStress(real_t E,
                           real_t nu)
{
   PlaneStrain(E,nu);
   _lambda = 2*_lambda*_G/(_lambda+2*_G);
   _E1 = _young/(1-_poisson*_poisson);
   _E2 = _E1*_poisson;
   _E3 = _E1;
   _E6 = _G;
}


void Elas2DT3::PlaneStress()
{
   PlaneStrain();
   _lambda = 2*_lambda*_G/(_lambda+2*_G);
   _E1 = _young/(1-_poisson*_poisson);
   _E2 = _E1*_poisson;
   _E3 = _E1;
   _E6 = _G;
}


void Elas2DT3::LMass(real_t coef)
{
   real_t c=_rho*coef*_el_geo.area*OFELI_THIRD;
   for (size_t i=0; i<3; ++i) {
      eA2(2*i+1,2*i+1) += c;
      eA2(2*i+2,2*i+2) += c;
   }
}


void Elas2DT3::Mass(real_t coef)
{
   real_t c=0.5*OFELI_SIXTH*_el_geo.area*_rho*coef;
   for (size_t i=0; i<3; ++i) {
      for (size_t j=0; j<3; ++j) {
         eA2(2*i+1,2*j+1) += c;
         eA2(2*i+2,2*j+2) += c;
      }
      eA2(2*i+1,2*i+1) += c;
      eA2(2*i+2,2*i+2) += c;
   }
}


void Elas2DT3::Deviator(real_t coef)
{
   real_t c=_G*_el_geo.area*coef;
   for (size_t i=0; i<3; ++i) {
      Point<real_t> a=c*_dSh[i];
      for (size_t j=0; j<3; ++j) {
         eA0(2*i+1,2*j+1) += 2*a.x*_dSh[j].x + a.y*_dSh[j].y;
         eA0(2*i+1,2*j+2) += a.y*_dSh[j].x;
         eA0(2*i+2,2*j+1) += a.x*_dSh[j].y;
         eA0(2*i+2,2*j+2) += 2*a.y*_dSh[j].y + a.x*_dSh[j].x;
      }
   }
   eMat += eA0;
}


void Elas2DT3::Dilatation(real_t coef)
{
   real_t c=_lambda*_el_geo.area*coef;
   for (size_t i=0; i<3; ++i) {
      Point<real_t> a=c*_dSh[i];
      for (size_t j=0; j<3; ++j) {
         eA0(2*i+1,2*j+1) += a.x*_dSh[j].x;
         eA0(2*i+1,2*j+2) += a.x*_dSh[j].y;
         eA0(2*i+2,2*j+1) += a.y*_dSh[j].x;
         eA0(2*i+2,2*j+2) += a.y*_dSh[j].y;
      }
   }
   eMat += eA0;
}


void Elas2DT3::BodyRHS()
{
   for (size_t k=0; k<3; ++k) {
      eRHS(2*k+1) += OFELI_THIRD*_el_geo.area*_ebf(2*k+1);
      eRHS(2*k+2) += OFELI_THIRD*_el_geo.area*_ebf(2*k+2);
   }
}


void Elas2DT3::BodyRHS(const Vect<real_t>& f)
{
   for (size_t k=0; k<3; ++k) {
      eRHS(2*k+1) += OFELI_THIRD*_el_geo.area*f(k+1,1);
      eRHS(2*k+2) += OFELI_THIRD*_el_geo.area*f(k+1,2);
   }
}


void Elas2DT3::BoundaryRHS(const Vect<real_t>& f)
{
   if (_theSide->getCode(1)>0) {
      real_t z = 0.5*_el_geo.length;
      if (f.getDOFType()==NODE_DOF) {
         sRHS(1) += z*f(2*(*_theSide)(1)->n()-1);
         sRHS(3) += z*f(2*(*_theSide)(2)->n()-1);
      }
      else if (f.getDOFType()==SIDE_DOF) {
         sRHS(1) += z*f(2*_theSide->n()-1);
         sRHS(3) += z*f(2*_theSide->n()-1);
      }
   }
   if (_theSide->getCode(2)>0) {
      real_t z = 0.5*_el_geo.length;
      if (f.getDOFType()==NODE_DOF) {
         sRHS(2) += z*f(2*(*_theSide)(1)->n());
         sRHS(4) += z*f(2*(*_theSide)(2)->n());
      }
      else if (f.getDOFType()==SIDE_DOF) {
         sRHS(2) += z*f(2*_theSide->n());
         sRHS(4) += z*f(2*_theSide->n());
      }
   }
}


void Elas2DT3::BoundaryRHS()
{
   real_t c = 0.5*_el_geo.length;
   if (_theSide->getCode(1) != CONTACT_BC) {
      sRHS(1) += c*_ssf[0];
      sRHS(3) += c*_ssf[0];
   }
   if (_theSide->getCode(2) != CONTACT_BC) {
      sRHS(2) += c*_ssf[1];
      sRHS(4) += c*_ssf[1];
   }
}


void Elas2DT3::Periodic(real_t coef)
{
   for (size_t i=1; i<=2; ++i) {
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
   sMat += sA0;
}


int Elas2DT3::Contact(real_t coef)
{
   if (_sf==nullptr)
      throw OFELIException("In Elas2DT3::Contact(coef): No contact distance provided.");
   int ret = 0;
   if (_theSide->getGlobalCode() != CONTACT_BC)
      return ret;
   real_t c=0.5*coef*_el_geo.length;
   real_t d1 = (*_sf)(_theSide->n(),1), d2 = (*_sf)(_theSide->n(),2);
   if (d1 > (*_u)((*_theSide)(1)->n(),1) && (*_theSide)(1)->getCode(1)==0) {
      ret = 1;
      sA0(1,1) += c;
      sRHS(1) += c*d1;
   }
   if (d1 > (*_u)((*_theSide)(2)->n(),1) && (*_theSide)(2)->getCode(1)==0) {
      ret = 1;
      sA0(3,3) += c;
      sRHS(3) += c*d1;
   }
   if (d2 > (*_u)((*_theSide)(1)->n(),2) && (*_theSide)(1)->getCode(2)==0) {
      ret = 1;
      sA0(2,2) += c;
      sRHS(2) += c*d2;
   }
   if (d2 > (*_u)((*_theSide)(2)->n(),2) && (*_theSide)(2)->getCode(2)==0) {
      ret = 1;
      sA0(4,4) += c;
      sRHS(4) += c*d2;
   }
   sMat += sA0;
   return ret;
}


void Elas2DT3::Reaction(Vect<real_t>& r)
{
   r.setSize(_nb_sides,2);
   MESH_EL {
      set(the_element);
      Deviator();
      Dilatation();
      for (size_t s=1; s<=3; ++s) {
         Side *sd = _theElement->getPtrSide(s);
         size_t t = (s+1)%3;
         if (sd) {
            size_t n1 = (*sd)(1)->n(), n2 = (*sd)(2)->n();
            real_t u1 = (*_u)(n1,1), u2 = (*_u)(n2,1);
            real_t v1 = (*_u)(n1,2), v2 = (*_u)(n2,2);
            if (sd->getCode(1)==CONTACT_BC) {
               r(n1,1) += eA0(2*s-1,2*s-1)*u1 + eA0(2*s-1,2*t-1)*u2;
               r(n2,1) += eA0(2*t-1,2*s-1)*u1 + eA0(2*t-1,2*t-1)*u2;
            }
            if (sd->getCode(2)==CONTACT_BC) {
               r(n1,2) += eA0(2*s  ,2*s  )*v1 + eA0(2*s  ,2*t  )*v2;
               r(n2,2) += eA0(2*t  ,2*s  )*v1 + eA0(2*t  ,2*t  )*v2;
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
   if (_theSide->getCode(1)==CONTACT_BC) {
      if (ff>=_eu(1) && ff>=_eu(3))
         p.x = penal*(ff-0.5*(_eu(1)+_eu(3)));
   }
   ff = f(_theSide->n(),2);
   if (_theSide->getCode(2)==CONTACT_BC) {
     if (ff>=_su(2) && ff>=_su(4))
         p.y = penal*(ff-0.5*(_su(2)+_su(4)));
   }
}


void Elas2DT3::Strain(Vect<real_t>& eps)
{
   eps.setSize(_nb_el,3);
   MESH_EL {
      size_t ne = element_label;
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
   MESH_EL {
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
      s(element_label,1) = s1, s(element_label,2) = s2;
      vm(element_label) = sqrt(0.5*((s1-s2)*(s1-s2) + s2*s2 + s1*s2));
   }
}

} /* namespace OFELI */