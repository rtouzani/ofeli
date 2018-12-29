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

namespace OFELI {

Elas2DT3::Elas2DT3()
{
}


Elas2DT3::Elas2DT3(const Element* el)
{
   set(el);
}


Elas2DT3::Elas2DT3(const Side* sd)
{
   set(sd);
}


Elas2DT3::Elas2DT3(const Element*      el,
                   const Vect<real_t>& u,
                         real_t        time)
{
   _time  = time;
   set(el);
   ElementVector(u);
}


Elas2DT3::Elas2DT3(const Element*      el,
                   const Vect<real_t>& u,
                         real_t        time,
                         real_t        deltat,
                         int           scheme)
{
   _time  = time;
   set(el);
   ElementVector(u);
   _time_step = deltat;
   setTimeIntegration(scheme);
}


Elas2DT3::Elas2DT3(const Side*         sd,
                   const Vect<real_t>& u,
                         real_t        time)
{
   _time  = time;
   set(sd);
   SideVector(u);
}


Elas2DT3::Elas2DT3(const Side*         sd,
                   const Vect<real_t>& u,
                   real_t              time,
                   real_t              deltat,
                   int                 scheme)
{
   _time  = time;
   set(sd);
   SideVector(u);
   _time_step = deltat;
   setTimeIntegration(scheme);
}


Elas2DT3::Elas2DT3(Mesh& ms) 
         : Equation<real_t,3,6,2,4>(ms)
{
   _equation_name = "Linearized Elasticity";
   _finite_element = "2-D, 3-Node Triangles (P1)";
}


Elas2DT3::~Elas2DT3() { }


void Elas2DT3::set(const Element* el)
{
   _nb_dof = 2;
   Init(el);
   setMaterial();
   PlaneStrain();
   Triang3 tr(el);
   _area = tr.getArea();
   _dSh(1) = tr.DSh(1);
   _dSh(2) = tr.DSh(2);
   _dSh(3) = tr.DSh(3);
   ElementNodeCoordinates();
   eMat = 0; eA0 = 0; eA1 = 0; eA2 = 0;
   eRHS = 0;
}


void Elas2DT3::set(const Side* sd)
{
   _nb_dof = 2;
   Init(sd);
   Line2 ln(sd);
   _length = ln.getLength();
   SideNodeCoordinates();
   sMat = 0;
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


void Elas2DT3::LMassToLHS(real_t coef)
{
   real_t c=_rho*coef*_area*OFELI_THIRD;
   for (size_t i=1; i<=3; i++) {
      eA2(2*i-1,2*i-1) += c;
      eA2(2*i  ,2*i  ) += c;
   }
}


void Elas2DT3::LMassToRHS(real_t coef)
{
   real_t c=_rho*coef*_area*OFELI_THIRD;
   for (size_t i=1; i<=3; i++) {
      eRHS(2*i-1) += c*ePrev(2*i-1);
      eRHS(2*i  ) += c*ePrev(2*i  );
   }
}


void Elas2DT3::MassToLHS(real_t coef)
{
   real_t c=0.5*OFELI_SIXTH*_area*_rho*coef;
   for (size_t i=1; i<=3; i++) {
      for (size_t j=1; j<=3; j++) {
         eA2(2*i-1,2*j-1) += c;
         eA2(2*i  ,2*j  ) += c;
      }
      eA2(2*i-1,2*i-1) += c;
      eA2(2*i  ,2*i  ) += c;
   }
}


void Elas2DT3::MassToRHS(real_t coef)
{
   real_t c=_rho*coef*_area*OFELI_THIRD;
   real_t d=0.5*c;
   eRHS(1) += c*ePrev(1) + d*(ePrev(3) + ePrev(5));
   eRHS(2) += c*ePrev(2) + d*(ePrev(4) + ePrev(6));
   eRHS(3) += c*ePrev(3) + d*(ePrev(1) + ePrev(5));
   eRHS(4) += c*ePrev(4) + d*(ePrev(2) + ePrev(6));
   eRHS(5) += c*ePrev(5) + d*(ePrev(1) + ePrev(3));
   eRHS(6) += c*ePrev(6) + d*(ePrev(2) + ePrev(4));
}


void Elas2DT3::Deviator(real_t coef)
{
_G = 1.;
   real_t c=_G*_area*coef;
   for (size_t i=1; i<=3; i++) {
      Point<real_t> a=c*_dSh(i);
      for (size_t j=1; j<=3; j++) {
         eA0(2*i-1,2*j-1) += 2*a.x*_dSh(j).x + a.y*_dSh(j).y;
         eA0(2*i-1,2*j  ) += a.y*_dSh(j).x;
         eA0(2*i  ,2*j-1) += a.x*_dSh(j).y;
         eA0(2*i  ,2*j  ) += 2*a.y*_dSh(j).y + a.x*_dSh(j).x;
      }
   }
   /*
   for (size_t i=1; i<=3; i++) {
      for (size_t j=1; j<=3; j++) {
         eA0(2*i-1,2*j-1) += _area*_dSh(i)*_dSh(j);
         eA0(2*i  ,2*j  ) += _area*_dSh(i)*_dSh(j);
      }
      }*/
   eMat = eA0;
//cout<<eMat;
}


void Elas2DT3::DeviatorToRHS(real_t coef)
{
   real_t c = _G*_area*coef;
   for (size_t i=1; i<=3; i++) {
      Point<real_t> a=c*_dSh(i);
      for (size_t j=1; j<=3; j++) {
         eRHS(2*i-1) -= (2*a.x*_dSh(j).x + a.y*_dSh(j).y)*ePrev(2*j-1)
                      + (a.y*_dSh(j).x)*ePrev(2*j);
         eRHS(2*i  ) -= (a.x*_dSh(j).y)*ePrev(2*j-1)
                      + (2*a.y*_dSh(j).y + a.x*_dSh(j).x)*ePrev(2*j);
      }
   }
}


void Elas2DT3::Dilatation(real_t coef)
{
   real_t c=_lambda*_area*coef;
   for (size_t i=1; i<=3; i++) {
      Point<real_t> a=c*_dSh(i);
      for (size_t j=1; j<=3; j++) {
         eA0(2*i-1,2*j-1) += a.x*_dSh(j).x;
         eA0(2*i-1,2*j  ) += a.x*_dSh(j).y;
         eA0(2*i  ,2*j-1) += a.y*_dSh(j).x;
         eA0(2*i  ,2*j  ) += a.y*_dSh(j).y;
      }
   }
   eMat = eA0;
}


void Elas2DT3::DilatationToRHS(real_t coef)
{
   real_t c=_lambda*_area*coef;
   for (size_t i=1; i<=3; i++) {
      Point<real_t> a = c*_dSh(i);
      for (size_t j=1; j<=3; j++) {
         eRHS(2*i-1) -= a.x*_dSh(j).x * ePrev(2*j-1) + a.x*_dSh(j).y * ePrev(2*j);
         eRHS(2*i  ) -= a.y*_dSh(j).x * ePrev(2*j-1) + a.y*_dSh(j).y * ePrev(2*j);
      }
   }
}


void Elas2DT3::BodyRHS(UserData<real_t>& ud)
{
   for (size_t k=1; k<=3; k++) {
      real_t fx = ud.BodyForce(_x[k-1], _time, 1),
             fy = ud.BodyForce(_x[k-1], _time, 2);
      eRHS(2*k-1) += OFELI_THIRD*_area*fx;
      eRHS(2*k  ) += OFELI_THIRD*_area*fy;
   }
}


void Elas2DT3::BodyRHS(const Vect<real_t>& f,
                             int           opt)
{
   if (opt==LOCAL_ARRAY) {
     for (size_t k=1; k<=3; k++) {
        eRHS(2*k-1) += OFELI_THIRD*_area*f((*_theElement)(k)->n(),1);
        eRHS(2*k  ) += OFELI_THIRD*_area*f((*_theElement)(k)->n(),2);
     }
   }
   else {
     for (size_t k=1; k<=3; k++) {
        eRHS(2*k-1) += OFELI_THIRD*_area*f(k,1);
        eRHS(2*k  ) += OFELI_THIRD*_area*f(k,2);
     }
   }
}


void Elas2DT3::BoundaryRHS(UserData<real_t>& ud)
{
   real_t c = 0.5*_length;
   for (size_t k=1; k<=2; k++) {
      if (_theSide->getCode(1) != CONTACT)
         sRHS(2*k-1) += c*ud.SurfaceForce(_x[k-1],_theSide->getCode(1),_time,1);
      if (_theSide->getCode(2) != CONTACT)
         sRHS(2*k  ) += c*ud.SurfaceForce(_x[k-1],_theSide->getCode(2),_time,2);
   }
}


void Elas2DT3::BoundaryRHS(const Vect<real_t>& f)
{
   real_t c = 0.5*_length;
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
      real_t c = 0.5*_length*coef;
      if (_theSide->getCode(1) == PERIODIC_A)
         sMat(2*i-1,2*i-1) += c;
      else if (_theSide->getCode(1) == PERIODIC_B)
         sMat(2*i-1,2*i-1) -= c;
      if (_theSide->getCode(2) == PERIODIC_A)
         sMat(2*i  ,2*i  ) += c;
      else if (_theSide->getCode(2) == PERIODIC_B)
         sMat(2*i  ,2*i  ) -= c;
   }
}


int Elas2DT3::SignoriniContact(UserData<real_t>& ud,
                               real_t            coef)
{
    int ret = 0;
    real_t g, c=0.5*coef*_length;
    (*_theSide)(1)->setCode(1,0); (*_theSide)(1)->setCode(2,0);
    (*_theSide)(2)->setCode(1,0); (*_theSide)(2)->setCode(2,0);
    int c1=_theSide->getCode(1);
    if (c1<0) {
       g = ud.SurfaceForce(_x[0], c1, _time, 1);
       if (g>ePrev(1)) {
          ret = 1;
          sMat(1,1) += c;
          sRHS(1) += c*g;
          (*_theSide)(1)->setCode(1,-1);
      }
      g = ud.SurfaceForce(_x[1], c1, _time, 1);
      if (g > ePrev(3)) {
         ret = 1;
         (*_theSide)(2)->setCode(1,-1);
         sMat(3,3) += c;
         sRHS(3) += c*g;
      }
    }
    int c2=_theSide->getCode(2);
    if (c2<0) {
       g = ud.SurfaceForce(_x[0], c2, _time, 2);
       if (g>ePrev(2)) {
          ret = 1;
          (*_theSide)(1)->setCode(2,-1);
          sMat(2,2) += c;
          sRHS(2) += c*g;
      }
      g = ud.SurfaceForce(_x[1], c2, _time, 2);
      if (g>ePrev(4)) {
         ret = 1;
         (*_theSide)(2)->setCode(2,-1);
         sMat(4,4) += c;
         sRHS(4) += c*g;
      }
    }
    return ret;
}


int Elas2DT3::SignoriniContact(Vect<real_t>& f,
                               real_t        coef)
{
   int ret = 0;
   real_t c=0.5*coef*_length;
   real_t fx=f(_theSide->n(),1), fy=f(_theSide->n(),2);
   if (_theSide->getCode(1) == CONTACT) {
      if (fx > ePrev(1)) {
         ret = 1;
         sMat(1,1) += c;
         sRHS(1) += c*fx;
      }
      if (fx > ePrev(3)) {
         ret = 1;
         sMat(3,3) += c;
         sRHS(3) += c*fx;
      }
   }

   if (_theSide->getCode(2) == CONTACT) {
      if (fy > ePrev(2)) {
         ret = 1;
         sMat(2,2) += c;
         sRHS(2) += c*fy;
      }
      if (fy > ePrev(4)) {
         ret = 1;
         sMat(4,4) += c;
         sRHS(4) += c*fy;
     }
   }
   return ret;
}


void Elas2DT3::Reaction(Vect<real_t>& r)
{
   size_t n1, n2;
   Deviator();
   Dilatation();
   for (size_t s=1; s<=3; s++) {
      Side *sd = _theElement->getPtrSide(s);
      size_t t = (s+1)%3;
      if (sd) {
	 n1 = sd->getNodeLabel(1);
         n2 = sd->getNodeLabel(2);
         if (sd->getCode(1)==CONTACT) {
            r(n1,1) += eMat(2*s-1,2*s-1)*ePrev(2*s-1) +  eMat(2*s-1,2*t-1)*ePrev(2*t-1);
            r(n2,1) += eMat(2*t-1,2*s-1)*ePrev(2*s-1) +  eMat(2*t-1,2*t-1)*ePrev(2*t-1);
         }
         if (sd->getCode(2)==CONTACT) {
            r(n1,2) += eMat(2*s  ,2*s  )*ePrev(2*s  ) +  eMat(2*s  ,2*t  )*ePrev(2*t  );
            r(n2,2) += eMat(2*t  ,2*s  )*ePrev(2*s  ) +  eMat(2*t  ,2*t  )*ePrev(2*t  );
         }
      }
   }
}


void Elas2DT3::ContactPressure(const Vect<real_t>&  f,
                                     real_t         penal,
                                     Point<real_t>& p)
{
   p = 0;
   real_t ff = f(_theSide->n(),1);
   if (_theSide->getCode(1)==CONTACT) {
      if (ff>=ePrev(1) && ff>=ePrev(3))
         p.x = penal*(f(_theSide->n(),1)-0.5*(ePrev(1)+ePrev(3)));
   }
   ff = f(_theSide->n(),2);
   if (_theSide->getCode(2)==CONTACT) {
      if (ff>=ePrev(2) && ff>=ePrev(4))
         p.y = penal*(f(_theSide->n(),2)-0.5*(ePrev(2)+ePrev(4)));
   }
}


void Elas2DT3::Strain(Vect<real_t>& eps)
{
    eps[0] = ePrev[0]*_dSh(1).x + ePrev[2]*_dSh(2).x + ePrev[4]*_dSh(3).x;
    eps[1] = ePrev[1]*_dSh(1).y + ePrev[3]*_dSh(2).y + ePrev[5]*_dSh(3).y;
    eps[2] = ePrev[0]*_dSh(1).y + ePrev[2]*_dSh(2).y + ePrev[4]*_dSh(3).y +
             ePrev[1]*_dSh(1).x + ePrev[3]*_dSh(2).x + ePrev[5]*_dSh(3).x;
}


void Elas2DT3::Stress(Vect<real_t>& s,
                      real_t&       vm)
{
    real_t e1 = ePrev[0]*_dSh(1).x + ePrev[2]*_dSh(2).x + ePrev[4]*_dSh(3).x,
           e2 = ePrev[1]*_dSh(1).y + ePrev[3]*_dSh(2).y + ePrev[5]*_dSh(3).y,
           e3 = ePrev[0]*_dSh(1).y + ePrev[2]*_dSh(2).y + ePrev[4]*_dSh(3).y +
                ePrev[1]*_dSh(1).x + ePrev[3]*_dSh(2).x + ePrev[5]*_dSh(3).x;
    real_t sx=_E1*e1+_E2*e2, sy=_E2*e2+_E3*e3;
    real_t txy=_E6*e3;
    real_t delta=sqrt((sx+sy)*(sx+sy)+4*(txy*txy-sx*sy));
    s[0] = 0.5*(sx+sy-delta);
    s[1] = 0.5*(sx+sy+delta);
    vm = (s[0]-s[1])*(s[0]-s[1]) + (s[1]-s[2])*(s[1]-s[2]) + (s[2]-s[0])*(s[2]-s[0]);
    vm = sqrt(0.5*vm);
}

} /* namespace OFELI */
