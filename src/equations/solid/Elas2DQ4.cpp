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

                         Implementation of class ELAS2DQ4
             for 2-D Linear Elasticity Equations with plane deformations
                    using 4-node quadrilateral Q1 finite element

  ==============================================================================*/


#include "equations/solid/Elas2DQ4.h"

namespace OFELI {

Elas2DQ4::Elas2DQ4(const Element* el)
{
   set(el);
}


Elas2DQ4::Elas2DQ4(const Side* sd)
{
   set(sd);
}


Elas2DQ4::Elas2DQ4(const Element*      el,
                   const Vect<real_t>& u,
                   const real_t&       time)
{
   set(el);
   _time = time;
   ElementVector(u);
}


Elas2DQ4::Elas2DQ4(const Side*         sd,
                   const Vect<real_t>& u,
                   const real_t&       time)
{
   set(sd);
   _time = time;
   SideVector(u);
}


Elas2DQ4::~Elas2DQ4()
{
   if (_quad) { delete _quad; _quad = NULL; }
   if (_ln) { delete _ln; _ln = NULL; }
}


void Elas2DQ4::set(const Element* el)
{
   _nb_dof = 2;
   Init(el);
   setMaterial();
   _quad = new Quad4(el);
   _ln = NULL;
   Gauss g(2);
   _g[0] = g.x(1); _g[1] = g.x(2);
   _w[0] = g.w(1); _w[1] = g.w(2);
   _cg = 0.;
   _xl[0] = _yl[0] = _yl[1] = _xl[3] = -1.;
   _xl[1] = _xl[2] = _yl[2] = _yl[3] =  0.;
   PlaneStrain();
   ElementNodeCoordinates();
   eMat = 0; eA0 = 0; eA1 = 0; eA2 = 0;
   eRHS = 0;
}


void Elas2DQ4::set(const Side* sd)
{
   _nb_dof = 2;
   Init(sd);
   _ln = new Line2(sd);
   _quad = NULL;
   Gauss g(2);
   _g[0] = g.x(1); _g[1] = g.x(2);
   _w[0] = g.w(1); _w[1] = g.w(2);
   _cg = 0.;
   SideNodeCoordinates();
   sMat = 0;
   sRHS = 0;
}


void Elas2DQ4::PlaneStrain(real_t E,
                           real_t nu)
{
   _E = E; _nu = nu;
   _lambda = _nu*_E/((1+_nu)*(1-2*_nu));
   _G = 0.5*_E/(1+_nu);
   _E1 = _lambda*(1-_nu)/_nu;
   _E2 = _lambda;
   _E3 = _E1;
   _E6 = _G;
}


void Elas2DQ4::PlaneStrain()
{
   _lambda = _nu*_E/((1+_nu)*(1-2*_nu));
   _G = 0.5*_E/(1+_nu);
   _E1 = _lambda*(1-_nu)/_nu;
   _E2 = _lambda;
   _E3 = _E1;
   _E6 = _G;
}


void Elas2DQ4::PlaneStress(real_t E,
                           real_t nu)
{
   PlaneStrain(E,nu);
   _lambda = 2*_lambda*_G/(_lambda+2*_G);
   _E1 = _E/(1-_nu*_nu);
   _E2 = _E1*_nu;
   _E3 = _E1;
   _E6 = _G;
}


void Elas2DQ4::PlaneStress()
{
   PlaneStrain();
   _lambda = 2*_lambda*_G/(_lambda+2*_G);
   _E1 = _E/(1-_nu*_nu);
   _E2 = _E1*_nu;
   _E3 = _E1;
   _E6 = _G;
}


void Elas2DQ4::LMassToLHS(real_t coef)
{
   for (size_t i=1; i<=4; i++) {
      _quad->setLocal(Point<real_t> (_xl[i-1],_yl[i-1]));
      real_t c = _rho*coef*_quad->getDet();
      eA2(2*i-1,2*i-1) += c;
      eA2(2*i,2*i) += c;
   }
}


void Elas2DQ4::LMassToRHS(real_t coef)
{
   for (size_t i=1; i<=4; i++) {
      _quad->setLocal(Point<real_t>(_xl[i-1],_yl[i-1]));
      real_t c = _rho*coef*_quad->getDet();
      eRHS(2*i-1) += c*ePrev(2*i-1);
      eRHS(2*i  ) += c*ePrev(2*i  );
   }
}


void Elas2DQ4::Deviator(real_t coef)
{
   for (size_t k=0; k<2; k++) {
      for (size_t l=0; l<2; l++) {
         Point<real_t> g(_g[k],_g[l]);
         real_t w = _w[k]*_w[l];
         _quad->setLocal(g);
         real_t c = _G*w*_quad->getDet()*coef;
         for (size_t i=1; i<=4; i++) {
            Point<real_t>  a = c*_quad->DSh(i);
            for (size_t j=1; j<=4; j++) {
               eA0(2*i-1,2*j-1) += 2*a.x*_quad->DSh(j).x + a.y*_quad->DSh(j).y;
               eA0(2*i-1,2*j  ) += a.y*_quad->DSh(j).x;
               eA0(2*i  ,2*j-1) += a.x*_quad->DSh(j).y;
               eA0(2*i  ,2*j  ) += 2*a.y*_quad->DSh(j).y + a.x*_quad->DSh(j).x;
            }
         }
      }
   }
   eMat = eA0;
}


void Elas2DQ4::DeviatorToRHS(real_t coef)
{
   for (size_t k=0; k<2; k++) {
      for (size_t l=0; l<2; l++) {
         Point<real_t> g(_g[k],_g[l]);
         real_t w = _w[k]*_w[l];
         _quad->setLocal(g);
         real_t c = _G*w*_quad->getDet()*coef;
         for (size_t i=1; i<=4; i++) {
            Point<real_t> a = c*_quad->DSh(i);
            for (size_t j=1; j<=4; j++) {
               eRHS(2*i-1) -= (2*a.x*_quad->DSh(j).x + a.y*_quad->DSh(j).y)*ePrev(2*j-1)
                            + (a.y*_quad->DSh(j).x)*ePrev(2*j  );
               eRHS(2*i  ) -= (a.x*_quad->DSh(j).y)*ePrev(2*j-1)
                            + (2*a.y*_quad->DSh(j).y+a.x*_quad->DSh(j).x)*ePrev(2*j  );
            }
         }
      }
   }
}


void Elas2DQ4::Dilatation(real_t coef)
{
   for (size_t k=0; k<2; k++) {
      for (size_t l=0; l<2; l++) {
         Point<real_t> g(_g[k],_g[l]);
         real_t w = _w[k]*_w[l];
         _quad->setLocal(g);
         real_t c = _lambda*w*_quad->getDet()*coef;
         for (size_t i=1; i<=4; i++) {
            Point<real_t> a = c*_quad->DSh(i);
            for (size_t j=1; j<=4; j++) {
               eA0(2*i-1,2*j-1) += a.x*_quad->DSh(j).x;
               eA0(2*i-1,2*j  ) += a.x*_quad->DSh(j).y;
               eA0(2*i  ,2*j-1) += a.y*_quad->DSh(j).x;
               eA0(2*i  ,2*j  ) += a.y*_quad->DSh(j).y;
            }
         }
      }
   }
   eMat = eA0;
}


void Elas2DQ4::DilatationToRHS(real_t coef)
{
   for (size_t k=0; k<2; k++) {
      for (size_t l=0; l<2; l++) {
         Point<real_t> g(_g[k],_g[l]);
         real_t w=_w[k]*_w[l];
         _quad->setLocal(g);
         real_t c=_lambda*w*_quad->getDet()*coef;
         for (size_t i=1; i<=4; i++) {
            Point<real_t> a=c*_quad->DSh(i);
            for (size_t j=1; j<=4; j++) {
               eRHS(2*i-1) -= a.x*_quad->DSh(j).x * ePrev(2*j-1)
                            + a.x*_quad->DSh(j).y * ePrev(2*j  );
               eRHS(2*i  ) -= a.y*_quad->DSh(j).x * ePrev(2*j-1)
                            + a.y*_quad->DSh(j).y * ePrev(2*j  );
            }
         }
      }
   }
}


void Elas2DQ4::BodyRHS(UserData<real_t>& ud)
{
   for (size_t k=0; k<2; k++)
      for (size_t l=0; l<2; l++) {
         Point<real_t> g(_g[k],_g[l]);
         real_t w = _w[k]*_w[l];
         _quad->setLocal(g);
         Point<real_t> x=_quad->getLocalPoint();
         real_t fx=ud.BodyForce(x, _time, 1);
         real_t fy=ud.BodyForce(x, _time, 2);
         real_t c=w*_quad->getDet();
         for (size_t i=1; i<=2; i++) {
            eRHS(2*i-1) += c*fx*_quad->Sh(i);
            eRHS(2*i  ) += c*fy*_quad->Sh(i);
         }
      }
}


void Elas2DQ4::BodyRHS(const Vect<real_t>& bf,
                             int           opt)
{
   for (size_t k=0; k<2; k++)
      for (size_t l=0; l<2; l++) {
         Point<real_t> g(_g[k],_g[l]);
         real_t w=_w[k]*_w[l];
         _quad->setLocal(g);
         real_t c=w*_quad->getDet();
         real_t fx=0., fy=0.;
         if (opt==LOCAL_ARRAY) {
            for (size_t j=1; j<=4; j++) {
               fx += _quad->Sh(j)*bf(2*j-1);
               fy += _quad->Sh(j)*bf(2*j  );
            }
         }
         else {
            for (size_t j=1; j<=4; j++) {
               fx += _quad->Sh(j)*bf(2*(*_theElement)(j)->n()-1);
               fy += _quad->Sh(j)*bf(2*(*_theElement)(j)->n()  );
            }
         }
         for (size_t i=1; i<=4; i++) {
            eRHS(2*i-1) += c*fx*_quad->Sh(i);
            eRHS(2*i  ) += c*fy*_quad->Sh(i);
         }
      }
}


void Elas2DQ4::BoundaryRHS(UserData<real_t>& ud)
{
   for (size_t k=0; k<2; k++) {
      Point<real_t> x=_ln->getLocalPoint(_g[k]);
      real_t fx=ud.SurfaceForce(x, _theSide->getCode(1), _time, 1);
      real_t fy=ud.SurfaceForce(x, _theSide->getCode(2), _time, 2);
      real_t c=0.5*_w[k]*_ln->getLength();
      for (size_t i=1; i<=2; i++) {
         if (_theSide->getCode(1) > 0)
            sRHS(2*i-1) += c*fx*_ln->Sh(i,_g[k]);
         if (_theSide->getCode(2) > 0)
            sRHS(2*i  ) += c*fy*_ln->Sh(i,_g[k]);
      }
   }
}


void Elas2DQ4::BoundaryRHS(const Vect<real_t>& sf)
{
   for (size_t k=0; k<2; k++) {
      real_t fx=0., fy=0.;
      for (size_t j=1; j<=2; j++) {
         fx += _ln->Sh(j,_g[k])*sf(2*_theSide->n()-1);
         fy += _ln->Sh(j,_g[k])*sf(2*_theSide->n()  );
      }
      real_t c=0.5*_w[k]*_ln->getLength();
      for (size_t i=1; i<=2; i++) {
         sRHS(2*i-1) += c*fx*_ln->Sh(i,_g[k]);
         sRHS(2*i  ) += c*fy*_ln->Sh(i,_g[k]);
      }
   }
}


int Elas2DQ4::SignoriniContact(UserData<real_t>& ud,
                               real_t            coef)
{
   int ret = 0;
   real_t g;
   real_t c=0.5*coef*_ln->getLength();
   (*_theSide)(1)->setCode(1,0); (*_theSide)(1)->setCode(2,0);
   (*_theSide)(2)->setCode(1,0); (*_theSide)(2)->setCode(2,0);
   int c1=_theSide->getCode(1), c2=_theSide->getCode(2);
   if (c1<0) {
      g = ud.SurfaceForce(_x[0], c1, _time, 1);
      if (g>ePrev(1)) {
         ret = 1;
         sMat(1,1) += c;
         sRHS(1) += c*g;
         (*_theSide)(1)->setCode(1,-1);
      }
      g = ud.SurfaceForce(_x[1], c1, _time, 1);
      if (g>ePrev(3)) {
         ret = 1;
         (*_theSide)(2)->setCode(1,-1);
         sMat(3,3) += c;
         sRHS(3) += c*g;
      }
    }
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


void Elas2DQ4::Strain(LocalVect<real_t,3>& eps)
{
    _quad->setLocal(Point<real_t> (0.,0.));
    eps[0] = ePrev[0]*_quad->DSh(1).x + ePrev[2]*_quad->DSh(2).x +
             ePrev[4]*_quad->DSh(3).x + ePrev[6]*_quad->DSh(4).x;
    eps[1] = ePrev[1]*_quad->DSh(1).y + ePrev[3]*_quad->DSh(2).y +
             ePrev[5]*_quad->DSh(3).y + ePrev[7]*_quad->DSh(4).y;
    eps[2] = ePrev[0]*_quad->DSh(1).y + ePrev[2]*_quad->DSh(2).y +
             ePrev[4]*_quad->DSh(3).y + ePrev[6]*_quad->DSh(4).y +
             ePrev[1]*_quad->DSh(1).x + ePrev[3]*_quad->DSh(2).x +
             ePrev[5]*_quad->DSh(3).x + ePrev[7]*_quad->DSh(4).x;
}


void Elas2DQ4::Stress(LocalVect<real_t,3>& s,
                      real_t&              vm)
{
    real_t dudx, dvdy, e, sx, sy, txy, delta;
    _quad->setLocal(Point<real_t> (0.,0.));
    dudx = ePrev[0]*_quad->DSh(1).x + ePrev[2]*_quad->DSh(2).x +
           ePrev[4]*_quad->DSh(3).x + ePrev[6]*_quad->DSh(4).x;
    dvdy = ePrev[1]*_quad->DSh(1).y + ePrev[3]*_quad->DSh(2).y +
           ePrev[5]*_quad->DSh(3).y + ePrev[7]*_quad->DSh(4).y;
    e = ePrev[0]*_quad->DSh(1).y + ePrev[2]*_quad->DSh(2).y +
        ePrev[4]*_quad->DSh(3).y + ePrev[6]*_quad->DSh(4).y +
        ePrev[1]*_quad->DSh(1).x + ePrev[3]*_quad->DSh(2).x +
        ePrev[5]*_quad->DSh(3).x + ePrev[7]*_quad->DSh(4).x;
    sx = _E1*dudx + _E2*dvdy;
    sy = _E2*dvdy + _E3*e;
    txy = _E6*e;
    delta = sqrt((sx+sy)*(sx+sy) + 4*(txy*txy-sx*sy));
    s[0] = 0.5*(sx+sy-delta);
    s[1] = 0.5*(sx+sy+delta);
    s[2] = _nu*(sx+sy);
    vm = (s[0]-s[1])*(s[0]-s[1]) + (s[1]-s[2])*(s[1]-s[2]) + (s[2]-s[0])*(s[2]-s[0]);
    vm = sqrt(0.5*vm);
}


void Elas2DQ4::Stress(LocalVect<real_t,3>& sigma,
                      LocalVect<real_t,3>& s,
                      real_t&              vm)
{
    real_t dudx, dvdy, e, sx, sy, txy, delta;
    _quad->setLocal(Point<real_t> (0.,0.));
    dudx = ePrev[0]*_quad->DSh(1).x + ePrev[2]*_quad->DSh(2).x +
           ePrev[4]*_quad->DSh(3).x + ePrev[6]*_quad->DSh(4).x;
    dvdy = ePrev[1]*_quad->DSh(1).y + ePrev[3]*_quad->DSh(2).y +
           ePrev[5]*_quad->DSh(3).y + ePrev[7]*_quad->DSh(4).y;
    e = ePrev[0]*_quad->DSh(1).y + ePrev[2]*_quad->DSh(2).y +
        ePrev[4]*_quad->DSh(3).y + ePrev[6]*_quad->DSh(4).y +
        ePrev[1]*_quad->DSh(1).x + ePrev[3]*_quad->DSh(2).x +
        ePrev[5]*_quad->DSh(3).x + ePrev[7]*_quad->DSh(4).x;
    sx = _E1*dudx + _E2*dvdy;
    sy = _E2*dvdy + _E3*e;
    txy = _E6*e;
    delta = sqrt((sx+sy)*(sx+sy) + 4*(txy*txy-sx*sy));
    s[0] = 0.5*(sx+sy-delta);
    s[1] = 0.5*(sx+sy+delta);
    s[2] = _nu*(sx+sy);
    vm = (s[0]-s[1])*(s[0]-s[1]) + (s[1]-s[2])*(s[1]-s[2]) + (s[2]-s[0])*(s[2]-s[0]);
    vm = sqrt(0.5*vm);
    real_t dudy = ePrev[0]*_quad->DSh(1).y + ePrev[2]*_quad->DSh(2).y +
                  ePrev[4]*_quad->DSh(3).y + ePrev[6]*_quad->DSh(4).y;
    real_t dvdx = ePrev[1]*_quad->DSh(1).x + ePrev[3]*_quad->DSh(2).x +
                  ePrev[5]*_quad->DSh(3).x + ePrev[7]*_quad->DSh(4).x;
    sigma[0] = _lambda*(dudx + dvdy) + _G*dudx;
    sigma[1] = _lambda*(dudx + dvdy) + _G*dvdy;
    sigma[2] = 0.5*_G*(dudy + dvdx);
}

} /* namespace OFELI */
