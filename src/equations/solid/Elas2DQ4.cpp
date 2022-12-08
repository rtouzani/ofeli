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

                         Implementation of class ELAS2DQ4
             for 2-D Linear Elasticity Equations with plane deformations
                    using 4-node quadrilateral Q1 finite element

  ==============================================================================*/


#include "equations/solid/Elas2DQ4.h"
#include "util/Gauss.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

Elas2DQ4::Elas2DQ4(Mesh& ms)
         : Equation<4,8,2,4>(ms), _quad(nullptr), _ln(nullptr)
{
   _equation_name = "Linearized elasticity";
   _finite_element = "2-D, 4-Node quadrilaterals (Q1)";
   setMatrixType(SPARSE|SYMMETRIC);
   setSolver(CG_SOLVER,DIAG_PREC);
   _terms = DEVIATORIC|DILATATION|BODY_FORCE|TRACTION;
}


Elas2DQ4::Elas2DQ4(Mesh&         ms,
                   Vect<real_t>& u)
         : Equation<4,8,2,4>(ms,u), _quad(nullptr), _ln(nullptr)
{
   _equation_name = "Linearized elasticity";
   _finite_element = "2-D, 4-Node quadrilaterals (Q1)";
   setMatrixType(SPARSE|SYMMETRIC);
   setSolver(CG_SOLVER,DIAG_PREC);
   _terms = DEVIATORIC|DILATATION|BODY_FORCE|TRACTION;
}


Elas2DQ4::~Elas2DQ4()
{
   if (_quad != nullptr)
      delete _quad;
   if (_ln != nullptr)
      delete _ln;
}


void Elas2DQ4::set(const Element* el)
{
   _theSide = nullptr, _theElement = el;
   setMaterial();
   if (_quad != nullptr)
      delete _quad, _quad = nullptr;
   if (_ln != nullptr)
      delete _ln, _ln = nullptr;
   _quad = new Quad4(el);
   _xl[0] = _yl[0] = _yl[1] = _xl[3] = -1.;
   _xl[1] = _xl[2] = _yl[2] = _yl[3] =  0.;
   _quad->atGauss(2,_sh,_dSh,_wg);
   _el_geo.center = _quad->getCenter();
   PlaneStrain();
   ElementNodeCoordinates();
   if (Equa::_u!=nullptr)
      ElementNodeVector(*_u,_eu);
   if (Equa::_bf!=nullptr)
      ElementNodeVector(*_bf,_ebf);
   if (_rho_set)
      _rho = _rho_fct(_el_geo.center,_TimeInt.time);
   if (_young_set)
      _young = _young_fct(_el_geo.center,_TimeInt.time);
   if (_poisson_set)
      _poisson = _poisson_fct(_el_geo.center,_TimeInt.time);
   eA0 = 0, eA1 = 0, eA2 = 0;
   eRHS = 0;
}


void Elas2DQ4::set(const Side* sd)
{
   _theSide = sd, _theElement = nullptr;
   if (_quad != nullptr)
      delete _quad, _quad = nullptr;
   if (_ln != nullptr)
      delete _ln, _ln = nullptr;
   _ln = new Line2(sd);
   Gauss g(2);
   _g[0] = g.x(1); _g[1] = g.x(2);
   _ww[0] = g.w(1); _ww[1] = g.w(2);
   SideNodeCoordinates();
   if (Equa::_u!=nullptr)
      SideNodeVector(*_u,_su);
   if (Equa::_sf!=nullptr)
      SideVector(*_sf,_ssf);
   sMat = 0;
   sRHS = 0;
}


void Elas2DQ4::PlaneStrain(real_t E,
                           real_t nu)
{
   _young = E; _poisson = nu;
   _lambda = _poisson*_young/((1+_poisson)*(1-2*_poisson));
   _G = 0.5*_young/(1+_poisson);
   _E1 = _lambda*(1-_poisson)/_poisson;
   _E2 = _lambda;
   _E3 = _E1;
   _E6 = _G;
}


void Elas2DQ4::PlaneStrain()
{
   _lambda = _poisson*_young/((1+_poisson)*(1-2*_poisson));
   _G = 0.5*_young/(1+_poisson);
   _E1 = _lambda*(1-_poisson)/_poisson;
   _E2 = _lambda;
   _E3 = _E1;
   _E6 = _G;
}


void Elas2DQ4::PlaneStress(real_t E,
                           real_t nu)
{
   PlaneStrain(E,nu);
   _lambda = 2*_lambda*_G/(_lambda+2*_G);
   _E1 = _young/(1-_poisson*_poisson);
   _E2 = _E1*_poisson;
   _E3 = _E1;
   _E6 = _G;
}


void Elas2DQ4::PlaneStress()
{
   PlaneStrain();
   _lambda = 2*_lambda*_G/(_lambda+2*_G);
   _E1 = _young/(1-_poisson*_poisson);
   _E2 = _E1*_poisson;
   _E3 = _E1;
   _E6 = _G;
}


void Elas2DQ4::LMass(real_t coef)
{
   for (size_t i=1; i<=4; i++) {
      _quad->setLocal(Point<real_t>(_xl[i-1],_yl[i-1]));
      real_t c = _rho*coef*_quad->getDet();
      eA2(2*i-1,2*i-1) += c;
      eA2(2*i  ,2*i  ) += c;
   }
}


void Elas2DQ4::Deviator(real_t coef)
{
   for (size_t k=0; k<4; ++k) {
      real_t c = _G*_wg[k]*coef;
      for (size_t i=0; i<4; i++) {
         Point<real_t> a = c*_dSh[4*i+k];
         for (size_t j=0; j<4; j++) {
            Point<real_t> b = _dSh[4*j+k];
            eA0(2*i+1,2*j+1) += 2*a.x*b.x + a.y*b.y;
            eA0(2*i+1,2*j+2) += a.y*b.x;
            eA0(2*i+2,2*j+1) += a.x*b.y;
            eA0(2*i+2,2*j+2) += 2*a.y*b.y + a.x*b.x;
         }
      }
   }
}


void Elas2DQ4::Dilatation(real_t coef)
{
   for (size_t k=1; k<=4; ++k) {
      real_t c = _lambda*_wg[k]*coef;
      for (size_t i=0; i<4; i++) {
         Point<real_t> a = c*_dSh[4*i+k];
         for (size_t j=0; j<4; j++) {
            Point<real_t> b = _dSh[4*j+k];
            eA0(2*i+1,2*j+1) += a.x*b.x;
            eA0(2*i+1,2*j+2) += a.x*b.y;
            eA0(2*i+2,2*j+1) += a.y*b.x;
            eA0(2*i+2,2*j+2) += a.y*b.y;
         }
      }
   }
}


void Elas2DQ4::BodyRHS()
{
   for (size_t k=0; k<4; ++k) {
      real_t fx=0., fy=0.;
      for (size_t j=1; j<=4; ++j) {
         fx += _sh[4*j+k]*_ebf(2*j-1);
         fy += _sh[4*j+k]*_ebf(2*j  );
      }
      for (size_t i=0; i<4; ++i) {
         eRHS(2*i+1) += _wg[k]*fx*_sh[4*i+k];
         eRHS(2*i+2) += _wg[k]*fy*_sh[4*i+k];
      }
   }
}


void Elas2DQ4::BodyRHS(const Vect<real_t>& f)
{
   for (size_t k=0; k<4; ++k) {
      real_t fx=0., fy=0.;
      for (size_t j=0; j<4; ++j) {
         fx += _sh[4*j+k]*f((*_theElement)(j+1)->n(),1);
         fy += _sh[4*j+k]*f((*_theElement)(j+1)->n(),2);
      }
      for (size_t i=0; i<4; ++i) {
         eRHS(2*i+1) += _wg[k]*fx*_sh[4*i+k];
         eRHS(2*i+2) += _wg[k]*fy*_sh[4*i+k];
      }
   }
}


void Elas2DQ4::BoundaryRHS()
{
   real_t c=0.5*_ln->getLength();
   for (size_t i=1; i<=2; i++) {
      sRHS(2*i-1) += c*_ssf[0];
      sRHS(2*i  ) += c*_ssf[1];
   }
}


void Elas2DQ4::BoundaryRHS(const Vect<real_t>& sf)
{
   if (_theSide->getCode(1)>0 || _theSide->getCode(2)>0) {
      if (sf.getDOFType()==NODE_DOF) {
         for (size_t k=0; k<2; k++) {
            real_t fx=0., fy=0.;
            for (size_t i=1; i<=2; i++) {
               fx += _ln->Sh(i,_g[k])*sf(2*(*_theSide)(i)->n()-1);
               fy += _ln->Sh(i,_g[k])*sf(2*(*_theSide)(i)->n()  );
            }
            real_t c=0.5*_ww[k]*_ln->getLength();
            for (size_t i=1; i<=2; i++) {
               sRHS(2*i-1) += c*fx*_ln->Sh(i,_g[k]);
               sRHS(2*i  ) += c*fy*_ln->Sh(i,_g[k]);
            }
         }
      }
      else if (sf.getDOFType()==SIDE_DOF) {
         for (size_t k=0; k<2; k++) {
            real_t fx=0., fy=0.;
            for (size_t i=1; i<=2; i++) {
               fx += _ln->Sh(i,_g[k])*sf(2*_theSide->n()-1);
               fy += _ln->Sh(i,_g[k])*sf(2*_theSide->n()  );
            }
            real_t c=0.5*_ww[k]*_ln->getLength();
            for (size_t i=1; i<=2; i++) {
               sRHS(2*i-1) += c*fx*_ln->Sh(i,_g[k]);
               sRHS(2*i  ) += c*fy*_ln->Sh(i,_g[k]);
            }
         }
      }
   }
}


void Elas2DQ4::Strain(Vect<real_t>& eps)
{
   eps.setSize(_nb_el,3);
   MESH_EL {
      size_t ne = The_element.n();
      size_t n1=The_element(1)->n(), n2=The_element(2)->n(),
             n3=The_element(3)->n(), n4=The_element(4)->n();
      Quad4 q(the_element);
      q.setLocal(Point<real_t>(0.,0.));
      q.atGauss(1,_sh,_dSh,_wg);
      eps(ne,1) = (*_u)(n1,1)*_dSh[0].x + (*_u)(n2,1)*_dSh[1].x +
                  (*_u)(n3,1)*_dSh[2].x + (*_u)(n4,1)*_dSh[3].x;
      eps(ne,2) = (*_u)(n1,2)*_dSh[0].y + (*_u)(n2,2)*_dSh[1].y +
                  (*_u)(n3,2)*_dSh[2].y + (*_u)(n4,2)*_dSh[3].y;
      eps(ne,3) = (*_u)(n1,1)*_dSh[0].y + (*_u)(n2,1)*_dSh[1].y +
                  (*_u)(n3,1)*_dSh[2].y + (*_u)(n4,1)*_dSh[3].y +
                  (*_u)(n1,2)*_dSh[0].x + (*_u)(n2,2)*_dSh[1].x +
                  (*_u)(n3,2)*_dSh[2].x + (*_u)(n4,2)*_dSh[3].x;
   }
}


void Elas2DQ4::Stress(Vect<real_t>& s,
                      Vect<real_t>& vm)
{
   vm.setSize(_nb_el);
   s.setSize(_nb_el,3);
   MESH_EL {
      size_t ne = The_element.n();
      size_t n1=The_element(1)->n(), n2=The_element(2)->n(),
             n3=The_element(3)->n(), n4=The_element(4)->n();
      Quad4 q(the_element);
      q.setLocal(Point<real_t>(0.,0.));
      q.atGauss(1,_sh,_dSh,_wg);
      real_t e, sx, sy, txy, delta;
      Point<real_t> Du = (*_u)(n1,1)*_dSh[0] + (*_u)(n2,1)*_dSh[1] + (*_u)(n3,1)*_dSh[2] + (*_u)(n4,1)*_dSh[3];
      Point<real_t> Dv = (*_u)(n1,2)*_dSh[0] + (*_u)(n2,2)*_dSh[1] + (*_u)(n3,2)*_dSh[2] + (*_u)(n4,2)*_dSh[3];
      e = (*_u)(n1,1)*_dSh[0].y + (*_u)(n2,1)*_dSh[0].y + (*_u)(n3,1)*_dSh[2].y + (*_u)(n4,1)*_dSh[3].y +
          (*_u)(n1,2)*_dSh[0].x + (*_u)(n2,2)*_dSh[1].x + (*_u)(n3,2)*_dSh[2].x + (*_u)(n4,2)*_dSh[3].x;
      sx = _E1*Du.x + _E2*Dv.y;
      sy = _E2*Dv.y + _E3*e;
      txy = _E6*e;
      delta = sqrt((sx+sy)*(sx+sy) + 4*(txy*txy-sx*sy));
      s(ne,1) = 0.5*(sx+sy-delta);
      s(ne,2) = 0.5*(sx+sy+delta);
      s(ne,3) = _poisson*(sx+sy);
      vm(ne) = sqrt(0.5*((s(ne,1)-s(ne,2))*(s(ne,1)-s(ne,2)) +
                         (s(ne,2)-s(ne,3))*(s(ne,2)-s(ne,3)) +
                         (s(ne,3)-s(ne,1))*(s(ne,3)-s(ne,1))));
   }
}


void Elas2DQ4::Stress(Vect<real_t>& sigma,
                      Vect<real_t>& s,
                      Vect<real_t>& st)
{
   sigma.setSize(_nb_el,3);
   s.setSize(_nb_el,3);
   st.setSize(_nb_el);
   MESH_EL {
      size_t ne = The_element.n();
      size_t n1=The_element(1)->n(), n2=The_element(2)->n(),
             n3=The_element(3)->n(), n4=The_element(4)->n();
      Quad4 q(the_element);
      q.setLocal(Point<real_t>(0.,0.));
      q.atGauss(1,_sh,_dSh,_wg);
      Point<real_t> Du = (*_u)(n1,1)*_dSh[0] + (*_u)(n2,1)*_dSh[1] + (*_u)(n3,1)*_dSh[2] + (*_u)(n4,1)*_dSh[3];
      Point<real_t> Dv = (*_u)(n1,2)*_dSh[0] + (*_u)(n2,2)*_dSh[1] + (*_u)(n3,2)*_dSh[2] + (*_u)(n4,2)*_dSh[3];
      real_t e = (*_u)(n1,1)*_dSh[0].y + (*_u)(n2,1)*_dSh[1].y +
                 (*_u)(n3,1)*_dSh[2].y + (*_u)(n4,1)*_dSh[3].y +
                 (*_u)(n1,2)*_dSh[0].x + (*_u)(n2,2)*_dSh[1].x +
                 (*_u)(n3,2)*_dSh[2].x + (*_u)(n4,2)*_dSh[3].x;
      real_t sx = _E1*Du.x+_E2*Dv.y, sy = _E2*Dv.y + _E3*e;
      real_t txy = _E6*e;
      real_t delta = sqrt((sx+sy)*(sx+sy) + 4*(txy*txy-sx*sy));
      s(ne,1) = 0.5*(sx+sy-delta);
      s(ne,2) = 0.5*(sx+sy+delta);
      s(ne,3) = _poisson*(sx+sy);
      st(ne) = sqrt(0.5*(s(ne,1)-s(ne,2))*(s(ne,1)-s(ne,2)) +
                        (s(ne,2)-s(ne,3))*(s(ne,2)-s(ne,3)) +
                        (s(ne,3)-s(ne,1))*(s(ne,3)-s(ne,1)));
      sigma(ne,1) = _lambda*(Du.x + Dv.y) + _G*Du.x;
      sigma(ne,2) = _lambda*(Du.x + Dv.y) + _G*Dv.y;
      sigma(ne,3) = 0.5*_G*(Du.y + Dv.x);
   }
}

} /* namespace OFELI */
