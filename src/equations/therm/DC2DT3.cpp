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

                     Class DC2DT3: Diffusion-Convection Element
               using 3-Node Triangular Finite element in two dimensions

  ==============================================================================*/


#include "equations/therm/DC2DT3.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"

namespace OFELI {


DC2DT3::DC2DT3()
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _stab = false;
}


DC2DT3::DC2DT3(const Element* el)
{
   set(el);
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _stab = false;
}


DC2DT3::DC2DT3(const Side* sd)
{
   set(sd);
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _stab = false;
}


DC2DT3::DC2DT3(const Element*      el,
               const Vect<real_t>& u,
               real_t              time)
       : Equation<real_t,3,3,2,2>(el,u,time)
{
   set(el);
   ElementVector(u);
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _stab = false;
}


DC2DT3::DC2DT3(const Element*      el,
               const Vect<real_t>& u,
               real_t              time,
               real_t              deltat,
               int                 scheme)
       : Equation<real_t,3,3,2,2>(el,u,time)
{
   set(el);
   _time_step = deltat;
   ElementVector(u);
   setTimeIntegration(scheme);
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _stab = false;
}


DC2DT3::DC2DT3(const Side*         sd,
               const Vect<real_t>& u,
               real_t              time)
       : Equation<real_t,3,3,2,2>(sd,u,time)
{
   set(sd);
   SideVector(u);
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _stab = false;
}


DC2DT3::DC2DT3(const Side*         sd,
               const Vect<real_t>& u,
               real_t              time,
               real_t              deltat,
               int                 scheme)
       : Equation<real_t,3,3,2,2>(sd,u,time)
{
   set(sd);
   SideVector(u);
   _time_step = deltat;
   setTimeIntegration(scheme);
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _stab = false;
}


DC2DT3::DC2DT3(Mesh& ms) 
       : Equation<real_t,3,3,2,2>(ms)
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _stab = false;
   setMatrixType(SPARSE);
   setSolver(BICG_STAB_SOLVER,DILU_PREC);
}


DC2DT3::DC2DT3(Mesh&         ms,
               Vect<real_t>& u)
       : Equation<real_t,3,3,2,2>(ms)
{
   setInput(SOLUTION,u);
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _stab = false;
}


DC2DT3::~DC2DT3() { }


void DC2DT3::set(const Element* el)
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
   eA0 = 0, eA1 = 0, eMat = 0;
   eRHS = 0;
}


void DC2DT3::set(const Side* sd)
{
   _nb_dof = 1;
   Init(sd);
   Line2 ln(sd);
   SideNodeCoordinates();
   _center = ln.getCenter();
   _length = ln.getLength();
   sMat = 0; sRHS = 0;
}


void DC2DT3::setInput(EqDataType    opt,
                      Vect<real_t>& u)
{
   AbsEqua<real_t>::setInput(opt,u);
   if (opt==VELOCITY_FIELD)
      _vel = &u;
}


void DC2DT3::LCapacityToLHS(real_t coef)
{
   real_t c=OFELI_THIRD*_area*_rhocp*coef;
   eA1(1,1) += c;
   eA1(2,2) += c;
   eA1(3,3) += c;
   eMat += eA1;
}


void DC2DT3::LCapacityToRHS(real_t coef)
{
   real_t c=OFELI_THIRD*_area*_rhocp*coef;
   eRHS(1) += c*ePrev(1);
   eRHS(2) += c*ePrev(2);
   eRHS(3) += c*ePrev(3);
}


void DC2DT3::CapacityToLHS(real_t coef)
{
   real_t c=OFELI_SIXTH*_area*_rhocp*coef;
   real_t d=0.5*c;
   eA1(1,1) += c; eA1(2,2) += c; eA1(3,3) += c;
   eA1(1,2) += d; eA1(2,1) += d; eA1(1,3) += d;
   eA1(3,1) += d; eA1(2,3) += d; eA1(3,2) += d;
}


void DC2DT3::CapacityToRHS(real_t coef)
{
   real_t c=OFELI_SIXTH*_area*_rhocp*coef;
   real_t d=0.5*c;
   eRHS(1) += c*ePrev(1) + d*(ePrev(2) + ePrev(3));
   eRHS(2) += c*ePrev(2) + d*(ePrev(1) + ePrev(3));
   eRHS(3) += c*ePrev(3) + d*(ePrev(1) + ePrev(2));
}


void DC2DT3::Diffusion(real_t coef)
{
   for (size_t i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eA0(i,j) += coef*_diff*_area*(_dSh[i-1],_dSh[j-1]);
   eMat += eA0;
}


void DC2DT3::Diffusion(const LocalMatrix<real_t,2,2>& diff,
                       real_t                         coef)
{
   real_t c=_area*coef;
   for (size_t i=0; i<3; i++)
      for (size_t j=0; j<3; j++)
         eA0(i+1,j+1) += c*(diff(1,1)*_dSh[i].x*_dSh[j].x
                          + diff(2,2)*_dSh[i].y*_dSh[j].y
                          + diff(1,2)*_dSh[i].y*_dSh[j].x
                          + diff(2,1)*_dSh[i].x*_dSh[j].y);
   eMat += eA0;
}


void DC2DT3::DiffusionToRHS(real_t coef)
{
   Point<real_t> u;
   for (size_t i=1; i<=3; i++)
      u += _dSh(i)*ePrev(i);
   for (size_t j=0; j<3; j++)
      eRHS(j+1) -= coef*_diff*_area*(_dSh[j]*u);
}


void DC2DT3::Convection(const Point<real_t>& v,
                        real_t               coef)
{
   LocalVect<real_t,3> dd;
   for (size_t i=1; i<=3; i++)
      dd(i) = v*_dSh(i);
   for (size_t i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eA0(i,j) += OFELI_THIRD*coef*_area*dd(j);
   if (_stab) {
      real_t c=coef*_area*_h/Abs(v), d;
      for (size_t i=1; i<=3; i++) {
         d = c*dd(i);
         for (size_t j=1; j<=3; j++)
            eA0(i,j) += d*dd(j);
      }
   }
   eMat += eA0;
}


void DC2DT3::Convection(real_t coef)
{
   LocalMatrix<real_t,3,2> ve;
   for (size_t i=1; i<=3; i++) {
      size_t n=_theElement->getNodeLabel(i);
      ve(i,1) = (*_vel)(2*n-1);
      ve(i,2) = (*_vel)(2*n  );
   }
   Point<real_t> v3(ve(1,1)+ve(2,1)+ve(3,1),ve(1,2)+ve(2,2)+ve(3,2));
   real_t c=OFELI_THIRD*OFELI_THIRD*_area*coef;
   for (size_t i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eA0(i,j) += c*(v3*_dSh(j));
   if (_stab) {
      real_t c=coef*_area*_h/Abs(v3), d;
      for (size_t i=1; i<=3; i++) {
         d = c*(v3*_dSh(i));
         for (size_t j=1; j<=3; j++)
            eA0(i,j) += d*(v3*_dSh(j));
      }
   }
   eMat += eA0;
}


void DC2DT3::Convection(const Vect<real_t>& v,
                        real_t              coef)
{
   size_t i;
   LocalMatrix<real_t,3,2> ve;
   for (i=1; i<=3; i++) {
      size_t n = _theElement->getNodeLabel(i);
      ve(i,1) = v(2*n-1);
      ve(i,2) = v(2*n  );
   }
   Point<real_t> vv(ve(1,1)+ve(2,1)+ve(3,1),ve(1,2)+ve(2,2)+ve(3,2));
   real_t c=_area*coef/9.0;
   for (i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eA0(i,j) += c*(vv*_dSh(j));
   eMat += eA0;
}


void DC2DT3::ConvectionToRHS(const Point<real_t>& v,
                             real_t               coef)
{
   Point<real_t> u;
   for (size_t i=1; i<=3; i++)
      u += _dSh(i) * ePrev(i);
   for (size_t j=1; j<=3; j++)
      eRHS(j) -= OFELI_THIRD*_area*coef*(v*u);
   if (_stab) {
      real_t c=coef*_h/Abs(v)*_area*(v*u);
      for (size_t i=1; i<=3; i++)
         eRHS(i) -= c*(v*_dSh(i));
   }
}


void DC2DT3::ConvectionToRHS(real_t coef)
{
   LocalVect<Point<real_t>,3> ve;
   for (size_t i=1; i<=3; i++) {
      size_t n=_theElement->getNodeLabel(i);
      ve(i).x = (*_vel)(2*n-1);
      ve(i).y = (*_vel)(2*n  );
   }
   Point<real_t> vv(OFELI_THIRD*(ve(1)+ve(2)+ve(3)));
   Point<real_t> du = _dSh(1)*ePrev(1) + _dSh(2)*ePrev(2) + _dSh(3)*ePrev(3);
   real_t c = OFELI_SIXTH*coef*_area;
   eRHS(1) -= c*(2*ve(1)+ve(2)+ve(3))*du;
   eRHS(2) -= c*(ve(1)+2*ve(2)+ve(3))*du;
   eRHS(3) -= c*(ve(1)+ve(2)+2*ve(3))*du;
   if (_stab) {
      c = 0.5*coef*_h/Abs(vv)*_area*(vv*du);
      eRHS(1) -= c*(vv*_dSh(1));
      eRHS(2) -= c*(vv*_dSh(2));
      eRHS(3) -= c*(vv*_dSh(3));
   }
}


void DC2DT3::LinearExchange(real_t coef,
                            real_t T)
{
   sMat(1,1) += 0.5*_length*coef;
   sMat(2,2) += 0.5*_length*coef;
   sRHS(1)   += 0.5*_length*coef*T;
   sRHS(2)   += 0.5*_length*coef*T;
}


void DC2DT3::BodyRHS(UserData<real_t>& ud,
                     real_t            coef)
{
   for (size_t i=1; i<=3; i++) {
      real_t f=ud.BodyForce(_x[i-1],_time);
      eRHS(i) += coef*OFELI_THIRD*_area*f;
   }
}


void DC2DT3::BodyRHS(real_t bf)
{
   real_t c=OFELI_THIRD*_area*bf;
   eRHS(1) += c;
   eRHS(2) += c;
   eRHS(3) += c;
}


void DC2DT3::BodyRHS(const Vect<real_t>& bf,
                     int                 opt)
{
   real_t c = OFELI_THIRD*_area;
   if (opt==LOCAL_ARRAY) {
      for (size_t i=1; i<=3; i++)
         eRHS(i) += c*bf(i);
   }
   else
      for (size_t i=1; i<=3; i++)
         eRHS(i) += c*bf((*_theElement)(i)->n());
}


void DC2DT3::BoundaryRHS(UserData<real_t>& ud,
                         real_t            coef)
{
   real_t c = 0.5*_length*coef;
   for (size_t i=1; i<=2; i++) {
      real_t f = ud.SurfaceForce((*_theSide)(i)->getCoord(),_theSide->getCode(1),_time);
      sRHS(i) += c*f;
   }
}


void DC2DT3::BoundaryRHS(real_t flux)
{
   real_t c = 0.5*_length*flux;
   for (size_t i=1; i<=2; i++)
       sRHS(i) += c;
}


void DC2DT3::BoundaryRHS(const Vect<real_t>& sf,
                               int           opt)
{
   real_t c = 0.5*_length;
   if (opt==LOCAL_ARRAY)
      for (size_t i=1; i<=2; i++)
         sRHS(i) += c*sf(i);
   else
      for (size_t i=1; i<=2; i++)
         sRHS(i) += c*sf((*_theSide)(i)->n());
}


Point<real_t> &DC2DT3::Flux() const
{
   _f = _diff*(ePrev(1)*_dSh[0] + ePrev(2)*_dSh[1] + ePrev(3)*_dSh[2]);
   return _f;
}


real_t DC2DT3::Energy(const Vect<real_t>& u)
{
   Point<real_t> du = u(1)*_dSh[0] + u(2)*_dSh[1] + u(3)*_dSh[2];
   return 0.5*_diff*_area*(du*du);
}


real_t DC2DT3::Energy(const Vect<real_t>& u,
                      UserData<real_t>&   ud)
{
   Point<real_t> du = u(1)*_dSh[0] + u(2)*_dSh[1] + u(3)*_dSh[2];
   real_t W = 0.5*_diff*_area*(du*du);
   for (size_t i=1; i<=3; i++)
      W -= OFELI_THIRD*_area*u(i)*ud.BodyForce((*_theElement)(i)->getCoord(),_time);
   return W;
}

void DC2DT3::EnerGrad(const Vect<real_t>& u,
                      UserData<real_t>&   ud)
{
   LocalMatrix<real_t,3,3> a;
   LocalVect<real_t,3> b;
   real_t c = _area*_diff;
   a(1,1) = c*(_dSh[0]*_dSh[0]);
   a(2,2) = c*(_dSh[1]*_dSh[1]);
   a(3,3) = c*(_dSh[2]*_dSh[2]);
   a(1,2) = c*(_dSh[0]*_dSh[1]);
   a(2,1) = a(1,2);
   a(1,3) = c*(_dSh[0]*_dSh[2]);
   a(3,1) = a(1,3);
   a(2,3) = c*(_dSh[1]*_dSh[2]);
   a(3,2) = a(2,3);
   b(1) = OFELI_THIRD*_area*ud.BodyForce(_x[0],_time);
   b(2) = OFELI_THIRD*_area*ud.BodyForce(_x[1],_time);
   b(3) = OFELI_THIRD*_area*ud.BodyForce(_x[2],_time);
   eRHS(1) = a(1,1)*u(1) + a(1,2)*u(2) + a(1,3)*u(3) - b(1);
   eRHS(2) = a(2,1)*u(1) + a(2,2)*u(2) + a(2,3)*u(3) - b(2);
   eRHS(3) = a(3,1)*u(1) + a(3,2)*u(2) + a(3,3)*u(3) - b(3);
}


void DC2DT3::EnerGrad(const Vect<real_t>& u,
                      const Vect<real_t>& f)
{
   LocalMatrix<real_t,3,3> a;
   LocalVect<real_t,3> b;
   real_t c = _area*_diff;
   a(1,1) = c*(_dSh[0]*_dSh[0]);
   a(2,2) = c*(_dSh[1]*_dSh[1]);
   a(3,3) = c*(_dSh[2]*_dSh[2]);
   a(1,2) = c*(_dSh[0]*_dSh[1]);
   a(2,1) = a(1,2);
   a(1,3) = c*(_dSh[0]*_dSh[2]);
   a(3,1) = a(1,3);
   a(2,3) = c*(_dSh[1]*_dSh[2]);
   a(3,2) = a(2,3);
   b(1) = OFELI_THIRD*_area*f(TheElement(1)->n());
   b(2) = OFELI_THIRD*_area*f(TheElement(2)->n());
   b(3) = OFELI_THIRD*_area*f(TheElement(3)->n());
   eRHS(1) = a(1,1)*u(1) + a(1,2)*u(2) + a(1,3)*u(3) - b(1);
   eRHS(2) = a(2,1)*u(1) + a(2,2)*u(2) + a(2,3)*u(3) - b(2);
   eRHS(3) = a(3,1)*u(1) + a(3,2)*u(2) + a(3,3)*u(3) - b(3);
}


void DC2DT3::EnerGrad(const Vect<real_t>& u)
{
   LocalMatrix<real_t,3,3> a;
   real_t c = _area*_diff;
   a(1,1) = c*(_dSh[0]*_dSh[0]);
   a(2,2) = c*(_dSh[1]*_dSh[1]);
   a(3,3) = c*(_dSh[2]*_dSh[2]);
   a(1,2) = c*(_dSh[0]*_dSh[1]);
   a(2,1) = a(1,2);
   a(1,3) = c*(_dSh[0]*_dSh[2]);
   a(3,1) = a(1,3);
   a(2,3) = c*(_dSh[1]*_dSh[2]);
   a(3,2) = a(2,3);
   eRHS(1) = a(1,1)*u(1) + a(1,2)*u(2) + a(1,3)*u(3);
   eRHS(2) = a(2,1)*u(1) + a(2,2)*u(2) + a(2,3)*u(3);
   eRHS(3) = a(3,1)*u(1) + a(3,2)*u(2) + a(3,3)*u(3);
}


Point<real_t>& DC2DT3::Grad(const LocalVect<real_t,3>& u) const
{
   _f = u[0]*_dSh[0] + u[1]*_dSh[1] + u[2]*_dSh[2];
   return _f;
}


Point<real_t>& DC2DT3::Grad(const Vect<real_t>& u) const
{
   LocalVect<real_t,3> ue(_theElement,u);
   _f = Grad(ue);
   return _f;
}


void DC2DT3::Periodic(real_t coef)
{
   for (size_t i=1; i<=2; i++) {
      real_t c=0.5*_length*coef;
      if (_theSide->getCode(1) == PERIODIC_A)
         eA0(i,i) += c;
      else if (_theSide->getCode(1) == PERIODIC_B)
         eA0(i,i) -= c;
   }
}

} /* namespace OFELI */
