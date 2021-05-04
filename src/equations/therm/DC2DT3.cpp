/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

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
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/LocalMatrix_impl.h"
#include "linear_algebra/LocalVect_impl.h"
#include "equations/Equa_impl.h"
#include "equations/Equation_impl.h"

namespace OFELI {


DC2DT3::DC2DT3()
{
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
   setSolver(CG_SOLVER,DILU_PREC);
}


DC2DT3::DC2DT3(Mesh&         ms,
               Vect<real_t>& u)
       : Equation<real_t,3,3,2,2>(ms,u)
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _stab = false;
   setMatrixType(SPARSE);
   setSolver(CG_SOLVER,DILU_PREC);
}


DC2DT3::~DC2DT3() { }


void DC2DT3::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   setMaterial();
   Triang3 tr(_theElement);
   _el_geo.area = tr.getArea();
   _el_geo.center = tr.getCenter();
   _dSh = tr.DSh();
   _el_geo.size = 2*tr.getCircumRadius();
   ElementNodeCoordinates();
   ElementNodeVector(*_u,_eu);
   if (_rho_set)
      _rho = _rho_fct(_el_geo.center,_TimeInt.time);
   if (_Cp_set)
      _cp = _Cp_fct(_el_geo.center,_TimeInt.time);
   if (_kappa_set)
      _diff = _kappa_fct(_el_geo.center,_TimeInt.time);
   eA0 = 0, eA1 = 0;
   eRHS = 0;
}


void DC2DT3::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Line2 ln(sd);
   SideNodeCoordinates();
   _el_geo.center = ln.getCenter();
   _el_geo.length = ln.getLength();
   sMat = 0; sRHS = 0;
}


void DC2DT3::setInput(EqDataType    opt,
                      Vect<real_t>& u)
{
   Equa<real_t>::setInput(opt,u);
   if (opt==VELOCITY_FIELD)
      _vel = &u;
}


void DC2DT3::LCapacity(real_t coef)
{
   real_t c=OFELI_THIRD*_el_geo.area*_rho*_cp*coef;
   eA1(1,1) += c;
   eA1(2,2) += c;
   eA1(3,3) += c;
   eMat += eA1;
}


void DC2DT3::Capacity(real_t coef)
{
   real_t c=OFELI_SIXTH*_el_geo.area*_rho*_cp*coef;
   real_t d=0.5*c;
   eA1(1,1) += c; eA1(2,2) += c; eA1(3,3) += c;
   eA1(1,2) += d; eA1(2,1) += d; eA1(1,3) += d;
   eA1(3,1) += d; eA1(2,3) += d; eA1(3,2) += d;
}


void DC2DT3::Diffusion(real_t coef)
{
   for (size_t i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eA0(i,j) += coef*_diff*_el_geo.area*(_dSh[i-1],_dSh[j-1]);
   eMat += eA0;
}


void DC2DT3::Diffusion(const LocalMatrix<real_t,2,2>& diff,
                       real_t                         coef)
{
   real_t c=_el_geo.area*coef;
   for (size_t i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eA0(i,j) += c*(diff(1,1)*_dSh[i-1].x*_dSh[j-1].x
                      + diff(2,2)*_dSh[i-1].y*_dSh[j-1].y
                      + diff(1,2)*_dSh[i-1].y*_dSh[j-1].x
                      + diff(2,1)*_dSh[i-1].x*_dSh[j-1].y);
   eMat += eA0;
}


void DC2DT3::Convection(const Point<real_t>& v,
                        real_t               coef)
{
   LocalVect<real_t,3> dd;
   for (size_t i=1; i<=3; i++)
      dd(i) = (v,_dSh[i-1]);
   for (size_t i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eA0(i,j) += OFELI_THIRD*coef*_el_geo.area*dd(j);
   if (_stab) {
      real_t c=coef*_el_geo.area*_el_geo.size/Abs(v), d;
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
      size_t n=(*_theElement)(i)->n();
      ve(i,1) = (*_vel)(2*n-1);
      ve(i,2) = (*_vel)(2*n  );
   }
   Point<real_t> v3(ve(1,1)+ve(2,1)+ve(3,1),ve(1,2)+ve(2,2)+ve(3,2));
   real_t c=OFELI_THIRD*OFELI_THIRD*_el_geo.area*coef;
   for (size_t i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eA0(i,j) += c*(v3,_dSh[j-1]);
   if (_stab) {
      real_t c=coef*_el_geo.area*_el_geo.size/Abs(v3), d;
      for (size_t i=1; i<=3; i++) {
         d = c*(v3,_dSh[i-1]);
         for (size_t j=1; j<=3; j++)
            eA0(i,j) += d*(v3,_dSh[j-1]);
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
      size_t n = (*_theElement)(i)->n();
      ve(i,1) = v(2*n-1);
      ve(i,2) = v(2*n  );
   }
   Point<real_t> vv(ve(1,1)+ve(2,1)+ve(3,1),ve(1,2)+ve(2,2)+ve(3,2));
   real_t c=_el_geo.area*coef/9.0;
   for (i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eA0(i,j) += c*(vv,_dSh[j-1]);
   eMat += eA0;
}


void DC2DT3::LinearExchange(real_t coef,
                            real_t T)
{
   sMat(1,1) += 0.5*_el_geo.length*coef;
   sMat(2,2) += 0.5*_el_geo.length*coef;
   sRHS(1)   += 0.5*_el_geo.length*coef*T;
   sRHS(2)   += 0.5*_el_geo.length*coef*T;
}


void DC2DT3::BodyRHS(real_t bf)
{
   real_t c=OFELI_THIRD*_el_geo.area*bf;
   eRHS(1) += c;
   eRHS(2) += c;
   eRHS(3) += c;
}


void DC2DT3::BodyRHS(const Vect<real_t>& f)
{
   real_t c = OFELI_THIRD*_el_geo.area;
   for (size_t i=1; i<=3; i++)
      eRHS(i) += c*f((*_theElement)(i)->n());
}


void DC2DT3::BoundaryRHS(real_t flux)
{
   real_t c = 0.5*_el_geo.length*flux;
   for (size_t i=1; i<=2; i++)
       sRHS(i) += c;
}


void DC2DT3::BoundaryRHS(const Vect<real_t>& f)
{
   real_t c = 0.5*_el_geo.length;
   for (size_t i=1; i<=2; i++)
      sRHS(i) += c*f((*_theSide)(i)->n());
}


Point<real_t> &DC2DT3::Flux() const
{
   _f = _diff*(_eu(1)*_dSh[0] + _eu(2)*_dSh[1] + _eu(3)*_dSh[2]);
   return _f;
}


real_t DC2DT3::Energy(Vect<real_t>& u)
{
   _u = &u;
   real_t E = 0.;
   MESH_EL {
      set(the_element);
      Point<real_t> du = _eu[0]*_dSh[0] + _eu[1]*_dSh[1] + _eu[2]*_dSh[2];
      E += 0.5*_diff*_el_geo.area*(du,du);
      if (_bf!=nullptr)
         E -= OFELI_THIRD*_el_geo.area*(_eu(1)*(*_bf)(The_element(1)->n()) +
                                        _eu(2)*(*_bf)(The_element(2)->n()) +
                                        _eu(3)*(*_bf)(The_element(3)->n())); 
   }
   return E;
}


void DC2DT3::EnergyGrad(Vect<real_t>& u,
                        Vect<real_t>& g)
{
   real_t f = 0.;
   g.clear();
   _u = &u;
   MESH_EL {
      set(the_element);
      Diffusion();
      for (size_t i=1; i<=3; ++i) {
         if (_bf!=nullptr)
            f = OFELI_THIRD*_el_geo.area*(*_bf)(The_element(i)->n());
         g(The_element(i)->n()) += eA0(i,1)*_eu(1) + eA0(i,2)*_eu(2) + eA0(i,3)*_eu(3) - f;
      }
   }
}


void DC2DT3::Grad(Vect<Point<real_t> >& g)
{
   MESH_EL {
      set(the_element);
      g(element_label) = _eu(1)*_dSh[0] + _eu(2)*_dSh[1] + _eu(3)*_dSh[2];
   }
}


void DC2DT3::Periodic(real_t coef)
{
   for (size_t i=1; i<=2; i++) {
      real_t c=0.5*_el_geo.length*coef;
      if (_theSide->getCode(1) == PERIODIC_A)
         eA0(i,i) += c;
      else if (_theSide->getCode(1) == PERIODIC_B)
         eA0(i,i) -= c;
   }
}

} /* namespace OFELI */
