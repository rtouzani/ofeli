/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2026 Rachid Touzani

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
#include "shape_functions/Tetra4.h"
#include "shape_functions/Triang3.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/LocalMatrix_impl.h"
#include "linear_algebra/LocalVect_impl.h"
#include "equations/Equation_impl.h"


namespace OFELI {

DC3DT4::DC3DT4()
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "3-D, 4-Node Tetrahedra (P1)";
}


DC3DT4::DC3DT4(Mesh& ms) 
       : Equation<4,4,3,3>(ms)
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "3-D, 4-Node Tetrahedra (P1)";
   setMatrixType(SPARSE);
   setSolver(CG_SOLVER,DILU_PREC);
}


DC3DT4::DC3DT4(Mesh&         ms,
               Vect<real_t>& u)
       : Equation<4,4,3,3>(ms,u)
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "3-D, 4-Node Tetrahedra (P1)";
   setMatrixType(SPARSE);
   setSolver(CG_SOLVER,DILU_PREC);
}


DC3DT4::~DC3DT4() { }


void DC3DT4::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   setMaterial();
   Tetra4 t(_theElement);
   _el_geo.volume = t.getVolume();
   _el_geo.det = t.getDet();
   _el_geo.center = t.getCenter();
   _dSh = t.DSh();
   ElementNodeCoordinates();
   ElementNodeVector(*_u,_eu);
   if (_rho_set)
      _rho = _rho_fct(_el_geo.center,_TimeInt.time);
   if (_Cp_set)
      _cp = _Cp_fct(_el_geo.center,_TimeInt.time);
   if (_kappa_set)
      _diff = _kappa_fct(_el_geo.center,_TimeInt.time);
   eA0 = 0, eA1 = 0, eMat = 0;
   eRHS = 0;
}


void DC3DT4::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Triang3 t(sd);
   _el_geo.area = t.getArea();
   SideNodeCoordinates();
   sA0 = 0;
   sRHS = 0;
}


void DC3DT4::LCapacity(real_t coef)
{
   real_t c = coef*0.25*_el_geo.volume*_rho*_cp;
   for (size_t i=1; i<=4; i++)
      eA1(i,i) += c;
}


void DC3DT4::Capacity(real_t coef)
{
   real_t c = 0.1*_el_geo.volume*_rho*_cp*coef;
   real_t d = 0.5*c;
   eA1(1,1) += c; eA1(2,2) += c; eA1(3,3) += c; eA1(4,4) += c;
   eA1(1,2) += d; eA1(2,1) += d; eA1(1,3) += d; eA1(1,4) += d;
   eA1(3,1) += d; eA1(2,3) += d; eA1(3,2) += d; eA1(2,4) += d;
   eA1(4,1) += d; eA1(4,2) += d; eA1(4,3) += d; eA1(3,4) += d;
}


void DC3DT4::Diffusion(real_t coef)
{
   real_t c = coef*_diff*_el_geo.volume;
   for (size_t i=1; i<=4; i++)
      for (size_t j=1; j<=4; j++)
         eA0(i,j) += c*(_dSh[i-1],_dSh[j-1]);
   eMat += eA0;
}


void DC3DT4::Diffusion(const DMatrix<real_t>& diff,
                       real_t                 coef)
{
   real_t c = coef*_el_geo.volume;
   for (size_t i=1; i<=4; i++)
       for (size_t j=1; j<=4; j++)
          eA0(i,j) += c*(diff(1,1)*_dSh[i-1].x*_dSh[j-1].x + diff(1,2)*_dSh[i-1].y*_dSh[j-1].x
                       + diff(1,3)*_dSh[i-1].z*_dSh[j-1].x + diff(2,1)*_dSh[i-1].x*_dSh[j-1].y
                       + diff(2,2)*_dSh[i-1].y*_dSh[j-1].y + diff(2,3)*_dSh[i-1].z*_dSh[j-1].y
                       + diff(3,1)*_dSh[i-1].x*_dSh[j-1].z + diff(3,2)*_dSh[i-1].y*_dSh[j-1].z
                       + diff(3,3)*_dSh[i-1].z*_dSh[j-1].z);
   eMat += eA0;
}


void DC3DT4::Convection(const Point<real_t>& v,
                        real_t               coef)
{
   for (size_t i=1; i<=4; i++)
      for (size_t j=1; j<=4; j++)
         eA0(i,j) += coef*0.25*_el_geo.volume*(v,_dSh[i-1]);
   eMat += eA0;
}


void DC3DT4::Convection(const Vect<Point<real_t> >& v,
                        real_t                      coef)
{
   size_t i;
   LocalMatrix<real_t,4,3> ve;
   for (i=1; i<=4; i++) {
      size_t n = (*_theElement)(i)->n();
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
         eA0(i,j) += coef*0.25*_el_geo.volume*(w,_dSh[i-1]);
   eMat += eA0;
}


void DC3DT4::Convection(real_t coef)
{
   LocalMatrix<real_t,4,3> ve;
   size_t i;
   for (i=1; i<=4; i++) {
      size_t n = (*_theElement)(i)->n();
      ve(i,1) = (*_vel)(3*n-2);
      ve(i,2) = (*_vel)(3*n-1);
      ve(i,3) = (*_vel)(3*n  );
   }
   Point<real_t> v;
   v.x = 0.25*(ve(1,1) + ve(2,1) + ve(3,1) + ve(4,1));
   v.y = 0.25*(ve(1,2) + ve(2,2) + ve(3,2) + ve(4,2));
   v.z = 0.25*(ve(1,3) + ve(2,3) + ve(3,3) + ve(4,3));
   real_t c = 0.25*_el_geo.volume*coef;
   for (i=1; i<=4; i++)
      for (size_t j=1; j<=4; j++)
         eA0(i,j) += c*(v,_dSh[j-1]);
   eMat += eA0;
}


void DC3DT4::BodyRHS(const Vect<real_t>& f)
{
   for (size_t i=1; i<=4; ++i)
      eRHS(i) += f((*_theElement)(i)->n())*0.25*_el_geo.volume;
}


void DC3DT4::BoundaryRHS(real_t flux)
{
   real_t c = flux*_el_geo.area*OFELI_THIRD;
   sRHS(1) += c;
   sRHS(2) += c;
   sRHS(3) += c;
}


void DC3DT4::BoundaryRHS(const Vect<real_t>& f)
{
   if (_theSide->getCode(1)>0) {
      real_t c = _el_geo.area*OFELI_THIRD;
      if (f.getDOFType()==NODE_DOF) {
         for (size_t i=1; i<=3; ++i)
            sRHS(i) += c*f((*_theSide)(i)->n());
      }
      else if (f.getDOFType()==SIDE_DOF) {
         real_t ff = c*f(_theSide->n());
         for (size_t i=1; i<=3; ++i)
            sRHS(i) += ff;
      }
   }
}


Point<real_t> DC3DT4::Flux() const
{
   Point<real_t> f;
   f = _diff*(_eu(1)*_dSh[0] + _eu(2)*_dSh[1] +
              _eu(3)*_dSh[2] + _eu(4)*_dSh[3]);
   return f;
}


void DC3DT4::Grad(Vect<Point<real_t> >& g)
{
   MESH_EL {
      set(the_element);
      g(element_label) = _eu(1)*_dSh[0] + _eu(2)*_dSh[1] + _eu(3)*_dSh[2] + _eu(4)*_dSh[3];
   }
}

void DC3DT4::Periodic(real_t coef)
{
   for (size_t i=1; i<=3; i++) {
      real_t c = OFELI_THIRD*_el_geo.area*coef;
      if (_theSide->getCode(1) == PERIODIC_A)
         eA0(i,i) += c;
      else if (_theSide->getCode(1) == PERIODIC_B)
         eA0(i,i) -= c;
   }
}

} /* namespace OFELI */
