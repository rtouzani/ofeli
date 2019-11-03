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

                         Implementation of class HelmholtzBT3
                       for Helmholtz Equation in a Bounded Domain
                         using 3-node triangular finite element

  ==============================================================================*/


#include "equations/electromagnetics/HelmholtzBT3.h"
#include "linear_algebra/Vect_impl.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"

namespace OFELI {

HelmholtzBT3::HelmholtzBT3(Mesh& ms)
             : Equation<complex_t,3,3,2,2>(ms)
{
   _wave_nb = 1.;
   _equation_name = "Helmholtz equation in a bounded domain";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _theMesh->removeImposedDOF();
   setMatrixType(SKYLINE);
   setSolver(DIRECT_SOLVER);
}


HelmholtzBT3::HelmholtzBT3(Mesh&            ms,
                           Vect<complex_t>& u)
             : Equation<complex_t,3,3,2,2>(ms,u)
{
   _wave_nb = 1.;
   _equation_name = "Helmholtz equation in a bounded domain";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _theMesh->removeImposedDOF();
   setMatrixType(SKYLINE);
   setSolver(DIRECT_SOLVER);
}


HelmholtzBT3::~HelmholtzBT3()
{
}


void HelmholtzBT3::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   setMaterial();
   Triang3 tr(_theElement);
   _el_geo.area = tr.getArea();
   _dSh = tr.DSh();
   ElementNodeCoordinates();
   ElementNodeVector(*_u,_eu);
   eMat = complex_t(0);
   eRHS = complex_t(0);
}


void HelmholtzBT3::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Line2 ln(sd);
   SideNodeCoordinates();
   _el_geo.length = ln.getLength();
   sMat = complex_t(0);
   sRHS = complex_t(0);
}


void HelmholtzBT3::LHS()
{
   real_t c = _wave_nb*_wave_nb*_el_geo.area/12.0;
   for (size_t i=1; i<=3; i++) {
      Point<real_t> a = _el_geo.area*_dSh[i-1];
      for (size_t j=1; j<=3; j++)
         eMat(i,j) += (a,_dSh[j-1]) - c;
      eMat(i,i) -= c;
   }
   eA0 = eMat;
}


void HelmholtzBT3::BodyRHS(Vect<complex_t>& f)
{
   real_t c = OFELI_THIRD*_el_geo.area;
   for (size_t i=1; i<=3; i++)
      eRHS(i) += c*f((*_theElement)(i)->n());
}


void HelmholtzBT3::BoundaryRHS(Vect<complex_t>& f)
{
   if (_theSide->getCode(1)>0) {
      sRHS(1) += 0.5*_el_geo.length*f((*_theSide)(1)->n());
      sRHS(2) += 0.5*_el_geo.length*f((*_theSide)(2)->n());
   }
}


void HelmholtzBT3::build()
{
   mesh_elements(*_theMesh) {
      set(the_element);
      LHS();
      AbsEqua<complex_t>::_A->Assembly(The_element,eMat.get());
      if (AbsEqua<complex_t>::_bf!=nullptr)
         BodyRHS(*AbsEqua<complex_t>::_bf);
      if (AbsEqua<complex_t>::_bc!=nullptr)
         this->updateBC(*_theElement,*AbsEqua<complex_t>::_bc);
      AbsEqua<complex_t>::_b->Assembly(The_element,eRHS.get());
   }
   mesh_sides(*_theMesh) {
      set(the_side);
      if (AbsEqua<complex_t>::_sf!=nullptr)
         BoundaryRHS(*_sf);
      AbsEqua<complex_t>::_b->Assembly(The_side,sRHS.get());
   }
}

} /* namespace OFELI */
