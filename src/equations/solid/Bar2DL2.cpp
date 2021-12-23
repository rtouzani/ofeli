/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani

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

                         Implementation of class Bar2DL2

  ==============================================================================*/


#include "equations/solid/Bar2DL2.h"
#include "shape_functions/Line2.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

Bar2DL2::Bar2DL2(Mesh& ms) 
        : Equation<2,4,1,2>(ms)
{
   _equation_name = "Truss";
   _finite_element = "2-D, 2-Node Bar (P1)";
   setMatrixType(SKYLINE|SYMMETRIC);
   setSolver(DIRECT_SOLVER);
   _section = 1.;
}


Bar2DL2::Bar2DL2(Mesh&         ms,
                 Vect<real_t>& u)
        : Equation<2,4,1,2>(ms,u)
{
   _equation_name = "Truss";
   _finite_element = "2-D, 2-Node Bar (P1)";
   setMatrixType(SKYLINE|SYMMETRIC);
   setSolver(DIRECT_SOLVER);
   _section = 1.;
}


void Bar2DL2::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   setMaterial();
   Line2 ln(_theElement);
   _el_geo.length = ln.getLength();
   _el_geo.center = ln.getCenter();
   ElementNodeCoordinates();
   ElementVector(*_u);
   real_t alpha = atan((_x[1].y-_x[0].y)/(_x[1].x-_x[0].x));
   real_t c = cos(alpha), s = sin(alpha);
   _cc = c*c;
   _ss = s*s;
   _sc = c*s;
   if (_young_set)
      _young = _young_fct(_el_geo.center,_TimeInt.time);
   eA0 = 0, eA1 = 0, eA2 = 0;
   eRHS = 0;
}


void Bar2DL2::setSection(real_t A)
{
   _section = A;
}


void Bar2DL2::LMass(real_t coef)
{
   real_t c = 0.5*_rho*coef*_el_geo.length;
   eA2(1,1) += c; eA2(2,2) += c;
   eA2(3,3) += c; eA2(4,4) += c;
}


void Bar2DL2::Mass(real_t coef)
{
   real_t c = coef*_rho*_el_geo.length*OFELI_SIXTH;
   real_t cc = c*_cc, ss = c*_ss, sc = c*_sc;
   eA2(1,1) += 2*cc; eA2(1,2) += 2*sc; eA2(1,3) +=   cc; eA2(1,4) +=   sc;
   eA2(2,1) += 2*sc; eA2(2,2) += 2*ss; eA2(2,3) +=   sc; eA2(2,4) +=   ss;
   eA2(3,1) +=   cc; eA2(3,2) +=   sc; eA2(3,3) += 2*cc; eA2(3,4) += 2*sc;
   eA2(4,1) +=   sc; eA2(4,2) +=   ss; eA2(4,3) += 2*sc; eA2(4,4) += 2*ss;
}


void Bar2DL2::Stiffness(real_t coef)
{
   real_t c = _young*_section*coef/_el_geo.length;
   real_t cc = c*_cc, ss = c*_ss, sc = c*_sc;
   eA0(1,1) += cc; eA0(1,2) += sc; eA0(1,3) -= cc; eA0(1,4) -= sc;
   eA0(2,1) += sc; eA0(2,2) += ss; eA0(2,3) -= sc; eA0(2,4) -= ss;
   eA0(3,1) -= cc; eA0(3,2) -= sc; eA0(3,3) += cc; eA0(3,4) += sc;
   eA0(4,1) -= sc; eA0(4,2) -= ss; eA0(4,3) += sc; eA0(4,4) += ss;
}


void Bar2DL2::BodyRHS(const Vect<real_t>& f)
{
   eRHS(1) += 0.5*_el_geo.length*f((*_theElement)(1)->n(),1);
   eRHS(2) += 0.5*_el_geo.length*f((*_theElement)(1)->n(),2);
   eRHS(3) += 0.5*_el_geo.length*f((*_theElement)(2)->n(),1);
   eRHS(4) += 0.5*_el_geo.length*f((*_theElement)(2)->n(),2);
}


void Bar2DL2::build()
{
   *Equa::_A = 0;
   MESH_EL {
      set(the_element);
      if (_analysis==TRANSIENT) {
         if (_terms&LUMPED_MASS)
            LMass();
         if (_terms&MASS)
            Mass();
      }
      Stiffness();
      Equa::_A->Assembly(The_element,eA0.get());
      Equa::_b->Assembly(The_element,eRHS.get());
   }
   if (_pf!=nullptr)
      Load();
}


void Bar2DL2::Load()
{
   for (size_t i=1; i<=_nb_nodes; ++i) {
      (*_b)(i,1) += (*_pf)(i,1);
      (*_b)(i,2) += (*_pf)(i,2);
   }
}

  
void Bar2DL2::getStresses(Vect<real_t>& s)
{
   s.setSize(_nb_el,2);
   MESH_EL {
      set(the_element);
      s(element_label,1) = _young*_section*(_eu(3)-_eu(1))/_el_geo.length;
      s(element_label,2) = _young*_section*(_eu(4)-_eu(2))/_el_geo.length;
   }
}

} /* namespace OFELI */
