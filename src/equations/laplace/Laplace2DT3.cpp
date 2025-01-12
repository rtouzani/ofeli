/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2025 Rachid Touzani

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

                      Class Laplace2DT3: Laplace Equation
              using 3-Node Triangular finite element in two dimensions

  ==============================================================================*/


#include "equations/laplace/Laplace2DT3.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"
#include "equations/Equation_impl.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

Laplace2DT3::Laplace2DT3()
            : Equation<3,3,2,2>()
{
   _equation_name = "Laplace";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   if (Verbosity>0)
      cout << "Solving the Laplace equation in 2D using P1 finite element triangle." << endl;
}


Laplace2DT3::Laplace2DT3(Mesh& ms)
            : Equation<3,3,2,2>(ms)
{
   _equation_name = "Laplace";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   setMatrixType(SKYLINE|SYMMETRIC);
   setSolver(DIRECT_SOLVER);
   if (Verbosity>0)
      cout << "Solving the Laplace equation in 2D using P1 finite element triangle." << endl;
}


Laplace2DT3::Laplace2DT3(Mesh&         ms,
                         Vect<real_t>& u)
            : Equation<3,3,2,2>(ms,u)
{
   setMatrixType(SKYLINE|SYMMETRIC);
   setSolver(DIRECT_SOLVER);
   _equation_name = "Laplace";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   if (Verbosity>0)
      cout << "Solving the Laplace equation in 2D using P1 finite element triangle." << endl;
}


Laplace2DT3::Laplace2DT3(Mesh&         ms,
                         Vect<real_t>& b,
                         Vect<real_t>& Dbc,
                         Vect<real_t>& Nbc,
                         Vect<real_t>& u)
            : Equation<3,3,2,2>(ms)
{
   _bf = &b;
   _u = &u;
   setMatrixType(SKYLINE|SYMMETRIC);
   setSolver(DIRECT_SOLVER);
   _bc = &Dbc;
   _sf = &Nbc;
}


void Laplace2DT3::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   Triang3 tr(_theElement);
   _el_geo.area = tr.getArea();
   _dSh = tr.DSh();
   eMat = 0, eA0 = 0, eA1 = 0;
   eRHS = 0;
}


void Laplace2DT3::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Line2 ln(_theSide);
   _el_geo.length = ln.getLength();
   sRHS = 0;
}


void Laplace2DT3::buildEigen(int opt)
 {
   for (size_t i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eA0(i,j) += _el_geo.area*(_dSh[i-1],_dSh[j-1]);
   if (opt==0) {
      real_t c = 0.5*OFELI_SIXTH*_el_geo.area;
      eA1(1,1) += 2*c; eA1(2,2) += 2*c; eA1(3,3) += 2*c;
      eA1(1,2) +=   c; eA1(2,1) +=   c; eA1(1,3) +=   c;
      eA1(3,1) +=   c; eA1(2,3) +=   c; eA1(3,2) +=   c;
   }
   else {
      real_t c = OFELI_THIRD*_el_geo.area;
      eA1(1,1) += c; eA1(2,2) += c; eA1(3,3) += c;
   }
}


void Laplace2DT3::LHS()
{
   for (size_t i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eMat(i,j) += _el_geo.area*(_dSh[i-1],_dSh[j-1]);
}


void Laplace2DT3::BodyRHS(const Vect<real_t>& f)
{
   for (size_t i=1; i<=3; i++)
      eRHS(i) += f((*_theElement)(i)->n())*_el_geo.area*OFELI_THIRD;
}


void Laplace2DT3::BoundaryRHS(const Vect<real_t>& h)
{
   if (_theSide->getCode(1)>0) {
      real_t z = 0.5*_el_geo.length;
      if (h.getDOFType()==NODE_DOF) {
         sRHS(1) += z*h((*_theSide)(1)->n());
         sRHS(2) += z*h((*_theSide)(2)->n());
      }
      else if (h.getDOFType()==SIDE_DOF) {
         sRHS(1) += z*h(_theSide->n());
         sRHS(2) += z*h(_theSide->n());
      }
   }
}

} /* namespace OFELI */
