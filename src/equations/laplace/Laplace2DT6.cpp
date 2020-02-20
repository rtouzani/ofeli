/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

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

                      Class Laplace2DT6: Laplace Equation
              using 6-Node Triangular finite element in two dimensions

  ==============================================================================*/


#include "equations/laplace/Laplace2DT6.h"
#include "shape_functions/Triang6S.h"
#include "shape_functions/Line3.h"
#include "equations/AbsEqua_impl.h"
#include "equations/Equation_impl.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

Laplace2DT6::Laplace2DT6()
            : Equation<real_t,6,6,3,3>()
{
   _equation_name = "Laplace";
   _finite_element = "2-D, 6-Node Triangles (P2)";
   if (Verbosity>0)
      cout << "Solving the Laplace equation in 2D using P2 finite element triangle." << endl;
}


Laplace2DT6::Laplace2DT6(Mesh& ms)
            : Equation<real_t,6,6,3,3>(ms)
{
   _equation_name = "Laplace";
   _finite_element = "2-D, 6-Node Triangles (P2)";
   setMatrixType(SPARSE);
   setSolver(CG_SOLVER,DILU_PREC);
   if (Verbosity>0)
      cout << "Solving the Laplace equation in 2D using P2 finite element triangle." << endl;
   if (Verbosity>1) {
      cout << "Matrix is stored in sparse format." << endl;
      cout << "Linear system is solved by Conjugate Gradient with DILU preconditioner." << endl;
   }
}


Laplace2DT6::Laplace2DT6(Mesh&         ms,
                         Vect<real_t>& u)
            : Equation<real_t,6,6,3,3>(ms,u)
{
   setMatrixType(SPARSE|SYMMETRIC);
   setSolver(CG_SOLVER,DILU_PREC);
   _equation_name = "Laplace";
   _finite_element = "2-D, 6-Node Triangles (P2)";
   if (Verbosity>0)
      cout << "Solving the Laplace equation in 2D using P2 finite element triangle." << endl;
   if (Verbosity>1) {
      cout << "Matrix is stored in sparse format." << endl;
      cout << "Linear system is solved by Conjugate Gradient with DILU preconditioner." << endl;
   }
}


void Laplace2DT6::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   Triang6S tr(el);
   _el_geo.area = tr.getArea();
   ElementNodeCoordinates();
   tr.atMidEdges(_dSh,_wg);
   eA0 = 0, eA1 = 0, eMat = 0;
   eRHS = 0;
}


void Laplace2DT6::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Line3 ln(sd);
   //   _el_geo.length = ln.getLength();
   SideNodeCoordinates();
   sRHS = 0;
}


void Laplace2DT6::buildEigen(int opt)
{
}


void Laplace2DT6::LHS()
{
   for (size_t i=0; i<6; ++i) {
      for (size_t j=0; j<6; ++j)
         eMat(i+1,j+1) = _wg[0]*((_dSh[3*i],_dSh[3*j]) + (_dSh[3*i+1],_dSh[3*j+1]) + (_dSh[3*i+2],_dSh[3*j+2]));
   }
}


void Laplace2DT6::BodyRHS(const Vect<real_t>& f)
{
   eRHS(4) += _el_geo.area*OFELI_THIRD*f((*_theElement)(4)->n());
   eRHS(5) += _el_geo.area*OFELI_THIRD*f((*_theElement)(5)->n());
   eRHS(6) += _el_geo.area*OFELI_THIRD*f((*_theElement)(6)->n());
}


void Laplace2DT6::BoundaryRHS(const Vect<real_t>& f)
{
   Line3 ln(_theSide);
   real_t c = ln.getDet()*OFELI_THIRD;
   ln.setLocal(-1.0);
   sRHS(1) += c*f((*_theSide)(1)->n());
   ln.setLocal(0.0);
   sRHS(2) += c*f((*_theSide)(2)->n());
   ln.setLocal(1.0);
   sRHS(3) += 4*c*f((*_theSide)(3)->n());
}

} /* namespace OFELI */
