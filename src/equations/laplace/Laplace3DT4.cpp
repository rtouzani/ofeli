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

                      Class Laplace3DT4 : Laplace Equation
          using 4-Node Tetrahedral finite element in three dimensions

  ==============================================================================*/


#include "equations/laplace/Laplace3DT4.h"
#include "shape_functions/Tetra4.h"
#include "shape_functions/Triang3.h"
#include "equations/AbsEqua_impl.h"
#include "equations/Equation_impl.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

Laplace3DT4::Laplace3DT4()
            : Equation<real_t,4,4,3,3>()
{
   _equation_name = "Laplace";
   _finite_element = "3-D, 4-Node Tetrahedra (P1)";
   if (Verbosity>0)
      cout << "Solving the Laplace equation in 3D using P1 finite element tetrahedron." << endl;
}


Laplace3DT4::Laplace3DT4(Mesh& ms)
            : Equation<real_t,4,4,3,3>(ms)
{
   _equation_name = "Laplace";
   _finite_element = "3-D, 4-Node Tetrahedra (P1)";
   setMatrixType(SPARSE);
   setSolver(CG_SOLVER,DILU_PREC);
   if (Verbosity>0)
      cout << "Solving the Laplace equation in 3D using P1 finite element tetrahedron." << endl;
   if (Verbosity>1) {
      cout << "Matrix stored in sparse format." << endl;
      cout << "Linear system is solved by Conjugate Gradient with DILU preconditioner." << endl;
   }
}

Laplace3DT4::Laplace3DT4(Mesh&         ms,
                         Vect<real_t>& u)
            : Equation<real_t,4,4,3,3>(ms,u)
{
   _equation_name = "Laplace";
   _finite_element = "3-D, 4-Node Tetrahedra (P1)";
   setMatrixType(SPARSE);
   setSolver(CG_SOLVER,DILU_PREC);
   if (Verbosity>0)
      cout << "Solving the Laplace equation in 3D using P1 finite element tetrahedron." << endl;
   if (Verbosity>1) {
      cout << "Matrix stored in sparse format." << endl;
      cout << "Linear system is solved by Conjugate Gradient with DILU preconditioner." << endl;
   }
}


void Laplace3DT4::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   Tetra4 t(_theElement);
   _el_geo.volume = t.getVolume();
   _el_geo.det = t.getDet();
   _dSh = t.DSh();
   eRHS = 0;
   eMat = 0, eA0 = 0, eA1 = 0;
}


void Laplace3DT4::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Triang3 t(_theSide);
   _el_geo.area = t.getArea();
   sRHS = 0;
}


void Laplace3DT4::buildEigen(int opt)
{
   for (size_t i=1; i<=4; ++i)
      for (size_t j=1; j<=4; ++j)
         eA0(i,j) += _el_geo.volume*(_dSh[i-1],_dSh[j-1]);
   if (opt==0) {
      real_t c = 0.1*_el_geo.volume;
      real_t d = 0.5*c;
      eA1(1,1) += c; eA1(2,2) += c; eA1(3,3) += c; eA1(4,4) += c;
      eA1(1,2) += d; eA1(2,1) += d; eA1(1,3) += d; eA1(1,4) += d;
      eA1(3,1) += d; eA1(2,3) += d; eA1(3,2) += d; eA1(2,4) += d;
      eA1(4,1) += d; eA1(4,2) += d; eA1(4,3) += d; eA1(3,4) += d;
   }
   else {
      real_t c = 0.25*_el_geo.volume;
      for (size_t i=1; i<=4; i++)
         eA1(i,i) += c;
   }
}


void Laplace3DT4::LHS()
{
   for (size_t i=1; i<=4; ++i)
      for (size_t j=1; j<=4; ++j)
         eMat(i,j) = _el_geo.volume*(_dSh[i-1],_dSh[j-1]);
}


void Laplace3DT4::BodyRHS(const Vect<real_t>& f)
{
   for (size_t i=1; i<=4; ++i)
      eRHS(i) += 0.25*_el_geo.volume*f((*_theElement)(i)->n());
}


void Laplace3DT4::BoundaryRHS(const Vect<real_t>& h)
{
   if (_theSide->getCode(1)>0) {
      real_t c = OFELI_THIRD*_el_geo.area;
      sRHS(1) += c*h((*_theSide)(1)->n());
      sRHS(2) += c*h((*_theSide)(2)->n());
      sRHS(3) += c*h((*_theSide)(3)->n());
   }
}

} /* namespace OFELI */
