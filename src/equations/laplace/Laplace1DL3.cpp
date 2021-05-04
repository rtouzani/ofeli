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

                    Class Laplace1DL3 : Laplace Equation
              using 3-Node Line finite element in one dimension

  ==============================================================================*/


#include "equations/laplace/Laplace1DL3.h"
#include "shape_functions/Line3.h"
#include "equations/Equa_impl.h"
#include "equations/Equation_impl.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

Laplace1DL3::Laplace1DL3()
            : Equation<real_t,3,3,1,1>()
{
   if (Verbosity>0)
      cout << "Solving the Laplace equation in 1D using P2 finite element." << endl;
}


Laplace1DL3::Laplace1DL3(Mesh& ms)
            : Equation<real_t,3,3,1,1>(ms)
{
   _lsf = _rsf = 0;
   setMatrixType(SPARSE|SYMMETRIC);
   setSolver(CG_SOLVER,DILU_PREC);
   if (Verbosity>0)
      cout << "Solving the Laplace equation in 1D using P2 finite element." << endl;
   if (Verbosity>1) {
      cout << "Matrix is stored in sparse format." << endl;
      cout << "Linear system is solved by Conjugate Gradient with DILU preconditioner." << endl;
   }
}


Laplace1DL3::Laplace1DL3(Mesh&         ms,
                         Vect<real_t>& u)
            : Equation<real_t,3,3,1,1>(ms,u)
{
   _u = &u;
   _lsf = _rsf = 0;
   setMatrixType(SPARSE|SYMMETRIC);
   setSolver(CG_SOLVER,DILU_PREC);
   if (Verbosity>0)
      cout << "Solving the Laplace equation in 1D using P1 finite element." << endl;
   if (Verbosity>1) {
      cout << "Matrix is stored in sparse format." << endl;
      cout << "Linear system is solved by Conjugate Gradient with DILU preconditioner." << endl;
   }
}


Laplace1DL3::~Laplace1DL3()
{
}


void Laplace1DL3::set(const Element* el)
{
   static real_t s[3] = {-1,0,1};
   _theElement = el, _theSide = nullptr;
   Line3 ln(_theElement);
   _el_geo.det = ln.getDet();
   for (size_t k=0; k<3; k++) {
      ln.setLocal(s[k]);
      _dSh[k] = ln.DSh();
   }
   eMat = 0, eA0 = 0, eA1 = 0;
   eRHS = 0;
}


void Laplace1DL3::setTraction(real_t f,
                              int    lr)
{
   if (lr==-1)
      _lsf = f;
   if (lr==1)
      _rsf = f;
}


void Laplace1DL3::LHS()
{
   real_t c = OFELI_THIRD*_el_geo.det;
   for (size_t i=1; i<=3; i++) {
      for (size_t j=1; j<=3; j++)
         eMat(i,j) += c*(_dSh[0](i).x*_dSh[0](j).x + 4*(_dSh[1](i).x*_dSh[1](j).x
                                                        + _dSh[2](i).x*_dSh[2](j).x));
   }
}


void Laplace1DL3::BodyRHS(const Vect<real_t>& f)
{
   real_t c = OFELI_THIRD*_el_geo.det;
   eRHS(1) += c*f((*_theElement)(1)->n());
   eRHS(2) += c*f((*_theElement)(2)->n());
   eRHS(3) += c*f((*_theElement)(3)->n());
}


void Laplace1DL3::BoundaryRHS(const Vect<real_t>& h)
{
   eRHS(1) -= h(1);
   eRHS(3) += h(_nb_nodes);
}

} /* namespace OFELI */
