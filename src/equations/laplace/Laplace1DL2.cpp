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

                       Class Laplace1DL2: Laplace Equation
                  using 2-Node Line finite element in one dimension

  ==============================================================================*/


#include "equations/laplace/Laplace1DL2.h"
#include "shape_functions/Line2.h"

namespace OFELI {

Laplace1DL2::Laplace1DL2()
            : _lsf(0), _rsf(0), _is_lbc(false), _is_rbc(false)
{
   _bc_given = _bf_given = false;
}


Laplace1DL2::Laplace1DL2(Mesh& ms)
            : Equation<real_t,2,2,1,1>(ms),
              _lsf(0), _rsf(0), _is_lbc(false), _is_rbc(false)
{
   _theMesh = &ms;
   setMatrixType(TRIDIAGONAL);
   setSolver(DIRECT_SOLVER);
   _bc_given = _bf_given = false;
}


Laplace1DL2::Laplace1DL2(Mesh&         ms,
                         Vect<real_t>& u)
            : Equation<real_t,2,2,1,1>(ms),
              _lsf(0), _rsf(0), _is_lbc(false), _is_rbc(false)
{
   _u = &u;
   _theMesh = &ms;
   setMatrixType(TRIDIAGONAL);
   setSolver(DIRECT_SOLVER);
   _bc_given = _bf_given = false;
}


Laplace1DL2::Laplace1DL2(Element* el)
{
   set(el);
}


void Laplace1DL2::set(const Element* el)
{
   _nb_dof = 1;
   _theElement = theElement;
   Line2 ln(_theElement);
   _length = ln.getLength();
   _dSh(1) = ln.DSh(1); _dSh(2) = ln.DSh(2);
   eMat = 0;
   eRHS = 0;
   eA0 = 0;
   eA1 = 0;
}


void Laplace1DL2::setBoundaryCondition(real_t f,
                                       int    lr)
{
   if (lr==-1) {
      _is_lbc = true;
      _lbc = f;
   }
   if (lr==1) {
      _is_rbc = true;
      _rbc = f;
   }
   _bc_given = true;
}


void Laplace1DL2::setTraction(real_t f,
                              int    lr)
{
   if (lr==-1)
      _lsf = f;
   if (lr==1)
      _rsf = f;
   _sf_given = true;
}


void Laplace1DL2::LHS(real_t coef)
{
   eMat(1,1) += coef/_length;
   eMat(2,2) += coef/_length;
   eMat(1,2) -= coef/_length;
   eMat(2,1) -= coef/_length;
}


void Laplace1DL2::BodyRHS(const Vect<real_t>& f)
{
   eRHS(1) += 0.5*_length*f((*_theElement)(1)->n());
   eRHS(2) += 0.5*_length*f((*_theElement)(2)->n());
}


void Laplace1DL2::BoundaryRHS(const Vect<real_t>& h)
{
   eRHS(1) -= h(1);
   eRHS(2) += h(_theMesh->getNbNodes());
}

} /* namespace OFELI */
