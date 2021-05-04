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

                       Class Laplace1DL2: Laplace Equation
                  using 2-Node Line finite element in one dimension

  ==============================================================================*/


#include "equations/laplace/Laplace1DL2.h"
#include "shape_functions/Line2.h"
#include "equations/Equa_impl.h"
#include "equations/Equation_impl.h"
#include "linear_algebra/Vect_impl.h"


namespace OFELI {

Laplace1DL2::Laplace1DL2()
            : _lsf(0), _rsf(0), _is_lbc(false), _is_rbc(false)
{
   _equation_name = "Laplace";
   _finite_element = "1-D, 2-Node Lines (P1)";
   if (Verbosity>0)
      cout << "Solving the Laplace equation in 1D using P1 finite element." << endl;
}


Laplace1DL2::Laplace1DL2(Mesh& ms)
            : Equation<real_t,2,2,1,1>(ms),
              _lsf(0), _rsf(0), _is_lbc(false), _is_rbc(false)
{
   _theMesh = &ms;
   setMatrixType(TRIDIAGONAL);
   setSolver(DIRECT_SOLVER);
   if (Verbosity>0)
      cout << "Solving the Laplace equation in 1D using P1 finite element." << endl;
   if (Verbosity>1) {
      cout << "Matrix is stored in tridiagonal storage." << endl;
      cout << "Linear system is solved by direct solver." << endl;
   }
}


Laplace1DL2::Laplace1DL2(Mesh&         ms,
                         Vect<real_t>& u)
            : Equation<real_t,2,2,1,1>(ms,u),
              _lsf(0), _rsf(0), _is_lbc(false), _is_rbc(false)
{
   _u = &u; 
   _theMesh = &ms;
   setMatrixType(TRIDIAGONAL);
   setSolver(DIRECT_SOLVER);
   if (Verbosity>0)
      cout << "Solving the Laplace equation in 1D using P1 finite element." << endl;
   if (Verbosity>1) {
      cout << "Matrix is stored in tridiagonal storage." << endl;
      cout << "Linear system is solved by direct solver." << endl;
   }
}


void Laplace1DL2::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   Line2 ln(_theElement);
   _el_geo.length = ln.getLength();
   _el_geo.det = ln.getDet();
   _dSh = ln.DSh();
   eRHS = 0;
   eMat = 0, eA0 = 0, eA1 = 0;
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
}


void Laplace1DL2::setTraction(real_t f,
                              int    lr)
{
   if (lr==-1)
      _lsf = f;
   if (lr==1)
      _rsf = f;
}


void Laplace1DL2::LHS()
{
   eMat(1,1) += 1./_el_geo.length;
   eMat(2,2) += 1./_el_geo.length;
   eMat(1,2) -= 1./_el_geo.length;
   eMat(2,1) -= 1./_el_geo.length;
}


void Laplace1DL2::buildEigen(int opt)
{
   eA0(1,1) += 1./_el_geo.length;
   eA0(2,2) += 1./_el_geo.length;
   eA0(1,2) -= 1./_el_geo.length;
   eA0(2,1) -= 1./_el_geo.length;
   if (opt==0) {
      real_t c = OFELI_SIXTH*_el_geo.length;
      eA1(1,1) += 2*c; eA1(2,2) += 2*c;
      eA1(1,2) +=   c; eA1(2,1) +=   c;
   }
   else {
      real_t c = 0.5*_el_geo.length;
      eA1(1,1) += c; eA1(2,2) += c;
   }
}


void Laplace1DL2::BodyRHS(const Vect<real_t>& f)
{
   eRHS(1) += 0.5*_el_geo.length*f((*_theElement)(1)->n());
   eRHS(2) += 0.5*_el_geo.length*f((*_theElement)(2)->n());
}


void Laplace1DL2::BoundaryRHS(const Vect<real_t>& h)
{
   eRHS(1) -= h(1);
   eRHS(2) += h(_nb_nodes);
}

} /* namespace OFELI */
