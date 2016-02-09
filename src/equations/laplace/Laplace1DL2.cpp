/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzani

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

                       Class Laplace1DL2 : Laplace Equation
                  using 2-Node Line finite element in one dimension

  ==============================================================================*/


#include "equations/laplace/Laplace1DL2.h"
#include "shape_functions/Line2.h"

namespace OFELI {

Laplace1DL2::Laplace1DL2(Mesh&         ms,
                         Vect<real_t>& u)
{
   _u = &u;
   _theMesh = &ms;
   _A.setSize(_theMesh->getNbDOF());
   _lsf = _rsf = 0;
   _is_lbc = _is_rbc = false;
}


Laplace1DL2::Laplace1DL2(Element* el)
{
   _theElement = el;
   Line2 ln(el);
   _length = ln.getLength();
   eMat = 0;
   eRHS = 0;
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


void Laplace1DL2::Matrix(real_t coef)
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


void Laplace1DL2::BoundaryRHS(int    n,
                              real_t p)
{
   if (n==-1)
      eRHS(1) -= p;
   else
      eRHS(2) += p;
}


int Laplace1DL2::run()
{
   *_u = 0;
   MESH_EL {
      set(theElement);
      Matrix();
      if (_bf_given)
         BodyRHS(*_bf);
      if (_sf_given) {
         if (theElementLabel==1)
            BoundaryRHS(-1,_lsf);
         if (theElementLabel==_theMesh->getNbElements())
            BoundaryRHS(1,_rsf);
      }
      ElementAssembly(_A);
      ElementAssembly(*_u);
   }
   if (_is_lbc) {
      _A.set(1,2,0);
      _u->set(1,_lbc*_A(1,1));
   }
   if (_is_rbc) {
      size_t n = _theMesh->getNbDOF();
      _A.set(n,n-1,0);
      _u->set(n,_rbc*_A(n,n));
   }
   return _A.Solve(*_u);
}

} /* namespace OFELI */
