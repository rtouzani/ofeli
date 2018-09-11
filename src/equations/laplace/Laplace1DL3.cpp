/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2018 Rachid Touzani

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

namespace OFELI {

Laplace1DL3::Laplace1DL3(Mesh& ms)
            : Equation<real_t,3,3,1,1>(ms)
{
   _lsf = _rsf = 0;
   _A.setSize(ms.getNbDOF(),2,2);
}


Laplace1DL3::Laplace1DL3(Mesh&         ms,
                         Vect<real_t>& u)
            : Equation<real_t,3,3,1,1>(ms)
{
   _u = &u;
   _lsf = _rsf = 0;
   _A.setSize(ms.getNbDOF(),2,2);
}


Laplace1DL3::Laplace1DL3(Element* el)
{
   set(el);
}


Laplace1DL3::~Laplace1DL3()
{
}


void Laplace1DL3::set(Element* el)
{
   static real_t s[3]={-1,0,1};
   _theElement = el;
   Line3 ln(_theElement);
   _det = ln.getDet();
   for (size_t k=0; k<3; k++) {
      ln.setLocal(s[k]);
      for (size_t i=1; i<=3; i++)
         _dSh(i,k+1) = ln.DSh(i);
   }
   eMat = 0;
   eRHS = 0;
}


void Laplace1DL3::setTraction(real_t f,
                              int    lr)
{
   if (lr==-1)
      _lsf = f;
   if (lr==1)
      _rsf = f;
   _sf_given = true;
}


void Laplace1DL3::Matrix(real_t coef)
{
   real_t c = coef*OFELI_THIRD*_det;
   for (size_t i=1; i<=3; i++) {
      for (size_t j=1; j<=3; j++)
         eMat(i,j) += c*(_dSh(i,1)*_dSh(j,1)+4*_dSh(i,2)*_dSh(j,2)+_dSh(i,3)*_dSh(j,3));
   }
}


void Laplace1DL3::BodyRHS(const Vect<real_t>& f)
{
   real_t c = OFELI_THIRD*_det;
   eRHS(1) +=   c*f(_theElement->getNodeLabel(1));
   eRHS(2) += 4*c*f(_theElement->getNodeLabel(2));
   eRHS(3) +=   c*f(_theElement->getNodeLabel(3));
}


void Laplace1DL3::BoundaryRHS(int    n,
                              real_t p)
{                  
   if (n==-1)
      eRHS(1) -= p;
   else
      eRHS(3) += p;
}


int Laplace1DL3::run()
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
   _A.Prescribe(*_u,*_bc);
   return _A.FactorAndSolve(*_u);
}

} /* namespace OFELI */
