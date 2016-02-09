/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   This file is part of the OFELI Library

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

                    Implementation of class SteklovPoincare2DBE

  ==============================================================================*/


#include "equations/laplace/SteklovPoincare2DBE.h"

namespace OFELI {

SteklovPoincare2DBE::SteklovPoincare2DBE(const Mesh& mesh,
                                               bool  ext)
{
   setMesh(mesh,ext);
}


SteklovPoincare2DBE::SteklovPoincare2DBE(const Mesh& mesh,
                                         const Vect<real_t>& g, 
                                               Vect<real_t>& b,
                                               bool          ext)
{
   setMesh(mesh,ext);
   Solve(b,g);
}


void SteklovPoincare2DBE::setMesh(const Mesh& mesh,
                                        bool  ext)
{
   _theMesh = &mesh;
   _nn.setSize(_theMesh->getNbSides());
   _length.setSize(_theMesh->getNbSides());
   _center.setSize(_theMesh->getNbSides());
   _ttg.setSize(_theMesh->getNbSides());
   _ext = -1;
   if (ext)
      _ext = 1; 
   _util();
   _nb_eq = _theMesh->getNbSides();
   _A.setSize(_nb_eq,_nb_eq);
}


int SteklovPoincare2DBE::Solve(      Vect<real_t>& b,
                               const Vect<real_t>& g)
{
   b = 0;
   _A = 0;
   for (size_t i=1; i<=_theMesh->getNbSides(); i++) {
      b(i) -= 0.5*g(i);
      for (size_t j=1; j<=_theMesh->getNbSides(); j++) {
         _h = _length(j);
         Point<real_t> z = _center(j) - _center(i);
         real_t s = single_layer(j,z);
         real_t d = double_layer(j,z);
         _A.add(i,j,_ext*s);
         b.add(i,_ext*g(j)*d);
      }
   }
   Vect<real_t> x(b.size());
   Prec<real_t> p(_A,IDENT_PREC);
   real_t toler = 1.e-8;
   int nb_it = GMRes(_A,p,b,x,_nb_eq/4,1000,toler,0);
   b = x;
   return nb_it;
}


void SteklovPoincare2DBE::_util()
{
   Point<real_t> N, x1, x2, T, c;
   for (size_t s=1; s<=_theMesh->getNbSides(); s++) {
      const Side *sd = _theMesh->getPtrSide(s);
      const Element *el = sd->getNeighborElement(1);
      try {
         if (el->getShape()!=TRIANGLE)
            THROW_RT("_util(): This class is valid for triangles only.");
      }
      CATCH("SteklovPoincare");
      Triang3 tr(el);
      x1 = sd->getPtrNode(1)->getCoord();
      x2 = sd->getPtrNode(2)->getCoord();
      T = x2 - x1;
      _center(s) = 0.5*(x1+x2);
      c = _center(s) - tr.getCenter();
      N.x =  T.y;
      N.y = -T.x;
      if (c*N<0) {
         N.x = -N.x;
         N.y = -N.y;
      }
      _nn(s) = N;
      _ttg(s) = T;
      _length(s) = T.Norm();
   }
}

} /* namespace OFELI */
