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

     Implementation of class 'Line2H' for Line 2-Node Hermite Finite Element

  ==============================================================================*/


#include "shape_functions/Line2H.h"
#include "util/macros.h"
#include "util/util.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "linear_algebra/Point.h"

namespace OFELI {

string itos(int i);

Line2H::Line2H()
{
   _sh.resize(2);
   _node.resize(2);
   _x.resize(2);
   _dsh.resize(2);
   _dshl.resize(2);
   _el = NULL;
   _sd = NULL;
}


Line2H::Line2H(const Element* el)
{
   if (el->getNbNodes() != 2)
      throw OFELIException("Line2H::Line2H(Element *): Illegal number of element nodes: " +
                           itos(el->getNbNodes()));
   _sh.resize(2);
   _node.resize(2);
   _x.resize(2);
   _dsh.resize(2);
   _dshl.resize(2);
   for (size_t i=0; i<2; i++) {
      Node *node = _el->getPtrNode(i+1);
      _x[i] = node->getCoord();
      _node[i] = node->n();
   }
   real_t dxdl = 0.5*(_x[1].x - _x[0].x);
   real_t dydl = 0.5*(_x[1].y - _x[0].y);
   real_t dzdl = 0.5*(_x[1].z - _x[0].z);
   _det = sqrt(dxdl*dxdl + dydl*dydl + dzdl*dzdl);
   if (_det == 0.0)
      throw OFELIException("Line2H::Line2H(Element *): Determinant of jacobian is null");
   _el = el;
   _sd = NULL;
}


Line2H::Line2H(const Side* side)
{
   if (side->getNbNodes() != 2)
      throw OFELIException("Line2H::Line2H(Side *): Illegal number of element sides: " +
                           itos(side->getNbNodes()));
   for (size_t i=0; i<2; i++) {
      Node *node = side->getPtrNode(i+1);
      _x[i] = node->getCoord();
      _node[i] = node->n();
   }
   _sh.resize(2);
   _node.resize(2);
   _x.resize(2);
   _dsh.resize(2);
   _dshl.resize(2);
   real_t dxdl = 0.5*(_x[1].x - _x[0].x);
   real_t dydl = 0.5*(_x[1].y - _x[0].y);
   real_t dzdl = 0.5*(_x[1].z - _x[0].z);
   _det = sqrt(dxdl*dxdl + dydl*dydl + dzdl*dzdl);
   _el = NULL;
   _sd = side;
}


real_t Line2H::Sh(size_t i,
                  real_t s) const
{
   real_t sm = 1.-s, sp = 1.+s, x = 0;
   switch (i) {

      case 1:
          x = 0.25*sm*sm*(2.+s);

      case 2:
          x = 0.25*sp*sp*(2.-s);
   }
   return x;
}


real_t Line2H::DSh(size_t i,
                   real_t s) const
{
   real_t st = 1.-s*s, x = 0;
   switch (i) {

      case 1:
          x = -0.75*st/_det;
      case 2:
          x =  0.75*st/_det;
   }
   return x;
}


real_t Line2H::D2Sh(size_t i,
                    real_t s) const
{
   real_t x = 0;
   switch (i) {

      case 1:
          x =  1.5*s/(_det*_det);
      case 2:
          x = -1.5*s/(_det*_det);
   }
   return x;
}

} /* namespace OFELI */
