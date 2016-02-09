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

  Implementation of class 'Triang6' for straight 6-Node Triangle Finite Element

  ==============================================================================*/

#include "util/macros.h"
#include "shape_functions/Triang6S.h"
#include "shape_functions/Line3.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "linear_algebra/Point.h"
#include "util/util.h"

namespace OFELI {

string itos(int i);

Triang6S::Triang6S()
{
   _sh.resize(6);
   _node.resize(6);
   _x.resize(6);
   _dsh.resize(6);
   _dshl.resize(6);
}


Triang6S::Triang6S(const Element* el)
{
   try {
      if (el->getNbNodes() != 6)
         THROW_RT("Triang6S(Element *): Illegal number of element nodes: " + itos(el->getNbNodes()));
   }
   CATCH("Triang6S");
   _sh.resize(6);
   _node.resize(6);
   _x.resize(6);
   _dsh.resize(6);
   _dshl.resize(6);
   for (size_t i=0; i<6; i++) {
      Node *node = (*el)(i+1);
      _x[i] = node->getCoord();
      _node[i] = node->n();
   }
   _x21 = _x[1] - _x[0];
   _x31 = _x[2] - _x[0];
   _x32 = _x[2] - _x[1];
   _det = _x21.x*_x31.y-_x21.y*_x31.x;
   try {
      if (_det < 0.0)
         THROW_RT("set(Element *): Negative determinant of jacobian");
   }
   CATCH("Triang6S");
   try {
      if (_det == 0.0)
         THROW_RT("set(Element *): Determinant of jacobian is null");
   }
   CATCH("Triang6S");
   _c =(_x[0] + _x[1] + _x[2])*OFELI_THIRD;
   _el = el;
}


real_t Triang6S::Sh(      size_t         i,
                    const Point<real_t>& s) const
{
   switch (i) {
      case 1: return (1.-s.x-s.y)*(2*(1.-s.x-s.y)-1.);
      case 2: return s.x*(2*s.x-1.);
      case 3: return s.y*(2*s.y-1.);
      case 4: return 4*s.x*(1.-s.x-s.y);
      case 5: return 4*s.x*s.y;
      case 6: return 4*s.y*(1.-s.x-s.y);
   }
   return 0;
}


Point<real_t> Triang6S::DSh(      size_t         i,
                            const Point<real_t>& s) const
{
   Point<real_t> dsh;
   switch (i) {
      case 1: dsh.x = (4*s.x+4*s.y-3.)*_x32.y/_det;
              dsh.y = (3.-4*s.x-4*s.y)*_x32.x/_det;
              return dsh;
      case 2: dsh.x = (4*s.x-1.)*_x31.y/_det;
              dsh.y = (1.-4*s.x)*_x31.x/_det;
              return dsh;
      case 3: dsh.x = (1.-4*s.y)*_x21.y/_det;
              dsh.y = (4*s.y-1.)*_x21.x/_det;
              return dsh;
      case 4: dsh.x = 4./_det*(_x31.y*(1-s.y-2*s.x)+_x21.y*s.x);
              dsh.y = 4./_det*(_x31.x*(2*s.x+s.y-1)-_x21.x*s.x);
              return dsh;
      case 5: dsh.x = 4./_det*(s.y*_x31.y-s.x*_x21.y);
              dsh.y = 4./_det*(s.x*_x21.x-s.y*_x31.x);
              return dsh;
      case 6: dsh.x = 4./_det*(_x21.y*(2*s.y+s.x-1)-_x31.y*s.y);
              dsh.y = 4./_det*(_x31.x*s.y+_x21.x*(1-s.x-2*s.y));
              return dsh;
   }
   return Point<real_t>(0.);
}


  Point<real_t> Triang6S::Grad(const LocalVect<real_t,6>& u,
                               const Point<real_t>&       s) const
{
   Point<real_t> g(0.,0.,0.);
   for (size_t i=1; i<=6; i++)
      g += u(i)*DSh(i,s);
   return g;
}


real_t Triang6S::getMaxEdgeLength() const
{
   real_t h = 0;
   for (size_t i=0; i<3; i++)
      h = std::max(h,Distance(_x[i],_x[(i+1)%3]));
   return h;
}


real_t Triang6S::getMinEdgeLength() const
{
   real_t h = Distance(_x[0],_x[1]);
   for (size_t i=1; i<3; i++)
      h =std:: min(h,Distance(_x[i],_x[(i+1)%3]));
   return h;
}

} /* namespace OFELI */
