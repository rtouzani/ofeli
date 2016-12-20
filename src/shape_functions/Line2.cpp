/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2017 Rachid Touzani

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

            Implementation of class 'Line2' for Line 3-D Finite Element

  ==============================================================================*/


#include "shape_functions/Line2.h"
#include "mesh/Element.h"
#include "linear_algebra/Point.h"
#include "linear_algebra/LocalVect.h"
#include "util/util.h"

namespace OFELI {

Line2::Line2()
{
   _sh.resize(2);
   _node.resize(2);
   _x.resize(2);
   _dsh.resize(2);
   _dshl.resize(2);
   _sd = NULL;
   _el = NULL;
}


Line2::Line2(const Element* el)
{
   try {
      if (el->getNbNodes() != 2)
         THROW_RT("Line2(Element *): Illegal number of element nodes: " + itos(el->getNbNodes()));
   }
   CATCH("Line2");
   _sh.resize(2);
   _node.resize(2);
   _x.resize(2);
   _dsh.resize(2);
   _dshl.resize(2);
   for (size_t i=0; i<2; i++) {
      Node *node = el->getPtrNode(i+1);
      _x[i] = node->getCoord();
      _node[i] = node->n();
   }
   _c = 0.5*(_x[0]+_x[1]);
   _el = el;
   _sd = NULL;
   Point<real_t> dl = 0.5*(_x[1] - _x[0]);
   _det = dl.Norm();
   try {
      if (_det == 0.0)
         THROW_RT("Line2(Element *): Determinant of jacobian is null");
   }
   CATCH("Line2");
   _length = 2*_det;
}


Line2::Line2(const Side* side)
{
   try {
      if (side->getNbNodes() != 2)
         THROW_RT("Line2(Side *): Illegal number of side nodes: " + itos(side->getNbNodes()));
   }
   CATCH("Line2");
   _sh.resize(2);
   _node.resize(2);
   _x.resize(2);
   _dsh.resize(2);
   _dshl.resize(2);
   for (size_t i=0; i<2; i++) {
      Node *node = side->getPtrNode(i+1);
      _x[i] = node->getCoord();
      _node[i] = node->n();
   }
   _c = 0.5*(_x[0] + _x[1]);
   _sd = side;
   _el = NULL;
   Point<real_t> dl = 0.5*(_x[1] - _x[0]);
   _det = dl.Norm();
   try {
      if (_det == 0.0)
         THROW_RT("Line2(Side *): Determinant of jacobian is null");
   }
   CATCH("Line2");
   _length = 2*_det;
}
   
   
Line2::Line2(const Edge* edge)
{
   _sh.resize(2);
   _node.resize(2);
   _x.resize(2);
   _dsh.resize(2);
   _dshl.resize(2);
   for (size_t i=0; i<2; i++) {
      Node *node = edge->getPtrNode(i+1);
      _x[i] = node->getCoord();
      _node[i] = node->n();
   }
   _c = 0.5*(_x[0] + _x[1]);
   _ed = edge;
   _sd = NULL;
   _el = NULL;
   Point<real_t> dl = 0.5*(_x[1] - _x[0]);
   _det = dl.Norm();
   try {
      if (_det == 0.0)
         THROW_RT("Line2(Edge *): Determinant of jacobian is null");
   }
   CATCH("Line2");
   _length = 2*_det;
}
   

Line2::~Line2() { }


Point<real_t> Line2::getNormal() const
{
   return Point<real_t>(_x[1].y-_x[0].y,_x[0].x-_x[1].x)/_length;
}


Point<real_t> Line2::getTangent() const
{
   return Point<real_t>(_x[1]-_x[0])/_length;
}


real_t Line2::Sh(size_t i,
                 real_t s) const
{
   switch (i) {
      case 1: return 0.5*(1.-s);
      case 2: return 0.5*(1.+s);
   }
   return 0;
}


real_t Line2::DSh(size_t i) const
{
   real_t sh=0;
   switch (i) {

      case 1:
          return (-0.5/_det);

      case 2:
          return ( 0.5/_det);

   }
   return sh;
}


Point<real_t> Line2::getRefCoord(const Point<real_t>& x)
{
   Point<real_t> a;
   return a+x;
}


bool Line2::isIn(const Point<real_t>& x)
{
   Point<real_t> s = getRefCoord(x);
   if (s.x>=-1 && s.x<=1)
      return true;
   else
      return false;
}


real_t Line2::getInterpolate(const Point<real_t>&       x,
                             const LocalVect<real_t,2>& v)
{
   Point<real_t> s = getRefCoord(x);
   return (Sh(1,s)*v(1) + Sh(2,s)*v(2));
}

} /* namespace OFELI */
