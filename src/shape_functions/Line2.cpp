/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2024 Rachid Touzani

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
#include "linear_algebra/LocalVect_impl.h"

using std::to_string;

namespace OFELI {

Line2::Line2()
{
   _sd = nullptr;
   _el = nullptr;
}


Line2::Line2(const Element* el)
{
   if (el->getNbNodes() != 2)
      throw OFELIException("Line2::Line2(Element *): Illegal number of element nodes: " +
                           to_string(el->getNbNodes()));
   for (size_t i=1; i<=2; i++) {
     Node *node = (*el)(i);
      _x.push_back(node->getCoord());
      _node.push_back(node->n());
   }
   _c = 0.5*(_x[0]+_x[1]);
   _el = el;
   _sd = nullptr;
   Point<real_t> dl = 0.5*(_x[1] - _x[0]);
   _det = dl.Norm();
   if (_det == 0.0)
      throw OFELIException("Line2::Line2(Element *): Determinant of jacobian is null");
   _sh.resize(2);
   _dsh.resize(2);
   _dsh[0].x = -0.5/_det;
   _dsh[1].x =  0.5/_det;
   _length = 2*_det;
}


Line2::Line2(const Side* side)
{
   if (side->getNbNodes() != 2)
      throw OFELIException("Line2::Line2(Side *): Illegal number of side nodes: " +
                           to_string(side->getNbNodes()));
   for (size_t i=1; i<=2; i++) {
      Node *node = (*side)(i);
      _x.push_back(node->getCoord());
      _node.push_back(node->n());
   }
   _c = 0.5*(_x[0] + _x[1]);
   _sd = side;
   _el = nullptr;
   Point<real_t> dl = 0.5*(_x[1] - _x[0]);
   _det = dl.Norm();
   if (_det == 0.0)
      throw OFELIException("Line2::Line2(Side *): Determinant of jacobian is null");
   _sh.resize(2);
   _dsh.resize(2);
   _dsh[0].x = -0.5/_det;
   _dsh[1].x =  0.5/_det;
   _length = 2*_det;
}


Line2::Line2(const Edge* edge)
{
   _sh.resize(2);
   _dsh.resize(2);
   for (size_t i=1; i<=2; i++) {
      Node *node = (*edge)(i);
      _x.push_back(node->getCoord());
      _node.push_back(node->n());
   }
   _c = 0.5*(_x[0] + _x[1]);
   _ed = edge;
   _sd = nullptr;
   _el = nullptr;
   Point<real_t> dl = 0.5*(_x[1] - _x[0]);
   _det = dl.Norm();
   if (_det == 0.0)
      throw OFELIException("Line2::Line2(Edge *): Determinant of jacobian is null");
   _dsh[0].x = -0.5/_det;
   _dsh[1].x =  0.5/_det;
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


std::vector<Point<real_t> > Line2::DSh() const
{
   return _dsh;
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
