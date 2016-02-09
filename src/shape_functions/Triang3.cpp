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

      Implementation of class 'Triang3' for Triangle 3-Node Finite Element

  ==============================================================================*/

#include "shape_functions/Triang3.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "linear_algebra/Point.h"
#include "linear_algebra/LocalVect.h"
#include "util/util.h"

namespace OFELI {

Triang3::Triang3()
{
   _sh.resize(3);
   _node.resize(3);
   _x.resize(3);
   _dsh.resize(3);
   _dshl.resize(3);
}

   
Triang3::Triang3(const Element* el)
{
   try {
      if (el->getNbNodes() != 3)
         THROW_RT("Triang3(Element *): Illegal number of element nodes: " + itos(el->getNbNodes()));
   }
   CATCH("Triang3");
   _el = el; _sd = NULL;
   set(el);
}


Triang3::Triang3(const Side* sd)
{
   try {
      if (sd->getNbNodes() != 3)
         THROW_RT("Triang3(Side *): Illegal number of side nodes: " + itos(sd->getNbNodes()));
   }
   CATCH("Triang3");
   _el = NULL; _sd = sd;
   set(sd);
}


Triang3::~Triang3() { }


void Triang3::set(const Element* el)
{
   _sh.resize(3);
   _node.resize(3);
   _x.resize(3);
   _dsh.resize(3);
   _dshl.resize(3);
   for (size_t i=0; i<3; i++) {
      Node *node = (*el)(i+1);
      _x[i] = node->getCoord();
      _node[i] = node->n();
   }
   _h1 = Distance(_x[0],_x[1]),
   _h2 = Distance(_x[1],_x[2]),
   _h3 = Distance(_x[2],_x[0]);
   _det = (_x[1].x-_x[0].x)*(_x[2].y-_x[0].y) - (_x[1].y-_x[0].y)*(_x[2].x-_x[0].x);
   try {
      if (_det<0.0)
         THROW_RT("set(Element *): Negative determinant of jacobian");
   }
   CATCH("Triang3");
   try {
      if (_det==0.0)
         THROW_RT("set(Element *): Determinant of jacobian is null");
   }
   CATCH("Triang3");
   _c = (_x[0]+_x[1]+_x[2])*OFELI_THIRD;
   _dsh[0].x = (_x[1].y - _x[2].y)/_det;
   _dsh[0].y = (_x[2].x - _x[1].x)/_det;
   _dsh[1].x = (_x[2].y - _x[0].y)/_det;
   _dsh[1].y = (_x[0].x - _x[2].x)/_det;
   _dsh[2].x = (_x[0].y - _x[1].y)/_det;
   _dsh[2].y = (_x[1].x - _x[0].x)/_det;
}


void Triang3::set(const Side* sd)
{
   _sh.resize(3);
   _node.resize(3);
   _x.resize(3);
   _dsh.resize(3);
   _dshl.resize(3);
   _x[0] = (*sd)(1)->getCoord();
   _x[1] = (*sd)(2)->getCoord();
   _x[2] = (*sd)(3)->getCoord();
   _node[0] = sd->getNodeLabel(1);
   _node[1] = sd->getNodeLabel(2);
   _node[2] = sd->getNodeLabel(3);
   _h1 = Distance(_x[0],_x[1]),
   _h2 = Distance(_x[1],_x[2]),
   _h3 = Distance(_x[2],_x[0]);
   real_t a = _x[0].y*(_x[1].z-_x[2].z) - _x[1].y*(_x[0].z-_x[2].z) + _x[2].y*(_x[0].z-_x[1].z);
   real_t b = _x[0].z*(_x[1].x-_x[2].x) - _x[1].z*(_x[0].x-_x[2].x) + _x[2].z*(_x[0].x-_x[1].x);
   real_t c = _x[0].x*(_x[1].y-_x[2].y) - _x[1].x*(_x[0].y-_x[2].y) + _x[2].x*(_x[0].y-_x[1].y);
   _det = sqrt(a*a + b*b + c*c);
   _c = (_x[0]+_x[1]+_x[2])*OFELI_THIRD;
   try {
      if (_det<0.0)
         THROW_RT("set(Side *): Negative determinant of jacobian");
   }
   CATCH("Triang3");
   try {
      if (_det==0.0)
         THROW_RT("set(Side *): Determinant of jacobian is null");
   }
   CATCH("Triang3");
}


double Triang3::getInterpolate(const Point<real_t>&       x,
                               const LocalVect<real_t,3>& v)
{
   Point<real_t> s = getRefCoord(x);
   return (Sh(1,s)*v(1)+Sh(2,s)*v(2)+Sh(3,s)*v(3));
}


Point<real_t> Triang3::Grad(const LocalVect<real_t,3>& u) const
{
   return u[0]*_dsh[0] + u[1]*_dsh[1] + u[2]*_dsh[2]; 
}


real_t Triang3::getMaxEdgeLength() const
{
   real_t h=0;
   for (size_t i=0; i<3; i++)
      h = std::max(h,Distance(_x[i],_x[(i+1)%3]));
   return h;
}


real_t Triang3::getMinEdgeLength() const
{
   real_t h = Distance(_x[0],_x[1]);
   for (size_t i=1; i<3; i++)
      h = std::min(h,Distance(_x[i],_x[(i+1)%3]));
   return h;
}

} /* namespace OFELI */
