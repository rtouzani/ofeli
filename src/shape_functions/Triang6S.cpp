/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani

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

                         Implementation of class Triang6S

  ==============================================================================*/

#include "util/macros.h"
#include "shape_functions/Triang6S.h"
#include "shape_functions/Line3.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "linear_algebra/Point.h"

namespace OFELI {

Triang6S::Triang6S()
{
   _sh.resize(6);
   _node.resize(6);
   _x.resize(6);
   _dshl.resize(6);
}


Triang6S::Triang6S(const Element* el)
{
   _el = el;
   if (_el->getNbNodes() != 6)
     throw OFELIException("Triang6S::Triang6S(Element *): Illegal number of element nodes: "
                          + std::to_string(el->getNbNodes()));
   _sh.resize(6);
   _dsh.resize(6);
   _node.resize(6);
   _x.resize(6);
   _dshl.resize(6);
   for (size_t i=0; i<6; i++) {
      Node *node = (*el)(i+1);
      _x[i] = node->getCoord();
      _node[i] = node->n();
   }
   _x21 = _x[1] - _x[0];
   _x31 = _x[2] - _x[0];
   _x32 = _x[2] - _x[1];
   _det = _x21.x*_x31.y - _x21.y*_x31.x;
   _area = 0.5*_det;
   if (_det < 0.0)
      throw OFELIException("Triang6S::Triang6S(Element *): Negative determinant of jacobian");
   if (_det == 0.0)
      throw OFELIException("Triang6S::Triang6S(Element *): Determinant of jacobian is null");
   _c =(_x[0] + _x[1] + _x[2])*OFELI_THIRD;
}


void Triang6S::Sh(real_t s,
                  real_t t,
                  real_t *sh) const
{
   sh[0] = (1.-s-t)*(2*(1.-s-t)-1.);
   sh[1] = s*(2*s-1.);
   sh[2] = t*(2*t-1.);
   sh[3] = 4*s*(1.-s-t);
   sh[4] = 4*s*t;
   sh[5] = 4*t*(1.-s-t);
}


void Triang6S::setLocal(real_t s,
                        real_t t)
{
   _dsh[0].x = (4*s+4*t-3.)*_x32.y/_det;
   _dsh[0].y = (3.-4*s-4*t)*_x32.x/_det;
   _dsh[1].x = (4*s-1.)*_x31.y/_det;
   _dsh[1].y = (1.-4*s)*_x31.x/_det;
   _dsh[2].x = (1.-4*t)*_x21.y/_det;
   _dsh[2].y = (4*t-1.)*_x21.x/_det;
   _dsh[3].x = 4*(_x31.y*(1-t-2*s)+_x21.y*s)/_det;
   _dsh[3].y = 4*(_x31.x*(2*s+t-1)-_x21.x*s)/_det;
   _dsh[4].x = 4*(t*_x31.y-s*_x21.y)/_det;
   _dsh[4].y = 4*(s*_x21.x-t*_x31.x)/_det;
   _dsh[5].x = 4*(_x21.y*(2*t+s-1)-_x31.y*t)/_det;
   _dsh[5].y = 4*(_x31.x*t+_x21.x*(1-s-2*t))/_det;
}


void Triang6S::atMidEdges(std::vector<Point<real_t> >& dsh,
                          std::vector<real_t>&         w)
{
   dsh.resize(18);
   w.resize(3);
   w[0] = w[1] = w[2] = OFELI_THIRD*_area;
   real_t s[3] = { 0.5,0.5,0.0 }, t[3] = { 0.0, 0.5,0.5 };
   for (size_t k=0; k<3; ++k) {
      setLocal(s[k],t[k]);
      for (size_t i=0; i<6; ++i)
         dsh[3*i+k] = _dsh[i];
   }
}


vector<Point<real_t> > Triang6S::DSh() const
{
  return _dsh;
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
