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

                           Implementation of class IPoint

  ==============================================================================*/

#include "equations/interface/IPoint.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace OFELI {

IPoint::IPoint(): _i(1), _j(1), _k(0)
{
   _value = 0;
   _sign = 1;
}


IPoint::IPoint(int i) : _i(i), _sign(Sgn(i))
{
}


IPoint::IPoint(size_t x,
               size_t y,
               real_t v)
       : _i(x), _j(y), _k(0)
{
   _value = Abs(v);
   _sign = Sgn(v);
}


IPoint::IPoint(size_t x,
               size_t y,
               size_t z,
               real_t v)
       : _i(x), _j(y), _k(z)
{
   _value = Abs(v);
   _sign = Sgn(v);
}


IPoint::IPoint(const IPoint& p)
{
   _i = p._i; 
   _j = p._j; 
   _k  = p._k;
   _value = p._value; 
   _sign = p._sign;
}


IPoint & IPoint::operator=(const IPoint& p)
{
   _i = p._i; 
   _j = p._j; 
   _k = p._k;
   _value = p._value; 
   _sign = p._sign;
   return *this;
}


IPoint &IPoint::operator=(const real_t& val)
{
   this->_value = val;
   return *this;
}


bool IPoint::operator==(const IPoint& p)
{
   return ((_i==p._i) && (_j==p._j) && (_k==p._k));
}


IPoint &IPoint::operator+(const IPoint& p)
{
   _i += p._i;
   _j += p._j;
   _k += p._k;
   return *this;
}


void IPoint::GenerateNeighbour(LocalVect<IPoint,6>& Neighbour)
{
   Neighbour[0]._i = _i - 1; Neighbour[0]._j = _j;
   Neighbour[1]._i = _i + 1; Neighbour[1]._j = _j;
   Neighbour[2]._i = _i;     Neighbour[2]._j = _j - 1;
   Neighbour[3]._i = _i;     Neighbour[3]._j = _j + 1;
   if (_k != 0) {
      Neighbour[0]._k = _k;
      Neighbour[1]._k = _k;
      Neighbour[2]._k = _k;
      Neighbour[3]._k = _k;
      Neighbour[4]._i = _i; 
      Neighbour[4]._j = _j;
      Neighbour[4]._k = _k - 1;
      Neighbour[5]._i = _i; 
      Neighbour[5]._j = _j;
      Neighbour[5]._k = _k + 1;
   }
}


void GenerateDisplacement(LocalVect<IPoint,6>& NXi,
                          bool                 three_D)
{
   NXi[0].setX(-1); NXi[0].setY( 0);
   NXi[1].setX( 1); NXi[1].setY( 0);
   NXi[2].setX( 0); NXi[2].setY(-1);
   NXi[3].setX( 0); NXi[3].setY( 1);

   if (three_D) {
      NXi[0].setZ( 0); NXi[1].setZ( 0);
      NXi[2].setZ( 0); NXi[3].setZ( 0);	
      NXi[4].setX( 0); NXi[4].setY( 0);
      NXi[4].setZ(-1); NXi[5].setX( 0);
      NXi[5].setY( 0); NXi[5].setZ( 1);
   }
}


IPoint operator*(const int&          i,
                       const IPoint& p)
{
   IPoint res(p);
   res.setX(i*res.getX());
   res.setY(i*res.getY());
   res.setZ(i*res.getZ());
   return res;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */
