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

                               Class IPoint 

  ==============================================================================*/

#ifndef _IPOINT_H__
#define _IPOINT_H__

#include "linear_algebra/Vect.h"
#include "linear_algebra/LocalVect.h"
#include "util/util.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace OFELI {

class IPoint
{
 private:
   int _i, _j, _k, _sign;
   real_t _value;

 public:

   IPoint();

   IPoint(int i);

   IPoint(size_t x,
          size_t y,
          real_t v=0.);

   IPoint(size_t x,
          size_t y,
          size_t z,
          real_t v=0.);

   IPoint(const IPoint& p);

   ~IPoint() {}

   bool operator>(IPoint& p) { return (_value > p._value); }

   IPoint &operator=(const IPoint& p);

   IPoint &operator=(const real_t& val);

   bool operator==(const IPoint& p);

   IPoint& operator+(const IPoint& p);

   void setValue(real_t v) { _value = Abs(v); }

   void setX(int ii) { _i = ii; }

   void setY(int jj)  { _j = jj; }

   void setZ(int kk) { _k = kk; }

   void setSgn(int sg) { _sign = sg; }

   real_t getValue() const { return _value; }

   int getX() const { return _i; }

   int getY() const { return _j; }

   int getZ() const { return _k; }		

   int getSgn() const { return _sign; }

   void GenerateNeighbour(LocalVect<IPoint,6>& Neighbour);
};

IPoint operator*(const int&    i,
                 const IPoint& p);

void GenerateDisplacement(LocalVect<IPoint,6>& NXi,
                          bool                 three_D=false);

} /* namespace OFELI */

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif
