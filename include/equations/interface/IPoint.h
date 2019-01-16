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

                               Class IPoint 

  ==============================================================================*/

#ifndef _IPOINT_H__
#define _IPOINT_H__

#include <vector>
#include "util/util.h"
#include <iostream>
using std::ostream;

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

class IPoint
{
 private:

 public:

   int i, j, k, sgn;
   real_t val;

   IPoint();

   IPoint(int ix);

   IPoint(int ix,
          int iy,
          real_t v=0.);

   IPoint(int ix,
          int iy,
          int iz,
          real_t v=0.);

   IPoint(const IPoint& p);

   ~IPoint() {}

   bool operator>(IPoint& p) { return (val > p.val); }

   IPoint &operator=(const IPoint& p);

   IPoint &operator=(const real_t& val);

   bool operator==(const IPoint& p);

   IPoint& operator+=(const IPoint& p);

   void set(int ix, int iy) { i = ix; j = iy; }

   void set(int ix, int iy, int iz) { i = ix; j = iy; k = iz; }

   void getNeighbour(vector<IPoint>& neig);
};

bool operator==(const IPoint& a,
                const IPoint& b);

IPoint operator*(const int&    ix,
                 const IPoint& p);

bool operator<(const IPoint& a,
               const IPoint& b);
 
void getDisplacement(vector<IPoint>& NXi,
                     bool            three_D=false);

std::ostream & operator<<(std::ostream& s,
                          const IPoint& a);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif
