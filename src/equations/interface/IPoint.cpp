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

                           Implementation of class IPoint

  ==============================================================================*/

#include "equations/interface/IPoint.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace OFELI {

IPoint::IPoint()
{
   i = j = 1, k = 0;
   val = 0;
   sgn = 1;
}


IPoint::IPoint(int ix)
{
   sgn = Sgn(ix);
   i = ix;
}


IPoint::IPoint(int    ix,
               int    iy,
               real_t v)
{
   i = ix, j = iy, k = 0;
   val = fabs(v);
   sgn = Sgn(v);
}


IPoint::IPoint(int    ix,
               int    iy,
               int    iz,
               real_t v)
{
   i = ix, j = iy, k = iz;
   val = fabs(v);
   sgn = Sgn(v);
}


IPoint::IPoint(const IPoint& p)
{
   i = p.i, j = p.j, k = p.k;
   val = p.val; 
   sgn = p.sgn;
}


IPoint & IPoint::operator=(const IPoint& p)
{
   i = p.i, j = p.j, k = p.k;
   val = p.val; 
   sgn = p.sgn;
   return *this;
}


IPoint &IPoint::operator=(const real_t& val)
{
   this->val = val;
   return *this;
}


bool IPoint::operator==(const IPoint& p)
{
   return (i==p.i && j==p.j && k==p.k);
}


IPoint &IPoint::operator+=(const IPoint& p)
{
   i += p.i;
   j += p.j; 
   k += p.k;
   return *this;
}


void IPoint::getNeighbour(vector<IPoint>& neig)
{
   neig[0].i = i-1; neig[0].j = j;
   neig[1].i = i+1; neig[1].j = j;
   neig[2].i = i;   neig[2].j = j-1;
   neig[3].i = i;   neig[3].j = j+1;
   if (k != 0) {
      neig[0].k = k; neig[1].k = k;
      neig[2].k = k; neig[3].k = k;
      neig[4].i = i; neig[4].j = j;
      neig[4].k = k-1; neig[5].i = i; 
      neig[5].j = j; neig[5].k = k+1;
   }
}


void getDisplacement(vector<IPoint>& NXi,
                     bool            three_D)
{
   NXi[0].i = -1; NXi[0].j =  0;
   NXi[1].i =  1; NXi[1].j =  0;
   NXi[2].i =  0; NXi[2].j = -1;
   NXi[3].i =  0; NXi[3].j =  1;

   if (three_D) {
      NXi[0].k =  0; NXi[1].k = 0;
      NXi[2].k =  0; NXi[3].k = 0;	
      NXi[4].i =  0; NXi[4].j = 0;
      NXi[4].k = -1; NXi[5].i = 0;
      NXi[5].j =  0; NXi[5].k = 1;
   }
}


bool operator==(const IPoint& a,
                const IPoint& b)
{ 
   return (a.i==b.i && a.j==b.j && a.k==b.k); 
}


bool operator<(const IPoint& a,
               const IPoint& b)
{
   if (a.i < b.i)
      return true;
   if (a.i==b.i && a.j<b.j)
      return true;
   return false;
}


IPoint operator*(const int&    ix,
                 const IPoint& p)
{
   IPoint res(p);
   res.i = ix*res.i;
   res.j = ix*res.j;
   res.k = ix*res.k;
   return res;
}


std::ostream & operator<<(std::ostream& s,
                          const IPoint& a)
{
   s << "( " << a.i << "," << a.j << "," << a.k << " - "
     << a.val << " - " << a.sgn << ")";
   return s;
}


#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */
