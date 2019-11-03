/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

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

                         Implementation of class 'Gauss'

 ==============================================================================*/

#include "util/Gauss.h"
#include "linear_algebra/LocalVect_impl.h"

namespace OFELI {

Gauss::Gauss(size_t np)
      : _np(np)
{
   _triang = false;
   real_t p, dp, r;
   for (size_t i=0; i<=(_np-1)/2; i++) {
      r = cos(OFELI_PI*(i+0.75)/(np+0.5));
      while (1) {
         legendre(r, p, dp);
         if (fabs(p) <= 2.0*DBL_EPSILON*fabs(r*dp))
            break;
         r -= p/dp;
      }
      _x[i] = -r;
      _x[_np-i-1] = r;
      _w[_np-i-1] = _w[i] = 2.0/((1.0-r*r)*dp*dp);
   }
}


void Gauss::legendre(real_t  y,
                     real_t& p,
                     real_t& dp)
{
   real_t p0, p1, p2=0;
   if (_np==0) {
      p = 1.0;
      dp = 0.0;
      return;
   }
   if (_np==1) {
      p = y;
      dp = 1.0;
      return;
   }

   p0 = 1.0;
   p1 = y;
   int nn = int(_np*(_np+1)/2);
   for (size_t j=2; j<=_np; j++, p0=p1, p1=p2)
      p2 = ((2*j-1)*y*p1 - (j-1)*p0)/j;
   p = p2;
   if (y==1.0)
      dp = nn;
   else if (y==-1.0)
      dp = _np%2==1 ? nn : -nn;
   else
      dp = _np*(y*p1 - p0) / (y*y-1);
}


void Gauss::setTriangle(LocalVect<real_t,7>&        w,
                        LocalVect<Point<real_t>,7>& x)
{
   x(1) = Point<real_t>(0.,0.);
   x(2) = Point<real_t>(1.,0.);
   x(3) = Point<real_t>(0.,1.);
   w(1) = w(2) = w(3) = 3./60.;

   x(4) = Point<real_t>(0.5,0. );
   x(5) = Point<real_t>(0.5,0.5);
   x(6) = Point<real_t>(0. ,0.5);
   w(4) = w(5) = w(6) = 8./60.;

   x(7) = Point<real_t>(OFELI_THIRD,OFELI_THIRD);
   w(7) = 27./60.;
}

} /* namespace OFELI */
