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

                         Implementation of class Muscl

       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/


#include "equations/cl/Muscl.h"
#include "mesh/Mesh.h"
#include "util/util.h"
#include <algorithm>

using std::min;
using std::max;

namespace OFELI {

Muscl::Muscl(Mesh& m)
{
   _theMesh = &m;
   _verbose = 1;
   _nb_nodes = _theMesh->getNbNodes();
   _nb_sides = _theMesh->getNbSides();
   _nb_elements = _theMesh->getNbElements();
   _solid_zone = false;
}


real_t Muscl::minmod(real_t val_plus,
                     real_t val_minus)
{
   if ((val_plus * val_minus<0.))
      return 0.;
//     ull signP = *((ull*)(&val_plus))  & 0x8000000000000000L;
//     ull signM = *((ull*)(&val_minus)) & 0x8000000000000000L;
//     if (signP ^ signM) return (0.0);
   real_t r = std::abs(val_plus)/(std::abs(val_minus) + OFELI_EPSMCH);
//     ull absPI = *((ull*)(&val_plus))  & 0x7FFFFFFFFFFFFFFFL;
//     ull absMI = *((ull*)(&val_moins)) & 0x7FFFFFFFFFFFFFFFL;
//     real_t absP = *((real_t * )(&absPI));
//     real_t absM = *((real_t * )(&absMI));
//     real_t r = absP / (absM + TOUT_PETIT);
   return min(r,1.);
//     ull rull = *((ull*)(&r));
//     return ( (rull & 0x7FF0000000000000L) >= 0x3FF0000000000000L)? r : 1.0 ;
}


real_t Muscl::superbee(real_t val_plus,
                       real_t val_minus)
{
   if (val_plus*val_minus<=0.)
      return 0.;
   real_t r = val_plus/val_minus;
   if (r<0.5)
      return (2.*r);
   if (r<1.)
      return 1.;
   if (r<2.)
      return r;
   return 2.;
}


real_t Muscl::vanleer(real_t val_plus,
                      real_t val_minus)
{
   if (val_plus*val_minus<=0.)
      return 0.;
   real_t r = val_plus/val_minus;
   return (2.*r/(1.+r));
}


real_t Muscl::vanalbada(real_t val_plus,
                        real_t val_minus)
{
   real_t r = val_plus/(std::abs(val_minus)+1.e-15)*Sgn(val_minus);
   if (r<-1)
      return 0.;
   else
      return (r*r+r)/(1+r*r);
}


real_t Muscl::m_limiter(real_t a,
                        real_t lim1,
                        real_t lim2)
{
   if (a<min(lim1,lim2))
      return min(lim1,lim2);
   else if (a>max(lim1,lim2))
      return max(lim1,lim2);
   else
      return a;
}


real_t Muscl::m_limiter2(real_t a,
                         real_t s,
                         real_t p,
                         real_t d)
{
   real_t coeff = 1;
   real_t Uq = s + p*coeff*d;
   if (a<min(s,Uq))
      return min(s,Uq);
   else if (a>max(s,Uq))
      return max(s,Uq);
   else
      return a;
}

} /* namespace OFELI */
