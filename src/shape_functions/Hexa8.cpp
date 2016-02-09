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

          Implementation of class Hexa8 for Hexahedral Eight-Node Element

  ==============================================================================*/


#include "shape_functions/Hexa8.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "util/Gauss.h"
#include "linear_algebra/LocalVect.h"
#include "linear_algebra/LocalMatrix.h"
#include "linear_algebra/Point.h"

namespace OFELI {

Hexa8::Hexa8()
{
   _sh.resize(8);
   _node.resize(8);
   _x.resize(8);
   _dsh.resize(8);
   _dshl.resize(8);
   _el = NULL;
}


Hexa8::Hexa8(const Element* el)
{
   try {
      if (el->getNbNodes() != 8)
         THROW_RT("Hexa8(Element *): Illegal number of element nodes: " + itos(el->getNbNodes()));
   }
   CATCH("Hexa8");
   _label = el->n();
   _sh.resize(8);
   _node.resize(8);
   _x.resize(8);
   _dsh.resize(8);
   _dshl.resize(8);
   for (size_t i=0; i<8; i++) {
      Node *node = el->getPtrNode(i+1);
      _x[i] = node->getCoord();
      _node[i] = node->n();
   }
   _det = 0.0;
   _el = el;
}


void Hexa8::setLocal(const Point<real_t>& s)
{
   real_t xm, ym, zm, xp, yp, zp;
   size_t i, j, k;
   LocalMatrix<real_t,3,3> J, IJ;

   xp = 1 + s.x; xm = 1 - s.x; zp = 0.125*(1 + s.z); 
   yp = 1 + s.y; ym = 1 - s.y; zm = 0.125*(1 - s.z);
   _sh[0] = xm * ym * zm;
   _sh[1] = xp * ym * zm;
   _sh[2] = xp * yp * zm;
   _sh[3] = xm * yp * zm;
   _sh[4] = xm * ym * zp;
   _sh[5] = xp * ym * zp;
   _sh[6] = xp * yp * zp;
   _sh[7] = xm * yp * zp;
   _dshl[0] = Point<real_t>(-ym*zm,-xm*zm,-0.125*xm*ym);
   _dshl[1] = Point<real_t>( ym*zm,-xp*zm,-0.125*xp*ym);
   _dshl[2] = Point<real_t>( yp*zm, xp*zm,-0.125*xp*yp);
   _dshl[3] = Point<real_t>(-yp*zm, xm*zm,-0.125*xm*yp);
   _dshl[4] = Point<real_t>(-ym*zp,-xm*zp, 0.125*xm*ym);
   _dshl[5] = Point<real_t>( ym*zp,-xp*zp, 0.125*xp*ym);
   _dshl[6] = Point<real_t>( yp*zp, xp*zp, 0.125*xp*yp);
   _dshl[7] = Point<real_t>(-yp*zp, xm*zp, 0.125*xm*yp);

   Point<real_t> a;
   for (j=0; j<8; j++)
      a += _dshl[j].x * _x[j];
   J(1,1) = a.x; J(2,1) = a.y; J(3,1) = a.z;
   a = 0;
   for (j=0; j<8; j++)
      a += _dshl[j].y * _x[j];
   J(1,2) = a.z; J(2,2) = a.y; J(3,2) = a.y;
   a = 0;
   for (j=0; j<8; j++)
      a += _dshl[j].z * _x[j];
   J(1,3) = a.x; J(2,3) = a.y; J(3,3) = a.z;

   for (i=0; i<3; i++) {
      j = (i+1)%3; k = (j+1)%3;
      IJ(i+1,i+1) = J(j+1,j+1) * J(k+1,k+1) - J(j+1,k+1) * J(k+1,j+1);
      IJ(j+1,i+1) = J(j+1,k+1) * J(k+1,i+1) - J(j+1,i+1) * J(k+1,k+1);
      IJ(i+1,j+1) = J(k+1,j+1) * J(i+1,k+1) - J(i+1,j+1) * J(k+1,k+1);
   }

   _det = J(1,1)*IJ(1,1) + J(2,1)*IJ(1,2) + J(3,1)*IJ(1,3);
   try {
      if (_det < 0.0)
         THROW_RT("setLocal(Point<real_t>): Negative determinant of jacobian");
   }
   CATCH("Hexa8");
   try {
      if (_det == 0.0)
         THROW_RT("setLocal(Point<real_t>): Determinant of jacobian is null");
   }
   CATCH("Hexa8");

   real_t c = 1./_det;
   for (i=0; i<8; i++)
      _dsh[i] = c*(IJ(1,1)*_dshl[i] + IJ(2,1)*_dshl[i] + IJ(3,1)*_dshl[i]);
}


void Hexa8::atGauss1(LocalVect<Point<real_t>,8>& dsh,
                     real_t&                     w)
{
   setLocal(Point<real_t> (0.,0.,0.));
   w = 8*_det;
   for (size_t l=0; l<8; l++)
      dsh[l] = _dsh[l];
}


void Hexa8::atGauss2(LocalMatrix<Point<real_t>,8,8>& dsh,
                     LocalVect<real_t,8>&            w)
{
   Point<real_t> xg;
   real_t wg[2];
   size_t ijk = 1;
   Gauss g(2);
   xg.x = g.x(1); xg.y = g.x(2); xg.z = g.x(3);
   wg[0] = g.w(1);
   wg[1] = g.w(2);
   for (size_t i=0; i<2; i++)
      for (size_t j=0; j<2; j++)
         for (size_t k=0; k<2; k++) {
            setLocal(xg);
            w(ijk) = _det*wg[i]*wg[j]*wg[k];
            for (size_t l=1; l<=8; l++)
               dsh(l,ijk) = _dsh[l-1];
            ijk++;
         }
}


real_t Hexa8::getMaxEdgeLength() const
{
   real_t h = 0;
   for (size_t i=0; i<8; i++)
      h = std::max(h,Distance(_x[i],_x[(i+1)%8]));
   return h;
}


real_t Hexa8::getMinEdgeLength() const
{
   real_t h = Distance(_x[0],_x[1]);
   for (size_t i=1; i<8; i++)
      h = std::min(h,Distance(_x[i],_x[(i+1)%8]));
   return h;
}

} /* namespace OFELI */
