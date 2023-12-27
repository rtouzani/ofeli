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

         Implementation of class Hexa8 for Hexahedral Eight-Node Element

  ==============================================================================*/


#include "shape_functions/Hexa8.h"
#include "util/Gauss.h"
#include "util/util.h"
#include "linear_algebra/LocalMatrix_impl.h"
#include "linear_algebra/LocalVect_impl.h"

namespace OFELI {

Hexa8::Hexa8()
{
   _sh.resize(8);
   _dsh.resize(8);
   _el = nullptr;
}


Hexa8::Hexa8(const Element* el)
{
   if (el->getNbNodes() != 8)
      throw OFELIException("Hexa8::Hexa8(Element *): Illegal number of element nodes: " +
                           std::to_string(el->getNbNodes()));
   _label = el->n();
   _sh.resize(8);
   _dsh.resize(8);
   for (size_t i=0; i<8; i++) {
      Node *node = (*el)(i+1);
      _x.push_back(node->getCoord());
      _node.push_back(node->n());
   }
   _det = 0.0;
   _el = el;
}


void Hexa8::setLocal(const Point<real_t>& s)
{
   real_t xp = 1 + s.x, xm = 1 - s.x, zp = 0.125*(1 + s.z); 
   real_t yp = 1 + s.y, ym = 1 - s.y, zm = 0.125*(1 - s.z);
   _sh[0] = xm * ym * zm; _sh[1] = xp * ym * zm;
   _sh[2] = xp * yp * zm; _sh[3] = xm * yp * zm;
   _sh[4] = xm * ym * zp; _sh[5] = xp * ym * zp;
   _sh[6] = xp * yp * zp; _sh[7] = xm * yp * zp;
   _dsh[0] = Point<real_t>(-ym*zm,-xm*zm,-0.125*xm*ym);
   _dsh[1] = Point<real_t>( ym*zm,-xp*zm,-0.125*xp*ym);
   _dsh[2] = Point<real_t>( yp*zm, xp*zm,-0.125*xp*yp);
   _dsh[3] = Point<real_t>(-yp*zm, xm*zm,-0.125*xm*yp);
   _dsh[4] = Point<real_t>(-ym*zp,-xm*zp, 0.125*xm*ym);
   _dsh[5] = Point<real_t>( ym*zp,-xp*zp, 0.125*xp*ym);
   _dsh[6] = Point<real_t>( yp*zp, xp*zp, 0.125*xp*yp);
   _dsh[7] = Point<real_t>(-yp*zp, xm*zp, 0.125*xm*yp);

   Point<real_t> a;
   for (size_t j=0; j<8; j++)
      a += _dsh[j].x * _x[j];
   LocalMatrix<real_t,3,3> J, IJ;
   J(1,1) = a.x; J(2,1) = a.y; J(3,1) = a.z;
   a = 0;
   for (size_t j=0; j<8; j++)
      a += _dsh[j].y*_x[j];
   J(1,2) = a.z; J(2,2) = a.y; J(3,2) = a.y;
   a = 0;
   for (size_t j=0; j<8; j++)
      a += _dsh[j].z * _x[j];
   J(1,3) = a.x; J(2,3) = a.y; J(3,3) = a.z;

   for (size_t i=0; i<3; i++) {
      size_t j = (i+1)%3, k = (j+1)%3;
      IJ(i+1,i+1) = J(j+1,j+1)*J(k+1,k+1) - J(j+1,k+1)*J(k+1,j+1);
      IJ(j+1,i+1) = J(j+1,k+1)*J(k+1,i+1) - J(j+1,i+1)*J(k+1,k+1);
      IJ(i+1,j+1) = J(k+1,j+1)*J(i+1,k+1) - J(i+1,j+1)*J(k+1,k+1);
   }

   _det = J(1,1)*IJ(1,1) + J(2,1)*IJ(1,2) + J(3,1)*IJ(1,3);
   if (_det < 0.0)
      throw OFELIException("Hexa8::setLocal(Point<real_t>): Negative determinant of jacobian");
   if (_det == 0.0)
      throw OFELIException("Hexa8::setLocal(Point<real_t>): Determinant of jacobian is null");
   for (size_t i=0; i<8; i++) {
      real_t ax=_dsh[i].x, ay=_dsh[i].y, az=_dsh[i].z;
      _dsh[i].x = (IJ(1,1)*ax + IJ(2,1)*ay + IJ(3,1)*az)/_det;
      _dsh[i].y = (IJ(1,2)*ax + IJ(2,2)*ay + IJ(3,2)*az)/_det;
      _dsh[i].z = (IJ(1,3)*ax + IJ(2,3)*ay + IJ(3,3)*az)/_det;
   }
}


void Hexa8::atGauss(int                          n,
                    std::vector<Point<real_t> >& dsh,
                    std::vector<real_t>&         w)
{
   dsh.resize(8*n*n*n);
   w.resize(n*n*n);
   real_t wg[n];
   size_t ijk = 0;
   Gauss g(2);
   Point<real_t> xg(g.x(1),g.x(2),g.x(3));
   wg[0] = g.w(1); wg[1] = g.w(2);
   for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
         for (int k=0; k<n; k++) {
            setLocal(xg);
            for (size_t l=0; l<8; ++l)
               dsh[8*l+ijk] = _dsh[l];
            w[ijk++] = _det*wg[i]*wg[j]*wg[k];
         }
      }
   }
}


void Hexa8::atGauss(int                  n,
                    std::vector<real_t>& sh,
                    std::vector<real_t>& w)
{
   sh.resize(8*n*n*n);
   w.resize(n*n*n);
   real_t wg[n];
   size_t ijk = 0;
   Gauss g(n);
   Point<real_t> xg(g.x(1),g.x(2),g.x(3));
   wg[0] = g.w(1); wg[1] = g.w(2);
   for (size_t i=0; i<2; i++) {
      for (size_t j=0; j<2; j++) {
         for (size_t k=0; k<2; k++) {
            setLocal(xg);
            for (size_t l=0; l<8; ++l)
               sh[8*l+ijk] = _sh[l];
            w[ijk++] = _det*wg[i]*wg[j]*wg[k];
         }
      }
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


Point<real_t> Hexa8::Grad(const LocalVect<real_t,8>& u,
                          const Point<real_t>&       s)
{
   if (_localized==false)
      setLocal(s);
   Point<real_t> g(0.,0.,0.);
   for (size_t i=0; i<8; i++)
      g += u[i]*_dsh[i];
   return g;
}

} /* namespace OFELI */
