/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2026 Rachid Touzani

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

         Implementation of class 'Tetra4' for Tetrahedral 4-Node Element

  ==============================================================================*/


#include "shape_functions/Tetra4.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "linear_algebra/LocalVect_impl.h"
#include "linear_algebra/LocalMatrix_impl.h"

namespace OFELI {

Tetra4::Tetra4()
{
   _sh.resize(4);
   _dsh.resize(4);
   _dshl.resize(4);
   _el = nullptr;
}

    
Tetra4::Tetra4(const Element* el)
{
   if (el->getNbNodes() != 4)
     throw OFELIException("Tetra4::Tetra4(Element *): Illegal number of element nodes: "
                          + std::to_string(el->getNbNodes()));
   set(el);
}


void Tetra4::set(const Element* el)
{
   _label = el->n();
   _sh.resize(4);
   _dsh.resize(4);
   _dshl.resize(4);
   for (size_t i=0; i<4; i++) {
      Node *node = (*el)(i+1);
      _x.push_back(node->getCoord());
      _node.push_back(node->n());
   }
   CalculateShape();
   _c = 0.25*(_x[0] + _x[1] + _x[2] + _x[3]);
   _el = el;
}


void Tetra4::CalculateShape()
{
   LocalMatrix<real_t,3,3> J, IJ;
   J = 0;
   _dshl[0].x = _dshl[0].y = _dshl[0].z = -1.0;
   _dshl[1].x = _dshl[2].y = _dshl[3].z =  1.0;
   _dshl[1].y = _dshl[1].z = _dshl[2].x =  0.0;
   _dshl[2].z = _dshl[3].x = _dshl[3].y =  0.0;
   for (size_t k=0; k<4; k++) {
      J(1,1) += _dshl[k].x*_x[k].x;
      J(1,2) += _dshl[k].y*_x[k].x;
      J(1,3) += _dshl[k].z*_x[k].x;
      J(2,1) += _dshl[k].x*_x[k].y;
      J(2,2) += _dshl[k].y*_x[k].y;
      J(2,3) += _dshl[k].z*_x[k].y;
      J(3,1) += _dshl[k].x*_x[k].z;
      J(3,2) += _dshl[k].y*_x[k].z;
      J(3,3) += _dshl[k].z*_x[k].z;
   }
   for (size_t i=1; i<=3; ++i) {
      size_t j = i%3 + 1;
      size_t k = j%3 + 1;
      IJ(i,i) = J(j,j)*J(k,k) - J(j,k)*J(k,j);
      IJ(j,i) = J(j,k)*J(k,i) - J(j,i)*J(k,k);
      IJ(i,j) = J(k,j)*J(i,k) - J(i,j)*J(k,k);
   }

   _det = J(1,1)*IJ(1,1) + J(2,1)*IJ(1,2) + J(3,1)*IJ(1,3);
   if (_det < 0.0)
      throw OFELIException("Tetra4::set(Element *): Negative determinant of jacobian");
   if (_det == 0.0)
      throw OFELIException("Tetra4::set(Element *): Determinant of jacobian is null");
   real_t c = 1./_det;
   for (size_t i=0; i<4; ++i) {
      _dsh[i].x = c*(IJ(1,1)*_dshl[i].x + IJ(2,1)*_dshl[i].y + IJ(3,1)*_dshl[i].z);
      _dsh[i].y = c*(IJ(1,2)*_dshl[i].x + IJ(2,2)*_dshl[i].y + IJ(3,2)*_dshl[i].z);
      _dsh[i].z = c*(IJ(1,3)*_dshl[i].x + IJ(2,3)*_dshl[i].y + IJ(3,3)*_dshl[i].z);
   }
}


vector<Point<real_t> > Tetra4::DSh() const
{
   return _dsh;
}


real_t Tetra4::Sh(size_t        i,
                  Point<real_t> s) const
{
   switch (i) {
      case 1:
         return (1.-s.x-s.y-s.z);

      case 2:
         return s.x;

      case 3:
         return s.y;

      case 4:
         return s.z;
   }
   return 0;
}


Point<real_t> Tetra4::getRefCoord(const Point<real_t>& x) const
{
   Point<real_t> s;
   s.x = (_x[2].y*_x[3].z-_x[2].y*_x[0].z-_x[0].y*_x[3].z-_x[3].y*_x[2].z+_x[3].y*_x[0].z+_x[0].y*_x[2].z) * (x.x-_x[0].x) +
         (-_x[2].x*_x[3].z+_x[2].x*_x[0].z+_x[0].x*_x[3].z+_x[3].x*_x[2].z-_x[3].x*_x[0].z-_x[0].x*_x[2].z) * (x.y-_x[0].y) +
         (_x[2].x*_x[3].y-_x[2].x*_x[0].y-_x[0].x*_x[3].y-_x[3].x*_x[2].y+_x[3].x*_x[0].y+_x[0].x*_x[2].y) * (x.z-_x[0].z);
   s.y = (-_x[1].y*_x[3].z+_x[1].y*_x[0].z+_x[0].y*_x[3].z+_x[3].y*_x[1].z-_x[3].y*_x[0].z-_x[0].y*_x[1].z) * (x.x-_x[0].x) +
         (_x[1].x*_x[3].z-_x[1].x*_x[0].z-_x[0].x*_x[3].z-_x[3].x*_x[1].z+_x[3].x*_x[0].z+_x[0].x*_x[1].z) * (x.y-_x[0].y) +
         (-_x[1].x*_x[3].y+_x[1].x*_x[0].y+_x[0].x*_x[3].y+_x[3].x*_x[1].y-_x[3].x*_x[0].y-_x[0].x*_x[1].y) * (x.z-_x[0].z);
   s.z = (_x[1].y*_x[2].z-_x[1].y*_x[0].z-_x[0].y*_x[2].z-_x[2].y*_x[1].z+_x[2].y*_x[0].z+_x[0].y*_x[1].z) * (x.x-_x[0].x) +
         (-_x[1].x*_x[2].z+_x[1].x*_x[0].z+_x[0].x*_x[2].z+_x[2].x*_x[1].z-_x[2].x*_x[0].z-_x[0].x*_x[1].z) * (x.y-_x[0].y) +
         (_x[1].x*_x[2].y-_x[1].x*_x[0].y-_x[0].x*_x[2].y-_x[2].x*_x[1].y+_x[2].x*_x[0].y+_x[0].x*_x[1].y) * (x.z-_x[0].z);
   return s/_det;
}


real_t Tetra4::getMaxEdgeLength() const
{
   real_t h = 0;
   for (size_t i=0; i<4; i++)
      h = std::max(h,Distance(_x[i],_x[(i+1)%4]));
   return h;
}


real_t Tetra4::getMinEdgeLength() const
{
   real_t h = Distance(_x[0],_x[1]);
   for (size_t i=1; i<4; i++)
      h = std::min(h,Distance(_x[i],_x[(i+1)%4]));
   return h;
}


bool Tetra4::isIn(const Point<real_t>& x)
{
   Point<real_t> s = getRefCoord(x);
   if (s.x>=0 && s.x<=1 && s.y>=0 && s.y<=1 && s.z>=0 && s.x+s.y+s.z<=1)
      return true;
   else
      return false;
}


real_t Tetra4::getInterpolate(const Point<real_t>&       x,
                              const LocalVect<real_t,4>& v)
{
   Point<real_t> s = getRefCoord(x);
   return (Sh(1,s)*v(1)+Sh(2,s)*v(2)+Sh(3,s)*v(3)+Sh(1,4)*v(4));
}


Point<real_t> Tetra4::EdgeSh(size_t        k,
                             Point<real_t> s)
{
   size_t i=k, j=k%4+1;
   return Sh(j,s)*_dsh[i-1]-Sh(i,s)*_dsh[j-1];
}


Point<real_t> Tetra4::CurlEdgeSh(size_t k)
{
   size_t i=k, j=k%3+1;
   real_t x = 2*(_dsh[j-1].y*_dsh[i-1].z - _dsh[j-1].z*_dsh[i-1].y);
   real_t y = 2*(_dsh[j-1].z*_dsh[i-1].x - _dsh[j-1].x*_dsh[i-1].z);
   real_t z = 2*(_dsh[j-1].x*_dsh[i-1].y - _dsh[j-1].y*_dsh[i-1].x);
   return Point<real_t>(x,y,z);
}

} /* namespace OFELI */
