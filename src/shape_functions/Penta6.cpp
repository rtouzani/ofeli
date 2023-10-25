/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani

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

         Implementation of class 'Penta6' for Pentahedral 6-Node Element

  ==============================================================================*/


#include "shape_functions/Penta6.h"
#include "util/util.h"
#include "linear_algebra/LocalMatrix_impl.h"

namespace OFELI {

Penta6::Penta6()
{
   _sh.resize(6);
   _dshl.resize(6);
   _dsh.resize(6);
   _localized = false;
}


Penta6::Penta6(const Element* el)
{
   if (el->getNbNodes() != 6)
      throw OFELIException("Penta6::Penta6(Element *): Illegal number of element nodes: " +
                           std::to_string(el->getNbNodes()));
   _sh.resize(6);
   _dshl.resize(6);
   _dsh.resize(6);
   for (size_t i=0; i<6; i++) {
      Node *node = (*el)(i+1);
      _x.push_back(node->getCoord());
      _node.push_back(node->n());
   }
   _det = 0.0;
   _label = el->n();
   _el = el;
   _localized = false;
}


void Penta6::set(const Element* el)
{
   for (size_t i=0; i<6; i++) {
      Node *node = (*el)(i+1);
      _x.push_back(node->getCoord());
      _node.push_back(node->n());
   }
   _det = 0.0;
   _label = el->n();
   _el = el;
}


void Penta6::setLocal(const Point<real_t>& s)
{
   _localized = true;

// Shape functions
   _sh[0] = s.x         *(1.-s.z);
   _sh[1] = s.y         *(1.-s.z);
   _sh[2] = (1.-s.x-s.y)*(1.-s.z);
   _sh[3] = s.x         *s.z;
   _sh[4] = s.y         *s.z;
   _sh[5] = (1.-s.x-s.y)*s.z;

// Local derivatives
   _dshl[0] = Point<real_t>(1.-s.z,0.,-s.x); 
   _dshl[1] = Point<real_t>(0.,1.-s.z,-s.y); 
   _dshl[2] = Point<real_t>(s.z-1.,s.z-1.,s.x+s.y-1.);
   _dshl[3] = Point<real_t>(s.z,0.,s.x); 
   _dshl[4] = Point<real_t>(0.,s.z,s.y); 
   _dshl[5] = Point<real_t>(-s.z,-s.z,1.-s.x-s.y);

// Jacobian matrix (dxds(i,j) := dx(i)/ds(j))
   LocalMatrix<real_t,3,3> dxds, dsdx;
   dxds = 0;
   for (size_t i=0; i<6; i++) {
      dxds(1,1) += _dshl[i].x*_x[i].x;
      dxds(1,2) += _dshl[i].x*_x[i].y; 
      dxds(1,3) += _dshl[i].x*_x[i].z; 
      dxds(2,1) += _dshl[i].y*_x[i].x;
      dxds(2,2) += _dshl[i].y*_x[i].y; 
      dxds(2,3) += _dshl[i].y*_x[i].z; 
      dxds(3,1) += _dshl[i].z*_x[i].x;
      dxds(3,2) += _dshl[i].z*_x[i].y; 
      dxds(3,3) += _dshl[i].z*_x[i].z;
   }

// Inverse of jacobian
   for (size_t i=1; i<=3; i++) {
      size_t j = i%3 + 1;
      size_t k = j%3 + 1;
      dsdx(i,i) = dxds(j,j)*dxds(k,k) - dxds(k,j)*dxds(j,k); 
      dsdx(i,j) = dxds(k,j)*dxds(i,k) - dxds(i,j)*dxds(k,k); 
      dsdx(j,i) = dxds(j,k)*dxds(k,i) - dxds(j,i)*dxds(k,k); 
   }
   _det = dxds(1,1)*dsdx(1,1) + dxds(1,2)*dsdx(2,1) + dxds(1,3)*dsdx(3,1);
   if (_det == 0.0)
      throw OFELIException("Penta6::setLocal(Point<real_t>): Determinant of jacobian is null");
   _det = fabs(_det);
   real_t c = 1./_det;
   for (size_t i=0; i<6; i++) {
      _dsh[i].x = c*(dsdx(1,1)*_dshl[i].x + dsdx(1,2)*_dshl[i].y + dsdx(1,3)*_dshl[i].z);
      _dsh[i].y = c*(dsdx(2,1)*_dshl[i].x + dsdx(2,2)*_dshl[i].y + dsdx(2,3)*_dshl[i].z);
      _dsh[i].z = c*(dsdx(3,1)*_dshl[i].x + dsdx(3,2)*_dshl[i].y + dsdx(3,3)*_dshl[i].z);
   }
}


vector<Point<real_t> > Penta6::DSh() const
{
   return _dsh;
}


real_t Penta6::getMaxEdgeLength() const
{
   real_t h = 0;
   for (size_t i=0; i<6; i++)
      h = std::max(h,Distance(_x[i],_x[(i+1)%6]));
   return h;
}


real_t Penta6::getMinEdgeLength() const
{
   real_t h = Distance(_x[0],_x[1]);
   for (size_t i=1; i<6; i++)
      h = std::min(h,Distance(_x[i],_x[(i+1)%6]));
   return h;
}

} /* namespace OFELI */
