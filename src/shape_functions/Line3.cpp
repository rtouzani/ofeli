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

            Implementation of class 'Line3' for Line 2-D Finite Element

  ==============================================================================*/


#include "shape_functions/Line3.h"
#include "util/macros.h"
#include "util/util.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "linear_algebra/Point.h"

namespace OFELI {

string itos(int i);

Line3::Line3()
{
   _sh.resize(3);
   _node.resize(3);
   _x.resize(3);
   _dsh.resize(3);
   _dshl.resize(3);
   _el = NULL;
   _sd = NULL;
}


Line3::Line3(const Element* el)
{
   if (el->getNbNodes() != 3)
      throw OFELIException("Line3::Line3(Element *): Illegal number of element nodes: " +
                           itos(el->getNbNodes()));
   _sh.resize(3);
   _node.resize(3);
   _x.resize(3);
   _dsh.resize(3);
   _dshl.resize(3);
   _label = el->n();
   for (size_t i=0; i<3; i++) {
      Node *node = (*el)(i+1);
      _x[i] = node->getCoord();
      _node[i] = node->n();
   }
   Point<real_t> dxdl = (_x[2]-_x[0]);
   _det = 0.5*dxdl.Norm();
   if (_det == 0.0)
      throw OFELIException("Line3::Line3(Element *): Determinant of jacobian is null");
   _el = el;
   _sd = NULL;
}


Line3::Line3(const Side *sd)
{
   if (sd->getNbNodes() != 3)
      throw OFELIException("Line3::Line3(Side *): Illegal number of side nodes: " +
                           itos(sd->getNbNodes()));
   _sh.resize(3);
   _node.resize(3);
   _x.resize(3);
   _dsh.resize(3);
   _dshl.resize(3);
   _label = sd->n();
   for (size_t i=0; i<3; i++) {
      Node *node = (*sd)(i+1);
      _x[i] = node->getCoord();
      _node[i] = node->n();
   }
   Point<real_t> dxdl = _x[2] - _x[0];
   _det = 0.5*dxdl.Norm();
   if (_det == 0.0)
      throw OFELIException("Line3::Line3(Side *): Determinant of jacobian is null");
}


void Line3::setLocal(real_t s)
{
   _sh[0] = -0.5*(1.-s)*s;
   _sh[1] = 1. - s*s;
   _sh[2] = 0.5*(1.+s)*s;
   _dshl[0].x = -0.5 + s;
   _dshl[1].x = -2.0*s;
   _dshl[2].x =  0.5 + s;
   Point<real_t> dxdl = (_x[2] - _x[0]);
   _det = 0.5*dxdl.Norm();
   for (size_t i=0; i<3; i++)
      _dsh[i] = _dshl[i]/_det;
}

} /* namespace OFELI */
