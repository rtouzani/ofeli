/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2018 Rachid Touzani

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

         Implementation of class 'Quad4' for Quadrilateral 4-Node Element

  ==============================================================================*/


#include "shape_functions/Quad4.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "linear_algebra/Point.h"
#include "util/util.h"

namespace OFELI {

Quad4::Quad4()
{
   _sh.resize(4);
   _node.resize(4);
   _x.resize(4);
   _dsh.resize(4);
   _dshl.resize(4);
   _el = NULL;
   _sd = NULL;
   _localized = false;
}


Quad4::Quad4(const Element* el)
{
   if (el->getNbNodes() != 4)
      throw OFELIException("Quad4::Quad4(Element *): Illegal number of element nodes: " +
                           itos(el->getNbNodes()));
   _sh.resize(4);
   _node.resize(4);
   _x.resize(4);
   _dsh.resize(4);
   _dshl.resize(4);
   for (size_t i=0; i<4; i++) {
      Node *node = (*el)(i+1);
      _x[i] = node->getCoord();
      _node[i] = node->n();
   }
   _det = 0.0;
   _label = el->n();
   _el = el;
   _sd = NULL;
   _localized = false;
}


Quad4::Quad4(const Side* side)
{
   if (side->getNbNodes() != 4)
      throw OFELIException("Quad4::Quad4(Side *): Illegal number of side nodes: " +
                           itos(side->getNbNodes()));
   _sh.resize(4);
   _node.resize(4);
   _x.resize(4);
   _dsh.resize(4);
   _dshl.resize(4);
   for (size_t i=0; i<4; i++) {
      Node *node = (*side)(i+1);
      _x[i] = node->getCoord();
      _node[i] = node->n();
   }
   _det = 0.0;
   _label = side->n();
   _el = NULL;
   _sd = side;
   _localized = false;
}


void Quad4::set(const Element* el)
{
   for (size_t i=0; i<4; i++) {
      Node *node = (*el)(i+1);
      _x[i] = node->getCoord();
      _node[i] = node->n();
   }
   _det = 0.0;
   _label = el->n();
   _el = el;
   _sd = NULL;
}


void Quad4::set(const Side* sd)
{
   for (size_t i=0; i<4; i++) {
      Node *node = (*sd)(i+1);
      _x[i] = node->getCoord();
      _node[i] = node->n();
   }
   _det = 0.0;
   _label = sd->n();
   _el = NULL;
   _sd = sd;
}


void Quad4::setLocal(const Point<real_t>& s)
{
   _localized = true;
   _sh[0] = 0.25*(1.-s.x)*(1.-s.y);
   _sh[1] = 0.25*(1.+s.x)*(1.-s.y);
   _sh[2] = 0.25*(1.+s.x)*(1.+s.y);
   _sh[3] = 0.25*(1.-s.x)*(1.+s.y);

   Point<real_t> dshl[4];
   dshl[0].x = -0.25*(1.-s.y); 
   dshl[1].x =  0.25*(1.-s.y);
   dshl[2].x =  0.25*(1.+s.y); 
   dshl[3].x = -0.25*(1.+s.y);
   dshl[0].y = -0.25*(1.-s.x); 
   dshl[3].y =  0.25*(1.-s.x);
   dshl[1].y = -0.25*(1.+s.x); 
   dshl[2].y =  0.25*(1.+s.x);
 
   real_t dxds=0., dyds=0., dxdt=0., dydt = 0.;
   for (size_t i=0; i<4; i++) {
      dxds += dshl[i].x*_x[i].x;
      dxdt += dshl[i].y*_x[i].x;
      dyds += dshl[i].x*_x[i].y;
      dydt += dshl[i].y*_x[i].y;
   }
   _det = dxds*dydt - dxdt*dyds;
   if (_det < 0.0)
      throw OFELIException("Quad4::setLocal(Point<real_t>): Negative determinant of jacobian");
   if (_det == 0.0)
      throw OFELIException("Quad4::setLocal(setLocal(Point<real_t>): Determinant of jacobian is null");
   real_t c = 1./_det;
   real_t dsdx =  c*dydt, dsdy = -c*dxdt, dtdx = -c*dyds, dtdy =  c*dxds;
   for (size_t i=0; i<4; i++) {
      _dsh[i].x = dsdx*dshl[i].x + dtdx*dshl[i].y;
      _dsh[i].y = dsdy*dshl[i].x + dtdy*dshl[i].y;
   }
}


real_t Quad4::getMaxEdgeLength() const
{
   real_t h = 0;
   for (size_t i=0; i<4; i++)
      h = std::max(h,Distance(_x[i],_x[(i+1)%4]));
   return h;
}


real_t Quad4::getMinEdgeLength() const
{
   real_t h = Distance(_x[0],_x[1]);
   for (size_t i=1; i<4; i++)
      h = std::min(h,Distance(_x[i],_x[(i+1)%4]));
   return h;
}


Point<real_t> Quad4::Grad(const LocalVect<real_t,4>& u,
                          const Point<real_t>&       s)
{
   if (_localized==false)
      setLocal(s);
   Point<real_t> g(0.,0.,0.);
   for (size_t i=1; i<=4; i++)
      g += u(i)*DSh(i);
   return g;
}

} /* namespace OFELI */
