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

                         Implementation of class 'Side'

 ==============================================================================*/

#include "mesh/Side.h"
#include "mesh/Mesh.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Quad4.h"
#include "shape_functions/Tetra4.h"
#include "shape_functions/Hexa8.h"
#include "util/util.h"
#include "OFELIException.h"

using std::to_string;

namespace OFELI {

Side::Side()
{
   _nb_nodes = _nb_eq = 0;
   _label = 0;
   _nb_nodes = _nb_edges = 0;
   _neig_el = 0;
   _el[0] = _el[1] = nullptr;
   _on_boundary = -1;
   _code.resize(MAX_NBDOF_SIDE);
   clear(_code);
   _dof.resize(MAX_NBDOF_SIDE);
   _nb_dof = 1;
   _nb_childs = 0;
   _parent = nullptr;
   _level = 0;
   _active = true;
}


Side::Side(size_t        label,
           const string& shape)
{
   if (label<1)
      throw OFELIException("Side::Side(size_t,string): Illegal side label "+to_string(label));
   shape_index(shape);
   _label = label;
   _nb_nodes = _nb_edges = _nb_eq = 0;
   _neig_el = 0;
   _el[0] = _el[1] = nullptr;
   _on_boundary = 0;
   _code.resize(MAX_NBDOF_SIDE);
   clear(_code);
   _dof.resize(MAX_NBDOF_SIDE);
   _nb_childs = 0;
   _level = 0;
   _parent = nullptr;
   _active = true;
}


Side::Side(size_t label,
           int    shape)
{
   if (label<1)
      throw OFELIException("Side::Side(size_t,int): Illegal side label "+to_string(label));
   _shape = shape;
   _label = label;
   _nb_nodes = _nb_edges = _nb_eq = 0;
   _neig_el = 0;
   _el[0] = _el[1] = nullptr;
   _on_boundary = 0;
   _code.resize(MAX_NBDOF_SIDE);
   clear(_code);
   _dof.resize(MAX_NBDOF_SIDE);
   _nb_childs = 0;
   _parent = nullptr;
   _level = 0;
   _active = true;
}


Side::Side(const Side& sd)
{
   _shape = sd._shape;
   _nb_nodes = sd._nb_nodes;
   _nb_edges = sd._nb_edges;
   _label = sd._label;
   _nb_eq = 0;
   _neig_el = sd._neig_el;
   for (size_t i=0; i<_nb_nodes; i++)
      _node[i] = sd._node[i];
   _on_boundary = sd._on_boundary;
   _nb_dof = sd._nb_dof;
   _code.resize(_nb_dof);
   _dof.resize(_nb_dof);
   for (size_t j=0; j<_nb_dof; j++) {
      _code[j] = sd._code[j];
      _dof[j] = sd._dof[j];
   }
   _el[0] = sd._el[0];
   _el[1] = sd._el[1];
   for (size_t i=0; i<_nb_edges; i++)
      _ed[i] = sd._ed[i];
   _nb_childs = sd._nb_childs;
   _level = sd._level;
   _active = sd._active;
   _parent = sd._parent;
}


void Side::setNode(size_t i,
                   Node*  node)
{
   _node[i-1] = node;
}


void Side::Add(Node* node)
{
   if (!node)
      throw OFELIException("Side::Add(Node *): Trying to add an undefined node");
   _node[_nb_nodes++] = node;
   _nb_eq += node->getNbDOF();
}
   
   
void Side::Add(Edge* edge)
{
   if (edge==nullptr)
      throw OFELIException("Side::Add(Edge): Trying to add an undefined edge");
   _ed[_nb_edges++] = edge;
}
   

void Side::Replace(size_t label,
                   Node*  node)
{
   if (!node)
      throw OFELIException("Side::Replace(size_t,Node *): Trying to replace "
                           "an undefined node to size "+to_string(label));
   _node[label-1] = node;
}


void Side::setOnBoundary()
{
   _on_boundary = 1;
}


void Side::setChild(Side* sd)
{
   _child[_nb_childs++] = sd;
   sd->_level = _level + 1;
   sd->_parent = this;
   sd->_nb_dof = _nb_dof;
   for (size_t j=0; j<_nb_dof; j++) {
      sd->_code[j] = _code[j];
      sd->_dof[j] = _dof[j];
   }
   _active = false;
}


Side *Side::getChild(size_t i) const
{
   if (i > _nb_childs)
      throw OFELIException("Side::getChild(size_t): Number of children is "+to_string(i));
   return _child[i-1];
}


int Side::isOnBoundary() const
{
   _on_boundary = -1;
   if ((_el[0] && !_el[1]) || (!_el[0] && _el[1]))
      _on_boundary = 1;
   if (_el[0] && _el[1])
      _on_boundary = 0;
   return _on_boundary;
}


size_t Side::getNbVertices() const
{
   size_t n=0;
   switch (_shape) {
      case LINE:          n = 2; break;
      case TRIANGLE:      n = 3; break;
      case QUADRILATERAL: n = 4; break;
      case TETRAHEDRON:   n = 4; break;
      case PENTAHEDRON:   n = 6; break;
      case HEXAHEDRON:    n = 8; break;
   }
   return n;
}


Point<real_t> Side::getCenter() const
{
   switch (_shape) {

      case LINE:
         return 0.5*(_node[0]->getCoord()+_node[1]->getCoord());

      case TRIANGLE:
         return OFELI_THIRD*(_node[0]->getCoord()+_node[1]->getCoord()+_node[2]->getCoord());

      case QUADRILATERAL:
         return 0.25*(_node[0]->getCoord() + _node[1]->getCoord() +
                      _node[2]->getCoord() + _node[3]->getCoord());

      default:
         return Point<real_t>(0.,0.,0.);
   }
}


real_t Side::getMeasure() const
{
   Point<real_t> x[4];
   real_t m=0;

   if (_shape==LINE) {
      x[0] = _node[0]->getCoord();
      x[1] = _node[1]->getCoord();
      m = sqrt((x[1].x-x[0].x)*(x[1].x-x[0].x) + (x[1].y-x[0].y)*(x[1].y-x[0].y));
      if (m==0)
         throw OFELIException("Side::getMeasure(): Length of line "+to_string(_label)+" is null.");
   }

   else if (_shape==TRIANGLE) {
      for (size_t i=0; i<3; i++)
         x[i] = _node[i]->getCoord();
      real_t a = x[0].y*(x[1].z-x[2].z) - x[1].y*(x[0].z-x[2].z) + x[2].y*(x[0].z-x[1].z);
      real_t b = x[0].z*(x[1].x-x[2].x) - x[1].z*(x[0].x-x[2].x) + x[2].z*(x[0].x-x[1].x);
      real_t c = x[0].x*(x[1].y-x[2].y) - x[1].x*(x[0].y-x[2].y) + x[2].x*(x[0].y-x[1].y);
      m = 0.5*sqrt(a*a + b*b + c*c);
      if (m==0)
         throw OFELIException("Side::getMeasure(): Area of triangle "+to_string(_label)+" is null.");
      if (m<0)
         throw OFELIException("Side::getMeasure(): Area of triangle " + to_string(_label) + "is negative.");
   }

   else if (_shape==QUADRILATERAL) {
      for (size_t i=0; i<4; i++)
         x[i] = _node[i]->getCoord();
      Point<real_t> dshl[4];
      dshl[0].x = dshl[3].x = dshl[0].y = dshl[1].y = -0.25; 
      dshl[1].x = dshl[2].x = dshl[3].y = dshl[2].y =  0.25;
      real_t dxds=0., dxdt=0., dyds=0., dydt=0.;
      for (size_t i=0; i<4; i++) {
         dxds += dshl[i].x*x[i].x;
         dxdt += dshl[i].y*x[i].x;
         dyds += dshl[i].x*x[i].y;
         dydt += dshl[i].y*x[i].y;
      }
      m = dxds*dydt - dxdt*dyds;
      if (m==0)
         throw OFELIException("Side::getMeasure(): Area of quadrilateral " + to_string(_label) + " is null.");
      if (m<0)
         throw OFELIException("Side::getMeasure(): Area of quadrilateral " + to_string(_label) + " is negative.");
   }

   return m;
}


void Side::setDOF(size_t& first_dof,
                  size_t  nb_dof)
{
   _nb_dof = nb_dof;
   _first_dof = first_dof;
   for (size_t i=0; i<_nb_dof; i++)
      _dof[i] = first_dof++;
}


Point<real_t> Side::getNormal() const
{
   Point<real_t> N, T, c;
   the_element = getNeighborElement(1);
   if (The_element.getShape()==TRIANGLE) {
      T = _node[0]->getCoord() - _node[1]->getCoord();
      N.x = -T.y; N.y = T.x;
      c = Triang3(the_element).getCenter();
      if (N*(_node[0]->getCoord()-c) < 0)
         N = -N;
   }
   else if (The_element.getShape()==QUADRILATERAL) {
      T = _node[0]->getCoord() - _node[1]->getCoord();
      N.x = -T.y; N.y = T.x;
      c = Quad4(the_element).getCenter();
      if (N*(_node[0]->getCoord()-c) < 0)
         N = -N;
   }
   else if (The_element.getShape()==TETRAHEDRON) {
      Point<real_t> AB = _node[1]->getCoord() - _node[0]->getCoord();
      Point<real_t> AC = _node[2]->getCoord() - _node[0]->getCoord();
      N = CrossProduct(AB,AC);
      c = Tetra4(the_element).getCenter();
      if (N*(_node[0]->getCoord()-c) < 0)
         N = -N;
   }
   else if (The_element.getShape()==HEXAHEDRON) {
      Point<real_t> AB = _node[1]->getCoord() - _node[0]->getCoord();
      Point<real_t> AC = _node[2]->getCoord() - _node[0]->getCoord();
      N = CrossProduct(AB,AC);
      c = Hexa8(the_element).getCenter();
      if (N*(_node[0]->getCoord()-c) < 0)
         N = -N;
   }
   else
      throw OFELIException("Side::getNormal(): Computation of normal to side is not"
                           " implemented for this element.");
   return N;
}


Point<real_t> Side::getUnitNormal() const
{
   return getNormal()/getMeasure();
}


Element *Side::getOtherNeighborElement(Element* el) const
{
   if (_el[0]==el)
      return _el[1];
   if (_el[1]==el)
      return _el[0];
   return nullptr; 
}


void Side::shape_index(const string& shape)
{
   if (shape=="line")
      _shape = LINE;
   else if (shape=="tria" || shape=="triangle")
      _shape = TRIANGLE;
   else if (shape=="quad" || shape=="quadrilateral")
      _shape = QUADRILATERAL;
}


size_t Side::Contains(const Node* nd) const
{
   for (size_t i=0; i<_nb_nodes; i++) {
     if (_node[i]==nd)
        return i+1;
   }
   return 0;
}


int Side::isReferenced()
{
   for (size_t i=0; i<_nb_dof; i++) {
      if (_code[i] > 0)
         return 1;
      else if (_code[i] < 0)
         return -1;
   }
   return 0;
}


int Side::getGlobalCode() const
{
   int code=0, p=1;
   for (int i=_nb_dof; i>0; --i) {
      code += p*_code[i-1];
      p *= 10;
   }
   return code;
}


ostream& operator<<(ostream&    s,
                    const Side& sd)
{
   s << "\n Side: " << setw(5) << sd.n();
   s << " **  " << sd.getNbDOF() << " d.o.f.  ** Codes ";
   for (size_t i=1; i<=sd.getNbDOF(); i++)
      s << setw(5) << sd.getCode(i);
   s << "  ** Nodes: ";
   for (size_t j=1; j<=sd.getNbNodes(); j++)
      s << sd.getNodeLabel(j) << "  ";
   s << endl;
   return s;
}

} /* namespace OFELI */
