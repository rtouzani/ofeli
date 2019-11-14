/*==============================================================================

                                   O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

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

                         Implementation of class 'Element'

  ==============================================================================*/

#include "mesh/Material.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "mesh/Node.h"
#include "linear_algebra/LocalMatrix_impl.h"
#include "io/fparser/fparser.h"
#include "util/util.h"
#include "OFELIException.h"

extern FunctionParser theParser;

namespace OFELI {

Element::Element()
        : _active(true), _nb_nodes(0), _nb_eq(0), _nb_sides(0), _label(0),
          _nbs(0), _nb_neig_el(0), _nb_childs(0), _nb_dof(1), _level(0),
          _code(1), _shape(0), _parent(nullptr)
{
   for (size_t i=0; i<MAX_NB_ELEMENT_SIDES; i++)
      _neig_el[i] = nullptr;
   _parent = nullptr;
}


Element::Element(size_t        label,
                 const string& shape)
        : _active(true), _nb_nodes(0), _nb_eq(0), _nb_sides(0), _label(label),
          _nbs(0), _nb_neig_el(0), _nb_childs(0), _nb_dof(1), _level(0),
          _code(1), _parent(nullptr)
{
   if (label<1)
      throw OFELIException("Element::Element(size_t,string): Illegal element label "+itos(label));
   shape_index(shape);
   calculate_nb_sides();
   for (size_t i=0; i<MAX_NB_ELEMENT_SIDES; i++)
      _neig_el[i] = nullptr;
}


Element::Element(size_t label,
                 int    shape)
        : _active(true), _nb_nodes(0), _nb_eq(0), _nb_sides(0), _label(label),
          _nbs(0), _nb_neig_el(0), _nb_childs(0), _nb_dof(1), _level(0),
          _code(1), _shape(shape), _parent(nullptr)
{
   if (label<1)
      throw OFELIException("Element::Element(size_t,int): Illegal element label "+itos(label));
   calculate_nb_sides();
   for (size_t i=0; i<MAX_NB_ELEMENT_SIDES; i++)
      _neig_el[i] = nullptr;
}


Element::Element(size_t        label,
                 const string& shape,
                 int           c)
        : _active(true), _nb_nodes(0), _nb_eq(0), _nb_sides(0), _label(label),
          _nbs(0), _nb_neig_el(0), _nb_childs(0), _nb_dof(1), _level(0),
          _code(c), _parent(nullptr)
{
   if (label<1)
      throw OFELIException("Element::Element(size_t,string,int): Illegal element label "+itos(label));
   shape_index(shape);
   calculate_nb_sides();
   for (size_t i=0; i<MAX_NB_ELEMENT_SIDES; i++)
      _neig_el[i] = nullptr;
}


Element::Element(size_t label,
                 int    shape,
                 int    c)
        : _active(true), _nb_nodes(0), _nb_eq(0), _nb_sides(0), _label(label),
          _nbs(0), _nb_neig_el(0), _nb_childs(0), _nb_dof(1), _level(0),
          _code(c), _shape(shape), _parent(nullptr)
{
   if (label<1)
      throw OFELIException("Element::Element(size_t,int,int): Illegal element label " + itos(label));
   calculate_nb_sides();
   for (size_t i=0; i<MAX_NB_ELEMENT_SIDES; i++)
      _neig_el[i] = nullptr;
}


Element::Element(const Element& el)
        : _active(el._active), _nb_nodes(el._nb_nodes), _nb_eq(el._nb_eq),
          _nb_sides(el._nb_sides), _label(el._label), _nbs(el._nbs), _nb_neig_el(el._nb_neig_el),
          _nb_childs(el._nb_childs), _nb_dof(el._nb_dof), _level(el._level),
          _code(el._code), _shape(el._shape), _parent(el._parent)
{
   for (size_t i=0; i<_nb_nodes; i++)
      _node[i] = el._node[i];
   for (size_t j=0; j<_nb_sides; j++)
      _side[j] = el._side[j];
}


size_t Element::getNbDOF() const
{
   size_t nb_dof = getPtrNode(1)->getNbDOF();
   for (size_t i=2; i<=_nb_nodes; i++) {
      if (getPtrNode(i)->getNbDOF() != nb_dof)
         nb_dof = 0;
   }
   return nb_dof;
}


void Element::calculate_nb_sides()
{
   if (_shape==LINE)
      _nbs = 0;
   else if (_shape==TRIANGLE)
     _nbs = 3;
   else if (_shape==QUADRILATERAL)
     _nbs = 4;
   else if (_shape==TETRAHEDRON)
     _nbs = 4;
   else if (_shape==HEXAHEDRON)
     _nbs = 6;
   else
     _nbs = 0;
}


int Element::setSide(size_t  n,
                     size_t* nd)
{
   int nb_nodes = 0;
   if (_shape==TRIANGLE) {
      nb_nodes = 2;
      _shape = LINE;
      nd[0] = n, nd[1] = n+1;
   }
   else if (_shape==QUADRILATERAL) {
      nb_nodes = 2;
      _shape = LINE;
      nd[0] = n, nd[1] = n+1;
   }
   else if (_shape==TETRAHEDRON) {
      nb_nodes = 3;
      _shape = TRIANGLE;
      if (n==1)
         nd[0] = 1, nd[1] = 2, nd[2] = 3;
      else if (n==2)
         nd[0] = 1, nd[1] = 3, nd[2] = 4;
      else if (n==3)
         nd[0] = 1, nd[1] = 4, nd[2] = 2;
      else if (n==4)
         nd[0] = 2, nd[1] = 4, nd[2] = 3;
   }
   else if (_shape==HEXAHEDRON) {
      nb_nodes = 4;
      _shape = QUADRILATERAL;
      if (n==1)
         nd[0] = 1, nd[1] = 2, nd[2] = 3, nd[3] = 4;
      else if (n==2)
         nd[0] = 5, nd[1] = 6, nd[2] = 7, nd[3] = 8;
      else if (n==3)
         nd[0] = 2, nd[1] = 6, nd[2] = 7, nd[3] = 3;
      else if (n==4)
         nd[0] = 1, nd[1] = 4, nd[2] = 8, nd[3] = 5;
      else if (n==5)
         nd[0] = 1, nd[1] = 2, nd[2] = 6, nd[3] = 5;
      else if (n==6)
         nd[0] = 3, nd[1] = 7, nd[2] = 8, nd[3] = 4;
   }
   return nb_nodes;
}


void Element::Add(Node* node)
{
   if (!node)
      throw OFELIException("Element::Add(Node *): Trying to add an undefined node");
   _node[_nb_nodes++] = node;
   _nb_eq += node->getNbDOF();
}


void Element::Add(Node* node,
                  int   n)
{
   if (n>0) {
      _node[n-1] = node;
      _nb_eq += node->getNbDOF();
      _nb_nodes++;
   }
}


void Element::Replace(size_t label,
                      Side*  side)
{
   if (!side)
      throw OFELIException("Element::Replace(size_t,Side *): Trying to replace "+itos(label) + 
                           "-th side by an undefined side.");
   _side[label-1] = side;
}


void Element::Replace(size_t label,
                      Node*  node)
{
   if (!node)
      throw OFELIException("Element::Replace(size_t,Node *): Trying to replace "+itos(label) +
                           "-th node by an undefined node.");
   _node[label-1] = node;
}


void Element::Add(Side* sd)
{
   Node *nd1=(*sd)(1), *nd2=(*sd)(2);
   if (!sd)
      throw OFELIException("Element::Add(Side *): Trying to add an undefined side.");

   switch (_shape) {

      case LINE:
         _nb_sides = 0;
         throw OFELIException("Element::Add(Side *): Shape of side is incompatible with element.");
         break;

      case TRIANGLE:
         if (sd->getShape()!=LINE)
            throw OFELIException("Element::Add(Side *): Shape of side is incompatible with element.");
         if ( (nd1==_node[0] && nd2==_node[1]) || (nd1==_node[1] && nd2==_node[0]) )
            _side[0] = sd;
         if ( (nd1==_node[1] && nd2==_node[2]) || (nd1==_node[2] && nd2==_node[1]) )
            _side[1] = sd;
         if ( (nd1==_node[2] && nd2==_node[0]) || (nd1==_node[0] && nd2==_node[2]) )
            _side[2] = sd;
         _nb_sides = 3;
         break;

      case QUADRILATERAL:
         if (sd->getShape()!=LINE)
            throw OFELIException("Element::Add(Side *): Shape of side is incompatible with element.");
         if ( (nd1==_node[0] && nd2==_node[1]) || (nd1==_node[1] && nd2==_node[0]) )
            _side[0] = sd;
         if ( (nd1==_node[1] && nd2==_node[2]) || (nd1==_node[2] && nd2==_node[1]) )
            _side[1] = sd;
         if ( (nd1==_node[2] && nd2==_node[3]) || (nd1==_node[3] && nd2==_node[2]) )
            _side[2] = sd;
         if ( (nd1==_node[3] && nd2==_node[0]) || (nd1==_node[0] && nd2==_node[3]) )
            _side[3] = sd;
         _nb_sides = 4;
         break;

      case TETRAHEDRON:
         if (sd->getShape()!=TRIANGLE)
            throw OFELIException("Element::Add(Side *): Shape of side is incompatible with element.");
         if ( (nd1==_node[2] && nd2==_node[0]) || (nd1==_node[0] && nd2==_node[2]) )
            _side[0] = sd;
         _nb_sides = 4;
         break;

      case HEXAHEDRON:
         if (sd->getShape()!=QUADRILATERAL)
            throw OFELIException("Element::Add(Side *): Shape of side is incompatible with element.");
         _nb_sides = 6;
         break;
   }
}


void Element::Add(Side* sd,
                  int   k)
{
   if (k>0) {
      if (!sd)
         throw OFELIException("Element::Add(Side *,int): Trying to add an undefined side.");
      _side[k-1] = sd;
   }
}


void Element::set(Element* el,
                  int      n)
{
   if (n>0) {
      _neig_el[n-1] = el;
      _nb_neig_el = n;
   }
}
   
   
int Element::Contains(const Node& nd) const
{
   for (size_t i=0; i<_nb_nodes; i++)
      if (_node[i]->n()==nd.n())
         return int(i+1);
   return 0;
}


int Element::Contains(const Node* nd) const
{
   for (size_t i=0; i<_nb_nodes; i++)
      if (_node[i]==nd)
         return int(i+1);
   return 0;
}
   
   
int Element::Contains(const Side& sd) const
{
   for (size_t i=0; i<_nb_sides; i++)
      if (_side[i]->n()==sd.n())
         return int(i+1);
   return 0;
}


int Element::Contains(const Side* sd) const
{
   for (size_t i=0; i<_nb_sides; i++)
      if (_side[i]==sd)
         return int(i+1);
   return 0;
}


void Element::setCode(const string& exp,
                            int     code)
{
   theParser.Parse(exp.c_str(),"x,y,z");
   real_t d[3];
   for (size_t i=0; i<_nb_nodes; i++) {
      d[0] = _node[i]->getCoord(1);
      d[1] = _node[i]->getCoord(2);
      d[2] = _node[i]->getCoord(3);
      if (!theParser.Eval(d))
         return;
   }
   _code = code;
}


size_t Element::getNbVertices() const
{
   size_t n=0;
   if (_shape==LINE)
      n = 2;
   else if (_shape==TRIANGLE)
      n = 3;
   else if (_shape==QUADRILATERAL || _shape==TETRAHEDRON)
      n = 4;
   else if (_shape==HEXAHEDRON)
      n = 8;
   return n;
}


void Element::setGlobalToLocal()
{
   for (size_t i=1; i<=_nb_nodes; i++)
      _g2l[_node[i-1]->n()] = i;
}


bool Element::isOnBoundary() const
{
   for (size_t i=0; i<_nb_sides; i++)
      if (_side[i]->isOnBoundary())
         return true;
   return false;
}


real_t Element::getMeasure() const
{
   Point<real_t> x[8], dshl[8];
   real_t m = 0;

   switch (_shape) {

      case LINE:
         {
            x[0] = _node[0]->getCoord();
            x[1] = _node[1]->getCoord();
            m = (x[1].x-x[0].x)*(x[1].x-x[0].x) - (x[1].y-x[0].y)*(x[1].y-x[0].y);
            if (m==0)
               throw OFELIException("Element::getMeasure(): Length of line "+itos(_label)+" is null.");
            break;
         }

      case TRIANGLE:
         {
            x[0] = _node[0]->getCoord();
            x[1] = _node[1]->getCoord();
            x[2] = _node[2]->getCoord();
            m = 0.5*((x[1].x-x[0].x)*(x[2].y-x[0].y) - (x[1].y-x[0].y)*(x[2].x-x[0].x));
            if (m==0)
               throw OFELIException("Element::getMeasure(): Area of triangle "+itos(_label)+" is null.");
            if (m<0)
               throw OFELIException("Element::getMeasure(): Area of triangle "+itos(_label) +
                                    " is negative. Triangle is incorrectly oriented.");
            break;
         }

      case TETRAHEDRON:
         {
            x[0] = _node[0]->getCoord();
            x[1] = _node[1]->getCoord();
            x[2] = _node[2]->getCoord();
            x[3] = _node[3]->getCoord();
            size_t i, j, k;
            LocalMatrix<real_t,3,3> J, IJ;
            dshl[0] = -1.0;
            dshl[1].x = dshl[2].y = dshl[3].z =  1.0;
            dshl[1].y = dshl[1].z = dshl[2].x =  0.0;
            dshl[2].z = dshl[3].x = dshl[3].y =  0.0;
            for (k=0; k<4; k++) {
               J(1,1) += dshl[k].x*x[k].x;
               J(1,2) += dshl[k].y*x[k].x;
               J(1,3) += dshl[k].z*x[k].x;
               J(2,1) += dshl[k].x*x[k].y;
               J(2,2) += dshl[k].y*x[k].y;
               J(2,3) += dshl[k].z*x[k].y;
               J(3,1) += dshl[k].x*x[k].z;
               J(3,2) += dshl[k].y*x[k].z;
               J(3,3) += dshl[k].z*x[k].z;
            }
            for (i=0; i<3; i++) {
               j = (i+1)%3; k = (j+1)%3;
               IJ(i+1,i+1) = J(j+1,j+1) * J(k+1,k+1) - J(j+1,k+1) * J(k+1,j+1);
               IJ(j+1,i+1) = J(j+1,k+1) * J(k+1,i+1) - J(j+1,i+1) * J(k+1,k+1);
               IJ(i+1,j+1) = J(k+1,j+1) * J(i+1,k+1) - J(i+1,j+1) * J(k+1,k+1);
            }
            m = OFELI_SIXTH*(J(1,1)*IJ(1,1) + J(2,1)*IJ(1,2) + J(3,1)*IJ(1,3));
            if (m==0)
               throw OFELIException("Element::getMeasure(): Area of tetrahedron "+itos(_label)+" is null.");
             if (m<0)
                throw OFELIException("Element::getMeasure(): Area of tetrahedron "+itos(_label) +
                                     " is negative. Tetrahedron is incorrectly oriented.");
            break;
         }

      case QUADRILATERAL:
         {
            x[0] = _node[0]->getCoord();
            x[1] = _node[1]->getCoord();
            x[2] = _node[2]->getCoord();
            x[3] = _node[3]->getCoord();
            dshl[0].x = dshl[3].x = dshl[0].y = dshl[1].y = -0.25;
            dshl[1].x = dshl[2].x = dshl[3].y = dshl[2].y =  0.25;
            real_t dxds=0, dxdt=0, dyds=0, dydt=0;
            for (size_t i=0; i<4; i++) {
               dxds += dshl[i].x*x[i].x;
               dxdt += dshl[i].y*x[i].x;
               dyds += dshl[i].x*x[i].y;
               dydt += dshl[i].y*x[i].y;
            }
            m = dxds*dydt - dxdt*dyds;
            if (m==0)
               throw OFELIException("Element::getMeasure(): Area of quadrilateral "+itos(_label)+" is null.");
            if (m<0)
                throw OFELIException("Element::getMeasure(): Area of quadrilateral " + itos(_label) +
                                     " is negative. Quadrilateral is incorrectly oriented.");
            break;
         }

      case HEXAHEDRON:
         {
            for (size_t j=0; j<8; j++)
               x[j] = _node[j]->getCoord();
            dshl[0] = Point<real_t>(-0.125,-0.125,-0.125);
            dshl[1] = Point<real_t>( 0.125,-0.125,-0.125);
            dshl[2] = Point<real_t>( 0.125, 0.125,-0.125);
            dshl[3] = Point<real_t>(-0.125, 0.125,-0.125);
            dshl[4] = Point<real_t>(-0.125,-0.125, 0.125);
            dshl[5] = Point<real_t>( 0.125,-0.125, 0.125);
            dshl[6] = Point<real_t>( 0.125, 0.125, 0.125);
            dshl[7] = Point<real_t>(-0.125, 0.125, 0.125);
            Point<real_t> a = 0;
            for (size_t j=0; j<8; j++)
               a += dshl[j].x * x[j];
            LocalMatrix<real_t,3,3> J, IJ;
            J(1,1) = a.x; J(2,1) = a.y; J(3,1) = a.z;
            a = 0;
            for (size_t j=0; j<8; j++)
               a += dshl[j].y * x[j];
            J(1,2) = a.z; J(2,2) = a.y; J(3,2) = a.y;
            J(1,3) = a.x; J(2,3) = a.y; J(3,3) = a.z;
            for (size_t i=0; i<3; i++) {
               size_t j = (i+1)%3; size_t k = (j+1)%3;
               IJ(i+1,i+1) = J(j+1,j+1) * J(k+1,k+1) - J(j+1,k+1) * J(k+1,j+1);
               IJ(j+1,i+1) = J(j+1,k+1) * J(k+1,i+1) - J(j+1,i+1) * J(k+1,k+1);
               IJ(i+1,j+1) = J(k+1,j+1) * J(i+1,k+1) - J(i+1,j+1) * J(k+1,k+1);
            }
            m = J(1,1)*IJ(1,1) + J(2,1)*IJ(1,2) + J(3,1)*IJ(1,3);
            if (m==0)
               throw OFELIException("Element::getMeasure(): Area of hexahedron "+itos(_label)+" is null.");
            if (m<0)
               throw OFELIException("Element::getMeasure(): Area of hexahedron "+itos(_label) +
                                    " is negative. Hexahedron is incorrectly oriented.");
            break;
         }
   }
   return fabs(m);
}


Point<real_t> Element::getUnitNormal(size_t i) const
{
   size_t j;
   Point<real_t> N, x[4];
   real_t s=0.;
   i--;
   if (_shape==TRIANGLE) {
      j = (i+1)%3;
      x[0] = _node[i]->getCoord();
      x[1] = _node[j]->getCoord();
      N = Point<real_t>(x[1].y-x[0].y,x[0].x-x[1].x);
      s = sqrt(N*N);
   }
   else if (_shape==QUADRILATERAL) {
      j = (i+1)%4;
      x[0] = _node[i]->getCoord();
      x[1] = _node[j]->getCoord();
      N = Point<real_t>(x[1].y-x[0].y,x[0].x-x[1].x);
      s = sqrt(N*N);
   }
   else
      throw OFELIException("Element::getUnitNormal(size_t): Calculation of outward normal "
                           "is not available for element shape " + itos(_shape));
   return N*(1./s);
}


void Element::shape_index(const string& shape)
{
   if (shape=="line")
      _shape = LINE;
   else if (shape=="tria" || shape=="triangle")
      _shape = TRIANGLE;
   else if (shape=="quad" || shape=="quadrilateral")
      _shape = QUADRILATERAL;
   else if (shape=="tetra" || shape=="tetrahedron")
      _shape = TETRAHEDRON;
   else if (shape=="hexa" || shape=="hexahedron")
      _shape = HEXAHEDRON;
   else if (shape=="penta" || shape=="pentahedron")
      _shape = PENTAHEDRON;
}


void Element::setChild(Element* el)
{
   _child[_nb_childs++] = el;
   el->_level = _level + 1;
   el->_parent = this;
   el->_code = _code;
   _active = false;
}


Element *Element::getChild(size_t i) const
{
   if (i>_nb_childs)
      throw OFELIException("Element::getChild(size_t): Number of children is too large: "+itos(i));
   return _child[i-1];
}


size_t Element::IsIn(const Node* nd)
{
   size_t n = 0;
   for (size_t i=0; i<_nb_nodes; i++) {
      if (_node[i]==nd) {
         n = i+1;
         break;
      }
   }
   return n;
}


bool Element::isActive() const
{
   if (_active)
      return true;
   else
      return false;
}

size_t Element::getNodeLabel(size_t n) const
{
   return _node[n-1]->n();
}


size_t Element::getSideLabel(size_t n) const
{
   return _side[n-1]->n();
}


ostream& operator<<(ostream&       s,
                    const Element& el)
{
   s << "\n Element: " << setw(6) << el.n();
   s << " ** Code: " << setw(8) << el.getCode();
   s << " ** Nodes: ";
   for (size_t i=1; i<=el.getNbNodes(); i++)
      s << " " << el.getNodeLabel(i);
   s << endl;
   return s;
}

} /* namespace OFELI */
