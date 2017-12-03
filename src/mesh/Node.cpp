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

                         Implementation of class 'Node'

  ==============================================================================*/

#include "util/macros.h"
#include "mesh/Node.h"
#include "io/output.h"
#include "linear_algebra/Point.h"
#include "io/fparser/fparser.h"
#include "OFELIException.h"

extern FunctionParser theParser;

namespace OFELI {

Node::Node()
{
   _nb_dof = 1;
   _label = _nb_neig_el = _neig_i = 0;
   _code[0] = 0;
   _x = 0;
   _on_boundary = false;
   _level = 0;
   for (size_t i=0; i<MAX_NBDOF_NODE; i++)
      _dof[i] = 0;
}


Node::Node(      size_t         label,
           const Point<double>& x)
{
   _nb_neig_el = _neig_i = 0;
   _x = x;
   _code[0] = 0;
   _nb_dof = 0;
   _label = label;
   _on_boundary = false;
   _level = 0;
   for (size_t i=0; i<MAX_NBDOF_NODE; i++)
      _dof[i] = 0;
}


Node::Node(const Node& node)
{
   _label = node._label;
   _nb_dof = node._nb_dof;
   _el.resize(_el.size());
   _el = node._el;
   _nb_neig_el = node._nb_neig_el;
   _neig_i = node._neig_i;
   _x = node._x;
   for (size_t i=0; i<MAX_NBDOF_NODE; i++) {
      _code[i] = node._code[i];
      _dof[i] = node._dof[i];
   }
   _first_dof = node._first_dof;
   _on_boundary = node._on_boundary;
   _level = node._level;
}


Node::~Node() { }


void Node::setNbDOF(size_t n)
{
   if (n<=0)
      throw OFELIException("Node::setNbDOF(size_t): Illegal argument.");
   _nb_dof = n;
   for (size_t i=0; i<_nb_dof; i++) {
      _dof[i] = _first_dof + i;
      _code[i] = 0;
   }
}


void Node::Add(Element* el)
{
   if (!el)
      throw OFELIException("Node::Add(Element *): Trying to add an illegal neighbour element.");
   _el.push_back(el);
   _neig_i++; _nb_neig_el++;
}


void Node::setCode(size_t dof,
                   int    code)
{
   _code[dof-1] = code;
}


void Node::setCode(const vector<int>& code)
{
   for (size_t i=0; i<_nb_dof; i++)
      _code[i] = code[i];
}


void Node::setCode(int* code)
{
   for (size_t i=0; i<_nb_dof; i++)
      _code[i] = code[i];
}


void Node::setCode(const string& exp,
                         int     code,
                         size_t  dof)
{
   PARSE(exp.c_str(),"x,y,z,t");
   double d[] = {_x.x,_x.y,_x.z};
   if (EVAL(d))
      _code[dof-1] = code;
}


void Node::setDOF(size_t& first_dof,
                  size_t  nb_dof)
{
   _nb_dof = nb_dof;
   _first_dof = first_dof;
   for (size_t i=0; i<_nb_dof; i++)
      _dof[i] = first_dof++;
}


ostream& operator<<(      ostream& s,
                    const Node&    nd)
{
   size_t i;
   s << "\n Node: " << setw(8) << nd.n() << "\n Coordinates: ";
   s.setf(ios::right|ios::scientific);
   for (i=1; i<=3; i++)
      s << std::setprecision(6) << setw(12) << nd.getCoord(i) << "  ";
   s << "\n " << setw(2) << nd.getNbDOF() << " d.o.f.: ";
   for (i=1; i<=nd.getNbDOF(); i++)
      s << setw(8) << nd.getDOF(i) << "  ";
   s << "\n Codes:     ";
   for (i=1; i<=nd.getNbDOF(); i++)
      s << setw(8) << nd.getCode(i);
   s << endl;
   if (nd.getNbNeigEl() > 0)
      s << "\n Nb. of neighbor elements: " << setw(6) << nd.getNbNeigEl() << endl;
   return s;
}

} /* namespace OFELI */
