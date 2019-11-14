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

                         Implementation of class 'Edge'

 ==============================================================================*/

#include "mesh/Edge.h"
#include "util/util.h"
#include "OFELIException.h"

namespace OFELI {

Edge::Edge()
{
   _nb_eq = _nb_nodes = 0;
   _label = 0;
   _neig_sd = 0;
   _sd[0] = _sd[1] = nullptr;
   _on_boundary = -1;
   _code.resize(MAX_NBDOF_EDGE);
   _dof.resize(MAX_NBDOF_EDGE);
}


Edge::Edge(size_t label)
{
   if (label<1)
      throw OFELIException("Edge::Edge(size_t): Illegal edge label "+itos(label));
   _label = label;
   _nb_eq = _nb_nodes = 0;
   _neig_sd = 0;
   _sd[0] = _sd[1] = nullptr;
   _on_boundary = 0;
   _code.resize(MAX_NBDOF_EDGE);
   _dof.resize(MAX_NBDOF_EDGE);
}


Edge::Edge(const Edge& ed)
{
   _label = ed._label;
   _nb_eq = 0;
   _neig_sd = ed._neig_sd;
   _node[0] = ed._node[0];
   _node[1] = ed._node[1];
   _on_boundary = ed._on_boundary;
   _nb_dof = ed._nb_dof;
   _code.resize(_nb_dof);
   _dof.resize(_nb_dof);
   for (size_t j=0; j<_nb_dof; j++) {
      _code[j] = ed._code[j];
      _dof[j] = ed._dof[j];
   }
   _sd[0] = ed._sd[0];
   _sd[1] = ed._sd[1];
}


Edge::~Edge() { }


void Edge::Add(Node* node)
{
   if (node==nullptr)
      throw OFELIException("Edge::Add(Node *):  Trying to add an undefined node.");
   _node[_nb_nodes++] = node;
   _nb_eq += node->getNbDOF();
}


void Edge::setOnBoundary()
{
   _on_boundary = 1;
}


int Edge::isOnBoundary() const
{
   _on_boundary = -1;
   if ((_sd[0] && !_sd[1]) || (!_sd[0] && _sd[1]))
      _on_boundary = 1;
   if (_sd[0] && _sd[1])
      _on_boundary = 0;
   return _on_boundary;
}


ostream& operator<<(ostream&    s,
                    const Edge& ed)
{
   size_t i;
   s << "\n Edge: " << setw(5) << ed.n();
   s << " **  " << ed.getNbDOF() << " d.o.f.  ** Codes ";
   for (i=1; i<=ed.getNbDOF(); i++)
      s << setw(5) << ed.getCode(i);
   s << "  ** Nodes: ";
   s << ed.getNodeLabel(1) << "  ";
   s << ed.getNodeLabel(2) << "  ";
   s << endl;
   return s;
}

} /* namespace OFELI */
