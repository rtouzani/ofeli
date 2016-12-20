/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2017 Rachid Touzani

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

            Implementation of class 'NodeList', 'ElementList', 'SideList'

  ==============================================================================*/

#include "mesh/MeshExtract.h"
#include "util/util.h"

namespace OFELI {

NodeList::NodeList(Mesh& ms)
{
   _ms = &ms;
   _nbn = _nb_nodes = 0;
}


void NodeList::selectCode(int code,
                          int dof)
{
   _nbn = _nb_nodes;
   mesh_nodes(*_ms) {
      if (the_node->getCode(dof)==code)
         _nb_nodes++;
   }
   _theList.resize(_nb_nodes);
   size_t i=_nbn;
   mesh_nodes(*_ms) {
      if (the_node->getCode(dof)==code)
         _theList[i++] = the_node;
   }
}


void NodeList::unselectCode(int code,
                            int dof)
{
   _nbn = _nb_nodes;
   mesh_nodes(*_ms) {
      if (the_node->getCode(dof)!=code)
         _nb_nodes++;
   }
   _theList.resize(_nb_nodes);
   size_t i=_nbn;
   mesh_nodes(*_ms) {
      if (the_node->getCode(dof)!=code)
         _theList[i++] = the_node;
   }
}


void NodeList::selectCoordinate(real_t x,
                                real_t y,
                                real_t z)
{
   _nbn = _nb_nodes;
   mesh_nodes(*_ms) {
//    Any
      if (x==ANY && y==ANY && z==ANY)
         _nb_nodes++;
//    x-coordinate only
      if (Equal(the_node->getX(),x))
         _nb_nodes++;
//    y-coordinate only
      else if (Equal(the_node->getY(),y))
         _nb_nodes++;
//    z-coordinate only
      else if (Equal(the_node->getZ(),z))
         _nb_nodes++;
//    xy-coordinate
      else if (Equal(the_node->getX(),x) && Equal(the_node->getY(),y))
         _nb_nodes++;
//    xz-coordinate
      else if (Equal(the_node->getX(),x) && Equal(the_node->getZ(),z)) 
         _nb_nodes++;
//    yz-coordinate
      else if (Equal(the_node->getY(),y) && Equal(the_node->getZ(),z))
         _nb_nodes++;
//    xyz-coordinate
      else if (Equal(the_node->getX(),x) && 
               Equal(the_node->getY(),y) && 
               Equal(the_node->getZ(),z))
         _nb_nodes++;
   }
   _theList.resize(_nb_nodes);
   size_t i=_nbn;
   mesh_nodes(*_ms) {
      if (x==ANY && y==ANY && z==ANY)
         _theList[i++] = the_node;
      if (Equal(the_node->getX(),x))
         _theList[i++] = the_node;
      else if (Equal(the_node->getY(),y))
         _theList[i++] = the_node;
      else if (Equal(the_node->getZ(),z))
         _theList[i++] = the_node;
      else if (Equal(the_node->getX(),x) && Equal(the_node->getY(),y))
         _theList[i++] = the_node;
      else if (Equal(the_node->getX(),x) && Equal(the_node->getZ(),z))
         _theList[i++] = the_node;
      else if (Equal(the_node->getY(),y) && Equal(the_node->getZ(),z))
         _theList[i++] = the_node;
      else if (Equal(the_node->getX(),x) && Equal(the_node->getY(),y) && Equal(the_node->getZ(),z))
         _theList[i++] = the_node;
   }
}


ElementList::ElementList(Mesh& ms)
{
   _ms = &ms;
   _nbn = _nb_elements = 0;
}


void ElementList::selectCode(int code)
{
   _nbn = _nb_elements;
   mesh_elements(*_ms) {
      if (The_element.getCode()==code)
         _nb_elements++;
   }
   _theList.resize(_nb_elements);
   size_t i=_nbn;
   mesh_elements(*_ms)
      if (The_element.getCode()==code)
         _theList[i++] = the_element;
}


void ElementList::selectLevel(int level)
{
   _nbn = _nb_elements;
   mesh_elements(*_ms) {
      if (The_element.getLevel()==level)
         _nb_elements++;
   }
   _theList.resize(_nb_elements);
   size_t i=_nbn;
   mesh_elements(*_ms)
      if (The_element.getLevel()==level)
         _theList[i++] = the_element;
}


void ElementList::selectActive()
{
   _nbn = _nb_elements;
   mesh_elements(*_ms) {
      if (The_element.isActive())
         _nb_elements++;
   }
   _theList.resize(_nb_elements);
   size_t i=_nbn;
   mesh_elements(*_ms)
      if (The_element.isActive())
         _theList[i++] = the_element;
}


void ElementList::selectMarked(const Vect<real_t>& e,
                                     real_t        threshold)
{
   _nbn = _nb_elements;
   mesh_elements(*_ms) {
      if (The_element.isActive() && e(element_label)>threshold)
         _nb_elements++;
   }
   _theList.resize(_nb_elements);
   size_t i=_nbn;
   mesh_elements(*_ms)
      if (The_element.isActive() && e(element_label)>threshold)
         _theList[i++] = the_element;
}


void ElementList::unselectCode(int code)
{
   _nbn = _nb_elements;
   mesh_elements(*_ms) {
      if (The_element.getCode() != code)
         _nb_elements++;
   }
   _theList.resize(_nb_elements);
   size_t i=_nbn;
   mesh_elements(*_ms)
      if (The_element.getCode() != code)
         _theList[i++] = the_element;
}


SideList::SideList(Mesh& ms)
{
   _ms = &ms;
   _nbn = _nb_sides = 0;
}


void SideList::selectCode(int code,
                          int dof)
{
   _nbn = _nb_sides;
   mesh_sides(*_ms)
      if (The_side.getCode(dof)==code)
         _nb_sides++;
   _theList.resize(_nb_sides);
   size_t i=_nbn;
   mesh_sides(*_ms)
      if (The_side.getCode(dof)==code)
         _theList[i++] = the_side;
}


void SideList::unselectCode(int code,
                            int dof)
{
   _nbn = _nb_sides;
   mesh_sides(*_ms)
      if (The_side.getCode(dof) != code)
         _nb_sides++;
   _theList.resize(_nb_sides);
   size_t i=_nbn;
   mesh_sides(*_ms)
      if (The_side.getCode(dof) != code)
         _theList[i++] = the_side;
}


EdgeList::EdgeList(Mesh& ms)
{
   _ms = &ms;
   _nbn = _nb_edges = 0;
}


void EdgeList::selectCode(int code,
                          int dof)
{
   _nbn = _nb_edges;
   mesh_edges(*_ms)
      if (The_edge.getCode(dof)==code)
         _nb_edges++;
   _theList.resize(_nb_edges);
   size_t i=_nbn;
   mesh_edges(*_ms)
      if (The_edge.getCode(dof)==code)
         _theList[i++] = the_edge;
}


void EdgeList::unselectCode(int code, int dof)
{
   _nbn = _nb_edges;
   mesh_edges(*_ms)
      if (The_edge.getCode(dof)!=code)
         _nb_edges++;
   _theList.resize(_nb_edges);
   size_t i=_nbn;
   mesh_edges(*_ms)
      if (The_edge.getCode(dof)!=code)
         _theList[i++] = the_edge;
}


ostream& operator<<(      ostream&  s,
                    const NodeList& nl)
{
   s << "\n\nLIST OF SELECTED NODES:\n";
   Node *nd;
   for (nl.top(); (nd=nl.get());)
      s << *nd;
   s << endl;
   return s;
}


ostream& operator<<(      ostream&     s,
                    const ElementList& el)
{
   s << "\n\nLIST OF SELECTED ELEMENTS:\n";
   Element *e;
   for (el.top(); (e=el.get());)
      s << *e;
   s << endl;
   return s;
}


ostream& operator<<(      ostream&  s,
                    const SideList& sl)
{
   s << "\n\nLIST OF SELECTED SIDES:\n";
   Side *sd;
   for (sl.top(); (sd=sl.get());)
      s << *sd;
   s << endl;
   return s;
}


ostream& operator<<(      ostream&  s,
                    const EdgeList& el)
{
   s << "\n\nLIST OF SELECTED EDGES:\n";
   Edge *ed;
   for (el.top(); (ed=el.get());)
      s << *ed;
   s << endl;
   return s;
}

} /* namespace OFELI */
