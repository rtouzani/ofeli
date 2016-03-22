/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzani

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

                         Implementation of class 'Mesh'

  ==============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <string>

#include <iostream>
using std::cout;
using std::ostream;
using std::endl;
using std::clog;

#include <fstream>
using std::ofstream;

#include <iomanip>
using std::setw;

#include <vector>
using std::vector;

#include <algorithm>
using std::string;

#include "mesh/Mesh.h"
#include "mesh/MeshUtil.h"
#include "linear_algebra/GraphOfMatrix.h"
#include "io/XMLParser.h"
#include "util/qksort.h"
#include "util/util.h"
#include "mesh/Grid.h"
#include "mesh/Material.h"
#include "mesh/saveMesh.h"
#include "mesh/getMesh.h"
#include <algorithm>

namespace OFELI {

Node *theNode;
Element *theElement;
Side *theSide;
Edge *theEdge;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
Node    *the_node;
Element *the_element;
Side    *the_side;
Edge    *the_edge;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

extern Material theMaterial;

Mesh::Mesh()
     : _nb_nodes(0), _nb_elements(0), _nb_sides(0), _nb_edges(0), _dim(2), _nb_dof(0),
       _nb_vertices(0), _first_dof(1), _nb_mat(0), _verb(0), _no_imposed_dof(false),
       _is_structured(false), _all_sides_created(false), _boundary_sides_created(false),
       _all_edges_created(false), _boundary_edges_created(false), _boundary_nodes_created(false),
       _node_neighbor_elements_created(false), _element_neighbor_elements_created(false),
       _n_view1(0), _n_view2(0), _e_view1(0), _e_view2(0), _s_view1(0), _s_view2(0),
       _ed_view1(0), _ed_view2(0)
{
   setNodesForDOF();
}


Mesh::Mesh(const string& file,
           bool          bc,
           int           opt)
     : _nb_nodes(0), _nb_elements(0), _nb_sides(0), _nb_edges(0), _dim(2), _nb_dof(0),
       _nb_vertices(0), _first_dof(1), _nb_mat(0), _verb(0), _no_imposed_dof(false),
       _is_structured(false), _all_sides_created(false), _boundary_sides_created(false),
       _all_edges_created(false), _boundary_edges_created(false), _boundary_nodes_created(false),
       _node_neighbor_elements_created(false), _element_neighbor_elements_created(false),
       _n_view1(0), _n_view2(0), _e_view1(0), _e_view2(0), _s_view1(0), _s_view2(0),
       _ed_view1(0), _ed_view2(0)
{
   setNodesForDOF();
   get(file);
   setDOFSupport(opt);
   if (bc)
      removeImposedDOF();
   string mat = theMaterial.getName(1);
}


Mesh::Mesh(real_t L,
           size_t nb_el,
           size_t p,
           size_t nb_dof)
     : _nb_nodes(0), _nb_elements(0), _nb_sides(0), _nb_edges(0), _dim(1), _nb_dof(0),
       _nb_vertices(0), _first_dof(1), _nb_mat(0), _nb_boundary_nodes(2), _verb(0), _no_imposed_dof(false),
       _is_structured(false), _all_sides_created(false), _boundary_sides_created(false),
       _all_edges_created(false), _boundary_edges_created(false), _boundary_nodes_created(false),
       _node_neighbor_elements_created(false), _element_neighbor_elements_created(false),
       _n_view1(0), _n_view2(0), _e_view1(0), _e_view2(0), _s_view1(0), _s_view2(0),
       _ed_view1(0), _ed_view2(0)
{
   size_t NbN = nb_el*p + 1;

// Insert nodes
   real_t xx=0., h=L/real_t(nb_el);
   size_t nn=1;
   for (size_t nnd=1; nnd<=NbN; nnd++) {
      Point<real_t> x(xx);
      Node *nd = new Node(nn++,x);
      nd->setNbDOF(nb_dof);
      for (size_t i=0; i<nb_dof; i++)
         _code[i] = 0;
      nd->setDOF(_first_dof,nb_dof);
      nd->setCode(_code);
      Add(nd);
      xx += h/p;
   }
   _nb_eq = _nb_dof;

// Insert elements
   nn = 0;
   for (size_t nne=1; nne<=nb_el; nne++) {
      Element *el = new Element(nne,LINE,1);
      for (size_t i=1; i<=p+1; i++)
         el->Add(_nodes[nn++]);
      nn--;
      el->setCode(1);
      el->getMeasure();
      Add(el);
   }
   setNodesForDOF();
   _boundary_nodes.push_back(_nodes[0]);
   _boundary_nodes.push_back(_nodes[_nb_nodes-1]);
}


Mesh::Mesh(const Grid& g,
           int         opt)
     : _nb_nodes(0), _nb_elements(0), _nb_sides(0), _nb_edges(0), _dim(3), _nb_dof(0),
       _nb_vertices(0), _first_dof(1), _nb_mat(1), _verb(0), _no_imposed_dof(false),
       _is_structured(true), _all_sides_created(false), _boundary_sides_created(false),
       _all_edges_created(false), _boundary_edges_created(false), _boundary_nodes_created(false),
       _node_neighbor_elements_created(false), _element_neighbor_elements_created(false),
       _n_view1(0), _n_view2(0), _e_view1(0), _e_view2(0), _s_view1(0), _s_view2(0),
       _ed_view1(0), _ed_view2(0)
{
   Node *nd;
   Element *el;
   size_t i, j, k;
   size_t nx=g.getNx(), ny=g.getNy(), nz=g.getNz();
   Vect<size_t> nn;
   setNodesForDOF();
   size_t n=0, ne=0, m=1;
   if (nz==0)
      _dim = 2;
   if (ny==0)
      _dim = 1;
   if (_dim==3)
      opt = HEXAHEDRON;

// Nodes
   switch (_dim) {

      case 1:
         nn.setSize(nx+1);
         nn = 0;
         for (i=1; i<=nx+1; i++) {
            nn(i) = ++n;
            nd = new Node(nn(i),g.getCoord(i));
            _code[0] = g.getCode(i);
            nd->setNbDOF(1);
            nd->setDOF(_first_dof,1);
            nd->setCode(_code);
            Add(nd);
         }
         break;

      case 2:
         nn.setSize(nx+1,ny+1);
         nn = 0;
         for (i=1; i<=nx; i++) {
            for (j=1; j<=ny; j++) {
               if (g.isActive(i,j))
                  nn(i,j) = nn(i+1,j) = nn(i+1,j+1) = nn(i,j+1) = 1;
            }
         }
         m = 0;
         for (i=1; i<=nx+1; i++) {
            for (j=1; j<=ny+1; j++) {
               if (nn(i,j))
                  nn(i,j) = ++m;
            }
         }
         for (i=1; i<=nx+1; i++) {
            for (j=1; j<=ny+1; j++) {
               if (nn(i,j)) {
                  nd = new Node(++n,g.getCoord(i,j));
                  _code[0] = g.getCode(i,j);
                  nd->setNbDOF(1);
                  nd->setDOF(_first_dof,1);
                  nd->setCode(_code);
                  Add(nd);
               }
            }
         }
         break;

      case 3:
         nn.setSize(nx+1,ny+1,nz+1);
         nn = 0;
         for (i=1; i<=nx; i++) {
            for (j=1; j<=ny; j++) {
               for (k=1; k<=nz; k++) {
                  if (g.isActive(i,j,k)) {
                     nn(i  ,j,k) = nn(i  ,j,k+1) = nn(i  ,j+1,k) = nn(i  ,j+1,k+1) = 1;
                     nn(i+1,j,k) = nn(i+1,j,k+1) = nn(i+1,j+1,k) = nn(i+1,j+1,k+1) = 1;
                  }
               }
            }
         }
         m = 0;
         for (i=1; i<=nx+1; i++) {
            for (j=1; j<=ny+1; j++) {
               for (k=1; k<=nz+1; k++) {
                  if (nn(i,j,k))
                     nn(i,j,k) = ++m;
               }
            }
         }
         for (i=1; i<=nx+1; i++) {
            for (j=1; j<=ny+1; j++) {
               for (k=1; k<=nz+1; k++) {
                  if (nn(i,j,k)) {
                     nd = new Node(++n,g.getCoord(i,j,k));
                     _code[0] = g.getCode(i,j,k);
                     nd->setNbDOF(1);
                     nd->setDOF(_first_dof,1);
                     nd->setCode(_code);
                     Add(nd);
                  }
               }
            }
         }
         break;
   }

// Elements
   n = ne = 0;
   m = 1;
   switch (_dim) {

      case 1:
         for (i=1; i<=nx; i++) {
            el = new Element(++ne,LINE,1);
            el->Add(_nodes[i-1]);
            el->Add(_nodes[i]);
            Add(el);
         }
         break;

      case 2:
         if (opt==QUADRILATERAL) {
            for (j=1; j<=ny; j++) {
               for (i=1; i<=nx; i++) {
                  if (g.isActive(i,j)) {
                     el = new Element(++ne,QUADRILATERAL,1);
                     el->Add(_nodes[nn(i  ,j  )-1]);
                     el->Add(_nodes[nn(i+1,j  )-1]);
                     el->Add(_nodes[nn(i+1,j+1)-1]);
                     el->Add(_nodes[nn(i  ,j+1)-1]);
                     Add(el);
                  }
               }
            }
         }
         else if (opt==TRIANGLE) {
           for (i=1; i<=nx; i++) {
               for (j=1; j<=ny; j++) {
                  if (g.isActive(i,j)) {
                     el = new Element(++ne,TRIANGLE,1);
                     el->Add(_nodes[nn(i  ,j  )-1]);
                     el->Add(_nodes[nn(i+1,j  )-1]);
                     el->Add(_nodes[nn(i+1,j+1)-1]);
                     Add(el);
                     el = new Element(++ne,TRIANGLE,1);
                     el->Add(_nodes[nn(i+1,j+1)-1]);
                     el->Add(_nodes[nn(i  ,j+1)-1]);
                     el->Add(_nodes[nn(i  ,j  )-1]);
                     Add(el);
                  }
               }
            }
         }
         break;

      case 3:
         if (opt==HEXAHEDRON) {
            for (k=1; k<=nz; k++) {
               for (j=1; j<=ny; j++) {
                  for (i=1; i<=nx; i++) {
                     if (g.isActive(i,j,k)) {
                        el = new Element(++ne,HEXAHEDRON,1);
                        el->Add(_nodes[nn(i  ,j  ,k  )-1]);
                        el->Add(_nodes[nn(i+1,j  ,k  )-1]);
                        el->Add(_nodes[nn(i+1,j+1,k  )-1]);
                        el->Add(_nodes[nn(i  ,j+1,k  )-1]);
                        el->Add(_nodes[nn(i  ,j  ,k+1)-1]);
                        el->Add(_nodes[nn(i  ,j+1,k+1)-1]);
                        el->Add(_nodes[nn(i+1,j+1,k+1)-1]);
                        el->Add(_nodes[nn(i  ,j+1,k+1)-1]);
                        Add(el);
                     }
                  }
               }
            }
         }
         break;
   }
   NumberEquations();
   mesh_elements(*this)
      theMaterial.check(the_element->getCode());
}


Mesh::Mesh(const Grid& g,
           int         shape,
           int         opt)
     : _nb_nodes(0), _nb_elements(0), _nb_sides(0), _nb_edges(0), _dim(3), _nb_dof(0),
       _nb_vertices(0), _first_dof(1), _nb_mat(1), _verb(0), _no_imposed_dof(false),
       _is_structured(true), _all_sides_created(false), _boundary_sides_created(false),
       _all_edges_created(false), _boundary_edges_created(false), _boundary_nodes_created(false),
       _node_neighbor_elements_created(false), _element_neighbor_elements_created(false),
       _n_view1(0), _n_view2(0), _e_view1(0), _e_view2(0),
       _s_view1(0), _s_view2(0), _ed_view1(0), _ed_view2(0)
{
   Node *nd;
   Element *el;
   size_t i, j, k, n1, n2;
   size_t nx=g.getNx(), ny=g.getNy(), nz=g.getNz();
   setNodesForDOF();
   size_t n=0, ne=0;
   if (nz==0)
      _dim = 2;
   if (ny==0)
      _dim = 1;
   _boundary_nodes_created = false;
   if (_dim==3)
      shape = HEXAHEDRON;

// Nodes
   switch (_dim) {

      case 1:
         for (i=1; i<=nx; i++) {
            Point<real_t> x = 0.5*(g.getCoord(i)+g.getCoord(i+1));
            nd = new Node(++n,x);
            _code[0] = g.getCode(i);
            nd->setNbDOF(1);
            nd->setDOF(_first_dof,1);
            nd->setCode(_code);
            Add(nd);
         }
         break;

      case 2:
         for (i=1; i<=nx; i++) {
            for (j=1; j<=ny; j++) {
              Point<real_t> x = 0.5*(g.getCoord(i,j)+g.getCoord(i+1,j+1));
               nd = new Node(++n,x);
               _code[0] = g.getCode(i,j);
               nd->setNbDOF(1);
               nd->setDOF(_first_dof,1);
               nd->setCode(_code);
               Add(nd);
            }
         }
         break;

      case 3:
         for (i=1; i<=nx; i++) {
            for (j=1; j<=ny; j++) {
               for (k=1; k<=nz; k++) {
                 Point<real_t> x = 0.5*(g.getCoord(i,j,k)+g.getCoord(i+1,j+1,k+1));
                  nd = new Node(++n,x);
                  _code[0] = g.getCode(i,j,k);
                  nd->setNbDOF(1);
                  nd->setDOF(_first_dof,1);
                  nd->setCode(_code);
                  Add(nd);
               }
            }
         }
         break;
   }

// Elements
   n = ne = 0;
   size_t nn = 1;
   switch (_dim) {

      case 1:
         for (i=1; i<=nx-1; i++) {
            el = new Element(++ne,LINE,1);
            el->Add(_nodes[i-1]);
            el->Add(_nodes[i]);
            Add(el);
         }
         break;

      case 2:
         if (shape==QUADRILATERAL) {
            for (j=1; j<=ny-1; j++) {
               for (i=1; i<=nx-1; i++) {
                  n1 = nn + i - 1;
                  el = new Element(++ne,QUADRILATERAL,1);
                  el->Add(_nodes[n1-1]);
                  el->Add(_nodes[n1]);
                  el->Add(_nodes[n1+nx+1]);
                  el->Add(_nodes[n1+nx]);
                  Add(el);
               }
               nn += nx;
            }
         }
         else if (shape==TRIANGLE) {
            for (i=1; i<=nx-1; i++) {
               for (j=1; j<=ny-1; j++) {
                  n1 = nn + i - 1;
                  el = new Element(++ne,TRIANGLE,1);
                  el->Add(_nodes[n1-1]);
                  el->Add(_nodes[n1+ny-1]);
                  el->Add(_nodes[n1]);
                  Add(el);
                  el = new Element(++ne,TRIANGLE,1);
                  el->Add(_nodes[n1+ny-1]);
                  el->Add(_nodes[n1]);
                  el->Add(_nodes[n1-1]);
                  Add(el);
               }
               nn += ny;
            }
         }
         break;

      case 3:
         if (shape==HEXAHEDRON) {
            for (k=1; k<=nz-1; k++) {
               for (j=1; j<=ny-1; j++) {
                  for (i=1; i<=nx-1; i++) {
                     n1 = nn + i - 1;
                     n2 = n1 + (nx+1)*(ny+1);
                     el = new Element(++ne,HEXAHEDRON,1);
                     el->Add(_nodes[n1-1]);
                     el->Add(_nodes[n1]);
                     el->Add(_nodes[n1+nx+1]);
                     el->Add(_nodes[n1+nx]);
                     el->Add(_nodes[n2-1]);
                     el->Add(_nodes[n2]);
                     el->Add(_nodes[n2+nx+1]);
                     el->Add(_nodes[n2+nx]);
                     Add(el);
                  }
                  nn += nx + 1;
               }
               nn += ny + 1;
            }
         }
         break;
   }
   NumberEquations();
   mesh_elements(*this)
      theMaterial.check(the_element->getCode());
}


Mesh::Mesh(real_t Lx,
           real_t Ly,
           size_t nx,
           size_t ny,
           int    c1,
           int    c2,
           int    c3,
           int    c4,
           int    opt)
     : _nb_nodes(0), _nb_elements(0), _nb_sides(0), _nb_edges(0), _dim(2), _nb_dof(0), _nb_vertices(0),
       _first_dof(1), _nb_mat(1), _verb(0), _no_imposed_dof(false),
       _is_structured(false), _all_sides_created(false), _boundary_sides_created(false),
       _all_edges_created(false), _boundary_edges_created(false), _boundary_nodes_created(false),
       _node_neighbor_elements_created(false), _element_neighbor_elements_created(false),
       _n_view1(0), _n_view2(0), _e_view1(0), _e_view2(0),
       _s_view1(0), _s_view2(0), _ed_view1(0), _ed_view2(0)
{
   Node *nd;
   Element *el;
   theMaterial.set(1,"Generic");
   setNodesForDOF();
   real_t hx, hy;
   size_t n=0, ne=0;
   Point<real_t> x;
   if (ny==0)
      _dim = 1;
   x.x = 0.;
   hx = Lx/nx;
   if (ny==0)
      hy = 0;
   else
      hy = Ly/ny;

   for (size_t i=0; i<=nx; i++) {
      x.y = 0.;
      for (size_t j=0; j<=ny; j++) {
         n++;
         nd = new Node(n,x);
         nd->setNbDOF(1);
         _code[0] = 0;
         if (j==0  && c1>0)
            _code[0] = c1;
         if (j==ny && c3>0)
            _code[0] = c3;
         if (i==0  && c4>0)
            _code[0] = c4;
         if (i==nx && c2>0)
            _code[0] = c2;
         nd->setDOF(_first_dof,1);
         nd->setCode(_code);
         Add(nd);
         x.y += hy;
         if (i>0 && j>0 && opt) {
            try {
               if (opt==QUADRILATERAL) {
                  ne++;
                  el = new Element(ne,QUADRILATERAL,1);
                  el->Add(_nodes[n-ny-3]);
                  el->Add(_nodes[n-2]);
                  el->Add(_nodes[n-1]);
                  el->Add(_nodes[n-ny-2]);
                  Add(el);
               }
               else if (opt==TRIANGLE) {
                  ne++;
                  el = new Element(ne,TRIANGLE,1);
                  el->Add(_nodes[n-ny-3]);
                  el->Add(_nodes[n-2]);
                  el->Add(_nodes[n-1]);
                  Add(el);
                  ne++;
                  el = new Element(ne,TRIANGLE,1);
                  el->Add(_nodes[n-1]);
                  el->Add(_nodes[n-ny-2]);
                  el->Add(_nodes[n-ny-3]);
                  Add(el);
               }
               else
                  THROW_RT("Mesh(real_t,real_t,size_t,size_t,int,int,int,int,int): Illegal option "+itos(opt));
            }
            CATCH_EXIT("Mesh");
         }
      }
      x.x += hx;
      _is_structured = true;
   }

   if (c1<0 || c2<0 || c3<0 || c4<0) {
      size_t is=1, n=1;
      for (size_t j=1; j<=ny; j++) {
         Side *sd = new Side(is++,LINE);
         sd->Add(getPtrNode(n));
         sd->Add(getPtrNode(n+1));
         n++;
         sd->setNbDOF(1);
         if (c4<0) {
            sd->setCode(1,-c4);
            Add(sd);
	 }
      }
      n = nx*(ny+1) + 1;
      for (size_t j=1; j<=ny; j++) {
         Side *sd = new Side(is++,LINE);
         sd->Add(getPtrNode(n));
         sd->Add(getPtrNode(n+1));
         n++;
         sd->setNbDOF(1);
         if (c2<0) {
            sd->setCode(1,-c2);
            Add(sd);
	 }
      }
      n = 1;
      for (size_t i=1; i<=nx; i++) {
         Side *sd = new Side(is++,LINE);
         sd->Add(getPtrNode(n));
         sd->Add(getPtrNode(n+ny+1));
         n += ny + 1;
         sd->setNbDOF(1);
         if (c1<0) {
            sd->setCode(1,-c1);
            Add(sd);
	 }
      }
      n = ny + 1;
      for (size_t i=1; i<=nx; i++) {
         Side *sd = new Side(is++,LINE);
         sd->Add(getPtrNode(n));
         sd->Add(getPtrNode(n+ny+1));
         n += ny + 1;
         sd->setNbDOF(1);
         if (c3<0) {
            sd->setCode(1,-c3);
            Add(sd);
	 }
      }
   }

   NumberEquations();
   mesh_elements(*this)
      theMaterial.check(the_element->getCode());
}


Mesh::Mesh(const Mesh& ms)
     : _nb_nodes(0), _nb_elements(0), _nb_sides(0),
       _nb_edges(0), _nb_side_nodes(ms._nb_side_nodes),
       _nb_element_nodes(ms._nb_element_nodes), _dim(ms._dim), _nb_dof(ms._nb_dof),
       _nb_vertices(ms._nb_vertices), _first_dof(ms._first_dof), _nb_mat(ms._nb_mat),
       _nb_eq(ms._nb_eq), _max_nb_nodes(ms._max_nb_nodes), _max_nb_elements(ms._max_nb_elements),
       _max_nb_sides(ms._max_nb_sides), _max_nb_edges(ms._max_nb_edges),
       _verb(ms._verb), _no_imposed_dof(ms._no_imposed_dof), _is_structured(ms._is_structured),
       _all_sides_created(ms._all_sides_created), _boundary_sides_created(ms._boundary_sides_created),
       _all_edges_created(ms._all_edges_created), _boundary_edges_created(ms._boundary_edges_created),
       _boundary_nodes_created(ms._boundary_nodes_created),
       _node_neighbor_elements_created(ms._node_neighbor_elements_created),
       _element_neighbor_elements_created(false),
       _n_view1(ms._n_view1), _n_view2(ms._n_view2), _e_view1(ms._e_view1), _e_view2(ms._e_view2),
       _s_view1(ms._s_view1), _s_view2(ms._s_view2), _ed_view1(ms._ed_view1), _ed_view2(ms._ed_view2)
{
// Insert nodes
   mesh_nodes(ms) {
      Add(new Node(The_node));
      _nb_eq = _nb_dof;
   }

// Insert elements
   mesh_elements(ms)
      Add(new Element(The_element));

// Insert sides
   mesh_sides(ms)
      Add(new Side(side_label,The_side.getShape()));
   
   if (ms._node_in_coarse_element.size()>0) {
      for (size_t i=0; i<_nb_nodes; i++)
         _node_in_coarse_element.push_back(ms._node_in_coarse_element[i]);
   }
   if (ms._node_in_fine_element.size()>0) {
      for (size_t i=0; i<_nb_nodes; i++)
         _node_in_fine_element.push_back(ms._node_in_fine_element[i]);
   }

   _nb_mat = ms._nb_mat;
   for (size_t k=0; k<_nb_mat; k++) {
      _code_mat[k] = ms._code_mat[k];
      _mat[k] = ms._mat[k];
   }
   _set_nodes = ms._set_nodes;
   _set_sides = ms._set_sides;
   _set_elements = ms._set_elements;
   _no_imposed_dof = ms._no_imposed_dof;
}


Mesh::Mesh(const Mesh&          m,
           const Point<real_t>& x_bl,
           const Point<real_t>& x_tr)
     : _nb_nodes(0), _nb_elements(0), _nb_sides(0), _nb_edges(0), _dim(2), _nb_dof(0),
       _nb_vertices(0), _first_dof(1), _nb_mat(1), _verb(0), _no_imposed_dof(false),
       _is_structured(false), _all_sides_created(false), _boundary_sides_created(false),
       _all_edges_created(false), _boundary_edges_created(false), _boundary_nodes_created(false),
       _node_neighbor_elements_created(false), _element_neighbor_elements_created(false),
       _n_view1(0), _n_view2(0), _e_view1(0), _e_view2(0),
       _s_view1(0), _s_view2(0), _ed_view1(0), _ed_view2(0)
{
   try {
      if (m.getDim()!=2)
         THROW_RT("Mesh(Mesh,Point<real_t>,Point<real_t>)\n"
                  "This constructor is valid for 2-D meshes only.");
   }
   CATCH_EXIT("Mesh");
   _max_nb_nodes = m.getNbNodes();
   _max_nb_elements = m.getNbElements();
   _max_nb_sides = m.getNbSides();
   setNodesForDOF();

   _node_old_label.resize(m.getNbNodes());
   _node_new_label.resize(m.getNbNodes());
   size_t label = 1;
   mesh_nodes(m) {
      Point<real_t> a = the_node->getCoord();
      if (a.x>=x_bl.x && a.y>=x_bl.y && a.x<=x_tr.x && a.y <= x_tr.y) {
         Node *nd = new Node(The_node);
         nd->setLabel(label++);
         _node_old_label[nd->n()-1] = node_label;
         _node_new_label[node_label-1] = nd->n();
         Add(nd);
      }
   }

   _element_old_label.resize(m.getNbElements());
   label = 1;
   mesh_elements(m) {
      size_t in=0;
      for (size_t i=1; i<=the_element->getNbNodes(); i++) {
         Point<real_t> a = (*the_element)(i)->getCoord();
         if (a.x>=x_bl.x && a.y>=x_bl.y && a.x<=x_tr.x && a.y <= x_tr.y)
            in++;
      }
      if (in==the_element->getNbNodes()) {
         Element *el = new Element(*the_element);
         el->setLabel(label++);
         _element_old_label[el->n()-1] = element_label;
         for (size_t i=1; i<=el->getNbNodes(); i++)
            el->Replace(i,getPtrNode(_node_new_label[el->getNodeLabel(i)-1]));
         Add(el);
      }
   }
}


Mesh::Mesh(const Mesh&  mesh,
           int          opt,
           size_t       dof1,
           size_t       dof2,
           bool         bc)
     : _n_view1(0), _n_view2(0), _e_view1(0), _e_view2(0),
       _s_view1(0), _s_view2(0), _ed_view1(0), _ed_view2(0)
{
   Node *nd;
   Side *sd;
   _first_dof = 1;
   size_t nb_dof = dof2 - dof1 + 1;
   _dim = mesh.getDim();
   if (opt==NODE_DOF)
      setNodesForDOF();
   else if (opt==ELEMENT_DOF)
      setElementsForDOF();
   else if (opt==SIDE_DOF)
      setSidesForDOF();

   mesh_elements(mesh)
      Add(new Element(The_element));

   if (opt==NODE_DOF) {
      mesh_nodes(mesh) {
         nd = new Node(The_node);
         nd->setNbDOF(nb_dof);
         for (size_t i=0; i<nb_dof; i++)
            nd->setCode(i+1,the_node->getCode(dof1+i));
         nd->setDOF(_first_dof,nb_dof);
         Add(nd);
      }
   }
   else {
      mesh_nodes(mesh)
         Add(new Node(The_node));
   }

   if (opt==SIDE_DOF) {
      mesh_sides(mesh) {
         sd = new Side(The_side);
         for (size_t i=0; i<nb_dof; i++)
            sd->setCode(i+1,the_side->getCode(dof1+i));
         sd->setDOF(_first_dof,nb_dof);
         Add(sd);
      }
   }
   else {
      mesh_sides(mesh)
         Add(new Side(The_side));
   }
   NumberEquations();
}


Mesh::~Mesh()
{
   if (_verb>2)
      cout << "Removing Mesh instance ..." << endl;
   for (size_t i=0; i<_nodes.size(); i++) {
      if (_nodes[i])
         delete _nodes[i];
   }
   for (size_t i=0; i<_elements.size(); i++) {
      if (_elements[i])
         delete _elements[i];
   }
   for (size_t i=0; i<_sides.size(); i++) {
      if (_sides[i])
         delete _sides[i];
   }
   for (size_t i=0; i<_edges.size(); i++) {
      if (_edges[i])
         delete _edges[i];
   }
}


Mesh &Mesh::operator*=(real_t a)
{
   mesh_nodes(*this) {
      Point<real_t> x(The_node.getCoord()*a);
      The_node.setCoord(1,x.x);
      The_node.setCoord(2,x.y);
      The_node.setCoord(3,x.z);
   }
   return *this;
}


int Mesh::getShape() const
{
   int t1=0, t2=0;
   size_t i=0;
   mesh_elements(*this) {
      t1 = t2;
      if (The_element.getShape()==LINE)
         t2 = LINE;
      else if (The_element.getShape()==TRIANGLE)
         t2 = TRIANGLE;
      else if (The_element.getShape()==TETRAHEDRON)
         t2 = TETRAHEDRON;
      else if (The_element.getShape()==HEXAHEDRON)
         t2 = HEXAHEDRON;
      else if (The_element.getShape()==PENTAHEDRON)
         t2 = PENTAHEDRON;
      else if (The_element.getShape()==QUADRILATERAL)
         t2 = QUADRILATERAL;
      else
         t1 = 0;
      if (t2!=t1 && i!=0)
         return 0;
      i = 1;
   }
   return t1;
}


void Mesh::RenumberNode(size_t n1,
                        size_t n2)
{
   Node *nd = getPtrNode(n1);
   try {
      if (nd)
         nd->setLabel(n2);
      else
         THROW_RT("RenumberNode(size_t,size_t): Node with label " + itos(n1) + " does not exist.");
   }
   CATCH_EXIT("Mesh");
}


void Mesh::RenumberElement(size_t n1,
                           size_t n2)
{
   Element *el = getPtrElement(n1);
   try {
      if (el)
         el->setLabel(n2);
      else
         THROW_RT("RenumberElement(size_t,size_t): Element with label " + itos(n1) + " does not exist.");
   }
   CATCH_EXIT("Mesh");
}


void Mesh::RenumberSide(size_t n1,
                        size_t n2)
{
   Side *sd = getPtrSide(n1);
   try {
      if (sd)
         sd->setLabel(n2);
      else
         THROW_RT("RenumberSide(size_t,size_t): Side with label " + itos(n1) + " does not exist.");
   }
   CATCH_EXIT("Mesh");
}


void Mesh::RenumberEdge(size_t n1,
                        size_t n2)
{
   Edge *ed = getPtrEdge(n1);
   try {
      if (ed)
         ed->setLabel(n2);
      else
         THROW_RT("RenumberEdge(size_t,size_t): Edge with label " + itos(n1) + " does not exist.");
   }
   CATCH_EXIT("Mesh");
}


void Mesh::setNodeView(size_t n1,
                       size_t n2)
{
   _n_view1 = n1, _n_view2 = n2;
   if (n1==0)
      _n_view1 = 1;

   if (n2>_nb_nodes)
      _n_view2 = _nb_nodes;
}


void Mesh::setElementView(size_t n1,
                          size_t n2)
{
   _e_view1 = n1, _e_view2 = n2;
   if (n1==0)
      _e_view1 = 1;
   if (n2>_nb_elements)
      _e_view2 = _nb_elements;
}


void Mesh::setSideView(size_t n1,
                       size_t n2)
{
   _s_view1 = n1, _s_view2 = n2;
   if (n1==0)
      _s_view1 = 1;
   if (n2>_nb_sides)
      _s_view2 = _nb_sides;
}


void Mesh::setEdgeView(size_t n1,
                       size_t n2)
{
   _ed_view1 = n1, _ed_view2 = n2;
   if (n1==0)
      _ed_view1 = 1;
   if (n2>_nb_edges)
      _ed_view2 = _nb_edges;
}


void Mesh::inCoarse(Mesh& ms,
                    bool  test_el)
{
   size_t i;
   real_t x[3], y[3];

   for (i=0; i<_nb_nodes; i++)
      _node_in_coarse_element.push_back(NULL);

   mesh_elements(ms) {
      try {
         if (The_element.getShape() != TRIANGLE)
            THROW_RT("inCoarse(Mesh,bool): Element " + itos(element_label) + " is not a triangle.");
      }
      CATCH_EXIT("Mesh");
      for (i=0; i<3; i++) {
         x[i] = (The_element)(i+1)->getCoord(1);
         y[i] = (The_element)(i+1)->getCoord(2);
      }
      real_t d = (x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]);

      mesh_nodes(*this) {
         Point<real_t> c = The_node.getCoord();
         real_t s = (y[2]-y[0])*(c.x-x[0])+(x[0]-x[2])*(c.y-y[0]);
         real_t t = (y[0]-y[1])*(c.x-x[0])+(x[1]-x[0])*(c.y-y[0]);
         if ((s>=0) && (s<=d) && (t>=0) && ((s+t)<=d))
            _node_in_coarse_element[node_label-1] = the_element;
      }

      if (test_el)
         mesh_elements(*this) {
            size_t host = 0;
            for (size_t k=0; k<3; k++) {
               if (the_element==_node_in_coarse_element[(The_element)(k+1)->n()-1])
                  host++;
            }
            if (host>=2)
               _node_in_coarse_element[The_element(1)->n()-1] = the_element;
         }
   }
}


void Mesh::inFine(Mesh& ms,
                  bool  test_el)
{
   size_t i;
   real_t x[3], y[3];

   for (i=0; i<_nb_nodes; i++)
      _node_in_fine_element.push_back(NULL);

   mesh_elements(ms) {
      try {
         if (the_element->getShape() != TRIANGLE)
            THROW_RT("inFine(Mesh,bool): This function is valid for triangles only.");
      }
      CATCH_EXIT("Mesh");
      for (i=0; i<3; i++) {
         x[i] = The_element(i+1)->getCoord(1);
         y[i] = The_element(i+1)->getCoord(2);
      }
      real_t d = (x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]);

      mesh_nodes(*this) {
         Point<real_t> c = the_node->getCoord();
         real_t s = (y[2]-y[0])*(c.x-x[0])+(x[0]-x[2])*(c.y-y[0]);
         real_t t = (y[0]-y[1])*(c.x-x[0])+(x[1]-x[0])*(c.y-y[0]);
         if ((s>=0) && (s<=d) && (t>=0) && ((s+t)<=d))
            _node_in_fine_element[node_label-1] = the_element;
      }

      if (test_el)
	mesh_elements(*this) {
            size_t host = 0;
            for (size_t k=0; k<3; k++) {
               if (the_element==_node_in_fine_element[The_element(k+1)->n()-1])
                  host++;
            }
            if (host>=2)
               _node_in_fine_element[The_element(1)->n()-1] = the_element;
         }
   }
}


void Mesh::Delete(Node* nd)
{
   if (_verb > 5)
      cout << "Deleting node " << nd->n() << endl;
   try {
      if (!nd)
         THROW_RT("Delete(Node *): Node does not exist.");
   }
   CATCH_EXIT("Mesh");
}


void Mesh::Delete(Element* el)
{
   if (_verb > 5)
      cout << "Deleting element " << el->n() << endl;
   try {
      if (!el)
         THROW_RT("Delete(Element *): Element does not exist.");
   }
   CATCH_EXIT("Mesh");
}


void Mesh::Delete(Side* sd)
{
   if (_verb > 5)
      cout << "Deleting side " << sd->n() << endl;
   try {
      if (!sd)
         THROW_RT("Delete(Side *): Side does not exist.");
   }
   CATCH_EXIT("Mesh");
}


Point<real_t> Mesh::getMaxCoord() const
{
   Point<real_t> a;
   real_t t;
   Node *nd;

   a.x = getPtrNode(1)->getCoord(1);
   if (_dim > 1)
      a.y = getPtrNode(1)->getCoord(2);
   if (_dim > 2)
      a.z = getPtrNode(1)->getCoord(3);

   for (size_t i=2; i<=_nb_nodes; i++) {
      nd = getPtrNode(i);
      t = nd->getCoord(1);
      if (t > a.x)
         a.x = t;
      if (_dim > 1) {
         t = nd->getCoord(2);
         if (t > a.y)
           a.y = t;
      }
      if (_dim > 2) {
         t = nd->getCoord(3);
         if (t > a.z)
            a.z = t;
      }
   }
   return a;
}


Point<real_t> Mesh::getMinCoord() const
{
   Point<real_t> a;
   real_t t;
   Node *nd;

   a.x = getPtrNode(1)->getCoord(1);
   if (_dim > 1)
      a.y = getPtrNode(1)->getCoord(2);
   if (_dim > 2)
      a.z = getPtrNode(1)->getCoord(3);

   for (size_t i=2; i<=_nb_nodes; i++) {
      nd = getPtrNode(i);
      t = nd->getCoord(1);
      if (t < a.x)
         a.x = t;
      if (_dim > 1) {
         t = nd->getCoord(2);
         if (t < a.y)
            a.y = t;
      }
      if (_dim > 2) {
         t = nd->getCoord(3);
         if (t < a.z)
            a.z = t;
      }
   }
   return a;
}


void Mesh::setNbDOFPerNode(size_t nb_dof)
{
   setDOFSupport(NODE_DOF);
   mesh_nodes(*this)
      The_node.setNbDOF(nb_dof);
   NumberEquations();
}


size_t Mesh::NumberEquations(size_t dof,
                             int    c)
{
   _nb_eq = 0;
   _first_dof = 1;
   size_t dof_type = 0;
   if (_set_nodes==true)
      dof_type = NODE_DOF;
   else if (_set_sides==true)
      dof_type = SIDE_DOF;
   else if (_set_elements==true)
      dof_type = ELEMENT_DOF;
   else if (_set_edges==true)
      dof_type = EDGE_DOF;

// Node supported d.o.f.
   if (dof_type == NODE_DOF) {
      if (_no_imposed_dof==true) {
         if (_verb > 1)
            cout << "Eliminating imposed d.o.f. from list of equations ..." << endl;
         mesh_nodes(*this) {
            for (size_t i=1; i<=The_node.getNbDOF(); i++) {
               if (The_node.getCode(i) != c)
                  The_node._dof[i-1] = ++_nb_eq;
               else
                  The_node._dof[i-1] = 0;
            }
         }
      }
      else {
         if (_verb > 1)
            cout << "Numbering node supported d.o.f. ..." << endl;
         _nb_eq = 0;
         mesh_nodes(*this) {
            if (dof)
               The_node._dof[0] = ++_nb_eq;
            for (size_t i=1; i<=The_node.getNbDOF(); i++)
               The_node._dof[i-1] = ++_nb_eq;
            _nb_dof = _nb_eq;
         }
      }
   }

// Side supported d.o.f.
   else if (dof_type==SIDE_DOF) {
      getAllSides();
      _nb_dof = 0;
      if (_no_imposed_dof==true) {
         if (_verb > 1)
            cout << "Numbering side supported d.o.f. ..." << endl;
         size_t l=0;
         mesh_sides(*this) {
            for (size_t k=1; k<=_nb_side_nodes; k++) {
               for (size_t i=1; i<=the_side->getNbDOF(); i++) {
                  if (the_side->getCode(i) != c)
                     the_side->_dof[l++] = ++_nb_eq;
                  else
                     the_side->_dof[l++] = 0;
                  _nb_dof++;
               }
            }
         }
      }
      else {
         if (_verb > 1)
            cout << "Eliminating imposed d.o.f. from list of equations ..." << endl;
         _nb_dof = 0;
         mesh_sides(*this) {
            the_side->setDOF(_first_dof,the_side->getNbDOF());
            _nb_dof += the_side->getNbDOF();
         }
         _nb_eq = _nb_dof;
      }
   }

// Edge supported d.o.f.
   else if (dof_type==EDGE_DOF) {
      if (_no_imposed_dof==false) {
         if (_verb > 1)
            cout << "Numbering edge supported d.o.f. ..." << endl;
         if (dof) {
            mesh_edges(*this)
               the_edge->DOF(dof,++_nb_eq);
         }
         else {
            mesh_edges(*this) {
               for (size_t i=1; i<=the_edge->getNbDOF(); i++)
                  if (the_edge->getCode(i) != c)
                     the_edge->DOF(i,++_nb_eq);
                  else
                     the_edge->DOF(i,0);
            }
         }
      }
      else {
         if (_verb > 1)
            cout << "Eliminating imposed d.o.f. from list of equations ..." << endl;
         if (dof) {
            mesh_edges(*this)
               for (size_t i=1; i<=the_edge->getNbDOF(); i++)
                  if (the_edge->getCode(i) != c)
                     the_edge->DOF(i,++_nb_eq);
                  else
                     the_edge->DOF(i,0);
         }
         else {
            mesh_edges(*this)
               if (the_edge->getCode(dof) != c)
                  the_edge->DOF(dof,++_nb_eq);
               else
                  the_edge->DOF(dof,0);
         }
      }
   }

// Element supported d.o.f.
   else if (dof_type==ELEMENT_DOF) {
      getAllSides();
      getElementNeighborElements();
      if (_no_imposed_dof==false) {
         if (_verb > 1)
            cout << "Numbering element supported d.o.f. ..." << endl;
         if (dof)
            mesh_elements(*this)
               the_element->setDOF(dof,++_nb_eq);
         else
            mesh_elements(*this)
               for (size_t i=1; i<=the_element->getNbDOF(); i++)
                  the_element->setDOF(i,++_nb_eq);
      }
   }
   return _nb_eq;
}


size_t Mesh::NumberEquations(size_t dof)
{
   _nb_eq = 0;
   _first_dof = 1;
   size_t dof_type = 0;
   if (_set_nodes==true)
      dof_type = NODE_DOF;
   else if (_set_sides==true)
      dof_type = SIDE_DOF;
   else if (_set_elements==true)
      dof_type = ELEMENT_DOF;
   else if (_set_edges==true)
      dof_type = EDGE_DOF;

// Node supported d.o.f.
   if (dof_type == NODE_DOF) {
      if (_no_imposed_dof==true) {
         if (_verb > 1)
            cout << "Eliminating imposed d.o.f. from list of equations ..." << endl;
         mesh_nodes(*this) {
            for (size_t i=1; i<=the_node->getNbDOF(); i++) {
               if (The_node.getCode(i) == 0)
                  The_node._dof[i-1] = ++_nb_eq;
               else
                  The_node._dof[i-1] = 0;
            }
         }
      }
      else {
         if (_verb > 1)
            cout << "Numbering node supported d.o.f. ..." << endl;
         _nb_eq = 0;
         mesh_nodes(*this) {
            if (dof)
               The_node._dof[0] = ++_nb_eq;
            else {
               for (size_t i=1; i<=The_node.getNbDOF(); i++)
                  The_node._dof[i-1] = ++_nb_eq;
            }
            _nb_dof = _nb_eq;
         }
      }
   }

// Side supported d.o.f.
   else if (dof_type == SIDE_DOF) {
      getAllSides();
      _nb_dof = 0;
      if (_no_imposed_dof==true) {
         if (_verb > 1)
            cout << "Numbering side supported d.o.f. ..." << endl;
         size_t l=0;
         mesh_sides(*this) {
            for (size_t k=1; k<=_nb_side_nodes; k++) {
               for (size_t i=1; i<=the_side->getNbDOF(); i++) {
                  if (the_side->getCode(i) <= 0)
                     the_side->_dof[l++] = ++_nb_eq;
                  else
                     the_side->_dof[l++] = 0;
                  _nb_dof++;
	       }
            }
         }
      }
      else {
         if (_verb > 1)
            cout << "Eliminating imposed d.o.f. from list of equations ..." << endl;
         _nb_dof = 0;
         mesh_sides(*this) {
            the_side->setDOF(_first_dof,the_side->getNbDOF());
            _nb_dof += the_side->getNbDOF();
         }
         _nb_eq = _nb_dof;
      }
   }

// Edge supported d.o.f.
   else if (dof_type == EDGE_DOF) {
      if (_no_imposed_dof==false) {
         if (_verb > 1)
            cout << "Numbering edge supported d.o.f. ..." << endl;
         if (dof) {
            mesh_edges(*this)
               the_edge->DOF(dof,++_nb_eq);
         }
         else {
            mesh_edges(*this) {
               for (size_t i=1; i<=the_edge->getNbDOF(); i++)
                  if (the_edge->getCode(i) == 0)
                     the_edge->DOF(i,++_nb_eq);
                  else
                     the_edge->DOF(i,0);
            }
         }
      }
      else {
         if (_verb > 1)
            cout << "Eliminating imposed d.o.f. from list of equations ..." << endl;
         Edge *ed;
         if (dof) {
            for (topEdge(); (ed=getEdge());)
               for (size_t i=1; i<=ed->getNbDOF(); i++)
                  if (ed->getCode(i) == 0)
                     ed->DOF(i,++_nb_eq);
                  else
                     ed->DOF(i,0);
         }
         else {
            for (topEdge(); (ed=getEdge());)
               if (ed->getCode(dof) == 0)
                  ed->DOF(dof,++_nb_eq);
               else
                  ed->DOF(dof,0);
         }
      }
   }

// Element supported d.o.f.
   else if (dof_type == ELEMENT_DOF) {
      getElementNeighborElements();
      if (_no_imposed_dof==false) {
         if (_verb > 1)
            cout << "Numbering element supported d.o.f. ..." << endl;
         if (dof)
            mesh_elements(*this)
               the_element->setDOF(dof,++_nb_eq);
         else
            mesh_elements(*this)
               for (size_t i=1; i<=the_element->getNbDOF(); i++)
                  the_element->setDOF(i,++_nb_eq);
         _nb_dof = _nb_eq;
      }
   }
   return _nb_eq;
}


void Mesh::selectDOF(int    dof_type,
                     size_t dof1,
                     size_t dof2,
                     bool   bc)
{
   _dof_nbeq[dof1-1] = 0;
   _nb_eq = 0;
   _first_dof = dof1;
   _no_imposed_dof = bc;

// Node supported d.o.f.
   if (dof_type==NODE_DOF) {
      if (_no_imposed_dof==true) {
         if (_verb > 1)
            cout << "Eliminating imposed d.o.f. from list of equations ..." << endl;
         mesh_nodes(*this) {
            for (size_t i=dof1; i<=dof2; i++) {
               if (the_node->getCode(i) == 0) {
                  the_node->_dof[i-1] = ++_nb_eq;
                  _dof_nbeq[dof1-1]++;
               }
               else
                  the_node->_dof[i-1] = 0;
            }
         }
      }
      else {
         if (_verb > 1)
            cout << "Numbering node supported d.o.f. ..." << endl;
         _dof_nbeq[dof1-1] = _nb_nodes*(dof2-dof1+1);
         mesh_nodes(*this)
            the_node->setDOF(_first_dof,dof2-dof1+1);
      }
   }

// Side supported d.o.f.
   else if (dof_type==SIDE_DOF) {
      getAllSides();
      _nb_dof = 0;
      if (_no_imposed_dof==true) {
         if (_verb > 1)
            cout << "Numbering side supported d.o.f. ..." << endl;
         size_t l=0;
         mesh_sides(*this) {
            for (size_t k=1; k<=_nb_side_nodes; k++) {
               for (size_t i=dof1; i<=dof2; i++) {
                  if (the_side->getCode(i) <= 0) {
                     the_side->_dof[l++] = ++_nb_eq;
                     _dof_nbeq[dof1-1]++;
                  }
                  else
                     the_side->_dof[l++] = 0;
                  _nb_dof++;
               }
            }
         }
      }
      else {
         if (_verb > 1)
            cout << "Eliminating imposed d.o.f. from list of equations ..." << endl;
         _nb_dof = 0;
         mesh_sides(*this) {
            the_side->setDOF(_first_dof,the_side->getNbDOF());
            _nb_dof += the_side->getNbDOF();
         }
         _dof_nbeq[dof1-1] = _nb_eq = _nb_dof;
      }
   }

// Edge supported d.o.f.
   else if (dof_type==EDGE_DOF) {
      if (_no_imposed_dof==false) {
         if (_verb > 1)
            cout << "Numbering edge supported d.o.f. ..." << endl;
         mesh_edges(*this) {
            for (size_t i=dof1; i<=dof2; i++)
               if (the_edge->getCode(i) == 0) {
                  the_edge->DOF(i,++_nb_eq);
                  _dof_nbeq[dof1-1]++;
               }
               else
                  the_edge->DOF(i,0);
         }
      }
      else {
         if (_verb > 1)
            cout << "Eliminating imposed d.o.f. from list of equations ..." << endl;
         mesh_edges(*this) {
            for (size_t i=dof1; i<=dof2; i++) {
               if (the_edge->getCode(i) == 0) {
                  the_edge->DOF(i,++_nb_eq);
                  _dof_nbeq[dof1-1]++;
               }
               else
                  the_edge->DOF(i,0);
            }
         }
      }
   }

// Element supported d.o.f.
   else if (dof_type==ELEMENT_DOF) {
      getElementNeighborElements();
      if (_no_imposed_dof==false) {
         if (_verb > 1)
            cout << "Numbering element supported d.o.f. ..." << endl;
         mesh_elements(*this)
            for (size_t i=dof1; i<=dof2; i++) {
               the_element->setDOF(i,++_nb_eq);
              _dof_nbeq[dof1-1]++;
            }
      }
   }
}


int Mesh::getAllSides(int opt)
{
   try {
      if (_dim==1)
         THROW_RT("getAllSides(): No sides can be created for 1-D meshes");
   }
   CATCH("Mesh");
   if (_all_sides_created==true) {
      if (_verb > 0) {
         cout << "Attempting to create all mesh sides:\n";
         cout << "All sides were created already. Skipping !" << endl;
      }
      return 0;
   }
   size_t k, ns=0, i;
   vector<vector<size_t> > nsd;
   int sh;
   if (_verb > 1)
      cout << "Creating all mesh sides ..." << endl;
   size_t nb_all_sides = 0;

// Store already defined sides in (*sides)
   vector<ND> sides;
   SideSet SL(_sides);
   size_t nbss=_nb_sides;

// Add all mesh sides
   nb_all_sides = 0;
   mesh_elements(*this) {
      try {
         if (_dim==3 && the_element->getShape()!=TETRAHEDRON)
            THROW_RT("getAllSides(): Member function is not valid for this element type.");
      }
      CATCH_EXIT("Mesh");
      ns = init_side_node_numbering(The_element.getShape(),nsd,sh);
      for (i=0; i<The_element.getNbSides(); i++) {
         size_t i1 = The_element(nsd[i][0])->n(),
                i2 = The_element(nsd[i][1])->n(),
                i3 = 0, i4 = 0;
         if (ns > 2)
            i3 = The_element(nsd[i][2])->n();
         if (ns > 3)
            i4 = The_element(nsd[i][3])->n();
         ND ssd(i1,i2,i3,i4);
         ssd.e1 = element_label, ssd.e2 = 0;
         sides.push_back(ssd);
         nb_all_sides++;
      }
   }

   sides.erase(sides.begin()+nb_all_sides,sides.end());
   for (i=0; i<nb_all_sides; i++)
      order_side_nodes(ns,sides[i]);
   std::vector<ND>::iterator it = sides.begin()+nb_all_sides;
   sort(sides.begin(),it,compare_sides);
   complete_sides(sides);
   std::vector<ND>::iterator new_end = std::unique(sides.begin(),sides.end());
   sides.erase(new_end,sides.end());
   _nb_sides = 0;
   _sides.clear();
   for (i=0; i<sides.size(); i++) {
      sh = LINE;
      if (sides[i].n==3)
         sh = TRIANGLE;
      Side *sd = new Side(i+1,sh);
      for (k=0; k<sides[i].n; k++)
         sd->Add(getPtrNode(sides[i].nd[k]));
      sd->setNbDOF(getPtrNode(sides[i].nd[0])->getNbDOF());

      if (opt==0) {
         if (sides[i].e1 && sides[i].e2==0)
            sd->Add(getPtrElement(sides[i].e1));
         if (sides[i].e2 && sides[i].e1==0)
            sd->Add(getPtrElement(sides[i].e2));
         if (sides[i].e1 && sides[i].e2) {
            Element *el1=getPtrElement(sides[i].e1),
                    *el2=getPtrElement(sides[i].e2);
            if (el1->n() < el2->n()) {
               sd->Add(el1);
               sd->Add(el2);
            }
            else {
               sd->Add(el2);
               sd->Add(el1);
            }
         }
          
      }
      else {
         if (sides[i].e1)
            sd->Add(getPtrElement(sides[i].e1));
         if (sides[i].e2)
            sd->Add(getPtrElement(sides[i].e2));
      }
      Add(sd);
   }

// Assign sides to elements
   if (opt==0) {
      mesh_elements(*this)
         The_element.setGlobalToLocal();
      mesh_sides(*this) {
         Element *el1=The_side.getNeighborElement(1),
                 *el2=The_side.getNeighborElement(2);
         if (el1)
            el1->Add(the_side);
         if (el2)
            el2->Add(the_side);
      }
   }
   else {
      mesh_sides(*this) {
         if ((The_side.getNeighborElement(1)))
            The_side.getNeighborElement(1)->Add(the_side);
         if ((The_side.getNeighborElement(2)))
            The_side.getNeighborElement(2)->Add(the_side);
      }
   }

// Assign DOFs to sides
   size_t n=1;
   mesh_sides(*this) {
      The_side.setFirstDOF(n);
      for (size_t i=1; i<=The_side.getNbDOF(); i++)
         The_side.DOF(i,n++);
   }

// Mark nodes on boundary
   mesh_sides(*this) {
      if (The_side.isOnBoundary())
         for (size_t j=1; j<=The_side.getNbNodes(); j++)
            getPtrNode(The_side.getNodeLabel(j))->setOnBoundary();
   }

   ND ss1, ss2;
   mesh_sides(*this) {
      ss1.nd[0] = The_side(1)->n();
      ss1.nd[1] = The_side(2)->n();
      if (The_side.getNbNodes()>2)
         ss1.nd[2] = The_side(3)->n();
      if (The_side.getNbNodes()>3)
         ss1.nd[3] = The_side(4)->n();
      for (i=0; i<nbss; i++) {
         Side *sd = SL[i];
         ss2.nd[0] = sd->getNodeLabel(1);
         ss2.nd[1] = sd->getNodeLabel(2);
         if (sd->getNbNodes()>2)
            ss2.nd[2] = sd->getNodeLabel(3);
         if (sd->getNbNodes()>3)
            ss2.nd[3] = sd->getNodeLabel(4);
         order_side_nodes(sd->getNbNodes(),ss2);
         if (ss1==ss2) {
            The_side.setNbDOF(sd->getNbDOF());
            for (size_t j=1; j<=The_side.getNbDOF(); j++)
               The_side.setCode(j,sd->getCode(j));
         }
      }
   }
   _all_sides_created = true;
   createBoundarySideList();
   return _nb_sides;
}


int Mesh::getAllEdges()
{
   try {
      if (_dim==1)
         THROW_RT("getAllEdges(): No edges can be created for 1-D meshes");
   }
   CATCH("Mesh");
   try {
      if (_dim==2)
         THROW_RT("getAllEdges(): No edges can be created for 2-D meshes");
   }
   CATCH("Mesh");
   if (_all_edges_created == true) {
      if (_verb > 0) {
         cout << "Attempting to create all mesh edges: \n";
         cout << "All edges were created already. Skipping !" << endl;
      }
      return 0;
   }
   getAllSides();
   size_t i;
   if (_verb > 1)
      cout << "Creating all mesh edges ..." << endl;
   size_t nb_all_edges = 0;

// Store already defined edges in (*edges)
   vector<ND> edges;
   EdgeSet SL(_edges);
   size_t nbss = _nb_edges;

// Add all mesh edges
   vector<vector<size_t> > nsd(3);
   nsd[0].push_back(1); nsd[0].push_back(2);
   nsd[1].push_back(2); nsd[1].push_back(3);
   nsd[2].push_back(3); nsd[2].push_back(1);
   ND eed;
   nb_all_edges = 0;
   mesh_sides(*this) {
      try {
         if (the_side->getShape()!=TRIANGLE)
            THROW_RT("getAllEdges(): Edges can be created for triangular faces only");
      }
      CATCH("Mesh");
      for (i=0; i<3; i++) {
         eed = ND(The_side(nsd[i][0])->n(),The_side(nsd[i][1])->n(),0,0);
         eed.e1 = side_label;
         eed.e2 = 0;
         edges.push_back(eed);
         nb_all_edges++;
      }
   }
   edges.erase(edges.begin()+nb_all_edges,edges.end());
   for (i=0; i<nb_all_edges; i++)
      order_edge_nodes(edges[i]);
   std::vector<ND>::iterator it = edges.begin()+nb_all_edges;
   sort(edges.begin(),it,compare_sides);
   complete_sides(edges);
   std::vector<ND>::iterator new_end = std::unique(edges.begin(),edges.end());
   edges.erase(new_end,edges.end());
   _nb_edges = 0;
   _edges.clear();
   for (i=0; i<edges.size(); i++) {
      Edge *ed = new Edge(i+1);
      ed->Add(getPtrNode(edges[i].nd[0]));
      ed->Add(getPtrNode(edges[i].nd[1]));
      ed->setNbDOF(getPtrNode(edges[i].nd[0])->getNbDOF());
      if (edges[i].e1)
         ed->AddNeighbor(getPtrSide(edges[i].e1));
      if (edges[i].e2)
         ed->AddNeighbor(getPtrSide(edges[i].e2));
      Add(ed);
   }

// Assign sides to elements
    mesh_edges(*this) {
      if ((the_edge->getNeighborSide(1)))
         the_edge->getNeighborSide(1)->Add(the_edge);
      if ((the_edge->getNeighborSide(2)))
         the_edge->getNeighborSide(2)->Add(the_edge);
   }

// Assign DOFs to edges
   size_t n = 1;
   mesh_edges(*this) {
      the_edge->setFirstDOF(n);
      for (size_t i=1; i<=the_edge->getNbDOF(); i++)
         the_edge->DOF(i,n++);
   }

// Mark nodes on boundary
   mesh_edges(*this) {
      if (the_edge->isOnBoundary()) {
         getPtrNode(the_edge->getNodeLabel(1))->setOnBoundary();
         getPtrNode(the_edge->getNodeLabel(2))->setOnBoundary();
      }
   }

//
   ND ss1, ss2;
   mesh_edges(*this) {
      ss1.nd[0] = the_edge->getNodeLabel(1);
      ss1.nd[1] = the_edge->getNodeLabel(2);
      for (i=0; i<nbss; i++) {
         Edge *ed = SL[i];
         ss2.nd[0] = ed->getNodeLabel(1);
         ss2.nd[1] = ed->getNodeLabel(2);
         order_side_nodes(2,ss2);
         if (ss1==ss2) {
            the_edge->setNbDOF(ed->getNbDOF());
            for (size_t j=1; j<=the_edge->getNbDOF(); j++)
               the_edge->setCode(j,ed->getCode(j));
         }
      }
   }
   _all_edges_created = true;
   return _nb_edges;
}


int Mesh::createBoundarySideList()
{
   _nb_boundary_sides = 0;
   mesh_sides(*this) {
      if (the_side->isOnBoundary()) {
         _boundary_sides.push_back(the_side);
         _nb_boundary_sides++;
      }
   }
   return _nb_boundary_sides;
}


int Mesh::createInternalSideList()
{
   _nb_internal_sides = 0;
   mesh_sides(*this) {
      if (the_side->isOnBoundary()==0) {
         _internal_sides.push_back(the_side);
         _nb_internal_sides++;
      }
   }
   return _nb_internal_sides;
}


int Mesh::getBoundarySides()
{
   try {
      if (_boundary_sides_created == true)
         THROW_RT("getBoundarySides(): Attempting to create boundary sides: \n"
                  "Boundary sides were created already. Skipping !");
   }
   CATCH("Mesh");
   size_t k, ns=0, i;
   vector<vector<size_t> > nsd;
   int sh;
   if (_verb > 1)
      cout << "Creating all mesh sides ..." << endl;
   size_t nb_all_sides = 0;

// Store already defined sides in (*sides)
   vector<ND> sides;
   SideSet SL(_sides);
   size_t nbss = _nb_sides;

// Add all mesh sides
   size_t i3, i4;
   nb_all_sides = 0;

   vector<size_t> s(2);
   mesh_elements(*this) {
     for (size_t i=1; i<=3; i++) {
        s[0] = The_element(i)->n(), s[1] = The_element(i%3+1)->n();
        mesh_sides(*this) {
           if (the_side->getNeighborElement(1)==NULL) {
              if (equal_sides(the_side,s)) {
                 The_side.set(the_element,1);
                 break;
              }
           }
        }
     }
   }

   mesh_elements(*this) {
      try {
         if (_dim==3 && the_element->getShape() != TETRAHEDRON)
            THROW_RT("getBoundarySides(): Member function is not valid for this element type.");
      }
      CATCH("Mesh");
      ns = init_side_node_numbering(the_element->getShape(), nsd, sh);
      for (i=0; i<the_element->getNbSides(); i++) {
         i3 = i4 = 0;
         if (ns>2)
            i3 = The_element(nsd[i][2])->n();
         if (ns > 3)
            i4 = The_element(nsd[i][3])->n();
         ND nnd(The_element(nsd[i][0])->n(),The_element(nsd[i][1])->n(),i3,i4);
         nnd.e1 = element_label, nnd.e2 = 0;
         sides.push_back(nnd);
         nb_all_sides++;
      }
   }

   for (i=0; i<nb_all_sides; i++)
      order_side_nodes(ns,sides[i]);
   std::vector<ND>::iterator it = sides.begin() + nb_all_sides;
   sort(sides.begin(),it,compare_sides);
   complete_sides(sides);
   sides.erase(sides.begin()+nb_all_sides,sides.end());
   size_t nns = remove_internal_sides(sides);

   _nb_sides = 0;
   for (i=0; i<nns; i++) {
      Side *sd = new Side(i+1,sh);
      for (k=0; k<sides[i].n; k++)
         sd->Add(getPtrNode(sides[i].nd[k]));
      sd->setNbDOF(getPtrNode(sides[i].nd[0])->getNbDOF());
      if (sides[i].e1)
         sd->Add(getPtrElement(sides[i].e1));
      if (sides[i].e2)
         sd->Add(getPtrElement(sides[i].e2));
      Add(sd);
   }

   mesh_sides(*this) {
      Element *el;
      if ((el=The_side.getNeighborElement(1)))
         The_side.getNeighborElement(1)->Add(the_side);
      if ((el=The_side.getNeighborElement(2)))
         The_side.getNeighborElement(2)->Add(the_side);
   }

   size_t n = 1;
   mesh_sides(*this) {
      The_side.setFirstDOF(n);
      for (size_t i=1; i<=The_side.getNbDOF(); i++)
         The_side.DOF(i,n++);
   }

   mesh_sides(*this) {
      if (The_side.isOnBoundary()) {
         for (size_t j=1; j<=The_side.getNbNodes(); j++)
            getPtrNode(The_side(j)->n())->setOnBoundary();
      }
   }

   ND ss1, ss2;
   mesh_sides(*this) {
      ss1.nd[0] = The_side(1)->n();
      ss1.nd[1] = The_side(2)->n();
      if (the_side->getNbNodes()>2)
         ss1.nd[2] = The_side(3)->n();
      if (the_side->getNbNodes()>3)
         ss1.nd[3] = The_side(4)->n();
      for (i=0; i<nbss; i++) {
         Side *sd = SL[i];
         ss2.nd[0] = (*sd)(1)->n();
         ss2.nd[1] = (*sd)(2)->n();
         if (sd->getNbNodes()>2)
            ss2.nd[2] = (*sd)(3)->n();
         if (sd->getNbNodes()>3)
            ss2.nd[3] = (*sd)(4)->n();
         order_side_nodes(sd->getNbNodes(),ss2);
         if (ss1==ss2) {
            The_side.setNbDOF(sd->getNbDOF());
            for (size_t j=1; j<=the_side->getNbDOF(); j++)
               The_side.setCode(j,sd->getCode(j));
         }
      }
   }

   _boundary_sides_created = true;
   _nb_boundary_sides = 0;
   mesh_sides(*this) {
      if (The_side.isOnBoundary()) {
         _boundary_sides.push_back(the_side);
         _nb_boundary_sides++;
      }
   }
   return _nb_sides;
}


int Mesh::getBoundaryNodes()
{
   if (_boundary_nodes_created)
      return _nb_boundary_nodes;
   getBoundarySides();
   _nb_boundary_nodes = 0;
   mesh_nodes(*this) {
      if (the_node->isOnBoundary()) {
         _boundary_nodes.push_back(the_node);
         _nb_boundary_nodes++;
      }
   }
   _boundary_nodes_created = true;
   return _nb_boundary_nodes;
}


void Mesh::getNodeNeighborElements()
{
   if (_node_neighbor_elements_created) {
      if (_verb > 1)
         cout << "List of node neighbor elements already created." << endl;
      return;
   }
   if (_verb > 1)
      cout << "Creating node neighbor elements ..." << endl;
/*   mesh_elements(*this)
      for (size_t i=1; i<=the_element->getNbNodes(); i++)
         The_element(i)->Neig();*/
   mesh_nodes(*this) {
      The_node.Add();
   }
   mesh_elements(*this)
      for (size_t i=1; i<=The_element.getNbNodes(); i++)
         The_element(i)->Add(the_element);
   _node_neighbor_elements_created = true;
}


void Mesh::getElementNeighborElements()
{
   if (_element_neighbor_elements_created) {
      if (_verb > 1)
         cout << "List of element neighbor elements already created." << endl;
      return;
   }
   if (_verb > 1)
      cout << "Creating element neighbor elements ..." << endl;
   getAllSides();
   mesh_sides(*this) {
      size_t n = side_label;
      Element *el1=The_side.getNeighborElement(1),
              *el2=The_side.getNeighborElement(2);
      if (el1 && el2) {
         for (size_t k=1; k<=el1->getNbSides(); k++)
            if (el1->getPtrSide(k))
               if (el1->getSideLabel(k)==n)
                 el1->set(el2,k);
         for (size_t k=1; k<=el2->getNbSides(); k++)
            if (el2->getPtrSide(k))
               if (el2->getSideLabel(k)==n)
                  el2->set(el1,k);
      }
   }
   _element_neighbor_elements_created = true;
}


void Mesh::setDiscontinuous(size_t p)
{
   p = 0;
   _nb_vertices = _nb_nodes;
   size_t k=1;
   _nb_nodes = 0;
   mesh_elements(*this) {
      for (size_t i=1; i<=The_element.getNbNodes(); i++) {
         The_element(i)->setLabel(k++);
         _nb_nodes++;
      }
   }
}


int Mesh::getDOFSupport() const
{
   if (_set_nodes)
      return NODE_DOF;
   else if (_set_sides)
      return SIDE_DOF;
   else if (_set_elements)
      return ELEMENT_DOF;
   else if (_set_edges)
      return EDGE_DOF;
   else
      return 0;
}


void Mesh::setDOFSupport(int opt,
                         int nb_nodes)
{
   _set_nodes = _set_sides = _set_edges = _set_elements = false;
   if (opt==NODE_DOF) {
      _set_nodes = true;
      _set_sides = _set_elements = _set_edges = false;
   }
   else if (opt==SIDE_DOF) {
      _set_sides = true;
      _nb_side_nodes = nb_nodes;
      _set_nodes = _set_elements = _set_edges = false;
   }
   else if (opt==EDGE_DOF) {
      _set_edges = true;
      _set_nodes = _set_sides = _set_elements = false;
   }
   else if (opt==ELEMENT_DOF) {
      _set_elements = true;
      _nb_element_nodes = nb_nodes;
      _set_nodes = _set_sides = _set_edges = false;
   }
   else
      ;
   NumberEquations();
}


void Mesh::AddMidNodes(int g)
{
   Point<real_t> x;
   size_t mid_label = _nb_nodes + 1;
   vector<int> code;
   getAllSides();

   mesh_sides(*this) {
      x.x = 0.5*(The_side(1)->getCoord(1) + The_side(2)->getCoord(1));
      x.y = 0.5*(The_side(1)->getCoord(2) + The_side(2)->getCoord(2));
      size_t nb_dof = std::max(The_side(1)->getNbDOF(),The_side(2)->getNbDOF());
      Node *nd = new Node(mid_label++,x);
      nd->setNbDOF(nb_dof);

      for (size_t i=0; i<The_side(1)->getNbDOF(); i++)
         code.push_back(std::min(The_side(1)->getCode(i+1),The_side(2)->getCode(i+1)));
      nd->setCode(code);
      nd->setDOF(_first_dof,nb_dof);
      Add(nd);
      the_side->Add(nd);
   }

   code.clear();
   mesh_elements(*this) {
      if (g) {
         x = 0;
         for (size_t j=1; j<=the_element->getNbNodes(); j++)
            x += The_element(j)->getCoord();
         x /= the_element->getNbNodes();
         Node *nd = new Node(mid_label++,x);
         nd->setNbDOF(The_element(1)->getNbDOF());
         for (size_t i=0; i<The_element(1)->getNbDOF(); i++)
            code.push_back(The_element(1)->getCode(i+1));
         nd->setCode(code);
         nd->setDOF(_first_dof,The_element(1)->getNbDOF());
         Add(nd);
         the_element->Add(nd);
      }
      size_t m1=The_element(1)->n(), m2=The_element(2)->n(), m3=The_element(3)->n();
      while (the_element->getNbNodes()!=6) {
         for (size_t j=1; j<=the_element->getNbSides(); j++) {
            Side *sd = the_element->getPtrSide(j);
            size_t n1=sd->getNodeLabel(1), n2=sd->getNodeLabel(2);
            if ((n1==m1 && n2==m2) || (n1==m2 && n2==m1))
               the_element->Add(the_element->getPtrSide(j)->getPtrNode(3),4);
            else if ((n1==m2 && n2==m3) || (n1==m3 && n2==m2))
               the_element->Add(the_element->getPtrSide(j)->getPtrNode(3),5);
            else if ( (n1==m3 && n2==m1) || (n1==m1 && n2==m3))
               the_element->Add(the_element->getPtrSide(j)->getPtrNode(3),6);
         }
      }
   }
   NumberEquations();
}


void Mesh::set(Node* nd)
{
   size_t n = nd->n();
   try {
      if (n>_nb_nodes)
         THROW_RT("set(Node *): The node label is beyond the total number of nodes.");
   }
   CATCH("Mesh");
   _nodes[n-1] = nd;
}


void Mesh::set(Element* el)
{
   size_t n = el->n();
   try {
      if (n>_nb_elements)
         THROW_RT("set(Element *): The element label is beyond the total number of elements.");
   }
   CATCH_EXIT("Mesh");
   _elements[n-1] = el;
}


void Mesh::set(Side* sd)
{
   size_t n = sd->n();
   try {
      if (n>_nb_sides)
         THROW_RT("set(Side *): The side label is beyond the total number of sides.");
   }
   CATCH_EXIT("Mesh");
   _sides[n-1] = sd;
}


void Mesh::Add(Node* nd)
{
   _nodes.push_back(nd);
   _nb_nodes++; _nb_vertices++;
   _nb_dof += nd->getNbDOF();
}


void Mesh::Add(Element* el)
{
   _elements.push_back(el);
   _nb_elements++;
}


void Mesh::Add(Side* sd)
{
   _sides.push_back(sd);
   _nb_sides++;
}


void Mesh::Add(Edge* ed)
{
   _edges.push_back(ed);
   _nb_edges++;
}


void Mesh::Rescale(real_t sx,
                   real_t sy,
                   real_t sz)
{
   if (sz==0.)
      sz = sx;
   if (sy==0.)
      sy = sx;
   mesh_nodes(*this) {
      The_node.setCoord(1,sx*The_node.getCoord(1));
      The_node.setCoord(2,sy*The_node.getCoord(2));
      The_node.setCoord(3,sz*The_node.getCoord(3));
   }
}


void Mesh::getList(vector<Node *>& nl) const
{
   for (size_t i=0; i<_nb_nodes; i++)
      nl.push_back(_nodes[i]);
}


void Mesh::getList(vector<Element *>& el) const
{
   for (size_t i=0; i<_nb_elements; i++)
      el.push_back(_elements[i]);
}


void Mesh::getList(vector<Side *>& sl) const
{
   for (size_t i=0; i<_nb_sides; i++)
      sl.push_back(_sides[i]);
}


void Mesh::setList(const std::vector<Node *>& nl)
{
   _nodes.clear();
   _nb_nodes = 0;
   for (size_t i=0; i<nl.size(); i++)
      Add(nl[i]);
}


void Mesh::setList(const std::vector<Element *>& el)
{
   _elements.clear();
   _nb_elements = 0;
   for (size_t i=0; i<el.size(); i++)
      Add(el[i]);
}


void Mesh::setList(const std::vector<Side *>& sl)
{
   _sides.clear();
   _nb_sides = 0;
   for (size_t i=0; i<sl.size(); i++)
      Add(sl[i]);
}


void Mesh::checkNodeLabels()
{
   for (size_t i=0; i<_nb_nodes; i++) {
      try {
         if (_nodes[i]->n()>_nb_nodes)
            THROW_RT("checkNodeLabels(): The node label: " + itos(_nodes[i]->n()) + 
                     " exceeds the total number of nodes.");
      }
      CATCH_EXIT("Mesh");
   }
   if (_nb_nodes>0)
      qksort(_nodes,0,_nb_nodes-1,_node_compare);
}


void Mesh::checkElementLabels()
{
   for (size_t i=0; i<_nb_elements; i++) {
      try {
         if (_elements[i]->n()>_nb_elements)
            THROW_RT("checkElementLabels(): The element label: " + itos(_elements[i]->n()) + 
                     " exceeds the total number of elements.");
      }
      CATCH("Mesh");
   }
   if (_nb_elements>0)
      qksort(_elements,0,_nb_elements-1,_element_compare);
}


void Mesh::checkSideLabels()
{
   for (size_t i=0; i<_nb_sides; i++) {
      try {
         if (_sides[i]->n()>_nb_sides)
            THROW_RT("checkSideLabels(): The side label: " + itos(_sides[i]->n()) + 
                     " exceeds the total number of sides.");
      }
      CATCH("Mesh");
   }
   if (_nb_sides>0)
      qksort(_sides,0,_nb_sides-1,_side_compare);
}


void Mesh::AddNodes(int p)
{
   mesh_elements(*this) {
      try {
         if (the_element->getShape() != TRIANGLE)
            THROW_RT("AddNodes(int): This function is valid for triangles only.");
      }
      CATCH("Mesh");
   }
   getAllSides();

   vector<int> code;
   size_t nb_dof, i;
   size_t mid_label = _nb_nodes + 1;
   Point<real_t> x;

   mesh_sides(*this) {
       for (int k=0; k<p-1; k++) {
          x = ((p-1)*The_side(1)->getCoord() + The_side(2)->getCoord())/real_t(p);
          nb_dof = std::max(The_side(1)->getNbDOF(),The_side(2)->getNbDOF());
          Node *nd = new Node(mid_label++,x);
          nd->setNbDOF(nb_dof);
          for (i=0; i<The_side(1)->getNbDOF(); i++)
             code.push_back(std::min(The_side(1)->getCode(i+1),The_side(2)->getCode(i+1)));
          nd->setCode(code);
          nd->setDOF(_first_dof,nb_dof);
          Add(nd);
          the_side->Add(nd);
       }
   }

   code.clear();
   mesh_elements(*this) {
      x = 0;
      for (size_t j=1; j<=the_element->getNbNodes(); j++)
         x += The_element(j)->getCoord();
      x /= the_element->getNbNodes();
      Node *nd = new Node(mid_label++,x);
      nd->setNbDOF(The_element(1)->getNbDOF());
      for (i=0; i<The_element(1)->getNbDOF(); i++)
         code.push_back(The_element(1)->getCode(i+1));
      nd->setCode(code);
      nd->setDOF(_first_dof,The_element(1)->getNbDOF());
      Add(nd);
      the_element->Add(nd);

      size_t m1=The_element(1)->n(), m2=The_element(2)->n(), m3=The_element(3)->n();
      while (the_element->getNbNodes()!=6) {
         for (size_t j=1; j<=the_element->getNbSides(); j++) {
            Side *sd = the_element->getPtrSide(j);
            size_t n1=(*sd)(1)->n(), n2=(*sd)(2)->n();
            if ((n1==m1 && n2==m2) || (n1==m2 && n2==m1))
               the_element->Add(the_element->getPtrSide(j)->getPtrNode(3),4);
            else if ( (n1==m2 && n2==m3) || (n1==m3 && n2==m2))
               the_element->Add(the_element->getPtrSide(j)->getPtrNode(3),5);
            else if ( (n1==m3 && n2==m1) || (n1==m1 && n2==m3))
               the_element->Add(the_element->getPtrSide(j)->getPtrNode(3),6);
         }
      }
   }
   NumberEquations();
}


void Mesh::Reset()
{
   _nb_elements = 0;
   ElementSet elements(_elements);
   _elements.clear();
   for (std::vector<Element *>::iterator it=elements.begin(); it!=elements.end(); it++) {
      if ((*it)->isActive()) {
         (*it)->setLabel(++_nb_elements);
         _elements.push_back(*it);
      }
   }
   _nb_sides = 0;
   SideSet sides(_sides);
   _sides.clear();
   for (std::vector<Side *>::iterator it=sides.begin(); it!=sides.end(); it++) {
      if ((*it)->isActive()) {
         (*it)->setLabel(++_nb_sides);
         _sides.push_back(*it);
      }
   }
}


void Mesh::RenumberNodes(size_t m)
/*---------------------------------------------------------------------------------
   xadj   : Array containing the adresses of the diagonal entries in the matrix
   adjncy : Array containing the graph of the matrix
  ---------------------------------------------------------------------------------*/
{
   size_t *nn = new size_t [_nb_nodes];
   _available_memory = m;
   if (_verb > 1)
     cout << "Reordering mesh nodes ..." << endl;

   {
     vector<long> xadj(_nb_nodes+1);
     vector<size_t> mask(_nb_nodes+1), xls(_nb_nodes+1), ls(_available_memory), adjncy(_available_memory);
     FindGraph(xadj,adjncy);
     GenRCM(xadj,adjncy,nn,mask,xls);
   }

   size_t n;
   Node *nd;
   NodeSet NL(_nodes);
   for (size_t j=0; j<_nb_nodes; j++) {
      nd = NL[j];
      n = nd->n();
      nd->setLabel(nn[n-1]);
      _nodes[nn[n-1]-1] = nd;
   }
   NumberEquations();
   delete [] nn;
}


void Mesh::RenumberNodes(vector<size_t>& perm,
                         size_t          m)
//---------------------------------------------------------------------------------
// xadj   : Array containing the adresses of the diagonal entries in the matrix
// adjncy : Array containing the graph of the matrix
//---------------------------------------------------------------------------------
{
   size_t *nn = new size_t [_nb_nodes];
   _available_memory = m;
   if (_verb > 1)
      cout << "Reordering mesh nodes ..." << endl;
   {
     vector<long> xadj(_nb_nodes+1);
     vector<size_t> mask(_nb_nodes+1), xls(_nb_nodes+1), ls(_available_memory);
     vector<size_t> adjncy(_available_memory);
     FindGraph(xadj,adjncy);
     GenRCM(xadj,adjncy,nn,mask,xls);
   }

   size_t n;
   Node *nd;
   NodeSet NL(_nodes);
   for (size_t j=0; j<_nb_nodes; j++) {
      nd = NL[j];
      n = nd->n();
      nd->setLabel(nn[n-1]);
      _nodes[nn[n-1]-1] = nd;
   }
   NumberEquations();
   for (size_t i=0; i<_nb_nodes; i++)
      perm[i] = nn[i];
}


void Mesh::get(const string &mesh_file)
{
   size_t i;
   string shape, mat_name, code_string;
   Point<real_t> x;
   if (_verb > 1)
      cout << "Getting mesh data from file: " << mesh_file << " ..." << endl;
   vector<string> bc_code_string;
   bc_code_string.push_back("PERIODIC_A");
   bc_code_string.push_back("PERIODIC_B");
   bc_code_string.push_back("CONTACT");
   bc_code_string.push_back("CONTACT_M");
   bc_code_string.push_back("CONTACT_S");
   bc_code_string.push_back("SLIP");
   vector<string> material_code_string(1);
   material_code_string[0] = "GENERIC";
   _dim = 0;
   XMLParser p(mesh_file);
   p.get(*this);
   checkElementLabels();
   if (_verb > 1)
      cout << "Mesh successfully read." << endl;
   NumberEquations();
   mesh_elements(*this)
      theMaterial.check(the_element->getCode());

// Fill list of marked nodes
   _nb_marked_nodes = 0;
   mesh_nodes(*this) {
      bool mark = false;
      for (i=1; i<=the_node->getNbDOF(); i++) {
         if (the_node->getCode(i)>0)
            mark = true;
      }
      if (mark) {
         _marked_nodes.push_back(the_node);
         _nb_marked_nodes++;
      }
   }
}


void Mesh::get(const string& mesh_file,
                     int     ff, 
                     int     nb_dof)
{
   try {
      if (ff==OFELI_FF)
         get(mesh_file);
      else if (ff==GMSH)
         getGmsh(mesh_file, *this, nb_dof);
      else if (ff==MATLAB)
         getMatlab(mesh_file, *this, nb_dof);
      else if (ff==EASYMESH)
         getEasymesh(mesh_file, *this, nb_dof);
      else if (ff==GAMBIT)
         getGambit(mesh_file, *this, nb_dof);
      else if (ff==BAMG)
         getBamg(mesh_file, *this, nb_dof);
      else if (ff==NETGEN)
         getNetgen(mesh_file, *this, nb_dof);
      else if (ff==TRIANGLE_FF)
         getTriangle(mesh_file, *this, nb_dof);
      else
         THROW_RT("get(string,int,int): Unknown file format "+itos(ff));
   }
   CATCH("Mesh");
}


void Mesh::save(const string& file) const
{
   string type = file.substr(file.rfind("."));
   if (type==".m")
      put(file);
   else if (type==".gpl")
      saveMesh(file,*this,GNUPLOT);
   else if (type==".msh" || type==".geo")
      saveMesh(file,*this,GMSH);
   else if (type==".vtk")
      saveMesh(file,*this,VTK);
}


void Mesh::put(const string& file) const
{
   ofstream fp;
   size_t i, k;
   int sign, m;

   string sh[10] = {"line",
                    "triangle",
                    "quadrilateral",
                    "tetrahedron",
                    "hexahedron",
                    "pentahedron"};
   fp.open(file.c_str(),ios::out);
   if (_verb > 1)
      cout << "Saving mesh data in XML file: " << file << " ..." << endl;
   fp << "<?xml version=\"1.0\"?>\n<OFELI_File>" << endl;
   fp << "<info>" << endl;
   fp << "   <title></title>" << endl;
   fp << "   <date></date>" << endl;
   fp << "   <author></author>" << endl;
   fp << "</info>" << endl;
   size_t n=0;
   if (getDOFSupport()==NODE_DOF)
      n = getPtrNode(1)->getNbDOF();
   else if (getDOFSupport()==ELEMENT_DOF)
      n = 1;
   else if (getDOFSupport()==SIDE_DOF)
      n = getPtrSide(1)->getNbDOF();
   fp << "<Mesh dim=\"" << _dim << "\" nb_dof=\"" << n << "\">" << endl;
   if (_nb_nodes>0) {
      k = 0;
      fp << "   <Nodes>" << endl;
      mesh_nodes(*this) {
         fp.setf(ios::right|ios::scientific);
         for (i=1; i<=_dim; i++)
            fp << "  " << setprecision(8) << setw(18) << the_node->getCoord(i);
         sign = 1;
         if (the_node->getCode(1)<0)
            sign = -1;
         m = 0;
         for (size_t j=1; j<=n; j++)
            m += abs(the_node->getCode(j))*size_t(pow(10.,real_t(n-j)));
         m *= sign;
         fp << setw(6) << m << "   ";
         k++;
         if (k%5==0)
            fp << endl;
      }
      if (k%5!=0)
         fp << endl;
      fp << "   </Nodes>" << endl;
   }

   if (_nb_elements>0) {
      k = 0;
      size_t nbn = getPtrElement(1)->getNbNodes();
      string shape = sh[getPtrElement(1)->getShape()];
      fp << "   <Elements shape=\"" << shape << "\"  nodes=\"" << nbn << "\">" << endl;
      mesh_elements(*this) {
         for (i=1; i<=nbn; i++)
            fp << setw(8) << The_element(i)->n();
         fp << setw(8) << the_element->getCode() << "   ";
         k++;
         if (k%5==0)
            fp << endl;
      }
      if (k%5!=0)
         fp << endl;
      fp << "   </Elements>" << endl;
   }

   if (_nb_sides>0) {
      k = 0;
      size_t nbn = getPtrSide(1)->getNbNodes();
      string shape = sh[getPtrSide(1)->getShape()];
      fp << "   <Sides shape=\"" << shape << "\"  nodes=\"" << nbn << "\">" << endl;
      mesh_sides(*this) {
         for (i=1; i<=nbn; i++)
            fp << setw(8) << the_side->getNodeLabel(i);
         m = 0;
         for (size_t j=1; j<=n; j++)
            m += the_side->getCode(j)*size_t(pow(10.,real_t(n-j)));
         fp << setw(6) << m << "   ";
         k++;
         if (k%5==0)
            fp << endl;
      }
      if (k%5!=0)
         fp << endl;
      fp << "   </Sides>" << endl;
   }
   if (_nb_mat>1 || theMaterial.getName(1)!="Generic") {
      fp << "   <Material>" << endl;
      for (size_t i=1; i<=_nb_mat; i++)
         fp << setw(9) << theMaterial.getCode(i) << "   " << theMaterial.getName(theMaterial.getCode(i)) << endl;
      fp << "   </Material>" << endl;
   }
   fp << "</Mesh>" << endl;
   fp << "</OFELI_File>" << endl;
   fp.close();
}


//=============================================================================
//                   MEMBER FUNCTIONS FOR INTERNAL USE
//=============================================================================


unsigned long Mesh::FindGraph(vector<long>&   xadj,
                              vector<size_t>& adjncy)
{
    long unsigned l, k, is, isn, m, memory;
    size_t i, j, ii, jj, in;
    int jtest;

    if (_verb > 1)
       cout << "Getting mesh graph ..." << endl;
    for (i=1; i<=_nb_nodes; ++i)
       xadj[i] = 0;
    xadj[0] = 1;
    for (i=0; i<_available_memory; ++i)
       adjncy[i] = 0;
    memory = 0;

    mesh_elements(*this) {
       for (i=1; i<=the_element->getNbNodes(); i++) {
          ii = the_element->getNodeLabel(i);
          for (j=1; j<=the_element->getNbNodes(); j++) {
             jj = The_element(j)->n();
             if (ii != jj) {
                k = 0;
                for (in=1; in<=ii; ++in)
                   k += xadj[in-1];
                m = k + xadj[ii];
                if (adjncy[k-1] == 0) {
                  adjncy[k-1] = jj;
                  ++memory;
                  ++xadj[ii];
                  goto L9;
                }
                for (l=k; l<=m-1; ++l) {
                   jtest = adjncy[l-1] - jj;
                   if (jtest < 0)
                      goto L5;
                   else if (jtest == 0)
                      goto L9;
                   else
                     goto L7;
L5:
                  ;
                }
                ++memory;
                try {
                   if (memory > _available_memory && _available_memory > 0)
                      THROW_RT("FindGraph(vector<long>,vector<size_t>): Insufficient Memory");
                }
                CATCH_EXIT("Mesh");
                ++xadj[ii];
                for (is=m; is<=memory; ++is) {
                   isn = m + memory - is;
                   adjncy[isn] = adjncy[isn-1];
                }
                adjncy[m-1] = jj;
                goto L9;
L7:
                ++memory;
                try {
                   if (memory > _available_memory && _available_memory > 0)
                      THROW_RT("FindGraph(vector<long>,vector<size_t>): Insufficient Memory");
                }
                CATCH_EXIT("Mesh");
                ++xadj[ii];
                for (is=l; is<=memory; ++is) {
                   isn = memory + l - is;
                   adjncy[isn] = adjncy[isn-1];
                }
                adjncy[l-1] = jj;
L9:
                ;
             }
          }
       }
    }

    for (i=1; i<=_nb_nodes; ++i)
       xadj[i] += xadj[i-1];
    return memory;
}


void Mesh::GenRCM(vector<long>&   xadj,
                  vector<size_t>& adjncy,
                  size_t*         perm,
                  vector<size_t>& mask,
                  vector<size_t>& xls)
/*-----------------------------------------------------------------------------
     Find the reverse Cuthill-McKee ordering for a general graph.
     For each connected component in the graph, GenRCM obtains the
     ordering by calling the function RCM.

     From : A. George and J.W. Liu
            Computer Solution of Large Sparse Positive Definite Systems
            Prentice-Hall Series in Computational Mathematics
            Prentice-Hall, Englewood Cliffs, 1981.
            (pages 66-75).
 ------------------------------------------------------------------------------*/
{
   size_t i, root, ccsize=0;
   long nlvl;
   for (i=0; i<_nb_nodes; ++i)
      mask[i] = 1;
   size_t num = 1;
   for (i=0; i<_nb_nodes; ++i) {
      if (mask[i] != 0) {
         root = i+1;
         FindRoot(root, xadj, adjncy, mask, nlvl, xls, perm+num-1);
         RCM(root, xadj, adjncy, mask, perm+num-1, ccsize, xls);
         num += ccsize;
         if (num > _nb_nodes)
            return;
      }
   }
}


Mesh & Mesh::operator=(Mesh& ms)
{
   _is_structured = ms._is_structured;
   _dim = ms._dim;
   _verb = ms._verb;
   _first_dof = 1;
   _nb_dof = 0;
   _nodes = ms._nodes;
   _elements = ms._elements;
   _sides = ms._sides;
   _nb_nodes = _nb_elements = _nb_sides = _nb_vertices = 0;
   _max_nb_nodes = ms._max_nb_nodes;
   _max_nb_elements = ms._max_nb_elements;
   _max_nb_sides = ms._max_nb_sides;
   _all_sides_created = ms._all_sides_created;
   _boundary_sides_created = ms._boundary_sides_created;
   _all_edges_created = ms._all_edges_created;
   _boundary_edges_created = ms._boundary_edges_created;
   _n_view1 = ms._n_view1; _n_view2 = ms._n_view2;
   _e_view1 = ms._e_view1; _e_view2 = ms._e_view2;
   _s_view1 = ms._s_view1; _s_view2 = ms._s_view2;
   _ed_view1 = ms._ed_view1; _ed_view2 = ms._ed_view2;

// Insert nodes
   mesh_nodes(ms) {
      Add(new Node(The_node));
      _nb_eq = _nb_dof;
   }

// Insert elements
   mesh_elements(ms)
      Add(new Element(The_element));

// Insert sides
   mesh_sides(ms)
      Add(new Side(side_label,the_side->getShape()));

   if (ms._node_in_coarse_element.size()>0) {
      for (size_t i=0; i<_nb_nodes; i++)
         _node_in_coarse_element.push_back(ms._node_in_coarse_element[i]);
   }
   if (ms._node_in_fine_element.size()>0) {
      for (size_t i=0; i<_nb_nodes; i++)
         _node_in_fine_element.push_back(ms._node_in_fine_element[i]);
   }

   _nb_mat = ms._nb_mat;
   for (size_t k=0; k<_nb_mat; k++) {
      _code_mat[k] = ms._code_mat[k];
      _mat[k] = ms._mat[k];
   }
   _set_nodes = ms._set_nodes;
   _set_sides = ms._set_sides;
   _set_elements = ms._set_elements;
   _no_imposed_dof = ms._no_imposed_dof;
   return *this;
}


ostream& operator<<(      ostream& s,
                    const Mesh&    ms)
{
   size_t n1=ms._n_view1, n2=ms._n_view2, e1=ms._e_view1, e2=ms._e_view2;
   size_t s1=ms._s_view1, s2=ms._s_view2, ed1=ms._ed_view1, ed2=ms._ed_view2;

   if (n1==0)
      n1 = 1;
   if (n2==0)
      n2 = ms.getNbNodes();
   if (e1==0)
      e1 = 1;
   if (e2==0)
      e2 = ms.getNbElements();
   if (s1==0)
      s1 = 1;
   if (s2==0)
      s2 = ms.getNbSides();
   if (ed1==0)
      ed1 = 1;
   if (ed2==0)
      ed2 = ms.getNbEdges();

   s << "\n\nM E S H     D A T A\n===================\n\n";
   s << "Space Dimension        : " << setw(6) << ms.getDim() << endl;
   s << "Number of nodes        : " << setw(6) << ms.getNbNodes() << endl;
   s << "Number of elements     : " << setw(6) << ms.getNbElements() << endl;
   if (ms.getNbSides()>0)
      s << "Number of sides        : " << setw(6) << ms.getNbSides() << endl;
   if (ms.getNbEdges()>0)
      s << "Number of edges        : " << setw(6) << ms.getNbEdges() << endl;

   size_t i;
   if (ms.getVerbose() > 1) {
      s << "\n\nLIST OF ELEMENTS:\n";
      for (i=e1; i<=e2; i++)
         s << *(ms.getPtrElement(i));
      s << endl;
   }

   if (ms.getVerbose() > 1) {
      s << "\n\nLIST OF NODES:\n";
      for (i=n1; i<=n2; i++)
         s << *(ms.getPtrNode(i));
      s << endl;
   }

   if (ms.getVerbose() > 2) {
      if (ms.getNbSides() > 0) {
         s << "\n\nLIST OF SIDES:\n";
         for (i=s1; i<=s2; i++)
            s << *(ms.getPtrSide(i));
         s << endl;
      }
   }
   s << endl;

   if (ms.getVerbose() > 3) {
      if (ms.getNbEdges() > 0) {
         s << "\n\nLIST OF EDGES:\n";
         for (i=ed1; i<=ed2; i++)
            s << *(ms.getPtrEdge(i));
         s << endl;
      }
   }
   s << endl;
   return s;
}

} /* namespace OFELI */
