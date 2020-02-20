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

                         Implementation of class 'Mesh'

  ==============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <string>

#include <iostream>
using std::cout;
using std::ostream;
using std::endl;

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
#include "linear_algebra/Vect_impl.h"
#include <algorithm>
#include "OFELIException.h"

namespace OFELI {

Node *theNode, *the_node;
Element *theElement, *the_element;
Side *theSide, *the_side;
Edge *theEdge, *the_edge;

extern Material theMaterial;
Mesh::Mesh()
     : _nb_nodes(0), _nb_elements(0), _nb_sides(0), _nb_boundary_sides(0), _nb_edges(0), _dim(2),
       _nb_dof(0), _nb_vertices(0), _first_dof(1), _nb_mat(0), _no_imposed_dof(false),
       _is_structured(false), _all_sides_created(false), _boundary_sides_created(false),
       _all_edges_created(false), _boundary_edges_created(false), _boundary_nodes_created(false),
       _node_neighbor_elements_created(false), _element_neighbor_elements_created(false)
{
   setDOFSupport(NODE_DOF);
}


Mesh::Mesh(const string& file,
           bool          bc,
           int           opt,
           int           nb_dof)
     : _nb_nodes(0), _nb_elements(0), _nb_sides(0), _nb_boundary_sides(0), _nb_edges(0), _dim(2),
       _nb_dof(0), _nb_vertices(0), _first_dof(1), _nb_mat(0), _no_imposed_dof(false),
       _is_structured(false), _all_sides_created(false), _boundary_sides_created(false),
       _all_edges_created(false), _boundary_edges_created(false), _boundary_nodes_created(false),
       _node_neighbor_elements_created(false), _element_neighbor_elements_created(false)
{
   string fs = file.substr(file.rfind(".")+1);
   int ff = 0;
   if (fs=="msh")
      ff = GMSH;
   else if (fs=="gpl")
      ff = GNUPLOT;
   else
      ff = OFELI_FF;
   _no_imposed_dof = bc;
   _set_nodes = true;
   _set_sides = _set_edges = _set_elements = false;
   get(file,ff,nb_dof);
   string mat = theMaterial.getName(1);
}


Mesh::Mesh(real_t xmax,
           size_t nb_el,
           size_t p,
           size_t nb_dof)
     : _nb_nodes(0), _nb_boundary_nodes(2), _nb_elements(0), _nb_sides(0), _nb_boundary_sides(0),
       _nb_edges(0), _dim(1), _nb_dof(0), _nb_vertices(0), _first_dof(1), _nb_mat(0),
       _no_imposed_dof(false), _is_structured(false), _all_sides_created(false),
       _boundary_sides_created(false), _all_edges_created(false), _boundary_edges_created(false),
       _boundary_nodes_created(false), _node_neighbor_elements_created(false),
       _element_neighbor_elements_created(false)
{
   set1D(0.,xmax,nb_el,p,nb_dof);
}


Mesh::Mesh(real_t xmin,
           real_t xmax,
           size_t nb_el,
           size_t p,
           size_t nb_dof)
     : _nb_nodes(0), _nb_boundary_nodes(2), _nb_elements(0), _nb_sides(0), _nb_boundary_sides(0),
       _nb_edges(0), _dim(1), _nb_dof(0), _nb_vertices(0), _first_dof(1), _nb_mat(0),
       _no_imposed_dof(false), _is_structured(false), _all_sides_created(false),
       _boundary_sides_created(false), _all_edges_created(false), _boundary_edges_created(false),
       _boundary_nodes_created(false), _node_neighbor_elements_created(false),
       _element_neighbor_elements_created(false)
{
   set1D(xmin,xmax,nb_el,p,nb_dof);
}


Mesh::Mesh(const Grid& g,
           int         opt)
     : _nb_nodes(0), _nb_elements(0), _nb_sides(0), _nb_boundary_sides(0), _nb_edges(0), _dim(3),
       _nb_dof(0), _nb_vertices(0), _first_dof(1), _nb_mat(1), _no_imposed_dof(false),
       _is_structured(true), _all_sides_created(false), _boundary_sides_created(false),
       _all_edges_created(false), _boundary_edges_created(false), _boundary_nodes_created(false),
       _node_neighbor_elements_created(false), _element_neighbor_elements_created(false)
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
   n = ne = 0, m = 1;
   switch (_dim) {

      case 1:
         for (i=1; i<=nx; i++) {
            el = new Element(++ne,LINE,1);
            el->Add(theNodes[i-1]);
            el->Add(theNodes[i]);
            Add(el);
         }
         break;

      case 2:
         if (opt==QUADRILATERAL) {
            for (j=1; j<=ny; j++) {
               for (i=1; i<=nx; i++) {
                  if (g.isActive(i,j)) {
                     el = new Element(++ne,QUADRILATERAL,1);
                     el->Add(theNodes[nn(i  ,j  )-1]);
                     el->Add(theNodes[nn(i+1,j  )-1]);
                     el->Add(theNodes[nn(i+1,j+1)-1]);
                     el->Add(theNodes[nn(i  ,j+1)-1]);
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
                     el->Add(theNodes[nn(i  ,j  )-1]);
                     el->Add(theNodes[nn(i+1,j  )-1]);
                     el->Add(theNodes[nn(i+1,j+1)-1]);
                     Add(el);
                     el = new Element(++ne,TRIANGLE,1);
                     el->Add(theNodes[nn(i+1,j+1)-1]);
                     el->Add(theNodes[nn(i  ,j+1)-1]);
                     el->Add(theNodes[nn(i  ,j  )-1]);
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
                        the_element = new Element(++ne,HEXAHEDRON,1);
                        The_element.Add(theNodes[nn(i  ,j  ,k  )-1]);
                        The_element.Add(theNodes[nn(i+1,j  ,k  )-1]);
                        The_element.Add(theNodes[nn(i+1,j+1,k  )-1]);
                        The_element.Add(theNodes[nn(i  ,j+1,k  )-1]);
                        The_element.Add(theNodes[nn(i  ,j  ,k+1)-1]);
                        The_element.Add(theNodes[nn(i  ,j+1,k+1)-1]);
                        The_element.Add(theNodes[nn(i+1,j+1,k+1)-1]);
                        The_element.Add(theNodes[nn(i  ,j+1,k+1)-1]);
                        Add(the_element);
                     }
                  }
               }
            }
         }
         break;
   }
   NumberEquations();
   element_loop(this)
      theMaterial.check(The_element.getCode());
}


Mesh::Mesh(const Grid& g,
           int         shape,
           int         opt)
     : _nb_nodes(0), _nb_elements(0), _nb_sides(0), _nb_boundary_sides(0), _nb_edges(0),
       _dim(3), _nb_dof(0), _nb_vertices(0), _first_dof(1), _nb_mat(1), _no_imposed_dof(false),
       _is_structured(true), _all_sides_created(false), _boundary_sides_created(false),
       _all_edges_created(false), _boundary_edges_created(false), _boundary_nodes_created(false),
       _node_neighbor_elements_created(false), _element_neighbor_elements_created(false)
{
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
            the_node = new Node(++n,x);
            _code[0] = g.getCode(i);
            The_node.setNbDOF(1);
            The_node.setDOF(_first_dof,1);
            The_node.setCode(_code);
            Add(the_node);
         }
         break;

      case 2:
         for (i=1; i<=nx; i++) {
            for (j=1; j<=ny; j++) {
              Point<real_t> x = 0.5*(g.getCoord(i,j)+g.getCoord(i+1,j+1));
               the_node = new Node(++n,x);
               _code[0] = g.getCode(i,j);
               The_node.setNbDOF(1);
               The_node.setDOF(_first_dof,1);
               The_node.setCode(_code);
               Add(the_node);
            }
         }
         break;

      case 3:
         for (i=1; i<=nx; i++) {
            for (j=1; j<=ny; j++) {
               for (k=1; k<=nz; k++) {
                 Point<real_t> x = 0.5*(g.getCoord(i,j,k)+g.getCoord(i+1,j+1,k+1));
                  the_node = new Node(++n,x);
                  _code[0] = g.getCode(i,j,k);
                  The_node.setNbDOF(1);
                  The_node.setDOF(_first_dof,1);
                  The_node.setCode(_code);
                  Add(the_node);
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
            the_element = new Element(++ne,LINE,1);
            The_element.Add(theNodes[i-1]);
            The_element.Add(theNodes[i]);
            Add(the_element);
         }
         break;

      case 2:
         if (shape==QUADRILATERAL) {
            for (j=1; j<=ny-1; j++) {
               for (i=1; i<=nx-1; i++) {
                  n1 = nn + i - 1;
                  the_element = new Element(++ne,QUADRILATERAL,1);
                  The_element.Add(theNodes[n1-1]);
                  The_element.Add(theNodes[n1]);
                  The_element.Add(theNodes[n1+nx+1]);
                  The_element.Add(theNodes[n1+nx]);
                  Add(the_element);
               }
               nn += nx;
            }
         }
         else if (shape==TRIANGLE) {
            for (i=1; i<=nx-1; i++) {
               for (j=1; j<=ny-1; j++) {
                  n1 = nn + i - 1;
                  the_element = new Element(++ne,TRIANGLE,1);
                  The_element.Add(theNodes[n1-1]);
                  The_element.Add(theNodes[n1+ny-1]);
                  The_element.Add(theNodes[n1]);
                  Add(the_element);
                  the_element = new Element(++ne,TRIANGLE,1);
                  The_element.Add(theNodes[n1+ny-1]);
                  The_element.Add(theNodes[n1]);
                  The_element.Add(theNodes[n1-1]);
                  Add(the_element);
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
                     the_element = new Element(++ne,HEXAHEDRON,1);
                     The_element.Add(theNodes[n1-1]);
                     The_element.Add(theNodes[n1]);
                     The_element.Add(theNodes[n1+nx+1]);
                     The_element.Add(theNodes[n1+nx]);
                     The_element.Add(theNodes[n2-1]);
                     The_element.Add(theNodes[n2]);
                     The_element.Add(theNodes[n2+nx+1]);
                     The_element.Add(theNodes[n2+nx]);
                     Add(the_element);
                  }
                  nn += nx + 1;
               }
               nn += ny + 1;
            }
         }
         break;
   }
   NumberEquations();
   element_loop(this)
      theMaterial.check(The_element.getCode());
}


Mesh::Mesh(real_t xmin,
           real_t xmax,
           size_t ne,
           int    c1,
           int    c2,
           int    p,
           size_t nb_dof)
     : _nb_nodes(0), _nb_elements(0), _nb_sides(0), _nb_boundary_sides(0), _nb_edges(0), _dim(1),
       _nb_dof(0), _nb_vertices(0), _first_dof(1), _nb_mat(1), _no_imposed_dof(false),
       _is_structured(true), _all_sides_created(false), _boundary_sides_created(false),
       _all_edges_created(false), _boundary_edges_created(false), _boundary_nodes_created(false),
       _node_neighbor_elements_created(false), _element_neighbor_elements_created(false)
{
   theMaterial.set(1,"Generic");
   setNodesForDOF();
   real_t x=xmin, hp=(xmax-xmin)/(ne*p);

   vector<Node *> nd;
   for (size_t n=1; n<=ne*p+1; n++) {
      the_node = new Node(n,Point<real_t>(x));
      The_node.setNbDOF(nb_dof);
      The_node.setDOF(_first_dof,nb_dof);
      The_node.setCode(1,0);
      if (n==1)
         The_node.setCode(1,c1);
      if (n==ne*p+1)
         The_node.setCode(1,c2);
      Add(the_node);
      nd.push_back(the_node);
      x += hp;
   }

   size_t n=0;
   for (size_t e=1; e<=ne; e++) {
      the_element = new Element(e,LINE);
      The_element.Add(nd[n]);
      n += p;
      The_element.Add(nd[n]);
      for (int m=1; m<p; m++)
         The_element.Add(nd[n-p+m]);
      Add(the_element);
   }

   the_side = new Side(1,POINT);
   The_side.Add(theNodes[0]);
   The_side.setCode(1,0);
   if (nd[0]->getCode(1)==0)
      The_side.setCode(1,1);
   Add(the_side);
   the_side = new Side(2,POINT);
   The_side.Add(nd[_nb_nodes-1]);
   The_side.setCode(1,0);
   if (nd[_nb_nodes-1]->getCode(1)==0)
      The_side.setCode(1,2);
   Add(the_side);

   NumberEquations();
   element_loop(this)
      theMaterial.check(The_element.getCode());
}


Mesh::Mesh(real_t xmin,
           real_t xmax,
           real_t ymin,
           real_t ymax,
           size_t nx,
           size_t ny,
           int    cx0,
           int    cxN,
           int    cy0,
           int    cyN,
           int    opt,
           size_t nb_dof)
     : _nb_nodes(0), _nb_elements(0), _nb_sides(0), _nb_boundary_sides(0), _nb_edges(0), _dim(2),
       _nb_dof(0), _nb_vertices(0), _first_dof(1), _nb_mat(1), _no_imposed_dof(false), _is_structured(true),
       _all_sides_created(false), _boundary_sides_created(false), _all_edges_created(false),
       _boundary_edges_created(false), _boundary_nodes_created(false),
       _node_neighbor_elements_created(false), _element_neighbor_elements_created(false)
{
   theMaterial.set(1,"Generic");
   setNodesForDOF();
   size_t n=0, ne=0;
   Point<real_t> x(xmin,ymin);
   real_t hx=(xmax-xmin)/nx, hy=(ymax-ymin)/ny;

   for (size_t i=0; i<=nx; i++) {
      x.y = ymin;
      for (size_t j=0; j<=ny; j++) {
         n++;
         the_node = new Node(n,x);
         The_node.setNbDOF(nb_dof);
         _code[0] = 0;
         if (j==0  && cy0>0)
            _code[0] = cy0;
         if (j==ny && cyN>0)
            _code[0] = cyN;
         if (i==0  && cx0>0)
            _code[0] = cx0;
         if (i==nx && cxN>0)
            _code[0] = cxN;
         The_node.setDOF(_first_dof,nb_dof);
         The_node.setCode(_code);
         Add(the_node);
         x.y += hy;
         if (i>0 && j>0 && opt) {
            if (opt==QUADRILATERAL) {
               ne++;
               the_element = new Element(ne,QUADRILATERAL,1);
               The_element.Add(theNodes[n-ny-3]);
               The_element.Add(theNodes[n-2]);
               The_element.Add(theNodes[n-1]);
               The_element.Add(theNodes[n-ny-2]);
               Add(the_element);
            }
            else if (opt==TRIANGLE) {
               ne++;
               the_element = new Element(ne,TRIANGLE,1);
               The_element.Add(theNodes[n-ny-3]);
               The_element.Add(theNodes[n-2]);
               The_element.Add(theNodes[n-1]);
               Add(the_element);
               ne++;
               the_element = new Element(ne,TRIANGLE,1);
               The_element.Add(theNodes[n-1]);
               The_element.Add(theNodes[n-ny-2]);
               The_element.Add(theNodes[n-ny-3]);
               Add(the_element);
            }
            else
               throw OFELIException("Mesh::Mesh(real_t,real_t,real_t,real_t,size_t,size_t,int,int,int,int,int): "
                                    "Illegal option "+itos(opt));
         }
      }
      x.x += hx;
   }

   size_t is=1;
   n =  1;
   if (cx0<0) {
      for (size_t j=1; j<=ny; j++) {
         the_side = new Side(is++,LINE);
         The_side.Add(getPtrNode(n));
         The_side.Add(getPtrNode(n+1));
         n++;
         The_side.setNbDOF(nb_dof);
         The_side.setCode(1,-cx0);
         Add(the_side);
      }
   }
   n = nx*(ny+1) + 1;
   if (cxN<0) {
      for (size_t j=1; j<=ny; j++) {
         the_side = new Side(is++,LINE);
         The_side.Add(getPtrNode(n));
         The_side.Add(getPtrNode(n+1));
         n++;
         The_side.setNbDOF(nb_dof);
         The_side.setCode(1,-cxN);
         Add(the_side);
      }
   }
   n = 1;
   if (cy0<0) {
      for (size_t i=1; i<=nx; i++) {
         the_side = new Side(is++,LINE);
         The_side.Add(getPtrNode(n));
         The_side.Add(getPtrNode(n+ny+1));
         n += ny + 1;
         The_side.setNbDOF(nb_dof);
         The_side.setCode(1,-cy0);
         Add(the_side);
      }
   }
   n = ny + 1;
   if (cyN<0) {
      for (size_t i=1; i<=nx; i++) {
         the_side = new Side(is++,LINE);
         The_side.Add(getPtrNode(n));
         The_side.Add(getPtrNode(n+ny+1));
         n += ny + 1;
         The_side.setNbDOF(nb_dof);
         The_side.setCode(1,-cyN);
         Add(the_side);
      }
   }
   NumberEquations();
   element_loop(this)
      theMaterial.check(The_element.getCode());
}


Mesh::Mesh(real_t xmin,
           real_t xmax,
           real_t ymin,
           real_t ymax,
           real_t zmin,
           real_t zmax,
           size_t nx,
           size_t ny,
           size_t nz,
           int    cx0,
           int    cxN,
           int    cy0,
           int    cyN,
           int    cz0,
           int    czN,
           int    opt,
           size_t nb_dof)
     : _nb_nodes(0), _nb_elements(0), _nb_sides(0), _nb_boundary_sides(0), _nb_edges(0), _dim(3),
       _nb_dof(0), _nb_vertices(0), _first_dof(1), _nb_mat(1), _no_imposed_dof(false),
       _is_structured(true), _all_sides_created(false), _boundary_sides_created(false),
       _all_edges_created(false), _boundary_edges_created(false), _boundary_nodes_created(false),
       _node_neighbor_elements_created(false), _element_neighbor_elements_created(false)
{
   theMaterial.set(1,"Generic");
   setNodesForDOF();
   Point<real_t> x(xmin,ymin,zmin);
   real_t hx=(xmax-xmin)/nx, hy=(ymax-ymin)/ny, hz=(zmax-zmin)/nz;

   size_t n=0;
   for (size_t k=0; k<=nz; k++) {
      x.y = ymin;
      for (size_t j=0; j<=ny; j++) {
         x.x = xmin;
         for (size_t i=0; i<=nx; i++) {
            the_node = new Node(++n,x);
            The_node.setNbDOF(nb_dof);
            _code[0] = 0;
            if (j==0  && cy0>0) _code[0] = cy0;
            if (j==ny && cyN>0) _code[0] = cyN;
            if (i==0  && cx0>0) _code[0] = cx0;
            if (i==nx && cxN>0) _code[0] = cxN;
            if (k==0  && cz0>0) _code[0] = cz0;
            if (k==nz && czN>0) _code[0] = czN;
            The_node.setDOF(_first_dof,nb_dof);
            The_node.setCode(_code);
            Add(the_node);
            x.x += hx;
         }
         x.y += hy;
      }
      x.z += hz;
   }

   n = 0;
   size_t ne=0, nn=(nx+1)*(ny+1);
   for (size_t k=1; k<=nz; k++) {
      for (size_t j=1; j<=ny; j++) {
         for (size_t i=1; i<=nx; i++) {
            Node *nd[8] = { theNodes[n], theNodes[n+1], theNodes[n+nx+2], theNodes[n+nx+1],
                            theNodes[n+nn], theNodes[n+nn+1], theNodes[n+nn+nx+2],
                            theNodes[n+nn+nx+1] };
            if (opt==HEXAHEDRON) {
               the_element = new Element(++ne,HEXAHEDRON,1);
               for (size_t i=0; i<8; i++)
                  The_element.Add(nd[i]);
               n++;
               Add(the_element);
            }
            else if (opt==TETRAHEDRON) {
               the_element = new Element(++ne,TETRAHEDRON,1);
               The_element.Add(nd[0]); The_element.Add(nd[1]);
               The_element.Add(nd[3]); The_element.Add(nd[7]);
               Add(the_element);
               the_element = new Element(++ne,TETRAHEDRON,1);
               The_element.Add(nd[0]); The_element.Add(nd[1]);
               The_element.Add(nd[4]); The_element.Add(nd[7]);
               Add(the_element);
               the_element = new Element(++ne,TETRAHEDRON,1);
               The_element.Add(nd[1]); The_element.Add(nd[2]);
               The_element.Add(nd[3]); The_element.Add(nd[7]);
               Add(the_element);
               the_element = new Element(++ne,TETRAHEDRON,1);
               The_element.Add(nd[1]); The_element.Add(nd[2]);
               The_element.Add(nd[6]); The_element.Add(nd[7]);
               Add(the_element);
               the_element = new Element(++ne,TETRAHEDRON,1);
               The_element.Add(nd[0]); The_element.Add(nd[4]);
               The_element.Add(nd[5]); The_element.Add(nd[7]);
               Add(the_element);
               the_element = new Element(++ne,TETRAHEDRON,1);
               The_element.Add(nd[1]); The_element.Add(nd[5]);
               The_element.Add(nd[6]); The_element.Add(nd[7]);
               Add(the_element);
            }
            else
               throw OFELIException("Mesh::Mesh(real_t,real_t,real_t,real_t,size_t,size_t,int,int,int,int,int):"
                                    " Illegal option "+itos(opt));
         }
         n++;
      }
      n += nx + 1;
   }

   size_t is=1;
   n = 1;
   for (size_t j=1; j<=ny; j++) {
      the_side = new Side(is++,LINE);
      The_side.Add(getPtrNode(n));
      The_side.Add(getPtrNode(n+1));
      n++;
      the_side->setNbDOF(nb_dof);
      if (cy0<0) {
         The_side.setCode(1,-cy0);
         Add(the_side);
      }
   }
   n = nx*(ny+1) + 1;
   for (size_t j=1; j<=ny; j++) {
      the_side = new Side(is++,LINE);
      The_side.Add(getPtrNode(n));
      The_side.Add(getPtrNode(n+1));
      n++;
      The_side.setNbDOF(nb_dof);
      if (cyN<0) {
         The_side.setCode(1,-cyN);
         Add(the_side);
      }
   }
   n = 1;
   for (size_t i=1; i<=nx; i++) {
      the_side = new Side(is++,LINE);
      The_side.Add(getPtrNode(n));
      The_side.Add(getPtrNode(n+ny+1));
      n += ny + 1;
      The_side.setNbDOF(nb_dof);
      if (cx0<0) {
         The_side.setCode(1,-cx0);
         Add(the_side);
      }
   }
   n = ny + 1;
   for (size_t i=1; i<=nx; i++) {
      the_side = new Side(is++,LINE);
      The_side.Add(getPtrNode(n));
      The_side.Add(getPtrNode(n+ny+1));
      n += ny + 1;
      The_side.setNbDOF(nb_dof);
      if (cxN<0) {
         The_side.setCode(1,-cxN);
         Add(the_side);
      }
   }

   NumberEquations();
   element_loop(this)
      theMaterial.check(The_element.getCode());
}


Mesh::Mesh(const Mesh& ms)
     : _nb_nodes(0), _nb_elements(0), _nb_sides(0), _nb_boundary_sides(0),
       _nb_edges(0), _nb_side_nodes(ms._nb_side_nodes),
       _nb_element_nodes(ms._nb_element_nodes), _dim(ms._dim), _nb_dof(0),
       _nb_vertices(ms._nb_vertices), _first_dof(ms._first_dof), _nb_eq(ms._nb_eq),
       _nb_mat(ms._nb_mat), _max_nb_nodes(ms._max_nb_nodes), _max_nb_elements(ms._max_nb_elements),
       _max_nb_sides(ms._max_nb_sides), _max_nb_edges(ms._max_nb_edges),
       _no_imposed_dof(ms._no_imposed_dof), _is_structured(ms._is_structured),
       _all_sides_created(ms._all_sides_created), _boundary_sides_created(ms._boundary_sides_created),
       _all_edges_created(ms._all_edges_created), _boundary_edges_created(ms._boundary_edges_created),
       _boundary_nodes_created(ms._boundary_nodes_created),
       _node_neighbor_elements_created(ms._node_neighbor_elements_created),
       _element_neighbor_elements_created(false)
{
// Insert nodes
   node_loop(&ms)
      Add(new Node(The_node));
   _nb_eq = _nb_dof;

// Insert elements
   for (size_t n=1; n<=ms.getNbElements(); ++n) {
      Element *el = ms(n);
      the_element = new Element(el->n(),el->getShape(),el->getCode());
      for (size_t i=1; i<=el->getNbNodes(); ++i)
         The_element.Add(getPtrNode((*el)(i)->n()));
      Add(the_element);
   }

// Insert sides
   for (size_t n=1; n<=ms.getNbSides(); ++n) {
      Side *sd = ms.getPtrSide(n);
      the_side = new Side(sd->n(),sd->getShape());
      for (size_t i=1; i<=sd->getNbNodes(); ++i)
         The_side.Add(getPtrNode((*sd)(i)->n()));
      The_side.setNbDOF(sd->getNbDOF());
      for (size_t i=1; i<=sd->getNbDOF(); ++i)
         The_side.setCode(i,sd->getCode(i));
      Add(the_side);
   }
   
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
     : _nb_nodes(0), _nb_elements(0), _nb_sides(0), _nb_boundary_sides(0), _nb_edges(0), _dim(2), _nb_dof(0),
       _nb_vertices(0), _first_dof(1), _nb_mat(1), _no_imposed_dof(false),
       _is_structured(false), _all_sides_created(false), _boundary_sides_created(false),
       _all_edges_created(false), _boundary_edges_created(false), _boundary_nodes_created(false),
       _node_neighbor_elements_created(false), _element_neighbor_elements_created(false)
{
   if (m.getDim()!=2)
      throw OFELIException("Mesh::Mesh(Mesh,Point<real_t>,Point<real_t>)\n"
                           "This constructor is valid for 2-D meshes only.");
   _max_nb_nodes = m.getNbNodes();
   _max_nb_elements = m.getNbElements();
   _max_nb_sides = m.getNbSides();
   setNodesForDOF();

   _node_old_label.resize(m.getNbNodes());
   _node_new_label.resize(m.getNbNodes());
   size_t label = 1;
   node_loop(&m) {
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
   element_loop(&m) {
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


Mesh::Mesh(const Mesh& mesh,
           int         opt,
           size_t      dof1,
           size_t      dof2,
           bool        bc)
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
   else if (opt==BOUNDARY_SIDE_DOF)
      setBoundarySidesForDOF();

   element_loop(&mesh)
      Add(new Element(The_element));

   if (opt==NODE_DOF) {
      node_loop(&mesh) {
         nd = new Node(The_node);
         nd->setNbDOF(nb_dof);
         for (size_t i=0; i<nb_dof; i++)
            nd->setCode(i+1,the_node->getCode(dof1+i));
         nd->setDOF(_first_dof,nb_dof);
         Add(nd);
      }
   }
   else {
      node_loop(&mesh)
         Add(new Node(The_node));
   }

   if (opt==SIDE_DOF) {
      side_loop(&mesh) {
         sd = new Side(The_side);
         for (size_t i=0; i<nb_dof; i++)
            sd->setCode(i+1,the_side->getCode(dof1+i));
         sd->setDOF(_first_dof,nb_dof);
         Add(sd);
      }
   }
   else {
      side_loop(&mesh)
         Add(new Side(The_side));
   }
   NumberEquations();
}


Mesh::~Mesh()
{
   if (Verbosity>5)
      cout << "Removing Mesh instance ..." << endl;
   /*   for (size_t i=0; i<_nodes.size(); i++) {
      if (_nodes[i])
         delete _nodes[i];
   }
   for (size_t i=0; i<_elements.size(); i++) {
      if (_elements[i])
         delete _elements[i];
   }
   for (size_t i=0; i<Sides.size(); i++) {
      if (Sides[i])
         delete Sides[i];
   }
   for (size_t i=0; i<_edges.size(); i++) {
      if (_edges[i])
         delete _edges[i];
         }*/
}


void Mesh::set1D(real_t xmin,
                 real_t xmax,
                 size_t nb_el,
                 size_t p,
                 size_t nb_dof)
{
   size_t NbN = nb_el*p + 1;

// Insert nodes
   real_t xx=xmin, h=(xmax-xmin)/real_t(nb_el);
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
         el->Add(theNodes[nn++]);
      nn--;
      el->setCode(1);
      el->getMeasure();
      Add(el);
   }
   setNodesForDOF();
   theBoundaryNodes.push_back(theNodes[0]);
   theBoundaryNodes.push_back(theNodes[_nb_nodes-1]);

// Insert sides
   the_side = new Side(1,POINT);
   The_side.Add(theNodes[0]);
   The_side.setCode(1,0);
   if (theNodes[0]->getCode(1)==0)
      The_side.setCode(1,1);
   Add(the_side);
   the_side = new Side(2,POINT);
   The_side.Add(theNodes[_nb_nodes-1]);
   The_side.setCode(1,0);
   if (theNodes[0]->getCode(1)==0)
      The_side.setCode(1,2);
   Add(the_side);
}


Mesh &Mesh::operator*=(real_t a)
{
   node_loop(this) {
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
   element_loop(this) {
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
   if (nd)
      nd->setLabel(n2);
   else
      throw OFELIException("Mesh::RenumberNode(size_t,size_t): Node with label "
                           + itos(n1) + " does not exist.");
}


void Mesh::RenumberElement(size_t n1,
                           size_t n2)
{
   Element *el = getPtrElement(n1);
   if (el)
      el->setLabel(n2);
   else
      throw OFELIException("Mesh::RenumberElement(size_t,size_t): Element with label " +
                           itos(n1)+" does not exist.");
}


void Mesh::RenumberSide(size_t n1,
                        size_t n2)
{
   Side *sd = getPtrSide(n1);
   if (sd)
      sd->setLabel(n2);
   else
      throw OFELIException("Mesh::RenumberSide(size_t,size_t): Side with label " +
                           itos(n1)+" does not exist.");
}


void Mesh::RenumberEdge(size_t n1,
                        size_t n2)
{
   Edge *ed = getPtrEdge(n1);
   if (ed)
      ed->setLabel(n2);
   else
      throw OFELIException("Mesh::RenumberEdge(size_t,size_t): Edge with label " +
                           itos(n1)+" does not exist.");
}


void Mesh::Deform(const Vect<real_t>& u,
                  real_t              rate)
{
   if (rate<0. || rate>1.)
      throw OFELIException("Mesh::Deform(Vect<real_t>,real_t): Illegal value of rate.");
   real_t a = 1.;
   node_loop(this) {
      for (size_t i=1; i<=_dim; ++i) {
         if (fabs(The_node.getCoord(i)/u(node_label,i))<OFELI_EPSMCH)
            a = std::max(a,The_node.getCoord(i)/u(node_label,i));
      }
   }

   node_loop(this) {
      Point<real_t> x = The_node.getCoord();
      for (size_t i=1; i<=_dim; ++i)
         The_node.setCoord(i,x(i)+a*rate*u(node_label,i));
   }
}


void Mesh::inCoarse(Mesh& ms,
                    bool  test_el)
{
   size_t i;
   real_t x[3], y[3];

   for (i=0; i<_nb_nodes; i++)
      _node_in_coarse_element.push_back(nullptr);

   element_loop(&ms) {
      if (The_element.getShape() != TRIANGLE)
         throw OFELIException("Mesh::inCoarse(Mesh,bool): Element "+itos(element_label) +
                              " is not a triangle.");
      for (i=0; i<3; i++) {
         x[i] = (The_element)(i+1)->getCoord(1);
         y[i] = (The_element)(i+1)->getCoord(2);
      }
      real_t d = (x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]);

      node_loop(this) {
         Point<real_t> c = The_node.getCoord();
         real_t s = (y[2]-y[0])*(c.x-x[0])+(x[0]-x[2])*(c.y-y[0]);
         real_t t = (y[0]-y[1])*(c.x-x[0])+(x[1]-x[0])*(c.y-y[0]);
         if ((s>=0) && (s<=d) && (t>=0) && ((s+t)<=d))
            _node_in_coarse_element[node_label-1] = the_element;
      }

      if (test_el)
         element_loop(this) {
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
      _node_in_fine_element.push_back(nullptr);

   element_loop(&ms) {
      if (The_element.getShape() != TRIANGLE)
         throw OFELIException("Mesh::inFine(Mesh,bool): This function is valid for "
                              "triangles only.");
      for (i=0; i<3; i++) {
         x[i] = The_element(i+1)->getCoord(1);
         y[i] = The_element(i+1)->getCoord(2);
      }
      real_t d = (x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]);

      node_loop(this) {
         Point<real_t> c = the_node->getCoord();
         real_t s = (y[2]-y[0])*(c.x-x[0])+(x[0]-x[2])*(c.y-y[0]);
         real_t t = (y[0]-y[1])*(c.x-x[0])+(x[1]-x[0])*(c.y-y[0]);
         if ((s>=0) && (s<=d) && (t>=0) && ((s+t)<=d))
            _node_in_fine_element[node_label-1] = the_element;
      }

      if (test_el)
         element_loop(this) {
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
   if (Verbosity > 5)
      cout << "Deleting node " << nd->n() << endl;
   if (!nd)
      throw OFELIException("Mesh::Delete(Node *): Node does not exist.");
}


void Mesh::Delete(Element* el)
{
   if (Verbosity > 5)
      cout << "Deleting element " << el->n() << endl;
   if (!el)
      throw OFELIException("Mesh::Delete(Element *): Element does not exist.");
}


void Mesh::Delete(Side* sd)
{
   if (Verbosity > 5)
      cout << "Deleting side " << sd->n() << endl;
   if (!sd)
      throw OFELIException("Mesh::Delete(Side *): Side does not exist.");
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
   node_loop(this)
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
         if (Verbosity > 6)
            cout << "Eliminating imposed d.o.f. from list of equations ..." << endl;
         node_loop(this) {
            for (size_t i=1; i<=The_node.getNbDOF(); i++) {
               if (The_node.getCode(i) != c)
                  The_node._dof[i-1] = ++_nb_eq;
               else
                  The_node._dof[i-1] = 0;
            }
         }
      }
      else {
         if (Verbosity > 6)
            cout << "Numbering node supported d.o.f. ..." << endl;
         _nb_eq = 0;
         node_loop(this) {
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
         if (Verbosity > 6)
            cout << "Numbering side supported d.o.f. ..." << endl;
         size_t l=0;
         side_loop(this) {
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
         if (Verbosity > 6)
            cout << "Eliminating imposed d.o.f. from list of equations ..." << endl;
         _nb_dof = 0;
         side_loop(this) {
            the_side->setDOF(_first_dof,the_side->getNbDOF());
            _nb_dof += the_side->getNbDOF();
         }
         _nb_eq = _nb_dof;
      }
   }

// Edge supported d.o.f.
   else if (dof_type==EDGE_DOF) {
      if (_no_imposed_dof==false) {
         if (Verbosity > 6)
            cout << "Numbering edge supported d.o.f. ..." << endl;
         if (dof) {
            edge_loop(this)
               The_edge.DOF(dof,++_nb_eq);
         }
      }
      else {
         if (Verbosity > 6)
            cout << "Eliminating imposed d.o.f. from list of equations ..." << endl;
         if (dof) {
            edge_loop(this)
               for (size_t i=1; i<=The_edge.getNbDOF(); i++)
                  if (The_edge.getCode(i) != c)
                     The_edge.DOF(i,++_nb_eq);
                  else
                     The_edge.DOF(i,0);
         }
         else {
            edge_loop(this)
               if (The_edge.getCode(dof) != c)
                  The_edge.DOF(dof,++_nb_eq);
               else
                  The_edge.DOF(dof,0);
         }
      }
   }

// Element supported d.o.f.
   else if (dof_type==ELEMENT_DOF) {
      getAllSides();
      getElementNeighborElements();
      if (_no_imposed_dof==false) {
         if (Verbosity > 6)
            cout << "Numbering element supported d.o.f. ..." << endl;
         if (dof)
            element_loop(this)
               The_element.setDOF(dof,++_nb_eq);
         else
            element_loop(this)
               for (size_t i=1; i<=The_element.getNbDOF(); i++)
                  The_element.setDOF(i,++_nb_eq);
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
         if (Verbosity > 6)
            cout << "Eliminating imposed d.o.f. from list of equations ..." << endl;
         node_loop(this) {
            for (size_t i=1; i<=the_node->getNbDOF(); i++) {
               if (The_node.getCode(i) == 0)
                  The_node._dof[i-1] = ++_nb_eq;
               else
                  The_node._dof[i-1] = 0;
            }
         }
      }
      else {
         if (Verbosity > 6)
            cout << "Numbering node supported d.o.f. ..." << endl;
         _nb_eq = 0;
         node_loop(this) {
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
         if (Verbosity > 6)
            cout << "Numbering side supported d.o.f. ..." << endl;
         size_t l=0;
         side_loop(this) {
            for (size_t k=1; k<=_nb_side_nodes; k++) {
               for (size_t i=1; i<=the_side->getNbDOF(); i++) {
                  if (The_side.getCode(i) <= 0)
                     The_side._dof[l++] = ++_nb_eq;
                  else
                     The_side._dof[l++] = 0;
                  _nb_dof++;
	       }
            }
         }
      }
      else {
         if (Verbosity > 6)
            cout << "Eliminating imposed d.o.f. from list of equations ..." << endl;
         _nb_dof = 0;
         side_loop(this) {
            The_side.setDOF(_first_dof,the_side->getNbDOF());
            _nb_dof += The_side.getNbDOF();
         }
         _nb_eq = _nb_dof;
      }
   }

// Edge supported d.o.f.
   else if (dof_type == EDGE_DOF) {
      if (_no_imposed_dof==false) {
         if (Verbosity > 6)
            cout << "Numbering edge supported d.o.f. ..." << endl;
         if (dof) {
            edge_loop(this)
               The_edge.DOF(dof,++_nb_eq);
         }
         else {
            edge_loop(this) {
               for (size_t i=1; i<=the_edge->getNbDOF(); i++)
                  if (The_edge.getCode(i) == 0)
                     The_edge.DOF(i,++_nb_eq);
                  else
                     The_edge.DOF(i,0);
            }
         }
      }
      else {
         if (Verbosity > 6)
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
         if (Verbosity > 6)
            cout << "Numbering element supported d.o.f. ..." << endl;
         if (dof)
            element_loop(this)
               The_element.setDOF(dof,++_nb_eq);
         else
            element_loop(this)
               for (size_t i=1; i<=the_element->getNbDOF(); i++)
                  The_element.setDOF(i,++_nb_eq);
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
         if (Verbosity > 6)
            cout << "Eliminating imposed d.o.f. from list of equations ..." << endl;
         node_loop(this) {
            for (size_t i=dof1; i<=dof2; i++) {
               if (The_node.getCode(i) == 0) {
                  The_node._dof[i-1] = ++_nb_eq;
                  _dof_nbeq[dof1-1]++;
               }
               else
                  The_node._dof[i-1] = 0;
            }
         }
      }
      else {
         if (Verbosity > 6)
            cout << "Numbering node supported d.o.f. ..." << endl;
         _dof_nbeq[dof1-1] = _nb_nodes*(dof2-dof1+1);
         node_loop(this)
            The_node.setDOF(_first_dof,dof2-dof1+1);
      }
   }

// Side supported d.o.f.
   else if (dof_type==SIDE_DOF) {
      getAllSides();
      _nb_dof = 0;
      if (_no_imposed_dof==true) {
         if (Verbosity > 6)
            cout << "Numbering side supported d.o.f. ..." << endl;
         size_t l=0;
         side_loop(this) {
            for (size_t k=1; k<=_nb_side_nodes; k++) {
               for (size_t i=dof1; i<=dof2; i++) {
                  if (The_side.getCode(i) <= 0) {
                     The_side._dof[l++] = ++_nb_eq;
                     _dof_nbeq[dof1-1]++;
                  }
                  else
                     The_side._dof[l++] = 0;
                  _nb_dof++;
               }
            }
         }
      }
      else {
         if (Verbosity > 6)
            cout << "Eliminating imposed d.o.f. from list of equations ..." << endl;
         _nb_dof = 0;
         side_loop(this) {
            The_side.setDOF(_first_dof,the_side->getNbDOF());
            _nb_dof += The_side.getNbDOF();
         }
         _dof_nbeq[dof1-1] = _nb_eq = _nb_dof;
      }
   }

// Edge supported d.o.f.
   else if (dof_type==EDGE_DOF) {
      if (_no_imposed_dof==false) {
         if (Verbosity > 6)
            cout << "Numbering edge supported d.o.f. ..." << endl;
         edge_loop(this) {
            for (size_t i=dof1; i<=dof2; i++)
               if (The_edge.getCode(i) == 0) {
                  The_edge.DOF(i,++_nb_eq);
                  _dof_nbeq[dof1-1]++;
               }
               else
                  The_edge.DOF(i,0);
         }
      }
      else {
         if (Verbosity > 6)
            cout << "Eliminating imposed d.o.f. from list of equations ..." << endl;
         edge_loop(this) {
            for (size_t i=dof1; i<=dof2; i++) {
               if (The_edge.getCode(i) == 0) {
                  The_edge.DOF(i,++_nb_eq);
                  _dof_nbeq[dof1-1]++;
               }
               else
                  The_edge.DOF(i,0);
            }
         }
      }
   }

// Element supported d.o.f.
   else if (dof_type==ELEMENT_DOF) {
      getElementNeighborElements();
      if (_no_imposed_dof==false) {
         if (Verbosity > 6)
            cout << "Numbering element supported d.o.f. ..." << endl;
         element_loop(this)
            for (size_t i=dof1; i<=dof2; i++) {
               The_element.setDOF(i,++_nb_eq);
              _dof_nbeq[dof1-1]++;
            }
      }
   }
}


int Mesh::getAllSides(int opt)
{
   if (_dim==1) {
      throw OFELIException("Mesh::getAllSides(): No sides can be created for 1-D meshes");
      return 0;
   }
   if (_all_sides_created==true)
      return _nb_sides;
   theBoundarySides.clear();

   size_t k, ns=0, i;
   vector<vector<size_t> > nsd;
   int sh;
   if (Verbosity > 1)
      cout << "Creating all mesh sides ..." << endl;
   size_t nb_all_sides = 0;

// Store already defined sides in (*sides)
   vector<ND> sides;
   vector<Side *> SL(theSides);
   size_t nbss=_nb_sides;

// Add all mesh sides
   nb_all_sides = 0;
   element_loop(this) {
      if (_dim==3 && The_element.getShape()!=TETRAHEDRON)
         throw OFELIException("Mesh::getAllSides(): Member function is not valid for this element type.");
      ns = init_side_node_numbering(The_element.getShape(),nsd,sh);
      for (i=0; i<The_element.getNbSides(); i++) {
         size_t i1 = The_element(nsd[i][0])->n(),
                i2 = The_element(nsd[i][1])->n(),
                i3 = 0, i4 = 0;
         if (ns>2)
            i3 = The_element(nsd[i][2])->n();
         if (ns>3)
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
   vector<ND>::iterator it = sides.begin()+nb_all_sides;
   sort(sides.begin(),it,compare_sides);
   complete_sides(sides);
   vector<ND>::iterator new_end = std::unique(sides.begin(),sides.end());
   sides.erase(new_end,sides.end());
   _nb_sides = 0;
   theSides.clear();
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
      element_loop(this)
         The_element.setGlobalToLocal();
      side_loop(this) {
         Element *el1=The_side.getNeighborElement(1),
                 *el2=The_side.getNeighborElement(2);
         if (el1)
            el1->Add(the_side);
         if (el2)
            el2->Add(the_side);
      }
   }
   else {
      side_loop(this) {
         if ((The_side.getNeighborElement(1)))
            The_side.getNeighborElement(1)->Add(the_side);
         if ((The_side.getNeighborElement(2)))
            The_side.getNeighborElement(2)->Add(the_side);
      }
   }

// Assign DOFs to sides
   size_t n=1;
   side_loop(this) {
      The_side.setFirstDOF(n);
      for (size_t i=1; i<=The_side.getNbDOF(); i++)
         The_side.DOF(i,n++);
   }

// Mark nodes on boundary
   side_loop(this) {
      if (The_side.isOnBoundary())
         for (size_t j=1; j<=The_side.getNbNodes(); j++)
            getPtrNode(The_side.getNodeLabel(j))->setOnBoundary();
   }

   ND ss1, ss2;
   side_loop(this) {
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
   _all_sides_created = _boundary_sides_created = true;
   createBoundarySideList();
   return _nb_sides;
}


int Mesh::getAllEdges()
{
   if (_dim==1)
      throw OFELIException("Mesh::getAllEdges(): No edges can be created for 1-D meshes");
   if (_dim==2)
      throw OFELIException("Mesh::getAllEdges(): No edges can be created for 2-D meshes");
   if (_all_edges_created == true) {
      if (Verbosity > 0) {
         cout << "All edges created already." << endl;
      }
      return 0;
   }
   getAllSides();
   size_t i;
   if (Verbosity > 1)
      cout << "Creating all mesh edges ..." << endl;
   size_t nb_all_edges = 0;

// Store already defined edges in (*edges)
   vector<ND> edges;
   vector<Edge *> SL(theEdges);
   size_t nbss = _nb_edges;

// Add all mesh edges
   vector<vector<size_t> > nsd(3);
   nsd[0].push_back(1); nsd[0].push_back(2);
   nsd[1].push_back(2); nsd[1].push_back(3);
   nsd[2].push_back(3); nsd[2].push_back(1);
   ND eed;
   nb_all_edges = 0;
   side_loop(this) {
      if (The_side.getShape()!=TRIANGLE)
         throw OFELIException("Mesh::getAllEdges(): Edges can be created for triangular faces only");
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
   vector<ND>::iterator it = edges.begin()+nb_all_edges;
   sort(edges.begin(),it,compare_sides);
   complete_sides(edges);
   vector<ND>::iterator new_end = std::unique(edges.begin(),edges.end());
   edges.erase(new_end,edges.end());
   _nb_edges = 0;
   theEdges.clear();
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
   edge_loop(this) {
      if ((The_edge.getNeighborSide(1)))
         The_edge.getNeighborSide(1)->Add(the_edge);
      if ((The_edge.getNeighborSide(2)))
         The_edge.getNeighborSide(2)->Add(the_edge);
   }

// Assign DOFs to edges
   size_t n = 1;
   edge_loop(this) {
      The_edge.setFirstDOF(n);
      for (size_t i=1; i<=the_edge->getNbDOF(); i++)
         The_edge.DOF(i,n++);
   }

// Mark nodes on boundary
   edge_loop(this) {
      if (The_edge.isOnBoundary()) {
         getPtrNode(the_edge->getNodeLabel(1))->setOnBoundary();
         getPtrNode(the_edge->getNodeLabel(2))->setOnBoundary();
      }
   }

//
   ND ss1, ss2;
   edge_loop(this) {
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
   side_loop(this) {
      if (The_side.isOnBoundary()) {
         theBoundarySides.push_back(the_side);
         _nb_boundary_sides++;
      }
   }
   return _nb_boundary_sides;
}


int Mesh::createInternalSideList()
{
   _nb_internal_sides = 0;
   side_loop(this) {
      if (The_side.isOnBoundary()==0) {
         theInternalSides.push_back(the_side);
         _nb_internal_sides++;
      }
   }
   return _nb_internal_sides;
}


int Mesh::getBoundarySides()
{
   if (_boundary_sides_created == true)
      return _nb_boundary_sides;
   size_t k, ns=0, i;
   vector<vector<size_t> > nsd;
   int sh;
   if (Verbosity > 1)
      cout << "Creating boundary sides ..." << endl;
   size_t nb_all_sides = 0;

// Store already defined sides in (*sides)
   vector<ND> sides;
   vector<Side *> SL(theSides);
   size_t nbss = _nb_sides;

// Add all mesh sides
   size_t i3, i4;
   nb_all_sides = 0;

   vector<size_t> s(2);
   element_loop(this) {
     for (size_t i=1; i<=3; i++) {
        s[0] = The_element(i)->n(), s[1] = The_element(i%3+1)->n();
        side_loop(this) {
           if (The_side.getNeighborElement(1)==nullptr) {
              if (equal_sides(the_side,s)) {
                 The_side.set(the_element,1);
                 break;
              }
           }
        }
     }
   }

   element_loop(this) {
      if (_dim==3 && the_element->getShape() != TETRAHEDRON)
         throw OFELIException("Mesh::getBoundarySides(): Member function not valid for this element type.");
      ns = init_side_node_numbering(the_element->getShape(), nsd, sh);
      for (i=0; i<The_element.getNbSides(); i++) {
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
   vector<ND>::iterator it = sides.begin() + nb_all_sides;
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

   side_loop(this) {
      Element *el;
      if ((el=The_side.getNeighborElement(1)))
         The_side.getNeighborElement(1)->Add(the_side);
      if ((el=The_side.getNeighborElement(2)))
         The_side.getNeighborElement(2)->Add(the_side);
   }

   size_t n = 1;
   side_loop(this) {
      The_side.setFirstDOF(n);
      for (size_t i=1; i<=The_side.getNbDOF(); i++)
         The_side.DOF(i,n++);
   }

   side_loop(this) {
      if (The_side.isOnBoundary()) {
         for (size_t j=1; j<=The_side.getNbNodes(); j++)
            getPtrNode(The_side(j)->n())->setOnBoundary();
      }
   }

   ND ss1, ss2;
   side_loop(this) {
      ss1.nd[0] = The_side(1)->n();
      ss1.nd[1] = The_side(2)->n();
      if (The_side.getNbNodes()>2)
         ss1.nd[2] = The_side(3)->n();
      if (The_side.getNbNodes()>3)
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
   side_loop(this) {
      theSides.push_back(the_side);
      if (The_side.isOnBoundary()) {
         theBoundarySides.push_back(the_side);
         _nb_boundary_sides++;
      }
   }
   _nb_sides = _nb_boundary_sides;
   return _nb_boundary_sides;
}


int Mesh::getBoundaryNodes()
{
   if (_boundary_nodes_created)
      return _nb_boundary_nodes;
   getBoundarySides();
   _nb_boundary_nodes = 0;
   node_loop(this) {
      if (The_node.isOnBoundary()) {
         theBoundaryNodes.push_back(the_node);
         _nb_boundary_nodes++;
      }
   }
   _boundary_nodes_created = true;
   return _nb_boundary_nodes;
}


void Mesh::getNodeNeighborElements()
{
   if (_node_neighbor_elements_created) {
      if (Verbosity > 1)
         cout << "List of node neighbor elements already created." << endl;
      return;
   }
   if (Verbosity > 1)
      cout << "Creating node neighbor elements ..." << endl;
/*   mesh_elements(*this)
      for (size_t i=1; i<=the_element->getNbNodes(); i++)
         The_element(i)->Neig();*/
   node_loop(this)
      The_node.Add();
   element_loop(this)
      for (size_t i=1; i<=The_element.getNbNodes(); i++)
         The_element(i)->Add(the_element);
   _node_neighbor_elements_created = true;
}


void Mesh::getElementNeighborElements()
{
   if (_element_neighbor_elements_created) {
      if (Verbosity > 1)
         cout << "List of element neighbor elements already created." << endl;
      return;
   }
   if (Verbosity > 1)
      cout << "Creating element neighbor elements ..." << endl;
   getAllSides();
   side_loop(this) {
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
   element_loop(this) {
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
   if (opt==NODE_DOF)
      _set_nodes = true;
   else if (opt==SIDE_DOF) {
      _set_sides = true;
      _nb_side_nodes = nb_nodes;
   }
   else if (opt==EDGE_DOF)
      _set_edges = true;
   else if (opt==ELEMENT_DOF) {
      _set_elements = true;
      _nb_element_nodes = nb_nodes;
   }
   else
      ;
   NumberEquations();
}


void Mesh::AddMidNodes(int g)
{
   Point<real_t> x;
   size_t mid_label = _nb_nodes + 1;
   int code[6];
   getAllSides();

   side_loop(this) {
      x.x = 0.5*(The_side(1)->getCoord(1) + The_side(2)->getCoord(1));
      x.y = 0.5*(The_side(1)->getCoord(2) + The_side(2)->getCoord(2));
      size_t nb_dof = std::max(The_side(1)->getNbDOF(),The_side(2)->getNbDOF());
      Node *nd = new Node(mid_label++,x);
      nd->setNbDOF(nb_dof);
      for (size_t i=0; i<The_side(1)->getNbDOF(); i++)
         code[i] = std::min(The_side(1)->getCode(i+1),The_side(2)->getCode(i+1));
      if (The_side.isOnBoundary())
         nd->setOnBoundary();
      else {
         for (size_t i=0; i<The_side(1)->getNbDOF(); i++)
            code[i] = 0;
      }
      nd->setCode(code);
      nd->setDOF(_first_dof,nb_dof);
      Add(nd);
      the_side->Add(nd);
   }

   element_loop(this) {
      if (g) {
         x = 0;
         for (size_t j=1; j<=the_element->getNbNodes(); j++)
            x += The_element(j)->getCoord();
         x /= the_element->getNbNodes();
         Node *nd = new Node(mid_label++,x);
         nd->setNbDOF(The_element(1)->getNbDOF());
         for (size_t i=0; i<The_element(1)->getNbDOF(); i++)
            code[i] = The_element(1)->getCode(i+1);
         nd->setCode(code);
         nd->setDOF(_first_dof,The_element(1)->getNbDOF());
         Add(nd);
         The_element.Add(nd);
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
   if (n>_nb_nodes)
     throw OFELIException("Mesh::set(Node *): Node label larger than total number of nodes.");
   theNodes[n-1] = nd;
}


void Mesh::set(Element* el)
{
   size_t n = el->n();
   if (n>_nb_elements)
      throw OFELIException("Mesh::set(Element *): Element label larger than total number of elements.");
   theElements[n-1] = el;
}


void Mesh::set(Side* sd)
{
   size_t n = sd->n();
   if (n>_nb_sides)
      throw OFELIException("Mesh::set(Side *): Side label larger than total number of sides.");
   theSides[n-1] = sd;
}


void Mesh::Add(Node* nd)
{
   theNodes.push_back(nd);
   _nb_nodes++, _nb_vertices++;
   _nb_dof += nd->getNbDOF();
}


void Mesh::Add(Element* el)
{
   theElements.push_back(el);
   _nb_elements++;
}


void Mesh::Add(Side* sd)
{
   theSides.push_back(sd);
   _nb_sides++;
}


void Mesh::Add(Edge* ed)
{
   theEdges.push_back(ed);
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
   node_loop(this) {
      The_node.setCoord(1,sx*The_node.getCoord(1));
      The_node.setCoord(2,sy*The_node.getCoord(2));
      The_node.setCoord(3,sz*The_node.getCoord(3));
   }
}


void Mesh::getList(vector<Node *>& nl) const
{
   for (size_t i=0; i<_nb_nodes; i++)
      nl.push_back(theNodes[i]);
}


void Mesh::getList(vector<Element *>& el) const
{
   for (size_t i=0; i<_nb_elements; i++)
      el.push_back(theElements[i]);
}


void Mesh::getList(vector<Side *>& sl) const
{
   for (size_t i=0; i<_nb_sides; i++)
      sl.push_back(theSides[i]);
}


void Mesh::setList(const vector<Node *>& nl)
{
   theNodes.clear();
   _nb_nodes = 0;
   for (size_t i=0; i<nl.size(); i++)
      Add(nl[i]);
}


void Mesh::setList(const vector<Element *>& el)
{
   theElements.clear();
   _nb_elements = 0;
   for (size_t i=0; i<el.size(); i++)
      Add(el[i]);
}


void Mesh::setList(const vector<Side *>& sl)
{
   theSides.clear();
   _nb_sides = 0;
   for (size_t i=0; i<sl.size(); i++)
      Add(sl[i]);
}


void Mesh::checkNodeLabels()
{
   for (size_t i=0; i<_nb_nodes; i++) {
      if (theNodes[i]->n()>_nb_nodes)
         throw OFELIException("Mesh::checkNodeLabels(): The node label: " + itos(theNodes[i]->n()) + 
                              " exceeds the total number of nodes.");
   }
   if (_nb_nodes>0)
      qksort(theNodes,0,_nb_nodes-1,_node_compare);
}


void Mesh::checkElementLabels()
{
   for (size_t i=0; i<_nb_elements; i++) {
      if (theElements[i]->n()>_nb_elements)
         throw OFELIException("Mesh::checkElementLabels(): The element label: " + itos(theElements[i]->n()) + 
                              " exceeds the total number of elements.");
   }
   if (_nb_elements>0)
      qksort(theElements,0,_nb_elements-1,_element_compare);
}


void Mesh::checkSideLabels()
{
   for (size_t i=0; i<_nb_sides; i++) {
      if (theSides[i]->n()>_nb_sides)
         throw OFELIException("Mesh::checkSideLabels(): The side label: " + itos(theSides[i]->n()) + 
                              " exceeds the total number of sides.");
   }
   if (_nb_sides>0)
      qksort(theSides,0,_nb_sides-1,_side_compare);
}


void Mesh::AddNodes(int p)
{
   element_loop(this) {
      if (The_element.getShape() != TRIANGLE)
         throw OFELIException("Mesh::AddNodes(int): This function is valid for triangles only.");
   }
   getAllSides();

   vector<int> code;
   size_t nb_dof, i, mid_label=_nb_nodes+1;
   Point<real_t> x;

   side_loop(this) {
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
          The_side.Add(nd);
       }
   }

   code.clear();
   element_loop(this) {
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
   vector<Element *> elements(theElements);
   theElements.clear();
   for (vector<Element *>::iterator it=elements.begin(); it!=elements.end(); it++) {
      if ((*it)->isActive()) {
         (*it)->setLabel(++_nb_elements);
         theElements.push_back(*it);
      }
   }
   _nb_sides = 0;
   vector<Side *> sides(theSides);
   theSides.clear();
   for (vector<Side *>::iterator it=sides.begin(); it!=sides.end(); it++) {
      if ((*it)->isActive()) {
         (*it)->setLabel(++_nb_sides);
         theSides.push_back(*it);
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
   if (Verbosity > 1)
     cout << "Reordering mesh nodes ..." << endl;

   {
     vector<long> xadj(_nb_nodes+1);
     vector<size_t> mask(_nb_nodes+1), xls(_nb_nodes+1), ls(_available_memory), adjncy(_available_memory);
     FindGraph(xadj,adjncy);
     GenRCM(xadj,adjncy,nn,mask,xls);
   }

   size_t n;
   Node *nd;
   vector<Node *> NL(theNodes);
   for (size_t j=0; j<_nb_nodes; j++) {
      nd = NL[j];
      n = nd->n();
      nd->setLabel(nn[n-1]);
      theNodes[nn[n-1]-1] = nd;
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
   if (Verbosity > 1)
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
   vector<Node *> NL(theNodes);
   for (size_t j=0; j<_nb_nodes; j++) {
      nd = NL[j];
      n = nd->n();
      nd->setLabel(nn[n-1]);
      theNodes[nn[n-1]-1] = nd;
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
   if (Verbosity > 1)
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
   if (Verbosity > 2)
      cout << "Mesh successfully read." << endl;
   NumberEquations();
   element_loop(this)
      theMaterial.check(The_element.getCode());

// Fill list of marked nodes
   _nb_marked_nodes = 0;
   node_loop(this) {
      bool mark = false;
      for (i=1; i<=The_node.getNbDOF(); i++) {
         if (The_node.getCode(i)>0)
            mark = true;
      }
      if (mark) {
         theMarkedNodes.push_back(the_node);
         _nb_marked_nodes++;
      }
   }
}


void Mesh::get(const string& mesh_file,
               int           ff, 
               int           nb_dof)
{
   if (ff==OFELI_FF)
      get(mesh_file);
   else if (ff==GMSH)
      getGmsh(mesh_file,*this,nb_dof);
   else if (ff==MATLAB)
      getMatlab(mesh_file,*this,nb_dof);
   else if (ff==EASYMESH)
      getEasymesh(mesh_file,*this,nb_dof);
   else if (ff==GAMBIT)
      getGambit(mesh_file,*this,nb_dof);
   else if (ff==BAMG)
      getBamg(mesh_file,*this,nb_dof);
   else if (ff==NETGEN)
      getNetgen(mesh_file,*this,nb_dof);
   else if (ff==TRIANGLE_FF)
      getTriangle(mesh_file,*this,nb_dof);
   else
      throw OFELIException("Mesh::get(string,int,int): Unknown file format "+itos(ff));
   NumberEquations();
   element_loop(this)
      theMaterial.check(The_element.getCode());
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
   size_t i;
   int sign, m;
   string sh[11] = {"none","point","line","triangle","quadrilateral","tetrahedron",
                    "hexahedron","pentahedron"};
   fp.open(file.c_str(),ios::out);
   if (Verbosity > 1)
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
      fp << "   <Nodes>" << endl;
      node_loop(this) {
         fp.setf(ios::right|ios::scientific);
         for (i=1; i<=_dim; i++)
            fp << "  " << setprecision(8) << setw(18) << The_node.getCoord(i);
         sign = 1;
         if (The_node.getCode(1)<0)
            sign = -1;
         m = 0;
         for (size_t j=1; j<=n; j++)
            m += abs(The_node.getCode(j))*size_t(pow(10.,real_t(n-j)));
         m *= sign;
         fp << setw(10) << m << endl;
      }
      fp << "   </Nodes>" << endl;
   }

   if (_nb_elements>0) {
      size_t nbn = getPtrElement(1)->getNbNodes();
      string shape = sh[getPtrElement(1)->getShape()];
      fp << "   <Elements shape=\"" << shape << "\"  nodes=\"" << nbn << "\">" << endl;
      element_loop(this) {
         for (i=1; i<=nbn; i++)
            fp << setw(8) << The_element(i)->n();
         fp << setw(8) << the_element->getCode() << endl;
      }
      fp << "   </Elements>" << endl;
   }

   if (_nb_sides>0) {
      size_t nbn = getPtrSide(1)->getNbNodes();
      string shape = sh[getPtrSide(1)->getShape()];
      fp << "   <Sides shape=\"" << shape << "\"  nodes=\"" << nbn << "\">" << endl;
      side_loop(this) {
         for (i=1; i<=nbn; i++)
            fp << setw(8) << The_side(i)->n();
         m = 0;
         for (size_t j=1; j<=n; j++)
            m += The_side.getCode(j)*size_t(pow(10.,real_t(n-j)));
         fp << setw(6) << m << endl;
      }
      fp << "   </Sides>" << endl;
   }
   if (_nb_mat>1 || theMaterial.getName(1)!="Generic") {
      fp << "   <Material>" << endl;
      for (size_t i=1; i<=_nb_mat; i++)
         fp << setw(9) << theMaterial.getCode(i) << "   "
            << theMaterial.getName(theMaterial.getCode(i)) << endl;
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
   if (Verbosity > 1)
      cout << "Getting mesh graph ..." << endl;
   for (size_t i=1; i<=_nb_nodes; ++i)
      xadj[i] = 0;
   xadj[0] = 1;
   for (size_t i=0; i<_available_memory; ++i)
      adjncy[i] = 0;
   size_t memory = 0;

   element_loop(this) {
      for (size_t i=1; i<=the_element->getNbNodes(); i++) {
         size_t ii = the_element->getNodeLabel(i);
         for (size_t j=1; j<=the_element->getNbNodes(); j++) {
            size_t jj = The_element(j)->n();
            if (ii != jj) {
               size_t k = 0;
               for (size_t in=1; in<=ii; ++in)
                  k += xadj[in-1];
               size_t m = k + xadj[ii];
               if (adjncy[k-1] == 0) {
                  adjncy[k-1] = jj;
                  ++memory;
                  ++xadj[ii];
                  goto L9;
               }
	       size_t l;
               for (l=k; l<=m-1; ++l) {
                  int jtest = adjncy[l-1] - jj;
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
               if (memory > _available_memory && _available_memory > 0)
                  throw OFELIException("Mesh::FindGraph(vector<long>,vector<size_t>): Insufficient Memory");
               ++xadj[ii];
               for (size_t is=m; is<=memory; ++is)
                  adjncy[m+memory-is] = adjncy[m+memory-is-1];
               adjncy[m-1] = jj;
               goto L9;
L7:
               ++memory;
               if (memory > _available_memory && _available_memory > 0)
                  throw OFELIException("Mesh::FindGraph(vector<long>,vector<size_t>): Insufficient Memory");
               ++xadj[ii];
               for (size_t is=l; is<=memory; ++is) {
                  size_t isn = memory + l - is;
                  adjncy[isn] = adjncy[isn-1];
               }
               adjncy[l-1] = jj;
L9:
               ;
            }
         }
      }
   }

   for (size_t i=1; i<=_nb_nodes; ++i)
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
   for (size_t i=0; i<_nb_nodes; ++i)
      mask[i] = 1;
   size_t num = 1;
   for (size_t i=0; i<_nb_nodes; ++i) {
      if (mask[i] != 0) {
         size_t root=i+1;
	 long nlvl;
         FindRoot(root,xadj,adjncy,mask,nlvl,xls,perm+num-1);
	 size_t ccsize = 0;
         RCM(root,xadj,adjncy,mask,perm+num-1,ccsize,xls);
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
   _first_dof = 1;
   _nb_dof = 0;
   theNodes = ms.theNodes;
   theElements = ms.theElements;
   theSides = ms.theSides;
   _nb_nodes = _nb_elements = _nb_sides = _nb_vertices = 0;
   _max_nb_nodes = ms._max_nb_nodes;
   _max_nb_elements = ms._max_nb_elements;
   _max_nb_sides = ms._max_nb_sides;
   _all_sides_created = ms._all_sides_created;
   _boundary_sides_created = ms._boundary_sides_created;
   _all_edges_created = ms._all_edges_created;
   _boundary_edges_created = ms._boundary_edges_created;

// Insert nodes
   node_loop(&ms) {
      Add(new Node(The_node));
      _nb_eq = _nb_dof;
   }

// Insert elements
   element_loop(&ms)
      Add(new Element(The_element));

// Insert sides
   side_loop(&ms)
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


void Mesh::outputNodes(ostream& s) const
{
   size_t n = _nb_nodes;
   if (Verbosity==2)
      n = !(10<n) ? n : 10;
   else if (Verbosity==3)
      n = !(50<n) ? n : 50;
   else if (Verbosity==4)
      n = !(100<n) ? n : 100;
   if (n>0)
      s << "\n\nLIST OF NODES:\n";
   for (size_t i=0; i<n; ++i)
      s << *(theNodes[i]);
   s << endl;
}


void Mesh::outputElements(ostream& s) const
{
   size_t n = _nb_elements;
   if (Verbosity==2)
      n = !(10<n) ? n : 10;
   else if (Verbosity==3)
      n = !(50<n) ? n : 50;
   else if (Verbosity==4)
      n = !(100<n) ? n : 100;
   if (n>0)
      s << "\n\nLIST OF ELEMENTS:\n";
   for (size_t i=0; i<n; ++i)
      s << *(theElements[i]);
   s << endl;
}


void Mesh::outputSides(ostream& s) const
{
   size_t n = _nb_sides;
   if (Verbosity==2)
      n = !(10<n) ? n : 10;
   else if (Verbosity==3)
      n = !(50<n) ? n : 50;
   else if (Verbosity==4)
      n = !(100<n) ? n : 100;
   if (n>0)
      s << "\n\nLIST OF SIDES:\n";
   for (size_t i=0; i<n; ++i)
      s << *(theSides[i]);
   s << endl;
}


void Mesh::outputBoundarySides(ostream& s) const
{
   size_t n = _nb_boundary_sides;
   if (Verbosity==2)
      n = !(10<n) ? n : 10;
   else if (Verbosity==3)
      n = !(50<n) ? n : 50;
   else if (Verbosity==4)
      n = !(100<n) ? n : 100;
   if (n>0)
      s << "\n\nLIST OF BOUNDARY SIDES:\n";
   for (size_t i=0; i<n; ++i)
      s << *(theBoundarySides[i]);
   s << endl;
}


void Mesh::outputEdges(ostream& s) const
{
   size_t n = _nb_edges;
   if (Verbosity==2)
      n = !(10<n) ? n : 10;
   else if (Verbosity==3)
      n = !(50<n) ? n : 50;
   else if (Verbosity==4)
      n = !(100<n) ? n : 100;
   if (n>0)
      s << "\n\nLIST OF EDGES:\n";
   for (size_t i=0; i<n; ++i)
      s << *(theEdges[i]);
   s << endl;
}


ostream& operator<<(ostream&    s,
                    const Mesh& ms)
{
   s << "\n\nM E S H     D A T A\n===================\n\n";
   s << "Space Dimension             : " << setw(6) << ms.getDim() << endl;
   s << "Number of nodes             : " << setw(6) << ms.getNbNodes() << endl;
   s << "Number of elements          : " << setw(6) << ms.getNbElements() << endl;
   if (ms.getNbSides()>0)
      s << "Number of sides          : " << setw(6) << ms.getNbSides() << endl;
   if (ms.getNbBoundarySides()>0)
      s << "Number of boundary sides : " << setw(6) << ms.getNbBoundarySides() << endl;
   if (ms.getNbEdges()>0)
      s << "Number of edges          : " << setw(6) << ms.getNbEdges() << endl;

   ms.outputNodes(s);
   ms.outputElements(s);
   ms.outputSides(s);
   ms.outputBoundarySides(s);
   ms.outputEdges(s);
   s << endl;
   return s;
}

} /* namespace OFELI */
