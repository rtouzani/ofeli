/*==============================================================================

                                    O  F  E  L  I

                          Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2025 Rachid Touzani
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

                         Implementation of class 'Domain'

  ==============================================================================*/

#include <algorithm>

#include "mesh/Domain.h"
#include "mesh/saveMesh.h"
#include "mesh/getMesh.h"
#include "io/XMLParser.h"
#include "linear_algebra/Vect_impl.h"
#include "mesh/bamg/Mesh2.h"
#include "mesh/bamg/Meshio.h"
#include "mesh/bamg/QuadTree.h"
#include "mesh/bamg/R2.h"
#include "OFELIException.h"

using std::to_string;

namespace OFELI {

Domain::Domain()
       : _ff(nullptr), _theMesh(nullptr), _output_file("out.m"), _nb_vertices(0), _nb_lines(0),
         _nb_contours(0), _nb_holes(0), _sub_domain(0), _nb_required_vertices(0),
         _nb_required_edges(0), _nb_sub_domains(0), _ret_cont(1), _ret_save(1), _ret_sd(1)
{
}


Domain::Domain(const string& file)
       : _ff(nullptr), _theMesh(nullptr), _domain_file(file), _output_file("out.m"),
         _nb_vertices(0), _nb_lines(0), _nb_contours(0), _nb_holes(0), _sub_domain(0),
         _nb_required_vertices(0), _nb_required_edges(0), _nb_sub_domains(0), _ret_cont(1),
         _ret_save(1), _ret_sd(1)
{
   get(file);
}


Domain::~Domain()
{
   if (_ff!=nullptr)
      delete _ff;
   if (_theMesh!=nullptr)
      delete _theMesh;
}


Point<real_t> Domain::getMinCoord() const
{
   Point<real_t> a(_v[0].x,_v[0].y);
   for (size_t i=1; i<_nb_vertices; i++) {
      a.x = std::min(a.x,_v[i].x);
      a.y = std::min(a.y,_v[i].y);
   }
   return a;
}


Point<real_t> Domain::getMaxCoord() const
{
   Point<real_t> a(_v[0].x,_v[0].y);
   for (size_t i=1; i<_nb_vertices; i++) {
      a.x = std::max(a.x,_v[i].x);
      a.y = std::max(a.y,_v[i].y);
   }
   return a;
}


real_t Domain::getMinh() const
{
   real_t h = _v[0].h;
   for (size_t i=1; i<_nb_vertices; i++)
      h = std::min(h,_v[i].h);
   return h;
}


void Domain::insertRequiredVertex(size_t v)
{
   _required_vertex.push_back(v);
   _nb_required_vertices++;
}


void Domain::insertRequiredEdge(size_t e)
{
   _required_edge.push_back(e);
   _nb_required_edges++;
}


void Domain::insertVertex(real_t x,
                          real_t y,
                          real_t h,
                          int    code)
{
   Vertex v;
   v.label = ++_nb_vertices;
   v.x = x; v.y = y; v.z = 0.; v.h = h;
   v.code = code;
   _v.push_back(v);
}


void Domain::insertVertex(real_t x,
                          real_t y,
                          real_t z,
                          real_t h,
                          int    code)
{
   Vertex v;
   v.label = ++_nb_vertices;
   v.x = x; v.y = y; v.z = z; v.h = h;
   v.code = code;
   _v.push_back(v);
}


void Domain::getVertex()
{
   Vertex v;
   v.label = ++_nb_vertices;
   v.x = _ff->getD("x-coordinate: ");
   v.y = _ff->getD("y-coordinate: ");
   v.h = _ff->getD("Spacing: ");
   v.code = _ff->getI("Code: ");
   _v.push_back(v);
}


int Domain::getCurve()
{
   size_t i=0;

// Label of curve

// First end point
   size_t n1 = _ff->getI("First End Point: ");
   Ln ll;
   ll.n1 = n1;

// Second end point
   size_t n2 = _ff->getI("Second End Point: ");
   ll.n2 = n2;

// Type of curve
   int type = _ff->getI("Curve's type: ");

// Case of a straight line
   Vertex v;
   v.label = _nb_vertices + 1;
   v.x = 0.5*(_v[n1-1].x + _v[n2-1].x);
   v.y = 0.5*(_v[n1-1].y + _v[n2-1].y);
   v.h = 0.5*(_v[n1-1].h + _v[n2-1].h);
   if (type == 0) {
      ll.nb = 3;
      ll.node.push_back(_v[n1-1]);
      ll.node.push_back(v);
      ll.node.push_back(_v[n2-1]);
   }

// Case of a curve given by equation
   else if (type == -1) {
      string regex = _ff->getE("Curve's equation: ");
      real_t x=0., y=0., z=0.;
      _theFct.set(regex);
      size_t nb = _ff->getI("nb. of discretization points: ");
      if (nb<3)
         nb = 3;
      ll.nb = nb;
      ll.node.resize(nb+1);
      real_t vv = _theFct(_v[n1-1],0.);
      if (fabs(vv) > 1.e-8)
         throw OFELIException("Domain::getCurve(): )");
      vv = _theFct(_v[n2-1],0.);
      if (fabs(vv) > 1.e-8)
         throw OFELIException("In Domain::getCurve(): ");
      ll.node[0] = _v[n1-1];
      ll.node[1] = _v[n2-1];
      x = 0.9*ll.node[0].x + 0.1*ll.node[1].x;
      y = 0.9*ll.node[0].y + 0.1*ll.node[1].y;
      z = 0.1*ll.node[1].z;
      for (i=1; i<nb-1; i++) {
         Position(i/real_t(nb-1),x,y,z);
         ll.node[i].x = x;
         ll.node[i].y = y;
         ll.node[i].z = z;
      }
   }

// Case of a predefined curve
   else if (type==1)
      string shape = _ff->getS("Curve's shape: ");

// Dirichlet and Neumann codes
   ll.Dcode = _ff->getI("Dirichlet Code: ");
   ll.Ncode = _ff->getI("Neumann Code: ");
   v.code = ll.Dcode;
   _v.push_back(v);
   _l.push_back(ll);
   _nb_lines++;
   return 0;
}


int Domain::Position(real_t  s,
                     real_t& x,
                     real_t& y,
                     real_t& z)
{
   real_t x0, y0, xx, yy, zz=0., dxf, dyf, ex, ey, ez;
   real_t det, g;
   int it=1;
   x0 = xx = x, y0 = yy = y;
   x -= _theFct(xx,yy,zz,0.);
   y -= _theFct(xx,yy,zz,0.);
   ex = x - xx; ey = y - yy; ez = z - zz;
   if (fabs(ex*ex+ey*ey+ez*ez) < 1.e-8)
      return 0;
   while (it<50) {
      dxf = (_theFct(x,y,z,0.)-_theFct(xx,yy,zz,0.))/(x-xx);
      dyf = (_theFct(x,y,z,0.)-_theFct(xx,yy,zz,0.))/(y-yy);
      g = 0.5*(s*s - (xx-x0)*(xx-x0) - (yy-y0)*(yy-y0));
      xx = x; yy = y;
      det = (xx-x0)*dyf - (yy-y0)*dxf;
      x += ((y-y0)*_theFct(x,y,z,0.) + dyf*g)/det;
      y -= ((x-x0)*_theFct(x,y,z,0.) + dxf*g)/det;
      ex = (x - xx)*(x - xx);
      ey = (y - yy)*(y - yy);
      ez = (z - zz)*(z - zz);
      it++;
      if (fabs(ex+ey+ez) < 1.e-8)
         return it;
   }
   return 0;
}


void Domain::insertLine(size_t n1,
                        size_t n2,
                        int    code)
{
   int dc=0, nc=0;
   (code<0) ? nc = -code : dc = code;
   insertLine(n1,n2,dc,nc);
}


void Domain::insertLine(size_t n1,
                        size_t n2,
                        int    dc,
                        int    nc)
{
   Ln ll;
   ll.n1 = n1; ll.n2 = n2;
   ll.Dcode = dc; ll.Ncode = nc;
   _mark.push_back(dc);
   ll.nb = 2;
   ll.node.push_back(_v[n1-1]);
   //   ll.node.push_back(0.5*(_v[n1-1]+_v[n2-1]));
   ll.node.push_back(_v[n2-1]);
   _l.push_back(ll);
   _nb_lines++;
}


int Domain::getLine()
{
   Ln ll;
   size_t n1 = _ff->getI("First End Point: ");
   if (n1 > _nb_vertices)
      throw OFELIException("In Domain::getLine(): Vertex "+to_string(n1)+" undefined.");
   size_t n2 = _ff->getI("Second End Point: ");
   if (n2 > _nb_vertices)
      throw OFELIException("In Domain::getLine(): Vertex "+to_string(n2)+" undefined.");
   ll.n1 = n1, ll.n2 = n2;
   ll.Dcode = _ff->getI("Dirichlet Code: ");
   ll.Ncode = _ff->getI("Neumann Code: ");
   _mark.push_back(ll.Dcode);
   ll.nb = 2;
   ll.node.push_back(_v[n1-1]);
   ll.node.push_back(0.5*(_v[n1-1]+_v[n2-1]));
   ll.node.push_back(_v[n2-1]);
   _l.push_back(ll);
   _nb_lines++;
   return 0;
}


void Domain::insertCircle(size_t n1,
                          size_t n2,
                          size_t n3,
                          int    code)
{
   int nc=0, dc=0;
   (code<0) ? nc = -code : dc = code;
   insertCircle(n1,n2,n3,dc,nc);
}


void Domain::insertCircle(size_t n1,
                          size_t n2,
                          size_t n3,
                          int    dc,
                          int    nc)
{
   size_t nb = 36;
   Ln ll;
   ll.n1 = n1; ll.n2 = n2; ll.nb = nb;
   ll.node.push_back(_v[n1-1]);
   ll.node.push_back(_v[n2-1]);
   ll.Dcode = dc; ll.Ncode = nc;
   Point<real_t> cc(_v[n3-1]);
   Point<real_t> a1=_v[n1-1]-cc, a2=_v[n2-1]-cc;
   real_t h=_v[n1-1].h, r=sqrt((a1,a1));
   real_t theta1=atan2(a1.y,a1.x), theta2=atan2(a2.y,a2.x);
   real_t theta = theta1;
   if (theta1 >= theta2)
      theta2 += 2*OFELI_PI;
   for (size_t i=0; i<nb; i++) {
      Vertex v;
      ll.n1 = n1; ll.n2 = n2; ll.nb = nb;
      ll.Dcode = dc; ll.Ncode = nc;
      theta += (theta2-theta1)/nb;
      if (i<nb-1) {
         ll.n2 = _nb_vertices + 1;
         v.x = cc.x + r*cos(theta);
         v.y = cc.y + r*sin(theta);
         v.h = h;
         v.code = dc;
         _v.push_back(v);
      }
      if (i>0)
         ll.n1 = _nb_vertices;
      _l.push_back(ll);
      _nb_lines++;
      _nb_vertices++;
   }
   _nb_vertices--;
}


void Domain::getCircle()
{
   size_t nb = 36;
   size_t n1 = _ff->getI("First end vertex: ");
   Ln ln;
   Vertex v;
   ln.n1 = n1;
   ln.node.push_back(_v[n1-1]);
   if (n1 > _nb_vertices)
      throw OFELIException("In Domain::getCircle(): Vertex "+to_string(n1)+" is not defined.");
   size_t n2 = _ff->getI("Second end vertex: ");
   if (n2 > _nb_vertices)
      throw OFELIException("In Domain::getCircle(): Vertex "+to_string(n2)+" is not defined.");
   ln.n2 = n2;
   ln.node.push_back(_v[n2-1]);
   size_t n3 = _ff->getI("Center vertex: ");
   if (n3 > _nb_vertices)
      throw OFELIException("In Domain::getCircle(): Vertex "+to_string(n3)+" is not defined.");
   ln.nb = nb;

   int dc = ln.Dcode = _ff->getI("Dirichlet Code: ");
   ln.Ncode = _ff->getI("Neumann Code: ");

   real_t h = _v[n1-1].h;
   Point<real_t> cc(_v[n3-1]), a1=_v[n1-1]-cc, a2=_v[n2-1]-cc;
   real_t r = sqrt((a1,a1));
   real_t theta1 = atan2(a1.y,a1.x), theta2 = atan2(a2.y,a2.x);
   real_t theta = theta1;
   if (theta1 >= theta2)
      theta2 += 2*OFELI_PI;
   for (size_t i=0; i<nb; i++) {
      theta += (theta2-theta1)/nb;
      if (i<nb-1) {
         ln.n2 = _nb_vertices + 1;
         v.x = cc.x + r*cos(theta);
         v.y = cc.y + r*sin(theta);
         v.h = h;
         v.code = dc;
         _v.push_back(v);
      }
      if (i>0)
         ln.n1 = _nb_vertices;
   }
   _l.push_back(ln);
   _nb_vertices = _v.size();
   _nb_lines = _l.size();
}


void Domain::insertContour(const vector<size_t>& c)
{
   size_t i, i1, i2, n1, n2;
   size_t nb = c.size();
   if (nb>_nb_lines)
      throw OFELIException("In Domain::insertContour(vector<size_t>): "
                           "Number of lines is larger than number of defined lines.");
   Cont cc;
   cc.nb = nb;
   for (i=0; i<nb; i++)
      cc.line.push_back(c[i]);
   for (i=0; i<nb-1; i++) {
      n1 = cc.line[i]; n2 = cc.line[i+1];
      i1 = _l[n1-1].n2; i2 = _l[n2-1].n1;
      cc.orientation.push_back(1);
      if (i1 != i2) {
         i1 = _l[n1-1].n1; i2 = _l[n2-1].n1;
         cc.orientation[i] = -1;
         if (i1 != i2)
            throw OFELIException("In Domain::insertContour(vector<size_t>): Lines "+to_string(n1)+" and "+to_string(n2) +
                                 " are not correctly connected.\nContour cancelled.");
      }
   }
   n1 = cc.line[0]; n2 = cc.line[nb-1];
   if (_l[n1-1].n1 != _l[n2-1].n2)
      throw OFELIException("In Domain::insertContour(vector<size_t>): Contour is not closed.");
   _c.push_back(cc);
   _nb_sub_domains++;
   _nb_contours++;
}


int Domain::getContour()
{
   size_t nb = _ff->getI("Number of lines: ");
   if (nb > _nb_lines)
      throw OFELIException("In Domain::getContour(): Number of contour lines is larger than number of lines.");
   Cont cc;
   cc.nb = nb;
   cout << "Give list of lines in direct order.\n";
   for (size_t i=0; i<nb; i++) {
      size_t n = _ff->getI("Line Label: ");
      if (n > _nb_lines)
         throw OFELIException("In Domain::getContour(): Label of line is larger than number of lines.");
      cc.line.push_back(n);
   }

   for (size_t i=0; i<nb-1; i++) {
      size_t n1 = cc.line[i], n2 = cc.line[i+1];
      size_t i1 = _l[n1-1].n2, i2 = _l[n2-1].n1;
      if (i1 == i2)
         cc.orientation.push_back(1);
      else {
         cc.orientation.push_back(-1);
         i1 = _l[n1-1].n1; i2 = _l[n2-1].n1;
         if (i1 != i2)
            throw OFELIException("In Domain::getContour(): Lines "+to_string(n1)+" and "+to_string(n2)+ 
                                 " are not correctly connected.\nContour cancelled.");
      }
   }

   if (_l[cc.line[0]-1].n1 != _l[cc.line[nb-1]-1].n2)
      throw OFELIException("In Domain::getContour(): Contour is not closed.");
   _c.push_back(cc);
   _nb_contours++;

   Sd sd;
   sd.code = 1;
   sd.contour = _nb_contours;
   sd.line = 1;
   sd.type = 1;
   sd.orient = 2;
   _sd.push_back(sd);
   _nb_sub_domains++;
   return 0;
}


void Domain::insertHole(const vector<size_t>& h)
{
   size_t nb=h.size();
   Cont hh;
   hh.nb = nb;
   for (size_t i=0; i<nb; i++)
      hh.line.push_back(h[i]);
   for (size_t i=0; i<nb-1; i++) {
      size_t n1=hh.line[i], n2=hh.line[i+1];
      size_t i1=_l[n1-1].n2, i2=_l[n2-1].n1;
      if (i1 != i2)
         throw OFELIException("In Domain::insertHole(vector<size_t>): Lines "+to_string(n1) +
                              " and "+to_string(n2)+" are not correctly connected.\nHole cancelled.");
   }
   _h.push_back(hh);
   _nb_holes++;
}


int Domain::getHole()
{
   size_t nb = _ff->getI("Nb. of lines: ");
   Cont hh;
   hh.nb = nb;
   hh.line.resize(nb);
   cout << "Please give list of lines in direct order.\n";
   for (size_t i=0; i<size_t(nb); i++)
      hh.line[i] = _ff->getI("Line Label: ");

   for (size_t i=0; i<size_t(nb-1); i++) {
      size_t n1=hh.line[i], n2=hh.line[i+1];
      size_t i1=_l[n1-1].n2, i2=_l[n2-1].n1;
      if (i1 != i2)
         throw OFELIException("In Domain::getHole(): Lines "+to_string(n1)+" and " +
                              to_string(n2) + " are not correctly connected.");
   }
   _h.push_back(hh);
   _nb_holes++;
   return 0;
}


void Domain::insertSubDomain(size_t ln,
                             int    orient,
                             int    code)
{
   if (ln==0 || ln>_nb_lines)
      throw OFELIException("In Domain::insertSubDomain(size_t,int,int): Illegal line label");
   if (orient!=1 && orient!=-1)
      throw OFELIException("In Domain::insertSubDomain(size_t,int,int): Orientation "
                           "must be equal to 1 or -1");
   Sd sd;
   sd.line = int(ln);
   sd.orient = orient;
   sd.code = code;
   _sd.push_back(sd);
   _nb_sub_domains++;
}


void Domain::insertSubDomain(size_t n,
                             int    code)
{
   size_t k=0;
   for (size_t i=0; i<_nb_contours; i++) {
      _c[i].first_line = k + 1;
      for (size_t j=0; j<_c[i].nb; j++) {
         size_t n = _c[i].line[j];
         for (size_t l=1; l<_l[n-1].nb; l++)
            k++;
         k++;
      }
   }
   if (_sub_domain+1 > _nb_sub_domains)
      throw OFELIException("In Domain::insertSubDomain(size_t,int): Number of read "
                           "subdomains is larger than number of contours");
   Sd sd;
   sd.contour = int(n);
   sd.code = code;
   _sd.push_back(sd);
   _sub_domain++;
}


int Domain::getSubDomain()
{
   Sd sd;
   if (_sub_domain+1 > _nb_sub_domains)
      throw OFELIException("In Domain::getSubDomain(): Number of read subdomains "
                           "is larger than number of contours");
   sd.contour = _ff->getI("Contour label: ");
   int code = _ff->getI("Material code for subdomain: ");
   sd.code = code;
   _sd.push_back(sd);
   _sub_domain++;
   return 0;
}


int Domain::Rectangle()
{
   size_t n1, n2, n3, n4, m1, m2, m3;
   vector<int> cd;

   real_t ax = _ff->getD("x min: ");
   real_t ay = _ff->getD("y min: ");
   real_t bx = _ff->getD("x max: ");
   real_t by = _ff->getD("y max: ");
   real_t lx = bx - ax, ly = by - ay;
   int type = 1;
   size_t ne1 = _ff->getI("Subdivision in the x-direction: ");
   size_t ne2 = _ff->getI("Subdivision in the y-direction: ");

   int c1 = _ff->getI("Node code for y=ymin: ");
   int c2 = _ff->getI("Node code for x=xmax: ");
   int c3 = _ff->getI("Node code for y=ymax: ");
   int c4 = _ff->getI("Node code for x=xmin: ");
   int cs1 = _ff->getI("Side code for y=ymin: ");
   if (cs1 && c1) {
      cout << "Error in data: You cannot give nonzeo codes for both nodes and sides." << endl;
      c1 = _ff->getI("Node code for y=ymin: ");
      cs1 = _ff->getI("Side code for y=ymin: ");
   }
   int cs2 = _ff->getI("Side code for x=xmax: ");
   if (cs2 && c2) {
      cout << "Error in data: You cannot give nonzeo codes for both nodes and sides." << endl;
      c2 = _ff->getI("Node code for x=xmax: ");
      cs2 = _ff->getI("Side code for x=xmax: ");
   }
   int cs3 = _ff->getI("Side code for y=ymax: ");
   if (cs3 && c3) {
      cout << "Error in data: You cannot give nonzeo codes for both nodes and sides." << endl;
      c3 = _ff->getI("Node code for y=ymax: ");
      cs3 = _ff->getI("Side code for y=ymac: ");
   }
   int cs4 = _ff->getI("Side code for x=xmin: ");
   if (cs4 && c4) {
      cout << "Error in data: You cannot give nonzeo codes for both nodes and sides." << endl;
      c4 = _ff->getI("Node code for x=xmin: ");
      cs4 = _ff->getI("Side code for x=xmin: ");
   }

   _theMesh = new Mesh;
   _theMesh->setDim(2);
   _nb_dof = 1;

// Nodes
   real_t hx, hy = ly/real_t(ne2);
   int c;
   Point<real_t> xx;
   real_t y = ay;
   size_t first_dof = 1, k = 0;
   for (size_t i=0; i<=ne2; i++) {
      hx = lx/real_t(ne1);
      real_t x = ax;

      for (size_t j=0; j<=ne1; j++) {
         the_node = new Node(++k,Point<real_t>(x,y,0));
         The_node.setNbDOF(_nb_dof);
         c = setCode(ne1,ne2,i,j,c1,c2,c3,c4);
         dof_code(c,cd);
         The_node.setCode(&cd[0]);
         The_node.setDOF(first_dof,_nb_dof);
         _theMesh->Add(the_node);
         x += hx;
      }
      y += hy;
   }

// Elements
   k = 0;
   size_t el=0;
   if (type==0) {
      for (size_t i=0; i<ne2; i++) {
         for (size_t j=0; j<ne1; j++) {
            n1 = ++k; n2 = n1 + 1; n3 = n2 + ne1 + 1; n4 = n3 - 1;
            the_element = new Element(el+1,QUADRILATERAL,1);
            The_element.Add((*_theMesh)[n1]);
            The_element.Add((*_theMesh)[n2]);
            The_element.Add((*_theMesh)[n3]);
            The_element.Add((*_theMesh)[n4]);
            _theMesh->Add(the_element);
            el++;
         }
         k++;
      }
   }
   else {
      for (size_t i=0; i<ne2; i++) {
         for (size_t j=0; j<ne1; j++) {
            n1 = ++k; n2 = n1 + 1; n3 = n2 + ne1 + 1;
            the_element = new Element(el+1,TRIANGLE,1);
            The_element.Add((*_theMesh)[n1]);
            The_element.Add((*_theMesh)[n2]);
            The_element.Add((*_theMesh)[n3]);
            _theMesh->Add(the_element);
            el++; m1 = n3; m2 = n3 - 1; m3 = n1;
            the_element = new Element(el+1,TRIANGLE,1);
            The_element.Add((*_theMesh)[m1]);
            The_element.Add((*_theMesh)[m2]);
            The_element.Add((*_theMesh)[m3]);
            _theMesh->Add(the_element);
            el++;
         }
         k++;
      }
   }

// Sides
   k = 0; m1 = 1;
   for (size_t i=0; i<ne1; i++) {
      m2 = m1 + 1;
      if (cs1!=0) {
         the_side = new Side(++k,LINE);
         The_side.Add((*_theMesh)[m1]);
         The_side.Add((*_theMesh)[m2]);
         dof_code(cs1,cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
         _theMesh->Add(the_side);
      }
      m1 = m2;
   }
   for (size_t i=0; i<ne2; i++) {
      m2 = m1 + ne1 + 1;
      if (cs2!=0) {
         the_side = new Side(++k,LINE);
         The_side.Add((*_theMesh)[m1]);
         The_side.Add((*_theMesh)[m2]);
         dof_code(cs2,cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
         _theMesh->Add(the_side);
      }
      m1 = m2;
   }
   for (size_t i=0; i<ne1; i++) {
      m2 = m1 - 1;
      if (cs3!=0) {
         the_side = new Side(++k,LINE);
         The_side.Add((*_theMesh)[m1]);
         The_side.Add((*_theMesh)[m2]);
         dof_code(cs3,cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
         _theMesh->Add(the_side);
      }
      m1 = m2;
   }
   for (size_t i=0; i<ne2; i++) {
      m2 = m1 - ne1 - 1;
      if (cs4!=0) {
         the_side = new Side(++k,LINE);
         The_side.Add((*_theMesh)[m1]);
         The_side.Add((*_theMesh)[m2]);
         dof_code(cs4,cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
         _theMesh->Add(the_side);
      }
      m1 = m2;
   }
   return 0;
}


int Domain::Rectangle(real_t* x,
                      size_t  n1, 
                      size_t  n2,
                      real_t  r,
                      int     c1,
                      int     c2,
                      int     c3,
                      int     c4)
{
   size_t m1, m2, m3, nn1, nn2, nn3;
   vector<int> cd;
   int cs1=0, cs2=0, cs3=0, cs4=0;
   if (c1<0)
      cs1 = -c1, c1 = 0;
   if (c2<0)
      cs2 = -c2, c2 = 0;
   if (c3<0)
      cs3 = -c3, c3 = 0;
   if (c4<0)
      cs4 = -c4, c4 = 0;

// Nodes
   _theMesh = new Mesh;
   _theMesh->setDim(2);
   size_t k = 0;
   real_t lx=x[2]-x[0], ly=x[3]-x[1];
   real_t hx, hy=ly/real_t(n2);
   if (r!=1 && r>0)
      hy = ly*(r-1.)/(pow(r,real_t(n2))-1.);
   real_t yy=x[1], cs=1.75;
   size_t first_dof = 1;
   _nb_dof = 1;
   for (size_t i=0; i<=n2; i++) {
      hx = lx/real_t(n1);
      if (r != 1)
         hx = lx*(r-1.)/(pow(r,int(n1))-1.);
      real_t xx=x[0];
      if (r < 0)
         yy = x[1] + 0.5*ly*(1.+tanh(cs*((2.*i)/n1-1.0))/tanh(cs));

      for (size_t j=0; j<=n1; j++) {
         if (r<0)
            xx = x[0] + 0.5*lx*(1.+tanh(cs*((2.*j)/n2-1.0))/tanh(cs));
         the_node = new Node(++k,Point<real_t>(xx,yy,0));
         The_node.setNbDOF(_nb_dof);
         dof_code(setCode(n1,n2,i,j,c1,c2,c3,c4),cd);
         The_node.setCode(&cd[0]);
         The_node.setDOF(first_dof,_nb_dof);
         _theMesh->Add(the_node);
         xx += hx; hx *= r;
      }
      yy += hy; hy *= r;
   }

// Elements
   k = 0;
   size_t el = 0;
   for (size_t i=0; i<n2; i++) {
      for (size_t j=0; j<n1; j++) {
         nn1 = ++k; nn2 = nn1 + 1; nn3 = nn2 + n1 + 1;
         the_element = new Element(el+1,TRIANGLE,1);
         The_element.Add((*_theMesh)[nn1]);
         The_element.Add((*_theMesh)[nn2]);
         The_element.Add((*_theMesh)[nn3]);
         _theMesh->Add(the_element);
         el++; m1 = nn3; m2 = nn3 - 1; m3 = nn1;
         the_element = new Element(el+1,TRIANGLE,1);
         The_element.Add((*_theMesh)[m1]);
         The_element.Add((*_theMesh)[m2]);
         The_element.Add((*_theMesh)[m3]);
         _theMesh->Add(the_element);
         el++;
      }
      k++;
   }

// Sides
   k = 0;
   m1 = 1;
   for (size_t i=0; i<n1; i++) {
      m2 = m1 + 1;
      if (cs1!=0) {
         the_side = new Side(++k,LINE);
         The_side.Add((*_theMesh)[m1]);
         The_side.Add((*_theMesh)[m2]);
         dof_code(cs1,cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
         _theMesh->Add(the_side);
      }
      m1 = m2;
   }
   for (size_t i=0; i<n2; i++) {
      m2 = m1 + n1 + 1;
      if (cs2!=0) {
         the_side = new Side(++k,LINE);
         The_side.Add((*_theMesh)[m1]);
         The_side.Add((*_theMesh)[m2]);
         dof_code(cs2,cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
         _theMesh->Add(the_side);
      }
      m1 = m2;
   }
   for (size_t i=0; i<n1; i++) {
      m2 = m1 - 1;
      if (cs3!=0) {
         the_side = new Side(++k,LINE);
         The_side.Add((*_theMesh)[m1]);
         The_side.Add((*_theMesh)[m2]);
         dof_code(cs3,cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
         _theMesh->Add(the_side);
      }
      m1 = m2;
   }
   for (size_t i=0; i<n2; i++) {
      m2 = m1 - n1 - 1;
      if (cs4!=0) {
         the_side = new Side(++k,LINE);
         The_side.Add((*_theMesh)[m1]);
         The_side.Add((*_theMesh)[m2]);
         dof_code(cs4,cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
        _theMesh->Add(the_side);
      }
      m1 = m2;
   }
   return 0;
}


int Domain::Rectangle(real_t* x,
                      size_t  n1, 
                      size_t  n2,
                      real_t  r,
                      int*    c,
                      int*    cv)
{
   size_t m1, m2, m3, nn1, nn2, nn3;
   vector<int> cd;
   int cs[4] = { 0, 0, 0, 0 };
   for (size_t i=0; i<4; ++i)
      if (c[i]<0)
         cs[i] = -c[i], c[i] = 0;

// Nodes
   _theMesh = new Mesh;
   _theMesh->setDim(2);
   size_t k=0;
   real_t lx=x[2]-x[0], ly=x[3]-x[1];
   real_t hx, hy=ly/real_t(n2);
   if (r!=1 && r>0)
      hy = ly*(r-1.)/(pow(r,real_t(n2))-1.);
   real_t yy=x[1], css=1.75;
   size_t first_dof = 1;
   _nb_dof = 1;
   for (size_t i=0; i<=n2; ++i) {
      hx = lx/real_t(n1);
      if (r != 1)
         hx = lx*(r-1.)/(pow(r,int(n1))-1.);
      real_t xx=x[0];
      if (r < 0)
         yy = x[1] + 0.5*ly*(1.+tanh(css*((2.*i)/n1-1.0))/tanh(css));

      for (size_t j=0; j<=n1; ++j) {
         if (r<0)
            xx = x[0] + 0.5*lx*(1.+tanh(css*((2.*j)/n2-1.0))/tanh(css));
         the_node = new Node(++k,Point<real_t>(xx,yy,0));
         The_node.setNbDOF(_nb_dof);
         dof_code(setCode(n1,n2,i,j,c[0],c[1],c[2],c[3]),cd);
         The_node.setCode(&cd[0]);
         The_node.setDOF(first_dof,_nb_dof);
         _theMesh->Add(the_node);
         xx += hx; hx *= r;
      }
      yy += hy; hy *= r;
   }
   (*_theMesh)[1]->setCode(1,cv[0]);
   (*_theMesh)[n1+1]->setCode(1,cv[1]);
   (*_theMesh)[(n1+1)*(n2+1)]->setCode(1,cv[2]);
   (*_theMesh)[n1*n2+n2+1]->setCode(1,cv[3]);

// Elements
   k = 0;
   size_t el=0;
   for (size_t i=0; i<n2; i++) {
      for (size_t j=0; j<n1; j++) {
         nn1 = ++k; nn2 = nn1 + 1; nn3 = nn2 + n1 + 1;
         the_element = new Element(el+1,TRIANGLE,1);
         The_element.Add((*_theMesh)[nn1]);
         The_element.Add((*_theMesh)[nn2]);
         The_element.Add((*_theMesh)[nn3]);
         _theMesh->Add(the_element);
         el++; m1 = nn3; m2 = nn3 - 1; m3 = nn1;
         the_element = new Element(el+1,TRIANGLE,1);
         The_element.Add((*_theMesh)[m1]);
         The_element.Add((*_theMesh)[m2]);
         The_element.Add((*_theMesh)[m3]);
         _theMesh->Add(the_element);
         el++;
      }
      k++;
   }

// Sides
   k = 0;
   m1 = 1;
   for (size_t i=0; i<n1; i++) {
      m2 = m1 + 1;
      if (cs[0]!=0) {
         the_side = new Side(++k,LINE);
         The_side.Add((*_theMesh)[m1]);
         The_side.Add((*_theMesh)[m2]);
         dof_code(cs[0],cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
         _theMesh->Add(the_side);
      }
      m1 = m2;
   }
   for (size_t i=0; i<n2; i++) {
      m2 = m1 + n1 + 1;
      if (cs[1]!=0) {
         the_side = new Side(++k,LINE);
         The_side.Add((*_theMesh)[m1]);
         The_side.Add((*_theMesh)[m2]);
         dof_code(cs[1],cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
         _theMesh->Add(the_side);
      }
      m1 = m2;
   }
   for (size_t i=0; i<n1; i++) {
      m2 = m1 - 1;
      if (cs[2]!=0) {
         the_side = new Side(++k,LINE);
         The_side.Add((*_theMesh)[m1]);
         The_side.Add((*_theMesh)[m2]);
         dof_code(cs[2],cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
         _theMesh->Add(the_side);
      }
      m1 = m2;
   }
   for (size_t i=0; i<n2; i++) {
      m2 = m1 - n1 - 1;
      if (cs[3]!=0) {
         the_side = new Side(++k,LINE);
         The_side.Add((*_theMesh)[m1]);
         The_side.Add((*_theMesh)[m2]);
         dof_code(cs[3],cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
        _theMesh->Add(the_side);
      }
      m1 = m2;
   }
   return 0;
}


int Domain::setCode(size_t ne1,
                    size_t ne2,
                    size_t i,
                    size_t j,
                    int    c1,
                    int    c2,
                    int    c3,
                    int    c4)
{
   int code;
   if (j==0)
      code = c4;
   else if (j==ne1)
      code = c2;
   else if (i==0)
      code = c1;
   else if (i==ne2)
      code = c3;
   else
      code = 0;
   return code;
}


int Domain::Disk()
{
   size_t      n, m, nbn, nbe, nbs;
   vector<int> dcc(MAX_NBDOF_NODE);
   Vertex c;
   real_t R = _ff->getI("Disk radius: ");
   c.x = _ff->getD("x-Coordinates of center: ");
   c.y = _ff->getD("y-Coordinates of center: ");

   size_t Nc = _ff->getI("Nb. of subdivisions of the circumference: ");
   size_t Nl = _ff->getI("Nb. of layers: ");
   int bs = _ff->getI("Store boundary sides (0/1) ? ");
   _nb_dof = _ff->getI("Nb. of degrees of freedom: ");
   int dc = _ff->getI("Dirichlet code: ");
   int nc = _ff->getI("Neumann code: ");
   string file = _ff->getS("Output mesh file name: ");

   ofstream fm(file.c_str());
   fm << "<?xml version=\"1.0\"?>\n<OFELI_File>" << endl;
   fm << "<info>" << endl;
   fm << "   <title></title>" << endl;
   fm << "   <date></date>" << endl;
   fm << "   <author></author>" << endl;
   fm << "</info>" << endl;
   fm << "<Mesh dim=\"" << _dim << "\" nb_dof=\"" << _nb_dof << "\">" << endl;

   nbn = Nc*Nl + 1;
   nbe = 2*Nc*(Nl-1) + Nc;
   nbs = Nc;
   if (!bs)
      nbs = 0;
   vector<Vertex> x(nbn);
   _el.resize(nbe);

   x[0] = c;
   x[0].code = 0;
   m = n = 1;
   for (size_t i=0; i<Nc; i++) {
      real_t theta = 2*i*OFELI_PI/real_t(Nc);
      _el[m-1].n1 = 1; _el[m-1].n2 = n+1; _el[m-1].n3 = n+1+Nl;
      if (i==Nc-1)
         _el[m-1].n3 = 2;
      size_t k = _el[m-1].n3;
      m++;
      for (size_t j=1; j<=Nl; j++) {
         real_t r = j*R/real_t(Nl);
         x[n].x = c.x + r*cos(theta);
         x[n].y = c.y + r*sin(theta);
         x[n].code = 0;
         if (j==Nl)
            x[n].code = dc;
         n++;
         if (j!=Nl) {
            _el[m-1].n1 = n; _el[m-1].n2 = n+1; _el[m-1].n3 = n+Nl;
           if (i==Nc-1)
              _el[m-1].n3 = k;
           m++;
           _el[m-1].n1 = n+1; _el[m-1].n2 = n+1+Nl; _el[m-1].n3 = n+Nl;
           if (i==Nc-1) {
              _el[m-1].n2 = k+1; _el[m-1].n3 = k++;
           }
           m++;
         }
      }
      Ln ln;
      ln.n1 = _el[m-2].n1;
      ln.n2 = _el[m-2].n2;
      ln.Ncode = nc;
      _l.push_back(ln);
   }

   int k=0;
   fm << "  <Nodes>" << endl;
   for (n=0; n<nbn; n++) {
      fm << "  " << setprecision(12) << setw(22) << x[n].x << setw(18) << x[n].y;
      dof_code(x[n].code,dcc);
      fm << endl;
      int sign = 1;
      if (dcc[0]<0)
         sign = -1;
      m = 0;
      for (size_t j=1; j<=_nb_dof; j++)
         m += abs(dcc[j-1])*size_t(pow(10.,real_t(_nb_dof-j)));
      m *= sign;
      fm << setw(6) << m << "   ";
      k++;
      if (k%5==0)
         fm << endl;
   }
   fm << "  </Nodes>" << endl;

   fm << "   <Elements shape=\"triangle\"  nodes=\"3\">" << endl;
   for (m=0; m<nbe; m++) {
      fm << setw(8) << _el[m].n1 << setw(8) << _el[m].n2 << setw(8) << _el[m].n3
         << setw(8) << 1 << "   ";
      k++;
      if (k%5==0)
         fm << endl;
   }
   if (k%5!=0)
      fm << endl;
   fm << "   </Elements>" << endl;

   if (bs) {
      k = 0;
      fm << "   <Sides shape=\"line\"  nodes=\"2\">" << endl;
      for (n=0; n<nbs; n++) {
         for (size_t i=0; i<_nb_dof; i++)
            fm << setw(4) << dcc[i];
         fm << endl;
         dof_code(_l[n].Dcode,dcc);
         fm << setw(8) << _l[n].n1 << setw(6) << _l[n].n2;
         m = 0;
         for (size_t j=1; j<=_nb_dof; j++)
            m += The_side.getCode(j)*size_t(pow(10.,real_t(_nb_dof-j)));
         fm << setw(6) << m << "   ";
         k++;
         if (k%5==0)
            fm << endl;
      }
      if (k%5!=0)
         fm << endl;
      fm << "   </Sides>" << endl;
   }
   fm << "</Mesh>\n</OFELI_File>" << endl;
   fm.close();

   cout << "Number of generated nodes:     " << nbn << endl;
   cout << "Number of generated triangles: " << nbe << endl;
   if (nbs)
      cout << "Number of generated sides:     " << nbs << endl;
   cout << "Mesh stored in file " << file << endl;
   return 0;
}


void Domain::list()
{
   cout.setf(ios::right|ios::scientific);
   if (!_nb_vertices)
      cout << "NO VERTICES.\n";
   else {
      cout << "\n****** List of vertices *****\n\n";
      cout << "   X               Y           Code\n";
      for (size_t i=0; i<_nb_vertices; i++) {
         cout << setw(3) << "    " << setprecision(5) << setw(10) << _v[i].x;
         cout << "    " << _v[i].y << setw(6) << _v[i].code << endl;
     }
     cout << endl;
   }

   if (!_nb_lines)
      cout << "NO LINES.\n";
   else {
      cout << "\n****** List of lines ******\n\n";
      cout << "   END 1  END 2   DCode   NCode\n";
      for (size_t i=0; i<_nb_lines; i++)
         cout << setw(7) << _l[i].n1 << setw(7) << _l[i].n2
              << setw(8) << _l[i].Dcode << setw(7) << _l[i].Ncode << endl;
      cout << endl;
   }
   if (!_nb_contours || _ret_cont)
      cout << "NO CONTOUR.\n";
   else {
      for (size_t i=0; i<_nb_contours; i++) {
         cout << "\n****** External Contour ******\n\n";
         cout << "NB. OF LINES     LIST OF LINES\n";
         cout << setw(7) << _c[i].nb << "          ";
         for (size_t j=0; j<_c[i].nb; j++)
            cout << setw(5) << _c[i].line[j];
         cout << endl;
      }
   }

   if (!_nb_holes)
      cout << "NO HOLES.\n";
   else {
      cout << "\n****** List of holes ******\n\n";
      cout << "LABEL   NB. OF LINES      LIST OF LINES\n";
      for (size_t i=0; i<_nb_holes; i++) {
         cout << setw(5) << _h[i].nb << "     ";
         for (size_t j=0; j<_h[i].nb; j++)
            cout << setw(5) << _h[i].line[j];
         cout << endl;
      }
   }
}


void Domain::deleteLine()
{
   size_t j = _ff->getI("Label of line to delete: ");
   int k = remove(j,_nb_lines,_l_label);
   if (k==-1)
      cout << "WARNING: This line was not present." << endl;
   else {
      for (size_t i=k; i<_nb_lines; i++)
         _l[i] = _l[i+1];
      cout << "LINE " << j << " REMOVED." << endl;
   }
}


void Domain::deleteVertex()
{
   size_t j = _ff->getI("Label of point to delete: ");
   int k = remove(j,_nb_vertices,_v_label);
   if (k==-1)
      cout << "WARNING: This point was not present.\n";
   else {
      for (size_t i=size_t(k); i<_nb_vertices; i++)
         _v[i] = _v[i+1];
      cout << "POINT " << j << " REMOVED.\n";
   }
}


int Domain::saveAsEasymesh()
{
   size_t  i, j, n, i1, nb, nbb, kl, k, kl0=0;
   int     dc, nc, c;
   real_t h;

   if (!_c[0].nb)
      throw OFELIException("In Domain::saveAsEasyMesh: No Contour is defined.");
   string mfile = _ff->getS("File Name (.d): ");
   string emfile = mfile + ".d";
   ofstream mf(emfile.c_str());

// Count nb. of points
   mf << "# POINTS #" << endl;
   nbb = 0;
   for (i=0; i<_c[0].nb; i++)
      nbb += _l[_c[0].line[i]-1].nb - 1;
   for (n=0; n<_nb_holes; n++)
      for (i=0; i<_h[n].nb; i++)
         nbb += _l[_h[n].line[i]-1].nb - 1;
   _nb_lines += nbb;
   _nb_vertices += nbb;
   mf << setw(5) << _nb_vertices << "# Nb. of points #" << endl;

// Output Points of the external contour
   mf <<  "# Points on External Contour #" << endl;
   _ln.resize(_c[0].nb+_nb_holes);
   kl = 0;
   for (i=0; i<_c[0].nb; i++) {
      n = _c[0].line[i];
      i1 = _l[n-1].n1;
      c = _v[i1-1].code;
      h = _v[i1-1].h;
      nb = _l[n-1].nb;
      dc = _l[n-1].Dcode;
      nc = _l[n-1].Ncode;
      _ln[kl].i = kl+1;
      if (i==0)
         kl0 = kl+1;
      _ln[kl].dc = dc; _ln[kl].nc = nc;
      mf << setw(5) << kl << ": " << _l[n-1].node[0].x << "  " << _l[n-1].node[0].y << "  "
         << h << "  " << c << endl;
      for (j=1; j<nb; j++) {
         mf << kl+1 << ":  " << _l[n-1].node[j].x << "  " << _l[n-1].node[j].y << "  "
            << h << "  " << dc << endl;
         _ln[kl+1].i = _ln[kl].j = kl + 2;
         _ln[kl+1].dc = dc; _ln[kl+1].nc = nc;
         kl++;
      }
      _ln[kl].j = kl + 2;
      kl++;
   }
   _ln[kl-1].j = kl0;

// Output Points of the holes
   for (k=0; k<_nb_holes; k++) {
      mf << "\n# Points on Hole Nb. " << k+1 << " #" << endl;
      for (i=0; i<_h[k].nb; i++) {
         n = _h[k].line[i];
         i1 = _l[n-1].n1;
         c = _v[i1-1].code;
         h = _v[i1-1].h;
         nb = _l[n-1].nb;
         dc = _l[n-1].Dcode; nc = _l[n-1].Ncode;
         if (i==0)
            kl0 = kl + 1;
         _ln[kl].i = kl + 1;
         _ln[kl].dc = dc; _ln[kl].nc = nc;
         mf << kl << ": " << _l[n-1].node[0].x << "  " << _l[n-1].node[0].y << "  "
            << h << "  " << dc << endl;
         for (j=1; j<nb; j++) {
            mf << kl+1 << ": " << _l[n-1].node[j].x << "  " << _l[n-1].node[j].y << h << "  "
               << c << endl;
            _ln[kl+1].i = _ln[kl].j = kl+2;
            _ln[kl+1].dc = dc;
            _ln[kl+1].nc = nc;
            kl++;
         }
         _ln[kl].j = kl+2;
         kl++;
      }
      _ln[kl-1].j = kl0;
   }

// Output Lines
   mf <<  "\n\n#  LINES #" << endl;
   mf << kl <<  "%5d  # Number of Lines #" << endl;
   for (i=0; i<kl; i++)
      mf << i << ": " << _ln[i].i-1 << "  " << _ln[i].j-1 << "  " << _ln[i].dc << "  "
         << _ln[i].nc << endl;
   cout << "Mesh stored in file: " << emfile << endl;
   cout << "You can generate mesh by typing\n easymesh " << mfile << endl;
   return 0;
}


int Domain::saveAsBamg()
{
   if (_c[0].nb == 0)
      throw OFELIException("In Domain::saveAsBamg: No Contour is defined.");
   string mfile = _ff->getS("File Name (.geo): ");
   string emfile = mfile + ".geo";
   ofstream mf(emfile.c_str());

// Count nb. of points
   mf << "MeshVersionFormatted 0" << endl;
   mf << "Dimension  2" << endl;
   size_t nbb = 0;
   for (size_t i=0; i<_c[0].nb; i++)
      nbb += _l[_c[0].line[i]-1].nb - 1;
   for (size_t n=0; n<_nb_holes; n++)
      for (size_t i=0; i<_h[n].nb; i++)
         nbb += _l[_h[n].line[i]-1].nb - 1;
   _nb_lines += nbb;
   _nb_vertices += nbb;
   mf << "\nVertices  " << setw(6) << _nb_vertices << endl;

// Output Points of the external contour
   size_t kl=0, kl0=0;
   for (size_t i=0; i<_c[0].nb; i++) {
      size_t n = _c[0].line[i];
      size_t i1 = _l[n-1].n1;
      int c = _v[i1-1].code;
      size_t nb = _l[n-1].nb;
      int dc=_l[n-1].Dcode, nc=_l[n-1].Ncode;
      _ln[kl].i = kl + 1;
      if (i==0)
         kl0 = kl + 1;
      _ln[kl].dc = dc;
	  _ln[kl].nc = nc;
      mf << _l[n-1].node[0].x << "  " << _l[n-1].node[0].y << setw(5) << c << endl;
      for (size_t j=1; j<nb; j++) {
         mf << _l[n-1].node[j].x << "  " << _l[n-1].node[j].y << setw(5) << dc << endl;
         _ln[kl+1].i = _ln[kl].j = kl + 2;
         _ln[kl+1].dc = dc;
         _ln[kl+1].nc = nc;
         kl++;
      }
      _ln[kl].j = kl + 2;
      kl++;
   }
   _ln[kl-1].j = kl0;

// Output Points of the holes
   for (size_t k=0; k<_nb_holes; k++) {
      for (size_t i=0; i<_h[k].nb; i++) {
         size_t n = _h[k].line[i];
         size_t i1 = _l[n-1].n1;
         int c = _v[i1-1].code;
         size_t nb = _l[n-1].nb;
         int dc=_l[n-1].Dcode, nc=_l[n-1].Ncode;
         if (i==0)
            kl0 = kl + 1;
         _ln[kl].i = kl + 1;
         _ln[kl].dc = dc, _ln[kl].nc = nc;
         mf << _l[n-1].node[0].x << "  " << _l[n-1].node[0].y << setw(5) << dc << endl;
         for (size_t j=1; j<nb; j++) {
            mf << _l[n-1].node[j].x << "  " << _l[n-1].node[j].y << setw(5) << c << endl;
            _ln[kl+1].i = _ln[kl].j = kl + 2;
            _ln[kl+1].dc = dc, _ln[kl+1].nc = nc;
            kl++;
         }
         _ln[kl].j = kl + 2;
         kl++;
      }
      _ln[kl-1].j = kl0;
   }

// Output Lines
   mf << "\nEdges    " << kl << endl;
   for (size_t i=0; i<kl; i++)
      mf << setw(6) << _ln[i].i << setw(6) << _ln[i].j << setw(6) << _ln[i].dc
         << setw(5) << _ln[i].nc << endl;

// Output mesh sizes
   mf << "\nhVertices  ";
   kl = 0;
   for (size_t i=0; i<_c[0].nb; i++) {
      size_t n = _c[0].line[i];
      size_t i1 = _l[n-1].n1;
      real_t h = _v[i1-1].h;
      size_t nb = _l[n-1].nb;
      mf << "  " << h;
      for (size_t j=1; j<nb; j++)
         mf << "  " << h;
      mf << endl;
   }
   for (size_t k=0; k<_nb_holes; k++) {
      for (size_t i=0; i<_h[k].nb; i++) {
         size_t n = _h[k].line[i];
         size_t i1 = _l[n-1].n1;
         real_t h = _v[i1-1].h;
         size_t nb = _l[n-1].nb;
         mf << "   " << h;
         for (size_t j=1; j<nb; j++)
            mf << "   " << h;
      }
      mf << endl;
   }

   mf << "\nSubDomain    1" << endl;
   mf << "    2    1    1    1" << endl;
   mf.close();

   cout << "Mesh stored in file: " << emfile << endl;
   cout << "You can generate mesh by typing\n bamg -g " << emfile
        << " -o " << mfile << ".bamg" << endl << endl;
   return 0;
}


int Domain::saveAsTriangle()
{
   size_t  i, j, n, i1, nb, nbb, kl, k, kl0=0, label;
   int     dc, nc;
   real_t  h, area=0;
   string mfile, emfile;

   if (_c[0].nb == 0)
      throw OFELIException("In Domain::saveAsTriangle(): No Contour is defined.");
   mfile = _ff->getS("File Name (.poly): ");
   emfile = mfile + ".poly";
   ofstream mf(emfile.c_str());

// Count nb. of points
   nbb = 0;
   for (i=0; i<_c[0].nb; i++)
      nbb += _l[_c[0].line[i]-1].nb - 1;
   for (n=0; n<_nb_holes; n++)
      for (i=0; i<_h[n].nb; i++)
         nbb += _l[_h[n].line[i]-1].nb - 1;
   _nb_lines += nbb;
   label = _nb_vertices;
   _nb_vertices += nbb;
   mf << "# Generated by OFELI" << endl;
   mf << _nb_vertices << "  2   0   1" << endl;

// Output Points of the external contour
   kl = 0;
   h = _v[_l[_c[0].line[0]-1].n1].h;
   area = 0.5*h*h;
   for (i=0; i<_c[0].nb; i++) {
      n = _c[0].line[i];
      i1 = _l[n-1].n1;
      nb = _l[n-1].nb;
      h = _v[i1-1].h;
      area = std::min(area,0.5*h*h);
      dc = _l[n-1].Dcode;
      nc = _l[n-1].Ncode;
      _ln[kl].i = kl+1;
      if (i==0)
        kl0 = kl + 1;
      _ln[kl].dc = dc;
      _ln[kl].nc = nc;
      mf << setw(6) << 1 << "  ";
      mf << _l[n-1].node[0].x << "  " << _l[n-1].node[0].y << setw(6) << dc << endl;
      for (j=1; j<nb; j++) {
         mf << setw(6) << label++ << "  ";
         mf << _l[n-1].node[j].x << "  " << _l[n-1].node[j].y << setw(6) << dc << endl;
         _ln[kl+1].i = _ln[kl].j = kl + 2;
         _ln[kl+1].dc = dc;
         _ln[kl+1].nc = nc;
         kl++;
      }
      _ln[kl].j = kl + 2;
      kl++;
   }
   _ln[kl-1].j = kl0;

// Output Points of the holes
   for (k=0; k<_nb_holes; k++) {
      for (i=0; i<_h[k].nb; i++) {
         n = _h[k].line[i];
         i1 = _l[n-1].n1;
         h = _v[i1-1].h;
         area = std::min(area,0.5*h*h);
         nb = _l[n-1].nb;
         dc = _l[n-1].Dcode;
         nc = _l[n-1].Ncode;
         if (i==0)
            kl0 = kl + 1;
         _ln[kl].i = kl + 1;
         _ln[kl].dc = dc;
         _ln[kl].nc = nc;
         mf << setw(6) << label++;
         mf << _l[n-1].node[0].x << "  " << _l[n-1].node[0].y << setw(6) << dc << endl;
         for (j=1; j<nb; j++) {
            label++;
            mf << setw(6) << label++;
            mf << _l[n-1].node[j].x << "  " << _l[n-1].node[j].y << setw(6) << dc << endl;
            _ln[kl+1].i = _ln[kl].j = kl + 2;
            _ln[kl+1].dc = dc;
            _ln[kl+1].nc = nc;
            kl++;
         }
         _ln[kl].j = kl + 2;
         kl++;
      }
      _ln[kl-1].j = kl0;
   }

// Output Lines
   mf << kl << setw(6) << 2 << endl;
   for (i=0; i<kl; i++) {
      mf << setw(6) << i+1 << setw(6) << _ln[i].i << setw(6) << _ln[i].j;
      mf << setw(6) << _ln[i].dc << setw(6) << _ln[i].nc << endl;
   }

// Output holes
   mf << _nb_holes << endl;
   real_t xx, yy;
   for (k=0; k<_nb_holes; k++) {
      cout << "Give x and y coordinates of one point in hole " << k+1;
      cin >> xx >> yy;
      mf << k+1 << "  " << xx << "  " << yy << endl;
   }
   cout << "You can now generate mesh by the typical command:" << endl;
   cout << "triangle -Dpea" << area << "  " << mfile << endl << endl;
   cout << "WARNING: Check triangle output file names (mesh iteration numbers) !" << endl;
   return 0;
}


Domain &Domain::operator*=(real_t a)
{
   for (size_t i=0; i<_nb_vertices; i++) {
      _v[i].x *= a;
      _v[i].y *= a;
      _v[i].h *= a;
   }
   return *this;
}


void Domain::genGeo(string file)
{
   ofstream mf(file.c_str());

   mf << "MeshVersionFormatted 0" << endl;
   mf << "Dimension  2" << endl;
   mf << "AngleOfCornerBound 46" << endl;
   if (_nb_vertices > 0) {
      mf << "\nVertices   " << _nb_vertices << endl;
      for (size_t i=0; i<_nb_vertices; ++i)
         mf << _v[i].x << "  " << _v[i].y << "  " << _v[i].code << endl;
   }
   if (_nb_lines > 0) {
      mf << "\nEdges    " << _nb_lines << endl;
      for (size_t i=0; i<_nb_lines; ++i) {
         int code = _l[i].Dcode;
         if (_l[i].Ncode)
            code = -_l[i].Ncode;
         mf << setw(5) << _l[i].n1 << "  " << _l[i].n2 << "  " << code << endl;
      }
   }
   if (_nb_vertices > 0) {
      mf << "\nhVertices\n";
      for (size_t i=0; i<_nb_vertices; i++)
         mf << "  " << _v[i].h;
      mf << endl;
   }
   if (_nb_required_vertices > 0) {
      mf << "\nRequiredVertices" << setw(5) << _nb_required_vertices << endl;
      for (size_t i=0; i<_nb_required_vertices; i++)
        mf << "  " << _required_vertex[i];
      mf << endl;
   }
   if (_nb_required_edges > 0) {
      mf << "\nRequiredEdges" << setw(5) << _nb_required_edges << endl;
      for (size_t i=0; i<_nb_required_edges; i++)
         mf << "  " << _required_edge[i];
      mf << endl;
   }
   if (_nb_sub_domains > 0) {
      mf << "\nSubDomain " << setw(6) << _nb_sub_domains << endl;
      for (size_t i=0; i<_nb_sub_domains; i++)
         mf << "  " << 2 << "  " << _sd[i].line << "  " << _sd[i].orient << "  " << _sd[i].code << endl;
   }
   mf << "\nEnd" << endl;
}


void Domain::genMesh(string geo_file,
                     string bamg_file,
                     string mesh_file)
{
   genGeo(geo_file);
   bamg::Geometry Gh(geo_file.c_str());
   int nbvx=100000;
   bamg::Triangles Th(nbvx,Gh);
   Th.MakeQuadrangles(2);
   Th.ReNumberingTheTriangleBySubDomain();
   Th.Write(bamg_file.c_str(),bamg::Triangles::BDMesh);
   _theMesh = new Mesh;
   getBamg(bamg_file,*_theMesh,_nb_dof);
   removeUnusedNodes();
   Vect<size_t> nd_ref(_theMesh->getNbNodes());
   nd_ref = 0;
   MESH_SD {
      for (size_t i=1; i<=The_side.getNbNodes(); i++)
         if (The_side.getCode(1))
            nd_ref[The_side(i)->n()-1]++;
   }
   MESH_ND {
      if (nd_ref[node_label-1]>1)
         The_node.setCode(1,0);
   }
   _theMesh->put(mesh_file);
}


void Domain::removeUnusedNodes()
{
   Mesh ms(*_theMesh);
   vector<int> used(_theMesh->getNbNodes());
   vector<Node *> nd(_theMesh->getNbNodes());
   clear(used);
   MESH_EL {
      for (size_t i=1; i<=The_element.getNbNodes(); i++)
         used[The_element(i)->n()-1]++;
   }

   delete _theMesh;
   _theMesh = new Mesh;
   _theMesh->setDim(ms.getDim());
   size_t n=0;
   for (size_t i=1; i<=ms.getNbNodes(); ++i) {
      size_t in=ms[i]->n()-1;
      nd[in] = nullptr;
      if (used[in]>0) {
         nd[in] = ms[in+1];
         nd[in]->setLabel(++n);
         _theMesh->Add(nd[in]);
      }
   }
   int ne = 0;
   element_loop(&ms) {
      Element *el = new Element(++ne,The_element.getShape(),The_element.getCode());
      for (size_t i=1; i<=The_element.getNbNodes(); i++)
          el->Add(The_element(i));
      _theMesh->Add(el);
   }
   side_loop(&ms) {
      Side *sd = new Side(side_label,The_side.getShape());
      for (size_t i=1; i<=The_side.getNbNodes(); i++)
         sd->Add(The_side(i));
      sd->setNbDOF(The_side.getNbDOF());
      for (size_t i=1; i<=sd->getNbDOF(); i++)
         sd->setCode(i,The_side.getCode(i));
      _theMesh->Add(sd);
   }
}


void Domain::genMesh()
{
   string prefix = _domain_file.substr(0,_domain_file.rfind("."));
   _geo_file = prefix + ".geo";
   _bamg_file = prefix + ".bamg";
   genMesh(_geo_file,_bamg_file,prefix+".m");
}


void Domain::genMesh(const string& file)
{
   string prefix = file.substr(0,file.rfind("."));
   _geo_file = prefix + ".geo";
   _bamg_file = prefix + ".bamg";
   genMesh(_geo_file,_bamg_file,file);
}


void Domain::generateMesh()
{
   if (_dim==1)
      gm1();
   else if (_dim==2)
      gm2();
   else if (_dim==3)
      gm3();
}


void Domain::gm1()
{
}


void Domain::gm2()
{
   if (_c[0].nb == 0)
      throw OFELIException("In Domain::gm2(): No Contour is defined.");
   string prefix = _output_file.substr(0,_output_file.rfind("."));
   string geo_file=prefix+".geo", bamg_file=prefix+".bamg";
   ofstream mf(geo_file.c_str());

// Count nb. of points
   mf << "MeshVersionFormatted 0\nDimension  2" << endl;
   size_t nbb = 0;
   for (size_t is=0; is<_nb_contours; is++) {
      for (size_t i=0; i<_c[is].nb; i++)
         nbb += _l[_c[is].line[i]-1].nb;
   }
   for (size_t n=0; n<_nb_holes; n++)
      for (size_t i=0; i<_h[n].nb; i++)
         nbb += _l[_h[n].line[i]-1].nb;
   _nb_lines += nbb;
   _nb_vertices = nbb;
   mf << "\nVertices  " << setw(5) << _nb_vertices << endl;

// Output Points of the external contour
   size_t kl=0, kl0=0;
   _ln.resize(_nb_lines);
   for (size_t is=0; is<_nb_contours; is++) {
      for (size_t i=0; i<_c[is].nb; i++) {
         size_t n=_c[is].line[i];
         real_t dc=_l[n-1].Dcode, nc=_l[n-1].Ncode;
         _ln[kl].i = kl + 1;
         if (i==0)
            kl0 = kl + 1;
         _ln[kl].dc = dc; _ln[kl].nc = nc;
         mf << "  " << _l[n-1].node[0].x << "  " << _l[n-1].node[0].y << "  "
            << _v[_l[n-1].n1-1].code << endl;
         for (size_t j=1; j<_l[n-1].nb; j++) {
            mf << "  " << _l[n-1].node[j].x << "  " << _l[n-1].node[j].y << "  " << dc << endl;
            _ln[kl+1].i = _ln[kl].j = kl + 2;
            _ln[kl+1].dc = dc, _ln[kl+1].nc = nc;
            kl++;
         }
         _ln[kl].j = kl + 2;
         kl++;
      }
      _ln[kl-1].j = kl0;
   }

// Output Points of holes
   for (size_t k=0; k<_nb_holes; k++) {
      for (size_t i=0; i<_h[k].nb; i++) {
         size_t n = _h[k].line[i];
         real_t dc = _l[n-1].Dcode, nc = _l[n-1].Ncode;
         if (i==0)
            kl0 = kl + 1;
         _ln[kl].i = kl + 1;
         _ln[kl].dc = dc, _ln[kl].nc = nc;
         mf << "  " << _l[n-1].node[0].x << "  " << _l[n-1].node[0].y << "  " << dc << endl;
         for (size_t j=1; j<_l[n-1].nb; j++) {
            mf << "  " << _l[n-1].node[j].x << "  " << _l[n-1].node[j].y << "  "
               << _v[_l[n-1].n1-1].code << endl;
            _ln[kl+1].i = _ln[kl].j = kl + 2;
            _ln[kl+1].dc = dc, _ln[kl+1].nc = nc;
            kl++;
         }
         _ln[kl].j = kl + 2;
         kl++;
      }
      _ln[kl-1].j = kl0;
   }

// Output Lines
   mf << "\nEdges  " << setw(5) << kl << endl;
   for (size_t i=0; i<kl; i++)
      mf << setw(6) << _ln[i].i << setw(6) << _ln[i].j
         << setw(6) << _ln[i].dc << setw(6) << _ln[i].nc << endl;

// Output mesh sizes
   mf << "\nhVertices\n";
   kl = 0;
   for (size_t k=0; k<_nb_contours; k++) {
      for (size_t i=0; i<_c[k].nb; i++) {
         size_t n = _c[k].line[i];
         size_t i1 = _l[n-1].n1;
         mf << setprecision(6) << setw(10) << _v[i1-1].h;
         for (size_t j=1; j<_l[n-1].nb; j++)
            mf << setprecision(6) << setw(10) << _v[i1-1].h;
         mf << endl;
      }
   }
   for (size_t k=0; k<_nb_holes; k++) {
      for (size_t i=0; i<_h[k].nb; i++) {
         size_t n = _h[k].line[i];
         size_t i1 = _l[n-1].n1;
         mf << setprecision(6) << setw(10) << _v[i1-1].h;
         for (size_t j=1; j<_l[n-1].nb; j++)
            mf << setprecision(6) << setw(10) << _v[i1-1].h;
      }
      mf << endl;
   }

   mf << "\nSubDomain " << setw(6) << _nb_sub_domains << endl;
   for (size_t k=0; k<_nb_sub_domains; k++)
      mf << setw(6) << 2 << setw(6) << _c[k].first_line << setw(6) << 1
         << setw(6) << _sd[k].code << endl;
   mf << "\nEnd" << endl;

   bamg1(geo_file,bamg_file);
   _theMesh = new Mesh;
   getBamg(bamg_file,*_theMesh,_nb_dof);
   ::remove(geo_file.c_str());
   ::remove(bamg_file.c_str());
}


void Domain::gm3()
{
}


void Domain::gm()
{
   cout << endl;
   if (_c[0].nb == 0)
      throw OFELIException("In Domain::gm(): No Contour is defined.");
   string emfile="out.geo", outfile="out.bamg", ofeli_file="out.m";
   ofstream mf(emfile.c_str());

// Count nb. of points
   mf << "MeshVersionFormatted 0\n\nDimension  2" << endl;
   size_t nbb = 0;
   for (size_t is=0; is<_nb_contours; is++) {
      for (size_t i=0; i<_c[is].nb; i++)
         nbb += _l[_c[is].line[i]-1].nb;
   }
   for (size_t n=0; n<_nb_holes; n++) {
      for (size_t i=0; i<_h[n].nb; i++)
         nbb += _l[_h[n].line[i]-1].nb;
   }
   _nb_lines += nbb;
   _nb_vertices = nbb;

// Output Points of the external contour
   mf << "\nVertices  " << setw(5) << _nb_vertices << endl;
   _ln.resize(_nb_lines);
   size_t kl=0, kl0=0;
   for (size_t is=0; is<_nb_contours; is++) {
      _c[is].first_line = kl + 1;
      for (size_t i=0; i<_c[is].nb; i++) {
         size_t n=_c[is].line[i];
         int dc=_l[n-1].Dcode, nc=_l[n-1].Ncode;
         _ln[kl].i = kl + 1;
         if (i==0)
            kl0 = kl + 1;
         _ln[kl].dc = dc, _ln[kl].nc = nc;
         mf << "  " << setprecision(10) << _l[n-1].node[0].x
            << "  " << setprecision(10) << _l[n-1].node[0].y
            << setw(6) << _v[_l[n-1].n1-1].code << endl;
         for (size_t j=1; j<_l[n-1].nb; j++) {
            mf << "  " << setprecision(10) << _l[n-1].node[j].x
               << "  " << setprecision(10) << _l[n-1].node[j].y
               << setw(6) << dc << endl;
            _ln[kl+1].i = _ln[kl].j = kl + 2;
            _ln[kl+1].dc = dc, _ln[kl+1].nc = nc;
            kl++;
         }
         _ln[kl].j = kl + 2;
         kl++;
      }
      _ln[kl-1].j = kl0;
   }

// Output Points of the holes
   for (size_t k=0; k<_nb_holes; k++) {
      for (size_t i=0; i<_h[k].nb; i++) {
         size_t n = _h[k].line[i];
         int dc = _l[n-1].Dcode, nc = _l[n-1].Ncode;
         if (i==0)
            kl0 = kl + 1;
         _ln[kl].i = kl + 1;
         _ln[kl].dc = dc; _ln[kl].nc = nc;
         mf << "  " << setprecision(10) << _l[n-1].node[0].x
            << "  " << setprecision(10) << _l[n-1].node[0].y
            << setw(6) << dc << endl;
         for (size_t j=1; j<_l[n-1].nb; j++) {
            mf << "  " << setprecision(10) << _l[n-1].node[j].x
               << "  " << setprecision(10) << _l[n-1].node[j].y
               << setw(6) << _v[_l[n-1].n1-1].code << endl;
            _ln[kl+1].i = _ln[kl].j = kl + 2;
            _ln[kl+1].dc = dc, _ln[kl+1].nc = nc;
            kl++;
         }
         _ln[kl].j = kl + 2;
         kl++;
      }
      _ln[kl-1].j = kl0;
   }

// Output Lines
   mf << "\nEdges  " << setw(5) << kl << endl;
   for (size_t i=0; i<kl; i++)
      mf << setw(6) << _ln[i].i << setw(6) << _ln[i].j
         << setw(6) << _ln[i].dc << endl;

// Output mesh sizes
   mf << "\nhVertices" << endl;
   kl = 0;
   for (size_t nct=0; nct<_nb_contours; nct++) {
      for (size_t i=0; i<_c[nct].nb; i++) {
         size_t n = _c[nct].line[i];
         real_t h = _v[_l[n-1].n1-1].h;
         mf << setprecision(6) << setw(10) << h;
         for (size_t j=1; j<_l[n-1].nb; j++)
            mf << setprecision(6) << setw(10) << h;
         mf << endl;
      }
   }
   for (size_t k=0; k<_nb_holes; k++) {
      for (size_t i=0; i<_h[k].nb; i++) {
         real_t h = _v[_l[_h[k].line[i]-1].n1-1].h;
         mf << setprecision(6) << h << "  ";
         for (size_t j=1; j<_l[_h[k].line[i]-1].nb; j++)
            mf << setprecision(6) << h << "  ";
      }
      mf << endl;
   }

// Output subdomains
   if (_nb_sub_domains > 0) {
      mf << "\nSubDomain " << setw(6) << _nb_sub_domains << endl;
      for (size_t k=0; k<_nb_sub_domains; k++)
         mf << setw(6) << 2 << setw(6) << _c[k].first_line << setw(6) << 1
            << setw(6) << 10 << endl;
   }
   mf << endl << "End" << endl;

   cout << "************************************************************************" << endl;
   cout << "                   Output provided by BAMG" << endl;
   cout << "************************************************************************" << endl;
   bamg1(emfile,outfile);
   cout << "************************************************************************" << endl;
   _theMesh = new Mesh;
   cout << "Converting output file to ofeli format ..." << endl;
   getBamg(outfile,*_theMesh,_nb_dof);
   ::remove(emfile.c_str());
   ::remove(outfile.c_str());
}


int Domain::get()
{
   int key;
   string mg;
   _ff = new FFI();
   _nb_dof = 1;
   _ff->setKeywords(_kw);
   int ret = 1;
   do {
      switch (key=_ff->getKW("> ")) {

         case -1:
            cout << "Unknown or Ambiguous Keyword !\n";
            break;

         case 0:
            return -1;

         case 1:
            cout << "\n";
            cout << "help:      display this help\n";
            cout << "dim:       Enter space dimension\n";
            cout << "dof:       Enter (const) nb. of degrees of freedom per node\n";
            cout << "vertex:    Enter a vertex\n";
            cout << "line:      Enter a line\n";
            //            cout << "curve     : Enter a curve\n";
            cout << "circle:    Enter a circular arc\n";
            cout << "contour:   Enter a (external) contour\n";
            cout << "hole:      Enter a hole\n";
            cout << "subdomain: Enter a subdomain\n";
            cout << "dof:       Nb. of degrees of freedom per node\n";
            cout << "rectangle: Mesh generation in a rectangle\n";
            //            cout << "disk      : Mesh generation in a disk\n";
            cout << "dv:        Remove a vertex\n";
            cout << "dl:        Remove a line\n";
            cout << "list:      List an item\n";
            cout << "mesh:      Generate mesh\n";
            cout << "quit:      end the session\n";
            break;

         case 2:
            getVertex();
            break;

         case 3:
            _ret_line = getLine();
            break;

         case 4:
            getCircle();
            break;

         case 5:
            getSubDomain();
            break;

         case 6:
            Rectangle();
            return -1;

         case 7:
            Disk();
            break;

         case 8:
            deleteVertex();
            break;

         case 9:
            deleteLine();
            break;

         case 10:
            list();
            break;

         case 11:
            mg = _ff->getS("Mesh Generator: ");
            if (mg == "easymesh")
               _ret_save = saveAsEasymesh();
            else if (mg == "bamg")
               _ret_save = saveAsBamg();
            else if (mg == "triangle")
               _ret_save = saveAsTriangle();
            else
               ;
            break;

         case 12:
            gm();
            break;

         case 13:
            _nb_dof = _ff->getI();
            break;

         case 14:
            _dim = _ff->getI();
            break;

         case 15:
            _ret_line = getCurve();
            break;

         case 16:
            _ret_cont = getContour();
            break;

         case 17:
            getHole();
            break;
       }
   }
   while ((key>=-1)&&(key<=17));
   return ret;
}


void Domain::get(const string& file)
{
   XMLParser p(file);
   p.get(*this);
}


int Domain::zero_code(const vector<int>& code)
{
   for (size_t i=0; i<_nb_dof; i++)
      if (code[i] != 0)
         return int(i+1);
   return 0;
}


int Domain::remove(size_t          item,
                   size_t&         length,
                   vector<size_t>& set)
{
   size_t i;
   int k = -1;
   for (i=0; i<length; i++)
      if (item == set[i])
         k = int(i);
   if (k!=-1) {
      for (i=k; i<length; i++)
         set[i] = set[i+1];
      length--;
   }
   return k;
}


void Domain::dof_code(int          mark,
                      vector<int>& code)
{
   int j = mark, kk;
   int sign = 1;
   if (mark < 0) {
      sign = -1;
      j *= sign;
   }
   code.clear();
   for (size_t k=0; k<_nb_dof; k++) {
      kk = int(pow(10.,real_t(_nb_dof-k-1)));
      int m = j/kk;
      code.push_back(m*sign);
      j -= m*kk;
   }
}


void Domain::dof_code(int  mark,
                      int* code)
{
   int j = mark, m, kk;
   int sign = 1;
   if (mark < 0) {
      sign = -1;
      j *= sign;
   }
   for (size_t k=0; k<_nb_dof; k++) {
      kk = int(pow(10.,real_t(_nb_dof-k-1)));
      code[k] = m = j/kk;
      code[k] *= sign;
      j -= m*kk;
   }
}

} /* namespace OFELI */
