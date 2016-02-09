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

                         Implementation of class 'Domain'

  ==============================================================================*/

#include "mesh/Domain.h"
#include "mesh/saveMesh.h"
#include "mesh/Mesh.h"
#include "io/XMLParser.h"
#include <algorithm>

namespace OFELI {

Domain::Domain()
{
   _nb_vertices = _nb_lines = _nb_holes = _nb_contours = 0;
   _nb_sub_domains = 0;
   _sub_domain = 0;
   _nb_required_vertices = _nb_required_edges = 0;
   _ret_cont = _ret_save = _ret_sd = 1;
   _v.resize(MAX_NB_VERTICES);
   _V.resize(MAX_NB_VERTICES);
   _l.resize(MAX_NB_LINES);
   _c.resize(MAX_NB_SUB_DOMAINS);
   _h.resize(MAX_NB_HOLES);
   _sd.resize(MAX_NB_SUB_DOMAINS);
   _ln.resize(MAX_NB_LINE_NODES*MAX_NB_LINES);
   init_kw();
   _ff = NULL;
}


Domain::Domain(const string& file)
{
   _nb_vertices = _nb_lines = _nb_holes = _nb_contours = 0;
   _nb_sub_domains = 0;
   _sub_domain = 0;
   _nb_required_vertices = _nb_required_edges = 0;
   _ret_cont = _ret_save = _ret_sd = 1;
   _v.resize(MAX_NB_VERTICES);
   _V.resize(MAX_NB_VERTICES);
   _l.resize(MAX_NB_LINES);
   _c.resize(MAX_NB_SUB_DOMAINS);
   _h.resize(MAX_NB_HOLES);
   _sd.resize(MAX_NB_SUB_DOMAINS);
   _ln.resize(MAX_NB_LINE_NODES*MAX_NB_LINES);
   init_kw();
   _domain_file = file;
   _ff = NULL;
   get(file);
}


Domain::~Domain()
{
   if (_ff)
      delete _ff;
}


void Domain::init_kw()
{
   _kw.resize(19);
   _kw[ 0] = "q$uit";
   _kw[ 1] = "he$lp";
   _kw[ 2] = "v$ertex";
   _kw[ 3] = "line$";
   _kw[ 4] = "circ$le";
   _kw[ 5] = "sub$domain";
   _kw[ 6] = "rectangle$";
   _kw[ 7] = "disk$";
   _kw[ 8] = "dv$";
   _kw[ 9] = "dl$";
   _kw[10] = "lis$t";
   _kw[11] = "save$";
   _kw[12] = "mesh$";
   _kw[13] = "dof$";
   _kw[14] = "dim$";
   _kw[15] = "curve$";
   _kw[16] = "con$tour";
   _kw[17] = "h$ole";
   _kw[18] = "EOF$";
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
   _required_vertex[_nb_required_vertices++] = v;
}


void Domain::insertRequiredEdge(size_t e)
{
   _required_edge[_nb_required_edges++] = e;
}


void Domain::insertVertex(real_t x,
                          real_t y,
                          real_t h,
                          int    code)
{
   _v[_nb_vertices].x = x;
   _v[_nb_vertices].y = y;
   _v[_nb_vertices].h = h;
   _v[_nb_vertices++].code = code;
}


void Domain::getVertex()
{
   size_t k = _ff->getI("Label of vertex: ");
   try {
      if (k < 1)
         THROW_RT("getVertex(): Illegal Number of vertices " + itos(k));
   }
   CATCH("Domain");
   try {
      if (k > MAX_NB_VERTICES)
         THROW_RT("getVertex(): Number of vertices " + itos(k) +
                  " exceeds its maximal value " + itos(MAX_NB_VERTICES));
   }
   CATCH("Domain");
   int l = insert(k,_nb_vertices,_v_label);
   try {
      if (l!=-1)
         THROW_RT("getVertex(): Vertex " + itos(k) + "is redefined.");
   }
   CATCH("Domain");

   _v[k-1].x = _ff->getD("x-coordinate: ");
   _v[k-1].y = _ff->getD("y-coordinate: ");
   _v[k-1].h = _ff->getD("Spacing: ");
   _v[k-1].code = _ff->getI("Code: ");
   cout << "VERTEX ENTERED.\n";
   _nb_vertices = std::max(k,_nb_vertices);
}


int Domain::getCurve()
{
   vector<real_t> s(MAX_NB_LINE_NODES);
   size_t i=0;

// Label of curve
   size_t label = _ff->getI("Label of the curve: ");
   try {
      if (label < 1)
         THROW_RT("getCurve(): Illegal label of curve " + itos(label));
   }
   CATCH("Domain");
   try {
      if (label > MAX_NB_LINES)
         THROW_RT("getCurve(): Number of curves " + itos(label) + 
                  " exceeds its maximal value " + itos(MAX_NB_LINES));
   }
   CATCH("Domain");
   int k = insert(label,_nb_lines,_l_label);
   try {
      if (k!=-1)
         THROW_RT("getCurve(): Curve " + itos(label) + " is redefined");
   }
   CATCH("Domain");

// First end point
   size_t n1 = _ff->getI("First End Point: ");
   _l[label-1].n1 = n1;
   k = insert(n1,_nb_vertices,_v_label);
   try {
      if (k==-1)
         THROW_RT("getCurve(): Curve " + itos(label) + " is redefined");
   }
   CATCH("Domain");

// Second end point
   size_t n2 = _ff->getI("Second End Point: ");
   _l[label-1].n2 = n2;
   k = insert(n2,_nb_vertices,_v_label);
   if (k==-1)
      cout << "WARNING: A new vertex is entered.\n";

// Type of curve
   int type = _ff->getI("Curve's type: ");

// Case of a straight line
   if (type == 0) {
      _l[label-1].nb = 3;
      _l[label-1].node[0] = _v[n1-1];
      _l[label-1].node[1] = 0.5*(_v[n1-1] + _v[n2-1]);
      _l[label-1].node[2] = _v[n2-1];
   }

// Case of a curve given by equation
   else if (type == -1) {
      real_t data[4];
      string regex = _ff->getE("Curve's equation: ");
      theParser.Parse(regex,"x,y,z,t");
      size_t nb = _ff->getI("nb. of discretization points: ");
      if (nb<3)
         nb = 3;
      _l[label-1].nb = nb;
      data[0] = _v[n1-1].x; data[1] = _v[n1-1].y; data[2] = _v[n1-1].z;
      real_t vv = theParser.Eval(data);
      try {
         if (fabs(vv) > 1.e-8)
            THROW_RT("getCurve(): )");
      }
      CATCH("Domain");
      data[0] = _v[n2-1].x; data[1] = _v[n2-1].y; data[2] = _v[n2-1].z;
      vv = theParser.Eval(data);
      try {
         if (fabs(vv) > 1.e-8)
            THROW_RT("getCurve(): ");
      }
     CATCH("Domain");
     _l[label-1].node[0] = _v[n1-1];
     _l[label-1].node[1] = _v[n2-1];
     data[0] = 0.9*_l[label-1].node[0].x + 0.1*_l[label-1].node[1].x;
     data[1] = 0.9*_l[label-1].node[0].y + 0.1*_l[label-1].node[1].y;
//    data[2] = 0.1*_l[label-1].node[1].z;
     for (i=1; i<nb-1; i++) {
        Position(i/real_t(nb-1),data);
        _l[label-1].node[i].x = data[0];
        _l[label-1].node[i].y = data[1];
//      _l[label-1].node[i].z = data[2];
     }
   }

// Case of a predefined curve
   else if (type == 1)
      string shape = _ff->getS("Curve's shape: ");

// Dirichlet and Neumann codes
   _l[label-1].Dcode = _ff->getI("Dirichlet Code: ");
   _l[label-1].Ncode = _ff->getI("Neumann Code: ");
   cout << "CURVE ENTERED." << endl;
   _nb_lines = std::max(label,_nb_lines);
   return 0;
}


int Domain::Position(real_t  s,
                     real_t* data)
{
   real_t d0[4], d[4], dxf, dyf, ex, ey, ez;
   real_t det, g;
   int it=1;
   d0[0] = d[0] = data[0];
   d0[1] = d[1] = data[1];
   d0[2] = d[2] = data[2];
   data[0] -= theParser.Eval(d);
   data[1] -= theParser.Eval(d);
   ex = data[0] - d[0]; ey = data[1] - d[1]; ez = data[3] - d[3];
   if (fabs(ex*ex+ey*ey+ez*ez) < 1.e-8)
      return 0;
//   data[2] -= theParser.Eval(d);
   while (it<50) {
      dxf = (EVAL(data)-EVAL(d))/(data[0]-d[0]);
      dyf = (EVAL(data)-EVAL(d))/(data[1]-d[1]);
      g = 0.5*(s*s - (d[0]-d0[0])*(d[0]-d0[0]) - (d[1]-d0[1])*(d[1]-d0[1]));
//      dzf = (theParser.Eval(data)-theParser.Eval(d))/(data[2]-d[2]);
      d[0] = data[0]; d[1] = data[1];// d[2] = data[2];
      det = (d[0]-d0[0])*dyf - (d[1]-d0[1])*dxf;
      data[0] += ((data[1]-d0[1])*EVAL(data) + dyf*g)/det;
      data[1] -= ((data[0]-d0[0])*EVAL(data) + dxf*g)/det;
      ex = (data[0] - d[0])*(data[0] - d[0]);
      ey = (data[1] - d[1])*(data[1] - d[1]);
      ez = (data[2] - d[2])*(data[2] - d[2]);
      it++;
      if (fabs(ex+ey+ez) < 1.e-8)
         return it;
   }
   return 0;
}


void Domain::insertLine(size_t n1,
                        size_t n2,
                        int    dc,
                        int    nc)
{
   _l[_nb_lines].n1 = n1;
   _l[_nb_lines].n2 = n2;
   _l[_nb_lines].Dcode = dc;
   _l[_nb_lines].Ncode = nc;
   _l[_nb_lines].nb = 2;
   _l[_nb_lines].node[0] = _v[n1-1];
   _l[_nb_lines].node[1] = 0.5*(_v[n1-1] + _v[n2-1]);
   _l[_nb_lines++].node[2] = _v[n2-1];
}


int Domain::getLine()
{
   vector<real_t> s(MAX_NB_LINE_NODES);
   size_t label = _ff->getI("Label of line: ");
   try {
      if (label < 1)
         THROW_RT("getLine(): Illegal Label of line " + itos(label));
   }
   CATCH("Domain");
   try {
      if (label > MAX_NB_VERTICES)
         THROW_RT("getLine(): Number of lines " + itos(label) +
                  " exceeds its maximal value " + itos(MAX_NB_LINES));
   }
   CATCH("Domain");

   int k = insert(label,_nb_lines,_l_label);
   try {
      if (k!=-1)
         THROW_RT("getLine(): Line " + itos(label) + " is redefined.");
   }
   CATCH("Domain");

   size_t n1 = _ff->getI("First End Point: ");
   _l[label-1].n1 = n1;
   k = insert(n1,_nb_vertices,_v_label);
   try {
      if (k==-1)
         THROW_RT("getLine(): Line " + itos(label) + " : A new vertex is entered.");
   }
   CATCH("Domain");

   size_t n2 = _ff->getI("Second End Point: ");
   _l[label-1].n2 = n2;
   k = insert(n2,_nb_vertices,_v_label);
   if (k==-1)
      cout << "WARNING: A new vertex is entered.\n";

   _l[label-1].Dcode = _ff->getI("Dirichlet Code: ");
   _l[label-1].Ncode = _ff->getI("Neumann Code: ");
   _l[label-1].nb = 2;
   _l[label-1].node[0] = _v[n1-1];
   _l[label-1].node[1] = 0.5*(_v[n1-1] + _v[n2-1]);
   _l[label-1].node[2] = _v[n2-1];
   cout << "LINE ENTERED." << endl;
   _nb_lines = std::max(label,_nb_lines);
   return 0;
}


void Domain::insertCircle(size_t n1,
                          size_t n2,
                          size_t n3,
                          int    dc,
                          int    nc)
{
   size_t nb = 36;
   _l[_nb_lines].n1 = n1;
   _l[_nb_lines].n2 = n2;
   _l[_nb_lines].nb = nb;
   _l[_nb_lines].node[0] = _v[n1-1];
   _l[_nb_lines].node[1] = _v[n2-1];
   _l[_nb_lines].Dcode = dc;
   _l[_nb_lines].Ncode = nc;
   real_t h = _v[n1-1].h;
   real_t cx = _v[n3-1].x, cy = _v[n3-1].y;
   real_t a1 = _v[n1-1].x - cx, a2 = _v[n2-1].x - cx;
   real_t b1 = _v[n1-1].y - cy, b2 = _v[n2-1].y - cy;
   real_t r = sqrt(a1*a1 + b1*b1);
   real_t theta1 = atan2(b1,a1), theta2 = atan2(b2,a2);
   real_t theta = theta1;
   if (theta1 >= theta2)
      theta2 += 2*OFELI_PI;
   for (size_t i=0; i<nb; i++) {
      _l[_nb_lines].n1 = n1;
      _l[_nb_lines].n2 = n2;
      _l[_nb_lines].nb = nb;
      _l[_nb_lines].Dcode = dc;
      _l[_nb_lines].Ncode = nc;
      theta += (theta2-theta1)/nb;
      if (i<nb-1) {
         _l[_nb_lines].n2 = _nb_vertices + 1;
         _v[_nb_vertices].x = cx + r*cos(theta);
         _v[_nb_vertices].y = cy + r*sin(theta);
         _v[_nb_vertices].h = h;
         _v[_nb_vertices].code = dc;
      }
      if (i>0)
         _l[_nb_lines].n1 = _nb_vertices;
      _nb_vertices++;
      _nb_lines++;
   }
   _nb_vertices--;
}


void Domain::getCircle()
{
   size_t nb=36;
   size_t n1 = _ff->getI("First end vertex: ");
   _l[_nb_lines].n1 = n1;
   int k = insert(n1,_nb_vertices,_v_label);
   _l[_nb_lines].node[0] = _v[n1-1];
   if (k==-1)
      cout << "WARNING: A new point is entered." << endl;

   size_t n2 = _ff->getI("Second end vertex: ");
   _l[_nb_lines].n2 = n2;
   k = insert(n2,_nb_vertices,_v_label);
   _l[_nb_lines].node[1] = _v[n2-1];
   if (k==-1)
      cout << "WARNING: A new point is entered." << endl;

   size_t n3 = _ff->getI("Center vertex: ");
   k = insert(n3,_nb_vertices,_v_label);
   _l[_nb_lines].nb = nb;
   if (k==-1)
      cout << "WARNING: A new point is entered." << endl;

   int dc = _l[_nb_lines].Dcode = _ff->getI("Dirichlet Code: ");
   _l[_nb_lines].Ncode = _ff->getI("Neumann Code: ");

   real_t h = _v[n1-1].h;
   real_t cx = _v[n3-1].x, cy = _v[n3-1].y;
   real_t a1 = _v[n1-1].x - cx, a2 = _v[n2-1].x - cx;
   real_t b1 = _v[n1-1].y - cy, b2 = _v[n2-1].y - cy;
   real_t r = sqrt(a1*a1 + b1*b1);
   real_t theta1 = atan2(b1,a1), theta2 = atan2(b2,a2);
   real_t theta = theta1;
   if (theta1 >= theta2)
      theta2 += 2*OFELI_PI;
   for (size_t i=0; i<nb; i++) {
      theta += (theta2-theta1)/nb;
      if (i<nb-1) {
         _l[_nb_lines].n2 = _nb_vertices + 1;
         _v[_nb_vertices].x = cx + r*cos(theta);
         _v[_nb_vertices].y = cy + r*sin(theta);
         _v[_nb_vertices].h = h;
         _v[_nb_vertices].code = dc;
      }
      if (i>0)
         _l[_nb_lines].n1 = _nb_vertices;
      _nb_vertices++;
      _nb_lines++;
   }
   _nb_vertices--;
   cout << "CIRCLE ENTERED.\n";
}


void Domain::insertContour(const vector<size_t>& c)
{
   size_t i, i1, i2, n1, n2;
   size_t nb = c.size();
   try {
      if (nb>_nb_lines)
         THROW_RT("insertContour(vector<size_t>): Number of lines is larger than number of defined lines.");
   }
   CATCH("Domain");
   _c[_nb_contours].nb = nb;
   for (i=0; i<nb; i++)
      _c[_nb_contours].line[i] = c[i];
   for (i=0; i<nb-1; i++) {
      n1 = _c[_nb_contours].line[i];
      n2 = _c[_nb_contours].line[i+1];
      i1 = _l[n1-1].n2;
      i2 = _l[n2-1].n1;
      _c[_nb_contours].orientation[i] = 1;
      if (i1 != i2) {
         i1 = _l[n1-1].n1;
         i2 = _l[n2-1].n1;
         _c[_nb_contours].orientation[i] = -1;
         try {
            if (i1 != i2)
               THROW_RT("insertContour(vector<size_t>): Lines " + itos(n1) + " and " + itos(n2) +
                        " are not correctly connected.\nContour cancelled.");
         }
         CATCH("Domain");
      }
   }
   n1 = _c[_nb_contours].line[0];
   n2 = _c[_nb_contours].line[nb-1];
   try {
      if (_l[n1-1].n1 != _l[n2-1].n2)
         THROW_RT("insertContour(vector<size_t>): Contour is not closed.");
   }
   CATCH("Domain");
   _nb_sub_domains++;
   _nb_contours++;
}


int Domain::getContour()
{
   try {
      if (_nb_contours+1 > MAX_NB_SUB_DOMAINS)
         THROW_RT("getContour(): Number of sub-domains " + itos(_nb_contours) +
                  " exceeds its maximum value " + itos(MAX_NB_HOLES+1));
   }
   CATCH("Domain");
   size_t nb = _ff->getI("Nb. of lines: ");
   try {
      if (nb > _nb_lines)
         THROW_RT("getContour(): Label of line is larger than number of defined lines.");
   }
   CATCH("Domain");
   _c[_nb_contours].nb = nb;
   cout << "Please give list of lines in direct order.\n";
   for (size_t i=0; i<nb; i++)
      _c[_nb_contours].line[i] = _ff->getI("Line Label: ");

   for (size_t i=0; i<nb-1; i++) {
      size_t n1 = _c[_nb_contours].line[i], n2 = _c[_nb_contours].line[i+1];
      size_t i1 = _l[n1-1].n2, i2 = _l[n2-1].n1;
      _c[_nb_contours].orientation[i] = 1;
      if (i1 != i2) {
         i1 = _l[n1-1].n1; i2 = _l[n2-1].n1;
         _c[_nb_contours].orientation[i] = -1;
         try {
            if (i1 != i2)
               THROW_RT("getContour(): Lines " + itos(n1) + " and " + itos(n2) + 
                        " are not correctly connected.\nContour cancelled.");
         }
         CATCH("Domain");
      }
   }

   size_t n1 = _c[_nb_contours].line[0], n2 = _c[_nb_contours].line[nb-1];
   try {
      if (_l[n1-1].n1 != _l[n2-1].n2)
         THROW_RT("getContour(): Contour is not closed.");
   }
   CATCH("Domain");
   _nb_sub_domains++;
   cout << "CONTOUR ENTERED." << endl;
   _nb_contours++;
   return 0;
}


void Domain::insertHole(const vector<size_t>& h)
{
   size_t nb=h.size();
   _h[_nb_holes].nb = nb;
   for (size_t i=0; i<nb; i++)
      _h[_nb_holes].line[i] = h[i];
   for (size_t i=0; i<nb-1; i++) {
      size_t n1 = _h[_nb_holes].line[i], n2 = _h[_nb_holes].line[i+1];
      size_t i1 = _l[n1-1].n2, i2 = _l[n2-1].n1;
      try {
         if (i1 != i2)
            THROW_RT("insertHole(vector<size_t>): Lines " + itos(n1) + " and " + itos(n2) +
                     " are not correctly connected.\nHole cancelled.");
      }
      CATCH("Domain");
   }
   _nb_holes++;
}


int Domain::getHole()
{
   size_t label = _ff->getI("Label of hole: ");
   while ((label<1) || (label>MAX_NB_HOLES)) {
      try {
         THROW_RT("getHole(): Label exceeds max. nb. of allowed holes (" + itos(MAX_NB_HOLES) + ")");
      }
      CATCH("Domain");
     label = _ff->getI("Label of hole: ");
   }
   int k = insert(label,_nb_holes,_h_label);
   try {
      if (k!=-1)
         THROW_RT("getHole(): Hole " + itos(label) + " is redefined.");
   }
   CATCH("Domain");
   size_t nb = _ff->getI("Nb. of lines: ");
   _h[label-1].nb = nb;
   cout << "Please give list of lines in direct order.\n";
   for (size_t i=0; i<size_t(nb); i++)
      _h[label-1].line[i] = _ff->getI("Line Label: ");

   for (size_t i=0; i<size_t(nb-1); i++) {
      size_t n1 = _h[label-1].line[i], n2 = _h[label-1].line[i+1];
      size_t i1 = _l[n1-1].n2, i2 = _l[n2-1].n1;
      try {
         if (i1 != i2)
            THROW_RT("getHole(): Lines " + itos(n1) + " and " + itos(n2) + " are not correctly connected.");
      }
      CATCH("Domain");
   }
   return 0;
}


void Domain::insertSubDomain(size_t ln,
                             int    orient,
                             int    code)
{
   try {
      if (ln==0 || ln>_nb_lines)
         THROW_RT("insertSubDomain(size_t,int,int): Illegal line label");
   }
   CATCH("Domain");
   try {
      if (orient!=1 && orient!=-1)
         THROW_RT("insertSubDomain(size_t,int,int): Orientation must be equal to 1 or -1");
   }
   CATCH("Domain");
   _sd[_nb_sub_domains].line = int(ln);
   _sd[_nb_sub_domains].orient = orient;
   _sd[_nb_sub_domains].code = code;
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
   try {
      if (_sub_domain+1 > _nb_sub_domains)
         THROW_RT("insertSubDomain(size_t,int): Number of read subdomains is larger than number of contours");
   }
   CATCH("Domain");
   _sd[_sub_domain].contour = int(n);
   _sd[_sub_domain].code = code;
   _sub_domain++;
}


int Domain::getSubDomain()
{
   if (_sub_domain+1 > _nb_sub_domains)
      cerr << "Nb of read subdomains is larger than number of contours" << endl;
   _sd[_sub_domain].contour = _ff->getI("Contour label: ");

   int code = _ff->getI("Material code for subdomain: ");
   _sd[_sub_domain++].code = code;
   cout << "SUBDOMAIN ENTERED." << endl;
   return 0;
}


int Domain::Rectangle()
{
   size_t n1, n2, n3, n4, m1, m2, m3;
   int cd[MAX_NBDOF_NODE];

   real_t ax = _ff->getD("x min: ");
   real_t ay = _ff->getD("y min: ");
   real_t bx = _ff->getD("x max: ");
   real_t by = _ff->getD("y max: ");
   real_t lx = bx - ax, ly = by - ay;
   int type = _ff->getI("Mesh Type (0: rectangles, 1: triangles): ");
   size_t ne1 = _ff->getI("Nb. of sub-intervals in the x-direction: ");
   size_t ne2 = _ff->getI("Nb. of sub-intervals in the y-direction: ");

   int c1 = _ff->getI("Node code for y=ymin: ");
   int c2 = _ff->getI("Node code for x=xmax: ");
   int c3 = _ff->getI("Node code for y=ymax: ");
   int c4 = _ff->getI("Node code for x=xmin: ");
   int cc1 = _ff->getI("Node code at x=xmin, y=ymin: ");
   int cc2 = _ff->getI("Node code at x=xmax, y=ymin: ");
   int cc3 = _ff->getI("Node code at x=xmax, y=ymax: ");
   int cc4 = _ff->getI("Node code at x=xmin, y=ymax: ");
   int cs1 = _ff->getI("Side code for y=ymin: ");
   if (cs1 && cc1) {
      cout << "Error in data: You cannot give nonzeo codes for both nodes and sides." << endl;
      c1 = _ff->getI("Node code for y=ymin: ");
      cs1 = _ff->getI("Side code for y=ymin: ");
   }
   int cs2 = _ff->getI("Side code for x=xmax: ");
   if (cs2 && cc2) {
      cout << "Error in data: You cannot give nonzeo codes for both nodes and sides." << endl;
      c2 = _ff->getI("Node code for x=xmax: ");
      cs2 = _ff->getI("Side code for x=xmax: ");
   }
   int cs3 = _ff->getI("Side code for y=ymax: ");
   if (cs3 && cc3) {
      cout << "Error in data: You cannot give nonzeo codes for both nodes and sides." << endl;
      c3 = _ff->getI("Node code for y=ymax: ");
      cs3 = _ff->getI("Side code for y=ymac: ");
   }
   int cs4 = _ff->getI("Side code for x=xmin: ");
   if (cs4 && cc4) {
      cout << "Error in data: You cannot give nonzeo codes for both nodes and sides." << endl;
      c4 = _ff->getI("Node code for x=xmin: ");
      cs4 = _ff->getI("Side code for x=xmin: ");
   }
   string file = _ff->getS("Output mesh file name: ");

   Mesh mesh;
   mesh.setDim(2);
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
         if (i==0 && j==0)
            c = cc1;
         else if (i==ne2 && j==0)
            c = cc4;
         else if (i==ne2 && j==ne1)
            c = cc3;
         else if (i==0   && j==ne1)
            c = cc2;
         else
            c = setCode(ne1,ne2,i,j,_nb_dof,c1,c2,c3,c4);
         dof_code(c,cd);
         The_node.setCode(cd);
         The_node.setDOF(first_dof,_nb_dof);
         mesh.Add(the_node);
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
            The_element.Add(mesh[n1]);
            The_element.Add(mesh[n2]);
            The_element.Add(mesh[n3]);
            The_element.Add(mesh[n4]);
            mesh.Add(the_element);
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
            The_element.Add(mesh[n1]);
            The_element.Add(mesh[n2]);
            The_element.Add(mesh[n3]);
            mesh.Add(the_element);
            el++; m1 = n3; m2 = n3 - 1; m3 = n1;
            the_element = new Element(el+1,TRIANGLE,1);
            The_element.Add(mesh[m1]);
            The_element.Add(mesh[m2]);
            The_element.Add(mesh[m3]);
            mesh.Add(the_element);
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
         The_side.Add(mesh[m1]);
         The_side.Add(mesh[m2]);
         dof_code(cs1,cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
         mesh.Add(the_side);
      }
      m1 = m2;
   }
   for (size_t i=0; i<ne2; i++) {
      m2 = m1 + ne1 + 1;
      if (cs2!=0) {
         the_side = new Side(++k,LINE);
         The_side.Add(mesh[m1]);
         The_side.Add(mesh[m2]);
         dof_code(cs2,cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
         mesh.Add(the_side);
      }
      m1 = m2;
   }
   for (size_t i=0; i<ne1; i++) {
      m2 = m1 - 1;
      if (cs3!=0) {
         the_side = new Side(++k,LINE);
         The_side.Add(mesh[m1]);
         The_side.Add(mesh[m2]);
         dof_code(cs3,cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
         mesh.Add(the_side);
      }
      m1 = m2;
   }
   for (size_t i=0; i<ne2; i++) {
      m2 = m1 - ne1 - 1;
      if (cs4!=0) {
         the_side = new Side(++k,LINE);
         The_side.Add(mesh[m1]);
         The_side.Add(mesh[m2]);
         dof_code(cs4,cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
         mesh.Add(the_side);
      }
      m1 = m2;
   }
   mesh.put(file);
   return 0;
}


int Domain::Rectangle(real_t* x,
                      size_t  n1, 
                      size_t  n2,
                      real_t  r,
                      int     c1,
                      int     c2,
                      int     c3,
                      int     c4,
                      int     cs1,
                      int     cs2,
                      int     cs3,
                      int     cs4,
                      string  file)
{
   size_t m1, m2, m3, nn1, nn2, nn3;
   int cd[MAX_NBDOF_NODE];
   int cc1=c1, cc2=c2, cc3=c3, cc4=c4;

// Nodes
   Mesh mesh;
   mesh.setDim(2);
   size_t k = 0;
   real_t lx=x[2]-x[0], ly=x[3]-x[1];
   real_t hx, hy=ly/real_t(n2);
   int c;
   if (r != 1 && r > 0)
      hy = ly*(r-1.)/(pow(r,real_t(n2))-1.);
   real_t yy=x[1], cs=1.75;
   size_t first_dof = 1;
   for (size_t i=0; i<=n2; i++) {
      hx = lx/real_t(n1);
      if (r != 1)
         hx = lx*(r-1.)/(pow(r,int(n1))-1.);
      real_t xx=x[0];
      if (r < 0)
         yy = x[1] + 0.5*ly*(1.+tanh(cs*((2.*i)/n1-1.0))/tanh(cs));

      for (size_t j=0; j<=n1; j++) {
         if (r < 0)
            xx = x[0] + 0.5*lx*(1.+tanh(cs*((2.*j)/n2-1.0))/tanh(cs));
         the_node = new Node(++k,Point<real_t>(xx,yy,0));
         The_node.setNbDOF(_nb_dof);
         if (i==0 && j==0)
            c = cc1;
         else if (i==n2 && j==0)
            c = cc4;
         else if (i==n2 && j==n1)
            c = cc3;
         else if (i==0 && j==n1)
            c = cc2;
         else
            c = setCode(n1,n2,i,j,_nb_dof,c1,c2,c3,c4);
         dof_code(c,cd);
         The_node.setCode(cd);
         The_node.setDOF(first_dof,_nb_dof);
         mesh.Add(the_node);
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
         The_element.Add(mesh[nn1]);
         The_element.Add(mesh[nn2]);
         The_element.Add(mesh[nn3]);
         mesh.Add(the_element);
         el++; m1 = nn3; m2 = nn3 - 1; m3 = nn1;
         the_element = new Element(el+1,TRIANGLE,1);
         The_element.Add(mesh[m1]);
         The_element.Add(mesh[m2]);
         The_element.Add(mesh[m3]);
         mesh.Add(the_element);
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
         The_side.Add(mesh[m1]);
         The_side.Add(mesh[m2]);
         dof_code(cs1,cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
         mesh.Add(the_side);
      }
      m1 = m2;
   }
   for (size_t i=0; i<n2; i++) {
      m2 = m1 + n1 + 1;
      if (cs2!=0) {
         the_side = new Side(++k,LINE);
         The_side.Add(mesh[m1]);
         The_side.Add(mesh[m2]);
         dof_code(cs2,cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
         mesh.Add(the_side);
      }
      m1 = m2;
   }
   for (size_t i=0; i<n1; i++) {
      m2 = m1 - 1;
      if (cs3!=0) {
         the_side = new Side(++k,LINE);
         The_side.Add(mesh[m1]);
         The_side.Add(mesh[m2]);
         dof_code(cs3,cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
         mesh.Add(the_side);
      }
      m1 = m2;
   }
   for (size_t i=0; i<n2; i++) {
      m2 = m1 - n1 - 1;
      if (cs4!=0) {
         the_side = new Side(++k,LINE);
         The_side.Add(mesh[m1]);
         The_side.Add(mesh[m2]);
         dof_code(cs4,cd);
         The_side.setNbDOF(_nb_dof);
         for (size_t i=0; i<_nb_dof; i++)
            The_side.setCode(i+1,cd[i]);
         mesh.Add(the_side);
      }
      m1 = m2;
   }
   mesh.put(file);
   return 0;
}


int Domain::setCode(size_t ne1,
                    size_t ne2,
                    size_t i,
                    size_t j,
                    size_t nb_dof,
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
   int         s;
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
   m = 1; n = 1; s = 0;
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
      _l[s].n1 = _el[m-2].n1;
      _l[s].n2 = _el[m-2].n2;
      _l[s].Ncode = nc; s++;
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
   size_t  i, j, k;
   cout.setf(ios::right|ios::scientific);
   if (!_nb_vertices)
      cout << "NO VERTICES.\n";
   else {
      cout << "\n****** List of vertices *****\n\n";
      cout << "LABEL      X               Y           Code\n";
      for (i=0; i<_nb_vertices; i++) {
         j = _v_label[i]-1;
         cout << setw(3) << j+1 << "    " << setprecision(5) << setw(10) << _v[j].x;
         cout << "    " << setprecision(5) << setw(12) << _v[j].y << setw(6);
         cout << _v[j].code << endl;
     }
     cout << endl;
   }

   if (!_nb_lines)
      cout << "NO LINES.\n";
   else {
      cout << "\n****** List of lines ******\n\n";
      cout << "LABEL   CODE   END 1  END 2\n";
      for (i=0; i<_nb_lines; i++) {
         j = _l_label[i]-1;
         cout << setw(3) << j+1 << setw(8) << _l[j].Dcode;
         cout << setw(7) << _l[j].n1 << setw(7) << _l[j].n2 << endl;
      }
      cout << endl;
   }
   if (!_nb_contours || _ret_cont)
      cout << "NO CONTOUR.\n";
   else {
      for (size_t k=0; k<_nb_contours; k++) {
         cout << "\n****** External Contour ******\n\n";
         cout << "NB. OF LINES     LIST OF LINES\n";
         cout << setw(7) << _c[k].nb << "          ";
         for (i=0; i<_c[k].nb; i++)
            cout << setw(5) << _c[k].line[i];
         cout << endl;
      }
   }

   if (!_nb_holes)
      cout << "NO HOLES.\n";
   else {
       cout << "\n****** List of holes ******\n\n";
       cout << "LABEL   NB. OF LINES      LIST OF LINES\n";
       for (j=0; j<_nb_holes; j++) {
          k = _h_label[j]-1;
          cout << setw(5) << _h[k].nb << "  " << setw(7) << k+1 << "        ";
          for (i=0; i<_h[k].nb; i++)
             cout << setw(5) << _h[k].line[i];
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
   size_t  i, j, n, i1, nb, nbb, kl, k, kl0=0, label;
   int     dc, nc, c;
   real_t h;

   try {
      if (!_c[0].nb)
         THROW_RT("saveAsEasyMesh: No Contour is defined.");
   }
   CATCH("Domain");
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
   label = _nb_vertices;
   _nb_vertices += nbb;
   mf << setw(5) << _nb_vertices << "# Nb. of points #" << endl;

// Output Points of the external contour
   mf <<  "# Points on External Contour #" << endl;
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
         label++;
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
            label++;
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
   size_t  i, j, n, i1, nb, nbb, kl, k, kl0=0, label;
   int     c, dc, nc;
   real_t  h;
   string mfile, emfile;

   try {
      if (_c[0].nb == 0)
         THROW_RT("saveAsBamg: No Contour is defined.");
   }
   CATCH("Domain");
   mfile = _ff->getS("File Name (.geo): ");
   emfile = mfile + ".geo";
   ofstream mf(emfile.c_str());

// Count nb. of points
   mf << "MeshVersionFormatted 0\nDimension  2" << endl;
   nbb = 0;
   for (i=0; i<_c[0].nb; i++)
      nbb += _l[_c[0].line[i]-1].nb - 1;
   for (n=0; n<_nb_holes; n++)
      for (i=0; i<_h[n].nb; i++)
         nbb += _l[_h[n].line[i]-1].nb - 1;
   _nb_lines += nbb;
   label = _nb_vertices;
   _nb_vertices += nbb;
   mf << "\nVertices  " << setw(6) << _nb_vertices << endl;

// Output Points of the external contour
   kl = 0;
   for (i=0; i<_c[0].nb; i++) {
      n = _c[0].line[i];
      i1 = _l[n-1].n1;
      c = _v[i1-1].code;
      nb = _l[n-1].nb;
      dc = _l[n-1].Dcode;
      nc = _l[n-1].Ncode;
      _ln[kl].i = kl+1;
      if (i==0)
         kl0 = kl + 1;
      _ln[kl].dc = dc;
	  _ln[kl].nc = nc;
      mf << _l[n-1].node[0].x << "  " << _l[n-1].node[0].y << setw(5) << c << endl;
      for (j=1; j<nb; j++) {
         label++;
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
   for (k=0; k<_nb_holes; k++) {
      for (i=0; i<_h[k].nb; i++) {
         n = _h[k].line[i];
         i1 = _l[n-1].n1;
         c = _v[i1-1].code;
         h = _v[i1-1].h;
         nb = _l[n-1].nb;
         dc = _l[n-1].Dcode;
         nc = _l[n-1].Ncode;
         if (i==0)
            kl0 = kl + 1;
         _ln[kl].i = kl + 1;
         _ln[kl].dc = dc;
		 _ln[kl].nc = nc;
         mf << _l[n-1].node[0].x << "  " << _l[n-1].node[0].y << setw(5) << dc << endl;
         for (j=1; j<nb; j++) {
            label++;
            mf << _l[n-1].node[j].x << "  " << _l[n-1].node[j].y << setw(5) << c << endl;
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
   mf << "\nEdges    " << kl << endl;
   for (i=0; i<kl; i++)
      mf << setw(6) << _ln[i].i << setw(6) << _ln[i].j << setw(6) << _ln[i].dc
         << setw(5) << _ln[i].nc << endl;

// Output mesh sizes
   mf << "\nhVertices  ";
   kl = 0;
   for (i=0; i<_c[0].nb; i++) {
      n = _c[0].line[i];
      i1 = _l[n-1].n1;
      h = _v[i1-1].h;
      nb = _l[n-1].nb;
      mf << "  " << h;
      for (j=1; j<nb; j++)
         mf << "  " << h;
      mf << endl;
   }
   for (k=0; k<_nb_holes; k++) {
      for (i=0; i<_h[k].nb; i++) {
         n = _h[k].line[i];
         i1 = _l[n-1].n1;
         c = _v[i1-1].code;
         h = _v[i1-1].h;
         nb = _l[n-1].nb;
         mf << "   " << h;
         for (j=1; j<nb; j++)
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

   try {
      if (_c[0].nb == 0)
         THROW_RT("saveAsTriangle(): No Contour is defined.");
   }
   CATCH("Domain");
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
      mf << "\nVertices  " << setw(5) << _nb_vertices << endl;
      for (size_t i=0; i<_nb_vertices; i++)
         mf << setprecision(8) << setw(12) << _v[i].x
            << setprecision(8) << setw(12) << _v[i].y
            << setw(6) << _v[i].code << endl;
   }
   if (_nb_lines > 0) {
      mf << "\nEdges  " << setw(5) << _nb_lines << endl;
      for (size_t i=0; i<_nb_lines; i++)
         mf << setw(5) << _l[i].n1 << setw(5) << _l[i].n2
            << setw(5) << _l[i].Dcode << setw(5) << endl;
   }
   if (_nb_vertices > 0) {
      mf << "\nhVertices\n";
      for (size_t i=0; i<_nb_vertices; i++)
         mf << "  " << setprecision(6) << setw(10) << _v[i].h;
      mf << endl;
   }
   if (_nb_required_vertices > 0) {
      mf << "\nRequiredVertices" << setw(5) << _nb_required_vertices << endl;
      for (size_t i=0; i<_nb_required_vertices; i++)
         mf << setw(5) << _required_vertex[i];
      mf << endl;
   }
   if (_nb_required_edges > 0) {
      mf << "\nRequiredEdges" << setw(5) << _nb_required_edges << endl;
      for (size_t i=0; i<_nb_required_edges; i++)
         mf << setw(5) << _required_edge[i];
      mf << endl;
   }
   if (_nb_sub_domains > 0) {
      mf << "\nSubDomain " << setw(6) << _nb_sub_domains << endl;
      for (size_t i=0; i<_nb_sub_domains; i++)
         mf << setw(6) << 2 << setw(6) << _sd[i].line << setw(6) << _sd[i].orient << setw(6)
            << _sd[i].code << endl;
   }
   mf << "\nEnd" << endl;
}


void Domain::genMesh(string geo_file,
                     string bamg_file,
                     string mesh_file)
{
// Generate geometry file
   genGeo(geo_file);

// Generate Bamg file
   main_bamg(geo_file,bamg_file);
   _theMesh = new Mesh;
   getBamg(bamg_file,*_theMesh,_nb_dof);
   _theMesh->put(mesh_file);
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


void Domain::generateMesh(const string& file)
{
   if (_dim==1)
      gm1(file);
   else if (_dim==2)
      gm2(file);
   else if (_dim==3)
      gm3(file);
}


void Domain::gm1(const string& file)
{
}


void Domain::gm2(const string& file)
{
   size_t is, i, j, n, i1, nb, nbb, kl, k, kl0=0, label;
   int c, dc, nc;
   real_t h;

   try {
      if (_c[0].nb == 0)
         THROW_RT("gm2(string): No Contour is defined.");
   }
   CATCH("Domain");
   string prefix = file.substr(0,file.rfind("."));
   string geo_file=prefix+".geo", bamg_file=prefix+".bamg";
   ofstream mf(geo_file.c_str());

// Count nb. of points
   mf << "MeshVersionFormatted 0" << endl;
   mf << "Dimension  2" << endl;
   nbb = 0;
   for (is=0; is<_nb_contours; is++) {
      for (i=0; i<_c[is].nb; i++)
         nbb += _l[_c[is].line[i]-1].nb;
   }
   for (n=0; n<_nb_holes; n++)
      for (i=0; i<_h[n].nb; i++)
         nbb += _l[_h[n].line[i]-1].nb;
   _nb_lines += nbb;
   label = _nb_vertices;
   _nb_vertices = nbb;
   mf << "\nVertices  " << setw(5) << _nb_vertices << endl;

// Output Points of the external contour
   kl = 0;
   for (is=0; is<_nb_contours; is++) {
      for (i=0; i<_c[is].nb; i++) {
         n = _c[is].line[i];
         i1 = _l[n-1].n1;
         c = _v[i1-1].code;
         nb = _l[n-1].nb;
         dc = _l[n-1].Dcode;
         nc = _l[n-1].Ncode;
         _ln[kl].i = kl+1;
         if (i==0)
            kl0 = kl + 1;
         _ln[kl].dc = dc;
         _ln[kl].nc = nc;
         mf << "  " << _l[n-1].node[0].x << "  " << _l[n-1].node[0].y << "  " << c << endl;
         for (j=1; j<nb; j++) {
            label++;
            mf << "  " << _l[n-1].node[j].x << "  " << _l[n-1].node[j].y << "  " << dc << endl;
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

// Output Points of holes
   for (k=0; k<_nb_holes; k++) {
      for (i=0; i<_h[k].nb; i++) {
         n = _h[k].line[i];
         i1 = _l[n-1].n1;
         c = _v[i1-1].code;
         h = _v[i1-1].h;
         nb = _l[n-1].nb;
         dc = _l[n-1].Dcode;
         nc = _l[n-1].Ncode;
         if (i==0)
            kl0 = kl + 1;
         _ln[kl].i = kl + 1;
         _ln[kl].dc = dc;
         _ln[kl].nc = nc;
         mf << "  " << _l[n-1].node[0].x << "  " << _l[n-1].node[0].y << "  " << dc << endl;
         for (j=1; j<nb; j++) {
            label++;
            mf << "  " << _l[n-1].node[j].x << "  " << _l[n-1].node[j].y << "  " << c << endl;
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
   mf << "\nEdges  " << setw(5) << kl << endl;
   for (i=0; i<kl; i++)
      mf << setw(6) << _ln[i].i << setw(6) << _ln[i].j
         << setw(6) << _ln[i].dc << setw(6) << _ln[i].nc << endl;

// Output mesh sizes
   mf << "\nhVertices\n";
   kl = 0;
   for (k=0; k<_nb_contours; k++) {
      for (i=0; i<_c[k].nb; i++) {
         n = _c[k].line[i];
         i1 = _l[n-1].n1;
         h = _v[i1-1].h;
         nb = _l[n-1].nb;
         mf << setprecision(6) << setw(10) << h;
         for (j=1; j<nb; j++)
            mf << setprecision(6) << setw(10) << h;
         mf << endl;
      }
   }
   for (k=0; k<_nb_holes; k++) {
      for (i=0; i<_h[k].nb; i++) {
         n = _h[k].line[i];
         i1 = _l[n-1].n1;
         c = _v[i1-1].code;
         h = _v[i1-1].h;
         nb = _l[n-1].nb;
         mf << setprecision(6) << setw(10) << h;
         for (j=1; j<nb; j++)
            mf << setprecision(6) << setw(10) << h;
      }
      mf << endl;
   }

   mf << "\nSubDomain " << setw(6) << _nb_sub_domains << endl;
   for (k=0; k<_nb_sub_domains; k++)
      mf << setw(6) << 2 << setw(6) << _c[k].first_line << setw(6) << 1 << setw(6)
         << _sd[k].code << endl;
   mf << "\nEnd" << endl;

   main_bamg(geo_file,bamg_file);
   Mesh mesh;
   getBamg(bamg_file,mesh,_nb_dof);
   mesh.put(file);
}


void Domain::gm3(const string &file)
{
}


void Domain::gm()
{
   size_t  is, i, j, n, i1, nb, nbb, kl, k, kl0=0, label;
   int     c, dc, nc;
   real_t h;

   cout << endl;
   try {
      if (_c[0].nb == 0)
         THROW_RT("gm(): No Contour is defined.");
   }
   CATCH("Domain");
   string mfile = _ff->getS("File Name (.geo): ");
   string emfile=mfile+".geo", outfile=mfile+".bamg", ofeli_file=mfile+".m";
   ofstream mf(emfile.c_str());

// Count nb. of points

   mf << "MeshVersionFormatted 0" << endl;
   mf << "Dimension  2" << endl;
   nbb = 0;
   for (is=0; is<_nb_contours; is++)
      for (i=0; i<_c[is].nb; i++)
         nbb += _l[_c[is].line[i]-1].nb;
   for (n=0; n<_nb_holes; n++)
      for (i=0; i<_h[n].nb; i++)
         nbb += _l[_h[n].line[i]-1].nb;
   _nb_lines += nbb;
   label = _nb_vertices;
   _nb_vertices = nbb;
   mf << "\nVertices  " << setw(5) << _nb_vertices << endl;

// Output Points of the external contour
   kl = 0;
   for (is=0; is<_nb_contours; is++) {
      _c[is].first_line = kl + 1;
      for (i=0; i<_c[is].nb; i++) {
         n = _c[is].line[i];
         i1 = _l[n-1].n1;
         c = _v[i1-1].code;
         nb = _l[n-1].nb;
         dc = _l[n-1].Dcode;
         nc = _l[n-1].Ncode;
         _ln[kl].i = kl + 1;
         if (i==0)
            kl0 = kl + 1;
         _ln[kl].dc = dc;
         _ln[kl].nc = nc;
         mf << "  " << setprecision(10) << setw(10) << _l[n-1].node[0].x
            << "  " << setprecision(10) << setw(10) << _l[n-1].node[0].y
            << setw(6) << c << endl;
         for (j=1; j<nb; j++) {
            label++;
            mf << "  " << setprecision(10) << setw(10) << _l[n-1].node[j].x
               << "  " << setprecision(10) << setw(10) << _l[n-1].node[j].y
               << setw(6) << dc << endl;
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

// Output Points of the holes
   for (k=0; k<_nb_holes; k++) {
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
         mf << "  " << setprecision(10) << setw(10) << _l[n-1].node[0].x
            << "  " << setprecision(10) << setw(10) << _l[n-1].node[0].y
            << setw(6) << dc << endl;
         for (j=1; j<nb; j++) {
            label++;
            mf << "  " << setprecision(10) << setw(10) << _l[n-1].node[j].x
               << "  " << setprecision(10) << setw(10) << _l[n-1].node[j].y
               << setw(6) << c << endl;
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
   mf << "\nEdges  " << setw(5) << kl << endl;
   for (i=0; i<kl; i++)
      mf << setw(6) << _ln[i].i << setw(6) << _ln[i].j
         << setw(6) << _ln[i].dc << setw(6) << _ln[i].nc << endl;

// Output mesh sizes
   mf << "\nhVertices\n";
   kl = 0;
   for (size_t nct=0; nct<_nb_contours; nct++) {
      for (i=0; i<_c[nct].nb; i++) {
         n = _c[nct].line[i];
         i1 = _l[n-1].n1;
         h = _v[i1-1].h;
         nb = _l[n-1].nb;
         mf << setprecision(6) << setw(10) << h;
         for (j=1; j<nb; j++)
            mf << setprecision(6) << setw(10) << h;
         mf << endl;
      }
   }
   for (k=0; k<_nb_holes; k++) {
      for (i=0; i<_h[k].nb; i++) {
         n = _h[k].line[i];
         i1 = _l[n-1].n1;
         c = _v[i1-1].code;
         h = _v[i1-1].h;
         nb = _l[n-1].nb;
         mf << setprecision(6) << setw(10) << h;
         for (j=1; j<nb; j++)
            mf << setprecision(6) << setw(10) << h;
      }
      mf << endl;
   }

// Output subdomains
   if (_nb_sub_domains > 0) {
     mf << "\nSubDomain " << setw(6) << _nb_sub_domains << endl;
     for (k=0; k<_nb_sub_domains; k++) {
        mf << setw(6) << 2 << setw(6) << _c[k].first_line << setw(6) << 1
           << setw(6) << _sd[k].code << endl;
     }
   }
   mf << endl << "End" << endl;

   cout << "************************************************************************" << endl;
   cout << "                   Output provided by BAMG" << endl;
   cout << "************************************************************************" << endl;
   main_bamg(emfile,outfile);
   cout << "************************************************************************" << endl;
   Mesh mesh;
   cout << "Converting output file to ofeli format ..." << endl;
   getBamg(outfile,mesh,_nb_dof);
   mesh.put(ofeli_file.c_str());
}


int Domain::get()
{
   int    key;
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
            return ret;

         case 1:
            cout << "\n";
            cout << "help:       display this help\n";
            cout << "dim:        Enter space dimension\n";
            cout << "dof:        Enter (const) nb. of degrees of freedom per node\n";
            cout << "vertex:     Enter a vertex\n";
            cout << "line:       Enter a line\n";
            //            cout << "curve     : Enter a curve\n";
            cout << "circle:     Enter a circular arc\n";
            //            cout << "contour   : Enter a (external) contour\n";
            //           cout << "hole      : Enter a hole\n";
            //           cout << "subdomain : Enter a subdomain\n";
            cout << "dof:        Nb. of degrees of freedom per node\n";
            cout << "rectangle:  Mesh generation in a rectangle\n";
            //            cout << "disk      : Mesh generation in a disk\n";
            cout << "dv:         Remove a vertex\n";
            cout << "dl:         Remove a line\n";
            cout << "list:       List an item\n";
            cout << "save:       Save data in file for a mesh generator\n";
            cout << "mesh:       Generate mesh\n";
            cout << "quit:       end the session\n";
            break;

         case 2:
            getVertex();
            break;

         case 3:
            _ret_line = getLine();
            break;

         case 15:
            _ret_line = getCurve();
            break;

         case 4:
            getCircle();
            break;

         case 16:
            _ret_cont = getContour();
            break;

         case 17:
            getHole();
            break;

         case 5:
            getSubDomain();
            break;

         case 6:
            Rectangle();
            ret = 0;
            break;

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
            cout << "This keyword will be deactivated in future versions\n";
            cout << "Use rather the member function GenerateMesh." << endl;
            gm();
            break;

         case 13:
            _nb_dof = _ff->getI();
            break;

         case 14:
            _dim = _ff->getI();
            break;
       }
   }
   while ((key>=-1)||(key<=17));
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


int Domain::insert(size_t  item,
                   size_t& length,
                   size_t* set)
{
   int k = -1;
   for (size_t i=0; i<length; i++)
      if (item == set[i])
         k = int(i);
   if (k==-1) {
      set[length] = item;
      length++;
   }
   return k;
}


int Domain::remove(size_t  item,
                   size_t& length,
                   size_t* set)
{
   size_t i;
   int k = -1;
   for (i=0; i<length; i++)
      if (item == set[i])
         k = int(i);
   if (k!=-1) {
      for (i=k;i<length;i++)
         set[i] = set[i+1];
      length--;
   }
   return k;
}


void Domain::dof_code(int          mark,
                      vector<int>& code)
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
