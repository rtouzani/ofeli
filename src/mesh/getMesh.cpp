/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2019 Rachid Touzani

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

                   Functions to read mesh files in various formats
                             and construct Mesh instance

  ==============================================================================*/

#include "mesh/getMesh.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::endl;

#include <iomanip>
using std::setw;
using std::ios;
using std::setprecision;


#include "OFELI_Config.h"
#include "io/FFI.h"
#include "mesh/Mesh.h"
#include "mesh/MeshUtil.h"
#include "linear_algebra/Point.h"
#include "util/util.h"
#include "linear_algebra/DMatrix.h"
#include "mesh/Material.h"
#include "OFELIException.h"

namespace OFELI {

extern Material theMaterial;

void getMesh(string             file,
             ExternalFileFormat form,
             Mesh&              mesh,
             size_t             nb_dof)
{
   switch (form) {

      case OFELI_FF:
         break;

      case GMSH:
         getGmsh(file,mesh,nb_dof);
         break;

      case GNUPLOT:
         break;

      case MATLAB:
         getMatlab(file,mesh,nb_dof);
         break;

      case VTK:
         break;

      case TECPLOT:
         break;

      case EASYMESH:
         getEasymesh(file,mesh,nb_dof);
         break;

      case GAMBIT:
         getGambit(file,mesh,nb_dof);
         break;

      case BAMG:
         getBamg(file,mesh,nb_dof);
         break;

      case NETGEN:
         getNetgen(file,mesh,nb_dof);
         break;

      case TETGEN:
         getTetgen(file,mesh,nb_dof);
         break;

      case TRIANGLE_FF:
         getTriangle(file,mesh,nb_dof);
         break;
   }
}


void getBamg(string file,
             Mesh&  mesh,
             size_t nb_dof)
{
   size_t  ii, jj, kk, k, nb_nodes, nb_elements, nb_sides, first_dof;
   int     key;
   Node    *nd;
   Element *el;
   Side    *sd;
   Point<real_t> x;
   string ww;
   static int mark;
   int code[MAX_NBDOF_NODE];
   vector<string> kw;
   kw.push_back("End$");
   kw.push_back("MeshVer$sionFormatted");
   kw.push_back("Dimension$");
   kw.push_back("Vertices$");
   kw.push_back("Edges$");
   kw.push_back("Triangles$");
   kw.push_back("SubDomain$FromMesh");
   kw.push_back("VertexOnGeometricVertex$");
   kw.push_back("VertexOnGeometricEdge$");
   kw.push_back("EdgeOnGeometricEdge$");
   kw.push_back("Identifier$");
   kw.push_back("Geometry$");

   mesh.setDim(2);
   FFI ff(file);
   ff.setKeywords(kw);
   do {
      switch (key=ff.getKW()) {

         case -1:
            break;

         case 0:
            ff.Close();
            break;

//       MeshVersionFormatted
         case 1:
            k = ff.getI();
            mesh.setDim(2);
            break;

//       Dimension
         case 2:
            k = ff.getI();
            break;

//       Vertices
         case 3:
            nb_nodes = ff.getI();
            first_dof = 1;
            for (size_t i=0; i<nb_nodes; ++i) {
               x.x = ff.getD();
               x.y = ff.getD();
               x.z = 0.;
               mark = ff.getI();
               if (mark<0)
                  mark = 0;
               nd = new Node(i+1,x);
               nd->setNbDOF(nb_dof);
               for (size_t j=0; j<nb_dof; ++j)
                  if (code[j]<0)
                     code[j] = 0;
               DOFCode(mark,nb_dof,code);
               nd->setDOF(first_dof,nb_dof);
               nd->setCode(code);
               mesh.Add(nd);
            }
            break;

//       Edges
         case 4:
            nb_sides = ff.getI();
            for (size_t n=0; n<nb_sides; ++n) {
               ii = ff.getI();
               jj = ff.getI();
               mark = ff.getI();
               sd = new Side(n+1,LINE);
               sd->Add(mesh[ii]);
               sd->Add(mesh[jj]);
               DOFCode(mark,nb_dof,code);
               sd->setNbDOF(nb_dof);
               for (size_t i=0; i<nb_dof; i++)
                  sd->setCode(i+1,code[i]);
               mesh.Add(sd);
            }
            break;

//       Triangles
         case 5:
            nb_elements = ff.getI();
            for (size_t n=0; n<nb_elements; ++n) {
               ii = ff.getI();
               jj = ff.getI();
               kk = ff.getI();
               code[0] = ff.getI();
               el = new Element(n+1,TRIANGLE,code[0]);
               el->Add(mesh[ii]);
               el->Add(mesh[jj]);
               el->Add(mesh[kk]);
               mesh.Add(el);
            }
            break;

//       SubDomain
         case 6:
            k = ff.getI();
            for (size_t n=0; n<k; ++n) {
               ii = ff.getI();
               jj = ff.getI();
               kk = ff.getI();
            }
            break;

//       VertexOnGeometricVertex
         case 7:
            k = ff.getI();
            for (size_t n=0; n<k; ++n) {
               ii = ff.getI();
               ff.getD();
            }
            break;

//       VertexOnGeometricEdge
         case 8:
            k = ff.getI();
            for (size_t n=0; n<k; ++n) {
               ii = ff.getI();
               jj = ff.getI();
               ff.getD();
            }
            break;

//       EdgeOnGeometricEdge
         case  9:
            k = ff.getI();
            for (size_t n=0; n<k; ++n) {
               ii = ff.getI();
               jj = ff.getI();
            }
            break;

//       Identifier
         case 10:
            ww = ff.getS();
            break;

//       Geometry
         case 11:
            ww = ff.getS();
            break;

     }
   }
   while (key);
}


void getEasymesh(string file,
                 Mesh&  mesh,
                 size_t nb_dof)
{
   size_t      i, n, nb_nodes, nb_elements, nb_sides, ii, jj, ei, ej;
   int         mark;
   size_t      kk, ek, si, sj, sk, first_dof;
   real_t      xx, yy;
   Point<real_t> x;
   char        cc[10];
   int         code[MAX_NBDOF_NODE];
   Node        *nd;
   Element     *el;
   Side        *sd;
   ifstream    ef, nf, sf;

   nf.open((file+".n").c_str());
   if (nf.fail())
      throw OFELIException("getEasyMesh(...): File " + file + ".n not found.");
   ef.open((file+".e").c_str());
   if (ef.fail())
      throw OFELIException("getEasyMesh(...): File " + file + ".e not found.");
   sf.open((file+".s").c_str());
   if (sf.fail())
      throw OFELIException("getEasyMesh(...): File " + file + ".s not found.");
   i = 0;
   sf >> nb_sides;
   for (n=0; n<nb_sides; n++) {
      sf >> cc >> ii >> jj >> ei >> ej >> mark;
      if (mark)
         i++;
   }
   nb_sides = i;
   sf.close();
   sf.open((file+".s").c_str());

   mesh.setDim(2);
   first_dof = 1;
   nf >> nb_nodes;
   for (i=0; i<nb_nodes; i++) {
      nf >> cc >> xx >> yy >> mark;
      x = Point<real_t>(xx,yy);
      if (mark == 999)
         mark = 0;
      nd = new Node(i+1,x);
      nd->setNbDOF(nb_dof);
      DOFCode(mark, nb_dof, code);
      nd->setDOF(first_dof,nb_dof);
      nd->setCode(code);
      mesh.Add(nd);
   }

   ef >> nb_elements;
   for (n=0; n<nb_elements; n++) {
      ef >> cc >> ii >> jj >> kk >> ei >> ej >> ek >> si >> sj >> sk;
      ef >> xx >> yy >> mark;
      el = new Element(n+1,TRIANGLE,mark);
      el->Add(mesh[ii+1]);
      el->Add(mesh[jj+1]);
      el->Add(mesh[kk+1]);
      mesh.Add(el);
   }

   sf >> i;
   i = 0;
   for (n=0; n<nb_sides; n++) {
      sf >> cc >> ii >> jj >> ei >> ej >> mark;
      if (mark == 999)
         mark = 0;
      if (mark) {
         sd = new Side(++i,LINE);
         sd->Add(mesh[ii+1]);
         sd->Add(mesh[jj+1]);
         code[0] = mark;
         sd->setNbDOF(nb_dof);
         sd->setCode(1,mark);
         mesh.Add(sd);
      }
   }
   nf.close(); ef.close(); sf.close();
}


void getGambit(string file,
               Mesh&  mesh,
               size_t nb_dof)
{
   int mark;
   size_t first_dof;
   real_t xx;
   string line;
   int code[MAX_NBDOF_NODE];
   string sh[10];
   sh[LINE] = "line";
   sh[TRIANGLE] = "tria";
   sh[QUADRILATERAL] = "quad";
   sh[TETRAHEDRON] = "tetra";
   sh[HEXAHEDRON] = "hexa";
   sh[PENTAHEDRON] = "penta";

   ifstream mf(file.c_str());
   if (mf.fail())
      throw OFELIException("getGambit(...): File " + file + " not found.");

   for (size_t i=0; i<6; i++)
      getline(mf,line);
   int nb_nodes, nb_elements, nb_el_gr, dim, n1, n2, n3, n4, n5;
   getline(mf,line);
   istringstream s(line);
   s >> nb_nodes >> nb_elements >> nb_el_gr >> n2 >> dim >> n3;
   mesh.setDim(dim);
   getline(mf,line);
   getline(mf,line);

// Nodes
   for (int i=0; i<nb_nodes; i++) {
      Point<real_t> x;
      getline(mf,line);
      istringstream s(line);
      s >> n1 >> xx;
      x.x = xx;
      if (n3>1) {
         s >> xx;
         x.y = xx;
      }
      if (n3>2) {
         s >> xx;
         x.z = xx;
      }
      theNode = new Node(n1,x);
      TheNode.setNbDOF(nb_dof);
      mesh.Add(theNode);
   }
   getline(mf,line);
   cout << "Number of nodes: " << mesh.getNbNodes() << endl;
   cout << "Number of dof per node: " << nb_dof << endl;

// Elements
   getline(mf,line);
   int shape=0;
   for (int i=0; i<nb_elements; i++) {
      mf >> n1 >> n2 >> n3;
      if (n2==1)
         shape = LINE;
      else if (n2==2)
         shape = QUADRILATERAL;
      else if (n2==3)
         shape = TRIANGLE;
      else if (n2==4)
         shape = HEXAHEDRON;
      else if (n2==5)
         shape = PENTAHEDRON;
      else if (n2==6)
         shape = TETRAHEDRON;
      theElement = new Element(n1,shape);
      for (int j=1; j<=n3; j++) {
         mf >> n4;
         TheElement.Add(mesh[n4]);
      }
      mesh.Add(theElement);
   }
   getline(mf,line);
   cout << "Number of elements: " << mesh.getNbElements() << endl;

// Material properties
   string str;
   size_t mc = 0;
   for (int i=0; i<nb_el_gr; i++) {
      getline(mf,line);
      s >> str >> n1; s >> str >> n2;
      mf >> str >> n3; mf >> str >> n4;
      mf >> str;
      ++mc;
      mf >> n4;
      for (int j=1; j<=n2; j++) {
         mf >> n1;
         mesh(n1)->setCode(mc);
      }
      cout << "Material '" << str << "' is assigned the material code " << mc << endl;
      getline(mf,line);
      getline(mf,line);
   }
   cout << "Number of materials: " << mc << endl;

// Boundary conditions
   size_t nn1=0, nn2=0, nn3=0, nn4=0;
   size_t sd_label = 0;
   while (getline(mf,line)) {
     mf >> mark >> n1 >> n2 >> n3 >> n4;
     DOFCode(mark,nb_dof,code);
     first_dof = 1;
     for (int i=1; i<=n2; i++) {
        mf >> n5;
        if (n1 == 1) {
           mf >> n3 >> n4;
           theElement = mesh(n5);
           if (TheElement.getShape()==TRIANGLE && TheElement.getNbNodes()==3) {
              theSide = new Side(++sd_label,LINE);
              TheSide.Add(mesh[TheElement(n4    )->n()]);
              TheSide.Add(mesh[TheElement(n4%3+1)->n()]);
           }
           else if (TheElement.getShape()==QUADRILATERAL && TheElement.getNbNodes()==4) {
              theSide = new Side(++sd_label,LINE);
              TheSide.Add(mesh[TheElement(n4    )->n()]);
              TheSide.Add(mesh[TheElement(n4%4+1)->n()]);
           }
           else if (TheElement.getShape()==TETRAHEDRON && TheElement.getNbNodes()==4) {
              if (n4==1)
                 nn1 = 2, nn2 = 1, nn3 = 3;
              else if (n4==2)
                 nn1 = 1, nn2 = 2, nn3 = 4;
              else if (n4==3)
                 nn1 = 2, nn2 = 3, nn3 = 4;
              else if (n4==4)
                 nn1 = 3, nn2 = 1, nn3 = 4;
              theSide = new Side(++sd_label,TRIANGLE);
              TheSide.Add(mesh[TheElement(nn1)->n()]);
              TheSide.Add(mesh[TheElement(nn2)->n()]);
              TheSide.Add(mesh[TheElement(nn3)->n()]);
           }
           else if (TheElement.getShape()==HEXAHEDRON && TheElement.getNbNodes()==8) {
              if (n4==1)
                 nn1=1, nn2=2, nn3=6, nn4=5;
              else if (n4==2)
                 nn1=2, nn2=4, nn3=8, nn4=6;
              else if (n4==3)
                 nn1=4, nn2=3, nn3=7, nn4=8;
              else if (n4==4)
                 nn1=3, nn2=1, nn3=5, nn4=7;
              else if (n4==5)
                 nn1=2, nn2=1, nn3=3, nn4=4;
              else if (n4==6)
                 nn1=5, nn2=6, nn3=8, nn4=7;
              theSide = new Side(++sd_label,QUADRILATERAL);
              TheSide.Add(mesh[TheElement(nn1)->n()]);
              TheSide.Add(mesh[TheElement(nn2)->n()]);
              TheSide.Add(mesh[TheElement(nn3)->n()]);
              TheSide.Add(mesh[TheElement(nn4)->n()]);
           }
           else
              throw OFELIException("getGambit(...): Element shape " + sh[TheElement.getShape()] +
                                   "incompatible with " + itos(TheElement.getNbNodes()) + "nodes.");
           TheSide.setNbDOF(nb_dof);
           for (size_t j=1; j<=nb_dof; j++)
              TheSide.setCode(j,code[j-1]);
           mesh.Add(theSide);
        }
        else {
           theNode = mesh[n5];
           TheNode.setDOF(first_dof,nb_dof);
           TheNode.setCode(code);
        }
     }
     getline(mf,line);
   }
   cout << "Number of marked sides: " << mesh.getNbSides() << endl;
   mf.close();
}

/*
void getGmsh(string file,
             Mesh&  mesh,
             size_t nb_dof)
{
   string         kw, w;
   size_t         nb_nodes, nb_elements, nb_el, n1, n2, n3, n4, n5, k, i, j, ln=0, dim=2;
   real_t         x1, x2, x3, a=0;
   Node           *nd;
   Element        *el;
   Side           *sd;
   Point<real_t>  x[10];
   int            code[MAX_NBDOF_NODE];
   int            format=1;
   vector<El>     elem;
   vector<int>    tag(5);
   vector<size_t> num, un;
   vector<Nd>     nod;
   vector<size_t> phys_d, phys_n;
   vector<string> phys_name;
   size_t nn[] = {0, 2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, 1, 8, 20, 15, 13};

// Read in file
   ifstream pf;
   pf.open(file.c_str());
   try {
      if (pf.fail())
         THROW_RT(" File " + file + " not found.");
   }
   CATCH("getGmsh(...):");

   while (pf.eof()==0) {
      pf >> kw, ln++;

      if (kw=="$MeshFormat") {
         pf >> w >> n1 >> n2, ln++;
         pf >> kw, ln++;
         format = 2;
      }

      else if (kw=="$PhysicalNames" && format==2) {
         pf >> n1, ln++;
         string s;
         phys_d.resize(n1); phys_n.resize(n1);
         phys_name.resize(n1);
         for (size_t i=0; i<n1; i++) {
            pf >> n2 >> n3 >> s, ln++;
            phys_d[i] = n2;
            phys_n[i] = n3;
            phys_name[i] = s;
            ln++;
         }
         pf >> kw, ln++;
      }

      else if (kw=="$NOD" || kw=="$Nodes") {
         pf >> nb_nodes, ln++;
         nod.resize(nb_nodes);
         un.resize(nb_nodes);
         un = 0;
         num.resize(2*nb_nodes);
         for (i=0; i<nb_nodes; i++) {
            pf >> j >> x1 >> x2 >> x3, ln++;
            if (i==0)
               a = x3;
            num[j-1] = i+1;
            nod[i].label = j;
            nod[i].x[0] = x1; nod[i].x[1] = x2; nod[i].x[2] = x3;
            if (a != x3)
               dim = 3;
         }
         pf >> kw, ln++;
      }

      else if (kw=="$ELEM" && format==1) {
         pf >> nb_elements, ln++;
         elem.resize(nb_elements);
         for (j=0; j<nb_elements; j++) {
            pf >> n1 >> n2 >> n3 >> n4 >> n5, ln++;
            elem[j].label = n1;
            elem[j].type = n2;
            elem[j].region = n3;
            elem[j].nb_nodes = n5;
            for (k=0; k<n5; k++) {
               pf >> elem[j].node[k], ln++;
               if (elem[j].type != 15)
                  un[elem[j].node[k]-1]++;
            }
         }
      }

      else if (kw=="$Elements" && format==2) {
         pf >> nb_el, ln++;
         nb_elements = 0;
         elem.resize(nb_el);
         for (j=0; j<nb_el; j++) {
            bool ok = false;
            pf >> n1 >> n2 >> n3, ln++;
            for (k=0; k<n3; k++)
               pf >> tag[k];
            n4 = 0;
            if (n2==1)
               n3 = 2;
            else if (n2==2)
               n3 = 3;
            else if (n3==2)
               n3 = 3;
            if (dim==1) {
               if (n2==1)
                  n3=2, nb_elements++, ok=true;
            }
            else if (dim==2) {
               if (n2==2)
                  n3=3, nb_elements++, ok=true;
               else if (n2==3)
                  n3=4, nb_elements++, ok=true;
               else
                  ;
            }
            else if (dim==3) {
               if (n2==4)
                  n3=4, nb_elements++, ok=true;
               else if (n2==5)
                  n3=8, nb_elements++, ok=true;
               else if (n2==6)
                  n3=6, nb_elements++, ok=true;
               else
                  ;
            }
            elem[j].label = n1; elem[j].type = n2; elem[j].region = tag[0];
            for (k=0; k<nn[n2]; k++) {
               pf >> n5, ln++;
               elem[j].node[k] = n5;
               if (n2 != 15)
                 //                  un[n5-1]++;
                  elem[j].node[k] = n5;
            }
         }
         pf >> kw, ln++;
      }
   }

   vector<size_t> NN(nb_nodes);
   NN = 0;
   size_t NbN = 0;
   for (size_t i=0; i<nb_nodes; i++) {
      size_t n = num[i];
      if (un[n-1])
         NN[n-1] = ++NbN;
   }

   vector<size_t> code_n;
   for (i=0; i<phys_d.size(); i++) {
      if (phys_name[i]=="\"NEUMANN\"")
         code_n.push_back(phys_n[i]);
   }

   mesh.setDim(dim);
   size_t first_dof = 1;
   for (j=0; j<nb_nodes; j++) {
      size_t n = NN[num[nod[j].label-1]-1];
      if (n>0) {
         nd = new Node(n,nod[j].x);
         nd->setNbDOF(nb_dof);
         nd->setDOF(first_dof,nb_dof);
         for (i=0; i<nb_dof; i++)
            code[i] = 0;
         nd->setCode(code);
         mesh.Add(nd);
      }
   }
   size_t ks=1, ke=1;
   for (j=0; j<nb_el; j++) {
      int cd = elem[j].region;
      for (i=0; i<code_n.size(); i++) {
         if (cd==int(code_n[i]))
            cd = 0;
      }
      if (dim==2) {
         if (elem[j].type==1) {
            sd = new Side(ks,LINE);
            sd->Add(mesh[NN[num[elem[j].node[0]-1]-1]]);
            sd->Add(mesh[NN[num[elem[j].node[1]-1]-1]]);
            //            if (cd) {
            //               for (i=0; i<sd->getNbNodes(); i++) {
            //                  nd = (*sd)(i+1);
            //                  DOFCode(cd,nb_dof,code);
            //                  for (k=0; k<nb_dof; k++)
            //                     nd->setCode(k+1,code[k]);
            //               }
//               }
            if (cd < 1000) {
               for (i=1; i<=sd->getNbNodes(); i++) {
                  nd = sd->getPtrNode(i);
                  DOFCode(cd,nb_dof,code);
                  for (k=0; k<nb_dof; k++)
                     nd->setCode(k+1,code[k]);
               }
               delete sd; sd = NULL;
            }
            else {
               ks++;
               sd->setNbDOF(nb_dof);
               DOFCode(cd-1000,nb_dof,code);
               for (k=0; k<nb_dof; k++)
                  sd->setCode(k+1,code[k]);
               mesh.Add(sd);
            }
            for (i=0; i<phys_d.size(); i++) {
               if (phys_d[i]==1 && phys_name[i]=="\"NEUMANN\"" && phys_n[i]==elem[j].region) {
                  ks++;
                  sd->setNbDOF(nb_dof);
                  DOFCode(phys_n[i],nb_dof,code);
                  for (k=0; k<nb_dof; k++)
                     sd->setCode(k+1,code[k]);
                  mesh.Add(sd);
               }
               if (phys_d[i]==1 && phys_name[i]=="\"DIRICHLET\"" && phys_n[i]==elem[j].region) {
                  DOFCode(phys_n[i],nb_dof,code);
                  nd = (*sd)(1);
                  for (k=0; k<nb_dof; k++) {
                     (*sd)(1)->setCode(k+1,code[k]);
                     (*sd)(2)->setCode(k+1,code[k]);
                  }
               }
            }
         }
         else if (elem[j].type==2) {
            el = new Element(ke++,TRIANGLE);
            x[0] = mesh[NN[num[elem[j].node[0]-1]-1]]->getCoord();
            x[1] = mesh[NN[num[elem[j].node[1]-1]-1]]->getCoord();
            x[2] = mesh[NN[num[elem[j].node[2]-1]-1]]->getCoord();
            el->Add(mesh[NN[num[elem[j].node[0]-1]-1]]);
            if ((x[1].x-x[0].x)*(x[2].y-x[0].y) - (x[1].y-x[0].y)*(x[2].x-x[0].x) < 0) {
               el->Add(mesh[NN[num[elem[j].node[2]-1]-1]]);
               el->Add(mesh[NN[num[elem[j].node[1]-1]-1]]);
            }
            else {
               el->Add(mesh[NN[num[elem[j].node[1]-1]-1]]);
               el->Add(mesh[NN[num[elem[j].node[2]-1]-1]]);
            }
            el->setCode(elem[j].region);
            mesh.Add(el);
         }
         else if (elem[j].type==3) {
            el = new Element(ke++,QUADRILATERAL);
            for (i=0; i<4; i++)
               x[i] = mesh[NN[num[elem[j].node[i]-1]-1]]->getCoord();
            for (i=0; i<4; i++)
               el->Add(mesh[NN[num[elem[j].node[i]-1]-1]]);
            for (i=0; i<phys_d.size(); i++)
               if (phys_d[i]==2 && phys_n[i]==elem[j].region)
                  el->setCode(elem[j].region);
            mesh.Add(el);
         }
      }
      else if (dim==3) {
         if (elem[j].type==2) {
            sd = new Side(ks++,TRIANGLE);
            for (i=0; i<3; i++)
               sd->Add(mesh[NN[num[elem[j].node[i]-1]-1]]);
            if (cd) {
               for (i=0; i<sd->getNbNodes(); i++) {
                  nd = (*sd)(i+1);
                  DOFCode(cd,nb_dof,code);
                  for (k=0; k<nb_dof; k++)
                     nd->setCode(k+1,code[k]);
               }
            }
         }
         if (elem[j].type==3) {
            sd = new Side(ks++,QUADRILATERAL);
            for (i=0; i<4; i++)
               sd->Add(mesh.getPtrNode(NN[num[elem[j].node[i]-1]-1]));
            if (cd) {
               for (i=0; i<sd->getNbNodes(); i++) {
                  nd = (*sd)(i+1);
                  DOFCode(cd,nb_dof,code);
                  for (k=0; k<nb_dof; k++)
                     nd->setCode(k+1,code[k]);
               }
            }
         }
         else if (elem[j].type==4) {
            el = new Element(ke++,TETRAHEDRON);
            for (i=0; i<4; i++)
               el->Add(mesh[NN[num[elem[j].node[i]-1]]]);
            for (i=0; i<phys_d.size(); i++)
               if (phys_d[i]==3 && phys_n[i]==elem[j].region)
                  el->setCode(elem[j].region);
            mesh.Add(el);
         }
         else if (elem[j].type==5) {
            el = new Element(ke++,HEXAHEDRON);
            for (i=0; i<8; i++)
               el->Add(mesh[NN[num[elem[j].node[i]-1]-1]]);
            for (i=0; i<phys_d.size(); i++)
               if (phys_d[i]==3 && phys_n[i]==elem[j].region)
                  el->setCode(elem[j].region);
            mesh.Add(el);
         }
         else if (elem[j].type==6) {
            el = new Element(ke++,PENTAHEDRON);
            for (i=0; i<6; i++)
               el->Add(mesh[NN[num[elem[j].node[i]-1]-1]]);
            for (i=0; i<phys_d.size(); i++)
               if (phys_d[i]==3 && phys_n[i]==elem[j].region)
                  el->setCode(elem[j].region);
            mesh.Add(el);
         }
         else if (elem[j].type==7) {
            el = new Element(ke++,PENTAHEDRON);
            for (i=0; i<5; i++)
               el->Add(mesh[NN[num[elem[j].node[i]-1]-1]]);
            for (i=0; i<phys_d.size(); i++)
               if (phys_d[i]==3 && phys_n[i]==elem[j].region)
                  el->setCode(elem[j].region);
            mesh.Add(el);
         }
      }
      if (elem[j].type==15) {
         size_t n = num[elem[j].label-1];
         if (NN[n-1]) {
            nd = mesh[NN[num[elem[j].label-1]-1]];
            for (size_t i=0; i<phys_d.size(); i++) {
               if (phys_d[i]==1 && phys_name[i]=="\"DIRICHLET\"" && phys_n[i]==elem[j].region)
                  DOFCode(phys_n[i],nb_dof,code);
            }
            for (k=0; k<nb_dof; k++)
               nd->setCode(k+1,code[k]);
         }
      }
   }
   pf.close();
}*/


void getGmsh(string file,
             Mesh&  mesh,
             size_t nb_dof)
{
   string        kw, w;
   int           nb_nodes, nb_elements, k;
   int           n1, n2, n3, n4, n5, n6, dim=0;
   real_t        x1, x2, x3, d1, d2, d3, d4, d5, d6;
   Point<real_t> x[3];
   Node          *nd;
   Element       *el;
   Side          *sd;
   int           code[MAX_NBDOF_NODE];
   Vect<El>      elem;
   Vect<int>     sdim;
   vector<int> nb_en(100), sh(100);
   nb_en[ 0] =  2; nb_en[ 1] =  3; nb_en[ 2] =  4; nb_en[ 3] =  4; nb_en[ 4] =  8; nb_en[ 5] =  6;
   nb_en[ 6] =  5; nb_en[ 7] =  3; nb_en[ 8] =  6; nb_en[ 9] =  9; nb_en[10] = 10; nb_en[11] = 27;
   nb_en[12] = 18; nb_en[13] = 14; nb_en[14] =  1; nb_en[15] =  8; nb_en[16] = 20; nb_en[17] = 15;
   sh[0] = sh[7] = LINE; sh[1] = sh[8] = TRIANGLE; sh[2] = sh[9] = sh[15] = QUADRILATERAL; 
   sh[3] = sh[10] = TETRAHEDRON; sh[4] = sh[11] = sh[16] = HEXAHEDRON; sh[5] = sh[12] = sh[17] = PRISM;
   sh[6] = sh[13] = PYRAMID; sh[14] = POINT;

// Read in file
   ifstream pf;
   pf.open(file.c_str());
   if (pf.fail())
      throw OFELIException("getGmsh(...): Can't open file "+file);
   
// Mesh Format
   pf >> kw;
   if (kw!="$MeshFormat")
      throw OFELIException("getGmsh(...): Keyword MeshFormat not found.");
   pf >> n1 >> n2 >> n3;
   pf >> kw;
   if (kw!="$EndMeshFormat")
      throw OFELIException("getGmsh(...): Keyword EndMeshFormat not found.");

// Physical name
   pf >> kw;
   if (kw=="$PhysicalName") {
      pf >> n1;
      for (int i=0; i<n1; ++i)
         pf >> n2 >> n3 >> w;
      pf >> kw;
      if (kw!="$EndPhysicalName")
         throw OFELIException("getGmsh(...): Keyword EndPhysicalName not found.");
   }

// Entities
   if (kw=="$Entities") {
      pf >> kw;
      while (kw!="$EndEntities")
         getline(pf,kw);
   }
/*   if (kw!="$Entities")
      throw OFELIException("getGmsh(...): Keyword Entities unfound.");
   pf >> n1 >> n2 >> n3 >> n4;
   for (int i=0; i<n1; ++i) {
      pf >> n5 >> d1 >> d2 >> d3 >> d4 >> d5 >> d6;
      pf >> n6;
      for (size_t j=0; j<n6; ++j)
         pf >> n6;
   }
   for (int i=0; i<n2; ++i) {
      pf >> n5 >> d1 >> d2 >> d3 >> d4 >> d5 >> d6;
      pf >> n6;
      for (int j=0; j<n6; ++j)
         pf >> n6;
      pf >> n6;
      for (int j=0; j<n6; ++j)
         pf >> n6;
   }
   for (int i=0; i<n3; ++i) {
      pf >> n5 >> d1 >> d2 >> d3 >> d4 >> d5 >> d6;
      pf >> n6;
      for (int j=0; j<n6; ++j)
         pf >> n6;
      pf >> n6;
      for (int j=0; j<n6; ++j)
         pf >> n6;
   }
   for (int i=0; i<n4; ++i) {
      pf >> n5 >> d1 >> d2 >> d3 >> d4 >> d5 >> d6;
      pf >> n6;
      for (int j=0; j<n6; ++j)
         pf >> n6;
      pf >> n6;
      for (int j=0; j<n6; ++j)
         pf >> n6;
   }
   pf >> kw;
   if (kw!="$EndEntities")
      throw OFELIException("getGmsh(...): Keyword EndEntities unfound.");*/

// Partitioned Entities
   pf >> kw;
   if (kw=="$PartitionedEntities") {
      pf >> n1;
      pf >> n2;
      for (int i=0; i<n2; ++i)
         pf >> n3 >> n4;
      pf >> n1 >> n2 >> n3 >> n4;
      for (int i=0; i<n1; ++i) {
         pf >> n5 >> d1 >> d2 >> d3 >> d4 >> d5 >> d6;
         pf >> n6;
         for (int j=0; j<n6; ++j)
            pf >> n6;
      }
      if (kw!="$EndPartitionedEntities")
         throw OFELIException("getGmsh(...): Keyword EndEntities not found.");
   }

   if (kw!="$Nodes")
      throw OFELIException("getGmsh(...): Keyword Nodes not found.");
   pf >> n1 >> nb_nodes;
   vector<Nd> nod(nb_nodes);
   for (int i=1; i<=n1; ++i) {
      pf >> n2 >> n3 >> n4 >> n5;
      if (n4)
         throw OFELIException("getGmsh(...): parametric not implemented.");
      for (int j=0; j<n5; ++j) {
         pf >> k >> x1 >> x2 >> x3;
         nod[k-1].n = k;
         nod[k-1].x[0] = x1, nod[k-1].x[1] = x2, nod[k-1].x[2] = x3;
      }
   }
   pf >> kw;
   if (kw!="$EndNodes")
      throw OFELIException("getGmsh(...): Keyword EndNodes not found.");

   pf >> kw;
   if (kw!="$Elements")
      throw OFELIException("getGmsh(...): Keyword Elements not found.");
   pf >> n1 >> nb_elements;
   vector<El> es, elements, sides;
   El ell;
   ell.region = 1;
   for (int i=0; i<n1; ++i) {
      pf >> n2 >> n3 >> n4 >> n5;
      dim = std::max(dim,n3);
      ell.dim = n3;
      for (int n=0; n<n5; ++n) {
         pf >> ell.n;
         ell.shape = sh[n4-1];
         ell.nb_nodes = nb_en[n4-1];
         for (int i=0; i<ell.nb_nodes; ++i)
             pf >> ell.node[i];
         es.push_back(ell);
      }
   }
   pf >> kw;
   if (kw!="$EndElements")
      throw OFELIException("getGmsh(...): Illegal file format");
   for (size_t n=0; n<nb_elements; ++n) {
      if (dim==es[n].dim)
         elements.push_back(es[n]);
      else if (dim==es[n].dim+1)
         sides.push_back(es[n]);
   }

// Build a mesh instance
   mesh.setDim(dim);
   size_t first_dof = 1;
   for (int j=0; j<nb_nodes; j++) {
      nd = new Node(nod[j].n,nod[j].x);
      nd->setNbDOF(nb_dof);
      nd->setDOF(first_dof,nb_dof);
      for (int i=0; i<nb_dof; i++)
         code[i] = 0;
      nd->setCode(code);
      mesh.Add(nd);
   }
   for (int j=0; j<elements.size(); j++) {
      ell = elements[j];
      el = new Element(ell.n,ell.shape,ell.region);
      for (int k=0; k<ell.nb_nodes; ++k)
         el->Add(mesh[ell.node[k]]);
      mesh.Add(el);    
   }
   for (int j=0; j<sides.size(); j++) {
      ell = sides[j];
      sd = new Side(ell.n,ell.shape);
      for (int k=0; k<ell.nb_nodes; ++k)
         sd->Add(mesh[ell.node[k]]);
      DOFCode(ell.region-1000,nb_dof,code);
      for (int k=0; k<nb_dof; k++)
         sd->setCode(k+1,code[k]);
      mesh.Add(sd);    
   }
   pf.close();
   cout << "Generated mesh:\n";
   cout << "Number of nodes:    " << nb_nodes << endl;
   cout << "Number of elements: " << elements.size() << endl;
   cout << "Number of sides:    " << sides.size() << endl;
}


void getMatlab(string file,
               Mesh&  mesh,
               size_t nb_dof)
{
   size_t i, j, k, nb_nodes, nb_elements, nb_sides, nbc;
   int code[MAX_NBDOF_NODE];
   ifstream mf;
   mf.open(file.c_str());
   if (mf.fail())
      throw OFELIException("getMatlab(...): File "+file+" not found.");

   mf >> nb_nodes;
   Vect<Point<real_t> > x(nb_nodes);
   Vect<int> c(nb_nodes);
   for (i=0; i<nb_nodes; i++)
      mf >> x[i].x >> x[i].y;

   mf >> nb_sides;
   DMatrix<float> e(7,nb_sides);
   Vect<size_t> cs(nb_sides);
   for (i=1; i<=nb_sides; i++)
      for (j=1; j<=7; j++)
         mf >> e(j,i);

   mf >> nb_elements;
   DMatrix<size_t> t(4,nb_elements);
   for (i=1; i<=nb_elements; i++)
      for (j=1; j<=4; j++)
         mf >> t(j,i);

   mf >> nbc;
   size_t N, M, l;
   string g;
   for (i=1; i<=nbc; i++) {
      mf >> N >> M;
      for (j=1; j<=N*N; j++)
         mf >> l;
      for (j=1; j<=N; j++)
         mf >> l;
      for (j=1; j<=M*N; j++)
         mf >> l;
      for (j=1; j<=M; j++)
         mf >> l;
      for (j=1; j<=N*N; j++)
         mf >> g;
      for (j=1; j<=N; j++) {
         mf >> g;
         for (k=1; k<=nb_sides; k++) {
            if (e(5,k)==i && M==0)
               cs[k-1] = stringTo<int>(g);
         }
      }
      for (j=1; j<=M*N; j++)
         mf >> g;
      for (j=1; j<=M; j++) {
         mf >> g;
         for (k=1; k<=nb_sides; k++) {
            if (e(5,k)==i)
               c[size_t(e(1,k))-1] = c[size_t(e(2,k))-1] = stringTo<int>(g);
         }
      }
   }

   mesh.setDim(2);
   size_t first_dof = 1;
   for (i=0; i<nb_nodes; i++) {
      theNode = new Node(i+1,x[i]);
      TheNode.setDOF(first_dof,nb_dof);
      DOFCode(c[i],nb_dof,code);
      TheNode.setCode(code);
      mesh.Add(theNode);
   }
   for (i=1; i<=nb_elements; i++) {
      theElement = new Element(i,TRIANGLE,t(4,i));
      for (j=1; j<=3; j++)
         TheElement.Add(mesh[t(j,i)]);
         mesh.Add(theElement);
   }
   for (i=1; i<=nb_sides; i++) {
      theSide = new Side(i,"line");
      TheSide.Add(mesh[size_t(e(1,i))]);
      TheSide.Add(mesh[size_t(e(2,i))]);
      DOFCode(cs[i-1],nb_dof,code);
      TheSide.setNbDOF(nb_dof);
      for (size_t j=0; j<nb_dof; j++)
         TheSide.setCode(j+1,code[j]);
      mesh.Add(theSide);
   }
}


void getNetgen(string file,
               Mesh&  mesh,
               size_t nb_dof)
{
   size_t      nb_nodes, nb_se, nb_edges, nb_elements, nb_sides;
   size_t      jj, i, k, l, m, n, dim, a, nn[10];
   int         j;
   Node        *nd;
   Element     *el;
   Side        *sd;
   int         code[MAX_NBDOF_NODE];
   ifstream    nf;

   nf.open((file+".vol").c_str());
   if (nf.fail())
      throw OFELIException("getNetgen(...): File " + file + ".vol not found.");

   string ch, cc;
   std::getline(nf,ch);
   std::getline(nf,ch);
   nf >> dim;
   do {
     std::getline(nf,ch);
     cc = string(ch,15);
   }
   while (cc!="surfaceelements");
   vector<Nd> nod(MAX_NB_NODES);
   Nd nnd;
   for (i=0; i<MAX_NB_NODES; i++)
      nod[i].cc = 0;
   nf >> nb_se;
   nb_sides = 0;
   vector<El> side;
   El ssd;
   for (a=0; a<nb_se; a++) {
      nf >> i >> j >> k >> l >> n;
      nf >> nn[0] >> nn[1] >> nn[2];
      nf >> i >> k >> l;
      if (j < 10000)
         nod[nn[0]-1].cc = nod[nn[1]-1].cc = nod[nn[2]-1].cc = j;
      else {
         j -= 10000;
         ssd.cc = j;
         side.push_back(ssd);
         nb_sides++;
      }
   }

   do {
      getline(nf,ch);
      cc = string(ch,14);
   }
   while (cc!="volumeelements");
   nf >> nb_elements;
   vector<El> elem(nb_elements);
   for (jj=0; jj<nb_elements; jj++) {
      nf >> i >> k;
      nf >> nn[0] >> nn[1] >> nn[2] >> nn[3];
      elem[jj].n = jj+1;
      elem[jj].shape = TETRAHEDRON;
      elem[jj].region = i;
      elem[jj].nb_nodes = 4;
      elem[jj].node[0] = nn[0];
      elem[jj].node[1] = nn[1];
      elem[jj].node[2] = nn[2];
      elem[jj].node[3] = nn[3];
   }

   do {
      getline(nf,ch);
      cc = string(ch,12);
   }
   while (cc!="edgesegments");
   nf >> nb_edges;
   for (a=0; a<nb_edges; a++) {
      nf >> i >> j >> k >> l >> m >> n;
      nf >> i >> j >> k >> l >> m >> n;
   }

   do {
      std::getline(nf,ch);
      cc = string(ch,MAX_NBDOF_NODE);
   }
   while (cc!="points");
   nf >> nb_nodes;
   for (i=0; i<nb_nodes; i++) {
      nf >> nod[i].x[0] >> nod[i].x[1] >> nod[i].x[2];
      nod[i].n = i+1;
      DOFCode(nod[i].cc,nb_dof,code);
      for (k=0; k<nb_dof; k++)
         nod[i].code[k] = code[k];
   }

   mesh.setDim(dim);
   size_t first_dof = 1;
   for (jj=0; jj<nb_nodes; jj++) {
      nd = new Node(nod[jj].n,nod[jj].x);
      nd->setDOF(first_dof,nb_dof);
      nd->setCode(nod[jj].code);
      mesh.Add(nd);
   }

   size_t ks=1, ke=1;
   for (jj=0; jj<nb_elements; jj++) {
      int cd = elem[jj].region;

      if (elem[jj].shape==POINT) {
        nd = mesh[elem[jj].n];
        DOFCode(cd,nb_dof,code);
        for (k=0; k<nb_dof; k++)
           nd->setCode(k+1,code[k]);
      }

      if (dim==2) {
         if (elem[jj].shape==LINE) {
            sd = new Side(ks,LINE);
            for (i=0; i<2; i++)
               sd->Add(mesh[elem[jj].node[i]]);
            if (cd%2 == 1) {
               for (i=0; i<sd->getNbNodes(); i++) {
	          nd = (*sd)(i+1);
                  DOFCode(cd,nb_dof,code);
                  for (k=0; k<nb_dof; k++)
                     nd->setCode(k+1,code[k]);
               }
               delete sd;
            }
            else {
               ks++;
               sd->setNbDOF(nb_dof);
               DOFCode(cd,nb_dof,code);
               for (k=0; k<nb_dof; k++)
                  sd->setCode(k+1,code[k]);
               mesh.Add(sd);
            }
         }
         else if (elem[jj].shape==TRIANGLE) {
            el = new Element(ke++,TRIANGLE,cd);
            for (i=0; i<3; i++)
               el->Add(mesh[elem[jj].node[i]]);
            mesh.Add(el);
         }
         else if (elem[jj].shape==QUADRILATERAL) {
            el = new Element(ke++,QUADRILATERAL,cd);
            for (i=0; i<4; i++)
               el->Add(mesh[elem[jj].node[i]]);
            mesh.Add(el);
         }
      }

      else if (dim==3) {
         if (elem[jj].shape==TRIANGLE) {
            sd = new Side(ks++,TRIANGLE);
            for (i=0; i<3; i++)
               sd->Add(mesh[elem[jj].node[i]]);
            if (cd%2 == 1) {
               for (i=0; i<sd->getNbNodes(); i++) {
                  nd = (*sd)(i+1);
                  DOFCode(cd,nb_dof,code);
                  for (k=0; k<nb_dof; k++)
                     nd->setCode(k+1,code[k]);
               }
            }
            else {
               sd->setNbDOF(nb_dof);
               DOFCode(cd,nb_dof,code);
               for (k=0; k<nb_dof; k++)
                  sd->setCode(k+1,code[k]);
               mesh.Add(sd);
            }
         }
         if (elem[jj].shape==QUADRILATERAL) {
            sd = new Side(ks++,QUADRILATERAL);
            for (i=0; i<4; i++)
               sd->Add(mesh[elem[jj].node[i]]);
            if (cd%2 == 1) {
               for (i=0; i<sd->getNbNodes(); i++) {
                  nd = (*sd)(i+1);
                  DOFCode(cd,nb_dof,code);
                  for (k=0; k<nb_dof; k++)
                     nd->setCode(k+1,code[k]);
               }
            }
            else {
               sd->setNbDOF(nb_dof);
               DOFCode(cd,nb_dof,code);
               for (k=0; k<nb_dof; k++)
                  sd->setCode(k+1,code[k]);
               mesh.Add(sd);
            }
         }
         else if (elem[jj].shape==TETRAHEDRON) {
            el = new Element(ke++,TETRAHEDRON,cd);
            for (i=0; i<4; i++)
               el->Add(mesh[elem[jj].node[i]]);
            mesh.Add(el);
         }
         else if (elem[jj].shape==HEXAHEDRON) {
            el = new Element(ke++,HEXAHEDRON,cd);
            for (i=0; i<8; i++)
               el->Add(mesh[elem[jj].node[i]]);
            mesh.Add(el);
         }
         else if (elem[jj].shape==PRISM) {
            el = new Element(ke++,PRISM,cd);
            for (i=0; i<6; i++)
               el->Add(mesh[elem[jj].node[i]]);
            mesh.Add(el);
         }
         else if (elem[jj].shape==7) {
            el = new Element(ke++,PYRAMID,cd);
            for (i=0; i<5; i++)
               el->Add(mesh[elem[jj].node[i]]);
            mesh.Add(el);
         }
      }
   }

   for (jj=0; jj<nb_sides; jj++) {
       sd = new Side(jj+1,TRIANGLE);
       for (i=0; i<3; i++)
          sd->Add(mesh[side[jj].node[i]]);
       mesh.Add(sd);
   }
   nf.close();
}


void getTetgen(string file,
               Mesh&  mesh,
               size_t nb_dof)
{
   size_t nb_nodes, nb_elements, n1, n2, n3, n4, k, i, j, dim;
   real_t x1, x2, x3;
   int    code[MAX_NBDOF_NODE];

// Read node file
   ifstream nf((file+".node").c_str());
   if (nf.fail())
      throw OFELIException("getTetgen(...): File "+file+".node not found.");
   nf >> nb_nodes >> dim >> i >> j;
   vector<Nd> nod(nb_nodes);
   vector<size_t> num(nb_nodes);
   for (i=0; i<nb_nodes; i++) {
      nf >> j >> x1 >> x2 >> x3;
      num[j-1] = j;
      nod[i].n = j;
      nod[i].x[0] = x1;
      nod[i].x[1] = x2;
      nod[i].x[2] = x3;
      for (size_t k=0; k<nb_dof; k++)
         nod[i].code[k] = 0;
   }
   nf.close();

// Read element file
   ifstream ef((file+".ele").c_str());
   if (ef.fail())
      throw OFELIException("getTetgen(...): File " + file + ".ele not found.");
   ef >> nb_elements >> k >> i;
   vector<El> elem(nb_elements);
   for (j=0; j<nb_elements; j++) {
      ef >> i >> n1 >> n2 >> n3 >> n4;
      elem[j].n = i;
      if (k==4)
         elem[j].shape = QUADRILATERAL;
      elem[j].region = 1; elem[j].nb_nodes = k;
      elem[j].node[0] = n1; elem[j].node[1] = n2;
      elem[j].node[2] = n3; elem[j].node[3] = n4;
   }

   mesh.setDim(dim);
   size_t first_dof = 1;
   for (j=0; j<nb_nodes; j++) {
      theNode = new Node(num[nod[j].n-1],nod[j].x);
      for (k=0; k<nb_dof; k++)
         code[k] = nod[j].code[k];
      TheNode.setDOF(first_dof,nb_dof);
      TheNode.setCode(code);
      mesh.Add(theNode);
   }

   size_t ks=1, ke=1;
   for (j=0; j<nb_elements; j++) {
      int cd = elem[j].region;

      if (elem[j].shape==POINT) {
         theNode = mesh[num[elem[j].n-1]];
         DOFCode(cd,nb_dof,code);
         for (k=0; k<nb_dof; k++)
            TheNode.setCode(k+1,code[k]);
      }

      if (dim==2) {
         if (elem[j].shape==LINE) {
            theSide = new Side(ks,LINE);
            for (i=0; i<2; i++)
               TheSide.Add(mesh[num[elem[j].node[i]-1]]);
            if (cd%2 == 1) {
               for (i=0; i<TheSide.getNbNodes(); i++) {
                  theNode = TheSide(i+1);
                  DOFCode(cd,nb_dof,code);
                  for (k=0; k<nb_dof; k++)
                     TheNode.setCode(k+1,code[k]);
               }
               delete theSide;
            }
            else {
               ks++;
               TheSide.setNbDOF(nb_dof);
               DOFCode(cd,nb_dof,code);
               for (k=0; k<nb_dof; k++)
                  TheSide.setCode(k+1,code[k]);
               mesh.Add(theSide);
            }
         }
         else if (elem[j].shape==TRIANGLE) {
            theElement = new Element(ke++,TRIANGLE,cd);
            for (i=0; i<3; i++)
               TheElement.Add(mesh[num[elem[j].node[i]-1]]);
            mesh.Add(theElement);
         }
         else if (elem[j].shape==QUADRILATERAL) {
            theElement = new Element(ke++,QUADRILATERAL,cd);
            for (i=0; i<4; i++)
               TheElement.Add(mesh[num[elem[j].node[i]-1]]);
            mesh.Add(theElement);
         }
      }

      else if (dim==3) {
         if (elem[j].shape==TRIANGLE) {
            theSide = new Side(ks++,TRIANGLE);
            for (i=0; i<3; i++)
               TheSide.Add(mesh[num[elem[j].node[i]-1]]);
            if (cd%2 == 1) {
               for (i=0; i<TheSide.getNbNodes(); i++) {
                  theNode = TheSide(i+1);
                  DOFCode(cd,nb_dof,code);
                  for (k=0; k<nb_dof; k++)
                     TheNode.setCode(k+1,code[k]);
               }
            }
            else {
               TheSide.setNbDOF(nb_dof);
               DOFCode(cd,nb_dof,code);
               for (k=0; k<nb_dof; k++)
                  TheSide.setCode(k+1,code[k]);
               mesh.Add(theSide);
            }
         }
         if (elem[j].shape==QUADRILATERAL) {
            theSide = new Side(ks++,QUADRILATERAL);
            for (i=0; i<4; i++)
               TheSide.Add(mesh[num[elem[j].node[i]-1]]);
            if (cd%2 == 1) {
               for (i=0; i<TheSide.getNbNodes(); i++) {
                  theNode = TheSide(i+1);
                  DOFCode(cd,nb_dof,code);
                  for (k=0; k<nb_dof; k++)
                     TheNode.setCode(k+1,code[k]);
               }
            }
            else {
               TheSide.setNbDOF(nb_dof);
               DOFCode(cd,nb_dof,code);
               for (k=0; k<nb_dof; k++)
                  TheSide.setCode(k+1,code[k]);
               mesh.Add(theSide);
            }
         }
         else if (elem[j].shape==TETRAHEDRON) {
            theElement = new Element(ke++,TETRAHEDRON,cd);
            for (i=0; i<4; i++)
               TheElement.Add(mesh[num[elem[j].node[i]-1]]);
            mesh.Add(theElement);
         }
         else if (elem[j].shape==HEXAHEDRON) {
            theElement = new Element(ke++,HEXAHEDRON,cd);
            for (i=0; i<8; i++)
               TheElement.Add(mesh[num[elem[j].node[i]-1]]);
            mesh.Add(theElement);
         }
         else if (elem[j].shape==PRISM) {
            theElement = new Element(ke++,PRISM,cd);
            for (i=0; i<6; i++)
               TheElement.Add(mesh[num[elem[j].node[i]-1]]);
            mesh.Add(theElement);
         }
         else if (elem[j].shape==PYRAMID) {
            theElement = new Element(ke++,PYRAMID,cd);
            for (i=0; i<5; i++)
               TheElement.Add(mesh[num[elem[j].node[i]-1]]);
            mesh.Add(theElement);
         }
      }
   }
   ef.close();
}


void getTriangle(string file,
                 Mesh&  mesh,
                 size_t nb_dof)
{
   size_t ne[10];
   int code[MAX_NBDOF_NODE];
   ifstream nf, ef, sf;

   nf.open((file+".node").c_str());
   if (nf.fail())
      throw OFELIException("getTriangle(...): File " + file + ".node not found.");
   ef.open((file+".ele").c_str());
   if (ef.fail())
      throw OFELIException("getTriangle(...): File " + file + ".ele not found.");

   mesh.setDim(2);
   size_t first_dof = 1;
   size_t nz=0, nb_nodes, n, na, nm;
   nf >> nb_nodes >> n >> na >> nm;
   for (size_t i=0; i<nb_nodes; i++) {
      nf >> n;
      if (i==0 && n==0)
         nz = 1;
      real_t xx, yy;
      nf >> xx >> yy;
      real_t attr;
      for (size_t j=0; j<na; j++)
         nf >> attr;
      int mark = 0;
      if (nm>0)
         nf >> mark;
      if (nz==1)
         n++;
      theNode = new Node(n,Point<real_t>(xx,yy));
      TheNode.setNbDOF(nb_dof);
      DOFCode(mark,nb_dof,code);
      TheNode.setDOF(first_dof,nb_dof);
      TheNode.setCode(code);
      mesh.Add(theNode);
   }
   size_t nb_elements, nbn;
   ef >> nb_elements >> nbn >> na;
   nz = 0;
   for (size_t i=0; i<nb_elements; i++) {
      ef >> n;
      if (i==0 && n==0)
         nz = 1;
      int mark = 0;
      for (size_t j=0; j<nbn; j++) {
         ef >> mark;
         ne[j] = mark;
      }
      real_t attr;
      for (size_t j=0; j<na; j++)
         ef >> attr;
      mark = int(attr);
      if (nz>0)
         n++;
      theElement = new Element(n,TRIANGLE,mark);
      for (size_t j=0; j<nbn; j++)
         TheElement.Add(mesh[ne[j]]);
      mesh.Add(theElement);
   }

   sf.open((file+".edge").c_str());
   if (!sf.fail()) {
      size_t nz=0, nb_sides;
      sf >> nb_sides >> nm;
      for (size_t i=0; i<nb_sides; i++) {
         sf >> n >> ne[0] >> ne[1];
         if (i==0 && n==0)
            nz = 1;
         if (nz>0)
            n++;
         real_t attr;
         for (size_t j=0; j<nm; j++)
            sf >> attr;
         int mark = int(attr);
         DOFCode(mark, nb_dof, code);
         theSide = new Side(n,LINE);
         TheSide.Add(mesh[ne[0]]);
         TheSide.Add(mesh[ne[1]]);
         TheSide.setNbDOF(nb_dof);
         for (size_t j=0; j<nb_dof; j++)
            TheSide.setCode(j+1,code[j]);
         mesh.Add(theSide);
      }
   }
   nf.close();
   ef.close();
   if (!sf.fail())
      sf.close();
}


} /* namespace OFELI */
