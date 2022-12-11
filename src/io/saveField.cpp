/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani

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

                Functions to convert XML field files to other formats

  ==============================================================================*/

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
using std::to_string;
using std::setw;
using std::ios;
using std::setprecision;

#include "OFELI_Config.h"
#include "mesh/MeshUtil.h"
#include "mesh/Mesh.h"
#include "io/saveField.h"
#include "io/IOField.h"
#include "linear_algebra/Vect_impl.h"
#include "util/util.h"
#include "OFELIException.h"

namespace OFELI {

void saveField(Vect<real_t>& v,
               string        output_file,
               int           opt)
{
   const Mesh *mesh = &(v.getMesh());
   if (mesh==nullptr)
      throw OFELIException("In saveField(v,*mesh,output_file,opt): "
                           " Vector does not contain mesh information.\n"
                           "Try the function saveField(Vect<real_t>,Mesh,string,int)");
   saveField(v,*mesh,output_file,opt);
}


#ifdef USE_PETSC
void saveField(PETScVect<real_t>& v,
               string             output_file,
               int                opt)
{
   const Mesh *mesh = &(v.getMesh());
   bool mts = false;
   static int ts = 0;
   size_t i, k, n, m=0;
   bool scalar = true;
   size_t nb_dof = v.getNbDOF();
   size_t dim = mesh->getDim();
   DOFSupport type = v.getDOFType();
   if (nb_dof>1)
      scalar = false;
   ofstream fp(output_file.c_str());
   fp.setf(ios::right|ios::scientific);
   fp << setprecision(16);

   switch (opt) {

      case GNUPLOT:
         if (dim==1) {
            node_loop(mesh)
              fp << The_node.getX() << " " << v(node_label) << endl;
            break;
         }
         element_loop(mesh) {
            m = 0;
            switch(The_element.getShape()) {
               case LINE:          m = 2; break;
               case TRIANGLE:      m = 3; break;
               case QUADRILATERAL: m = 4; break;
            }
            for (n=1; n<=m; n++) {
               Node *nd = The_element(n);
               fp << nd->getX() << " " << nd->getY() << " " << v(nd->n()) << endl;
            }
            Node *nd = The_element(1);
            fp << nd->getX() << " " << nd->getY() << " " << v(nd->n()) << endl << endl;
         }
         break;

      case VTK:
      {
         size_t size = 0;
         element_loop(mesh) {
            switch(The_element.getShape()) {
               case LINE:          m = 2; break;
               case TRIANGLE:      m = 3; break;
               case QUADRILATERAL: m = 4; break;
               case TETRAHEDRON:   m = 4; break;
               case HEXAHEDRON:    m = 8; break;
               case PENTAHEDRON:   m = 6; break;
            }
            size += m+1;
         }
         fp << "# vtk DataFile Version 2.0\nImported from OFELI files\nASCII" << endl;
         fp << "DATASET UNSTRUCTURED_GRID\nPOINTS " << mesh->getNbNodes() << " double" << endl;
         mesh_nodes(*mesh)
            fp << The_node.getX() << "  " << The_node.getY() << " " << The_node.getZ() << endl;
         fp << "\nCELLS " << mesh->getNbElements() << " " << size << endl;
         element_loop(mesh) {
            switch (The_element.getShape()) {
               case LINE:          m = 2; break;
               case TRIANGLE:      m = 3; break;
               case QUADRILATERAL: m = 4; break;
               case TETRAHEDRON:   m = 4; break;
               case HEXAHEDRON:    m = 8; break;
               case PENTAHEDRON:   m = 6; break;
            }
            fp << setw(9) << m;
            for (i=0; i<m; i++)
               fp << The_element(i+1)->n()-1;
            fp << endl;
         }
         fp << "\nCELL_TYPES  " << mesh->getNbElements() << endl;
         k = 0;
         element_loop(mesh) {
            switch(The_element.getShape()) {
               case LINE:          m =  3; break;
               case TRIANGLE:      m =  5; break;
               case QUADRILATERAL: m =  9; break;
               case TETRAHEDRON:   m = 10; break;
               case HEXAHEDRON:    m = 12; break;
               case PENTAHEDRON:   m = 13; break;
            }
            fp << setw(4) << m;
            if (++k%30 == 0)
               fp << endl;
         }
         if (type==NODE_DOF)
            fp << "\nPOINT_DATA  " << mesh->getNbNodes() << endl;
         else if (type==ELEMENT_DOF)
            fp << "\nCELL_DATA  " << mesh->getNbElements() << endl;
         if (scalar)
            fp << "SCALARS  u  double  1\nLOOKUP_TABLE  default" << endl;
         else
            fp << "VECTORS  u  double" << endl;

         if (type==NODE_DOF) {
            node_loop(mesh) {
               fp << v(node_label,1) << " ";
               if (!scalar) {
                  fp << v(node_label,2) << " ";
                  if (nb_dof > 2)
                     fp << v(node_label,3) << " ";
                  else
                     fp << 0. << " ";
               }
               fp << endl;
            }
         }
         else if (type==ELEMENT_DOF) {
            element_loop(mesh) {
               fp << v(element_label,1) << " ";
               if (!scalar) {
                  fp << v(element_label,2) << " ";
                  if (nb_dof > 2)
                     fp << v(element_label,3) << " ";
                  else
                     fp << 0. << " ";
               }
               fp << endl;
            }
         }
      }
      break;

      case GMSH:
      {
         char tt = 'S';
         if (nb_dof == dim)
            tt = 'V';
         fp << "View \"" << v.getName() << "\" {" << endl;
         switch (dim) {

            case 1:
               element_loop(mesh) {
                  fp << "SL(";
                  fp << The_element(1)->getX() <<  ", 0., 0., "
                     << The_element(2)->getX() <<  ", 0., 0. ) {\n";
                  fp << v(The_element(1)->n(),1) << ","
                     << v(The_element(2)->n(),1) << "};\n";
               }
               fp << "};" << endl;
               break;

            case 2:
               element_loop(mesh) {
                  size_t nb_en = The_element.getNbNodes();
                  fp << tt;
                  if (nb_en==3)
                     fp << "T";
                  else if (nb_en==4)
                     fp << "Q";
                  fp << "(";
                  for (k=1; k<nb_en; k++)
                     fp << The_element(k)->getX() << "," << 
                           The_element(k)->getY() << ",0.,";
                  fp << The_element(nb_en)->getX() << ","
                     << The_element(nb_en)->getY() << ",0.) {" << endl;
                  for (k=1; k<=nb_en; k++) {
                     fp << v(The_element(k)->n(),1);
                     if (nb_dof > 1)
                        fp << "," << v(The_element(k)->n(),2) << ",0.0";
                     if (k<nb_en)
                        fp << ",";
                  }
                  fp << "\n};\n";
               }
               fp << "};" << endl;
               break;

            case 3:
               element_loop(mesh) {
                  size_t nb_en = The_element.getNbNodes();
                  fp << tt;
                  if (nb_en==4)
                     fp << "S";
                  else if (nb_en==8)
                     fp << "H";
                  else if (nb_en==6)
                     fp << "I";
                  else if (nb_en==5)
                     fp << "Y";
                  fp << "(";
                  for (k=1; k<=nb_en-1; k++)
                     fp << The_element(k)->getX() << ","
                        << The_element(k)->getY() << ","
                        << The_element(k)->getZ() << ",";
                  fp << The_element(nb_en)->getX() << ","
                     << The_element(nb_en)->getY() << ","
                     << The_element(nb_en)->getZ() << ") {" << endl;
                  for (k=1; k<=nb_en; k++) {
                     fp << v(The_element(k)->n(),1);
                     if (nb_dof > 1)
                        fp << "," << v(The_element(k)->n(),2) << ","
                                  << v(The_element(k)->n(),3) << endl;
                     if (k<nb_en)
                        fp << ",";
                  }
                  fp << "};" << endl;
               }
               fp << "};" << endl;
               break;
         }
      }
      break;

      case TECPLOT:
      {
         string shape;
         switch(mesh->getShape()) {
            case LINE:           shape = "LINESEG";       break;
            case TRIANGLE:       shape = "TRIANGLE";      break;
            case QUADRILATERAL:  shape = "QUADRILATERAL"; break;
            case TETRAHEDRON:    shape = "TETRAHEDRON";   break;
            case HEXAHEDRON:     shape = "BRICK";         break;
         }
         fp << "TITLE = \" \"\n" << endl;
         fp << "VARIABLES = \"X\", \"Y\"";
         if (dim==3)
            fp << ", \"Z\"";
         if (nb_dof == 1)
            fp << ", \"T\"";
         else if (nb_dof == 2)
            fp << ", \"UX\", \"UY\"";
         else if (nb_dof == 3)
            fp << ", \"UX\", \"UY\", \"UZ\"";
         fp << endl;
         fp << "\nZONE T=\"" << "step-" << 1 << "\", N=" << mesh->getNbNodes() << ", E="
            << mesh->getNbElements() << ", F=FEPOINT, ET=" << shape
            << ", SOLUTIONTIME=" << 0.;
         if (mts) {
            fp << ", D=(1,";
            if (dim>1)
               fp << "2,";
            if (dim==3)
               fp << "3,";
            fp << "FECONNECT)";
         }
         fp << endl;
         node_loop(mesh) {
            for (size_t i=1; i<=dim; i++)
               fp << "  " << The_node.getCoord(i);
            for (size_t j=0; j<nb_dof; j++)
               fp << "  " << v(node_label,j);
            fp << endl;
         }
         element_loop(mesh) {
            for (size_t n=1; n<=The_element.getNbNodes(); n++)
               fp << "  " << The_element(n)->n();
            fp << endl;
         }
      }
      break;

      case MATLAB:
         break;
   }
   ts++;
   fp.close();
}


void saveField(PETScVect<real_t>& v,
               const Mesh&        mesh,
               string             output_file,
               int                opt)
{
   bool mts = false;
   static int ts = 0;
   size_t i, k, n, m=0;
   bool scalar = true;
   size_t nb_dof = v.getNbDOF();
   size_t dim = mesh.getDim();
   DOFSupport type = v.getDOFType();
   if (nb_dof>1)
      scalar = false;
   ofstream fp(output_file.c_str());
   fp.setf(ios::right|ios::scientific);
   fp << setprecision(16);

   switch (opt) {

      case GNUPLOT:
         if (dim==1) {
            node_loop(&mesh) {
               fp << setprecision(4) << setw(18) << The_node.getX() << " "
                  << setprecision(4) << setw(18) << v(node_label) << endl;
            }
            break;
         }
         element_loop(&mesh) {
            m = 0;
            switch(The_element.getShape()) {
               case LINE:          m = 2; break;
               case TRIANGLE:      m = 3; break;
               case QUADRILATERAL: m = 4; break;
            }
            for (n=1; n<=m; n++) {
               Node *nd = The_element(n);
               fp << setprecision(4) << setw(18) << nd->getX() << " "
                  << setprecision(4) << setw(18) << nd->getY()
                  << setprecision(4) << setw(18) << v(nd->n()) << endl;
            }
            Node *nd = The_element(1);
            fp << setprecision(4) << setw(18) << nd->getX() << " "
               << setprecision(4) << setw(18) << nd->getY()
               << setprecision(4) << setw(18) << v(nd->n()) << endl << endl;
         }
         break;

      case VTK:
      {
         size_t size = 0;
         element_loop(&mesh) {
            switch(The_element.getShape()) {
               case LINE:          m = 2; break;
               case TRIANGLE:      m = 3; break;
               case QUADRILATERAL: m = 4; break;
               case TETRAHEDRON:   m = 4; break;
               case HEXAHEDRON:    m = 8; break;
               case PENTAHEDRON:   m = 6; break;
            }
            size += m+1;
         }
         fp << "# vtk DataFile Version 2.0\nImported from OFELI files\nASCII" << endl;
         fp << "DATASET UNSTRUCTURED_GRID\nPOINTS " << mesh.getNbNodes() << " double" << endl;
         mesh_nodes(mesh) {
            fp << The_node.getX() << "  " << The_node.getY() << "  "
               << The_node.getZ() << endl;
         }
         fp << "\nCELLS " << mesh.getNbElements() << " " << size << endl;
         mesh_elements(mesh) {
            switch (The_element.getShape()) {
               case LINE:          m = 2; break;
               case TRIANGLE:      m = 3; break;
               case QUADRILATERAL: m = 4; break;
               case TETRAHEDRON:   m = 4; break;
               case HEXAHEDRON:    m = 8; break;
               case PENTAHEDRON:   m = 6; break;
            }
            fp << setw(9) << m;
            for (i=0; i<m; i++)
               fp << The_element(i+1)->n()-1;
            fp << endl;
         }
         fp << "\nCELL_TYPES  " << mesh.getNbElements() << endl;
         k = 0;
         element_loop(&mesh) {
            switch(The_element.getShape()) {
               case LINE:          m =  3; break;
               case TRIANGLE:      m =  5; break;
               case QUADRILATERAL: m =  9; break;
               case TETRAHEDRON:   m = 10; break;
               case HEXAHEDRON:    m = 12; break;
               case PENTAHEDRON:   m = 13; break;
            }
            fp << setw(4) << m;
            if (++k%30 == 0)
               fp << endl;
         }
         if (type==NODE_DOF)
            fp << "\nPOINT_DATA  " << mesh.getNbNodes() << endl;
         else if (type==ELEMENT_DOF)
            fp << "\nCELL_DATA  " << mesh.getNbElements() << endl;
         if (scalar)
            fp << "SCALARS  u  double  1\nLOOKUP_TABLE  default" << endl;
         else
            fp << "VECTORS  u  double" << endl;

         if (type==NODE_DOF) {
            node_loop(&mesh) {
               fp << v(node_label,1) << " ";
               if (!scalar) {
                  fp << v(node_label,2) << " ";
                  if (nb_dof > 2)
                     fp << v(node_label,3) << " ";
                  else
                     fp << 0. << " ";
               }
               fp << endl;
            }
         }
         else if (type==ELEMENT_DOF) {
            element_loop(&mesh) {
               fp << v(element_label,1) << " ";
               if (!scalar) {
                  fp << v(element_label,2) << " ";
                  if (nb_dof > 2)
                     fp << v(element_label,3) << " ";
                  else
                     fp << 0. << " ";
               }
               fp << endl;
            }
         }
      }
      break;

      case GMSH:
      {
         char tt = 'S';
         if (nb_dof == dim)
            tt = 'V';
         fp << "View \"" << v.getName() << "\" {" << endl;
         switch (dim) {

            case 1:
               element_loop(&mesh) {
                  fp << "SL(";
                  fp << The_element(1)->getX() <<  ", 0., 0., "
                     << The_element(2)->getX() <<  ", 0., 0. ) {\n";
                  fp << v(The_element(1)->n(),1) << ","
                     << v(The_element(2)->n(),1) << "};\n";
               }
               fp << "};" << endl;
               break;

            case 2:
               element_loop(&mesh) {
                  size_t nb_en = The_element.getNbNodes();
                  fp << tt;
                  if (nb_en==3)
                     fp << "T";
                  else if (nb_en==4)
                     fp << "Q";
                  fp << "(";
                  for (k=1; k<nb_en; k++)
                     fp << The_element(k)->getX() << "," << 
                           The_element(k)->getY() << ",0.,";
                  fp << The_element(nb_en)->getX() << ","
                     << The_element(nb_en)->getY() << ",0.) {" << endl;
                  for (k=1; k<=nb_en; k++) {
                     fp << v(The_element(k)->n(),1);
                     if (nb_dof > 1)
                        fp << "," << v(The_element(k)->n(),2) << ",0.0";
                     if (k<nb_en)
                        fp << ",";
                  }
                  fp << "\n};\n";
               }
               fp << "};" << endl;
               break;

            case 3:
               element_loop(&mesh) {
                  size_t nb_en = The_element.getNbNodes();
                  fp << tt;
                  if (nb_en==4)
                     fp << "S";
                  else if (nb_en==8)
                     fp << "H";
                  else if (nb_en==6)
                     fp << "I";
                  else if (nb_en==5)
                     fp << "Y";
                  fp << "(";
                  for (k=1; k<=nb_en-1; k++)
                     fp << The_element(k)->getX() << ","
                        << The_element(k)->getY() << ","
                        << The_element(k)->getZ() << ",";
                  fp << The_element(nb_en)->getX() << ","
                     << The_element(nb_en)->getY() << ","
                     << The_element(nb_en)->getZ() << ") {" << endl;
                  for (k=1; k<=nb_en; k++) {
                     fp << v(The_element(k)->n(),1);
                     if (nb_dof > 1)
                        fp << "," << v(The_element(k)->n(),2) << ","
                                  << v(The_element(k)->n(),3) << endl;
                     if (k<nb_en)
                        fp << ",";
                  }
                  fp << "};" << endl;
               }
               fp << "};" << endl;
               break;
         }
      }
      break;

      case TECPLOT:
      {
         string shape;
         switch(mesh.getShape()) {
            case LINE:           shape = "LINESEG";       break;
            case TRIANGLE:       shape = "TRIANGLE";      break;
            case QUADRILATERAL:  shape = "QUADRILATERAL"; break;
            case TETRAHEDRON:    shape = "TETRAHEDRON";   break;
            case HEXAHEDRON:     shape = "BRICK";         break;
         }
         fp << "TITLE = \" \"\n" << endl;
         fp << "VARIABLES = \"X\", \"Y\"";
         if (dim==3)
            fp << ", \"Z\"";
         if (nb_dof == 1)
            fp << ", \"T\"";
         else if (nb_dof == 2)
            fp << ", \"UX\", \"UY\"";
         else if (nb_dof == 3)
            fp << ", \"UX\", \"UY\", \"UZ\"";
         fp << endl;
         fp << "\nZONE T=\"" << "step-" << 1 << "\", N=" << mesh.getNbNodes() << ", E="
            << mesh.getNbElements() << ", F=FEPOINT, ET=" << shape
            << ", SOLUTIONTIME=" << 0.;
         if (mts) {
            fp << ", D=(1,";
            if (dim>1)
               fp << "2,";
            if (dim==3)
               fp << "3,";
            fp << "FECONNECT)";
         }
         fp << endl;
         node_loop(&mesh) {
            for (size_t i=1; i<=dim; i++)
               fp << "  " << The_node.getCoord(i);
            for (size_t j=0; j<nb_dof; j++)
               fp << "  " << v(node_label,j);
            fp << endl;
         }
         element_loop(&mesh) {
            for (size_t n=1; n<=The_element.getNbNodes(); n++)
               fp << "  " << The_element(n)->n();
            fp << endl;
         }
      }
      break;

      case MATLAB:
         break;
   }
   ts++;
   fp.close();
}
#endif


void saveField(const Vect<real_t>& v,
               const Mesh&         mesh,
               string              output_file,
               int                 opt)
{
   bool mts = false;
   static int ts = 0;
   size_t i, k, n, m=0;
   bool scalar = true;
   size_t nb_dof = v.getNbDOF();
   DOFSupport type = v.getDOFType();
   if (nb_dof>1)
      scalar = false;
   ofstream fp(output_file.c_str());
   map<int,int> sh = {{LINE,2},{TRIANGLE,3},{QUADRILATERAL,4},{TETRAHEDRON,4},
                      {HEXAHEDRON,8},{PENTAHEDRON,6}};
   map<int,int> sh1 = {{LINE,2},{TRIANGLE,5},{QUADRILATERAL,9},{TETRAHEDRON,10},
                      {HEXAHEDRON,12},{PENTAHEDRON,13}};
   fp.setf(ios::right|ios::scientific);

   switch (opt) {

      case GNUPLOT:
         if (mesh.getDim()==1) {
            node_loop(&mesh) {
               fp << setprecision(4) << setw(18) << The_node.getX() << " "
                  << setprecision(4) << setw(18) << v(node_label) << endl;
            }
            break;
         }
         element_loop(&mesh) {
            m = sh[The_element.getShape()];
            for (n=1; n<=m; n++) {
               Node *nd = The_element(n);
               fp << setprecision(4) << setw(18) << nd->getX() << " "
                  << setprecision(4) << setw(18) << nd->getY()
                  << setprecision(4) << setw(18) << v(nd->n()) << endl;
            }
            Node *nd = The_element(1);
            fp << setprecision(4) << setw(18) << nd->getX() << " "
               << setprecision(4) << setw(18) << nd->getY()
               << setprecision(4) << setw(18) << v(nd->n()) << endl << endl;
         }
         break;

      case VTK:
      {
         size_t size = 0;
         element_loop(&mesh)
            size += sh[The_element.getShape()]+1;
         fp << "# vtk DataFile Version 2.0\nImported from OFELI files\nASCII" << endl;
         fp << "DATASET UNSTRUCTURED_GRID\nPOINTS " << mesh.getNbNodes() << " double" << endl;
         node_loop(&mesh) {
            fp << The_node.getX() << "  " << The_node.getY() << "  "
               << The_node.getZ() << endl;
         }
         fp << "\nCELLS " << mesh.getNbElements() << " " << size << endl;
         element_loop(&mesh) {
            m = sh[The_element.getShape()];
            fp << setw(9) << m;
            for (i=0; i<m; i++)
               fp << setw(9) << The_element(i+1)->n()-1;
            fp << endl;
         }
         fp << "\nCELL_TYPES  " << mesh.getNbElements() << endl;
         k = 0;
         element_loop(&mesh) {
            fp << setw(4) << sh1[The_element.getShape()];
            if (++k%30 == 0)
               fp << endl;
         }
         if (type==NODE_DOF)
            fp << "\nPOINT_DATA  " << mesh.getNbNodes() << endl;
         else if (type==ELEMENT_DOF)
            fp << "\nCELL_DATA  " << mesh.getNbElements() << endl;
         if (scalar)
            fp << "SCALARS  u  double  1\nLOOKUP_TABLE  default" << endl;
         else
            fp << "VECTORS  u  double" << endl;

         if (type==NODE_DOF) {
            node_loop(&mesh) {
               fp << v(node_label,1) << " ";
               if (!scalar) {
                  fp << v(node_label,2) << " ";
                  if (nb_dof > 2)
                     fp << v(node_label,3) << " ";
                  else
                     fp << 0. << " ";
               }
               fp << endl;
            }
         }
         else if (type==ELEMENT_DOF) {
            element_loop(&mesh) {
               fp << v(element_label,1) << " ";
               if (!scalar) {
                  fp << v(element_label,2) << " ";
                  if (nb_dof > 2)
                     fp << v(element_label,3) << " ";
                  else
                     fp << 0. << " ";
               }
               fp << endl;
            }
         }
      }
      break;

      case GMSH:
      {
         char tt = 'S';
         if (nb_dof == mesh.getDim())
            tt = 'V';
         fp << "View \"" << v.getName() << "\" {" << endl;
         switch (mesh.getDim()) {

            case 1:
               element_loop(&mesh) {
                  fp << "SL(";
                  fp << The_element(1)->getX() <<  ", 0., 0., "
                     << The_element(2)->getX() <<  ", 0., 0. ) {\n";
                  fp << v(The_element(1)->n(),1) << ","
                     << v(The_element(2)->n(),1) << "};\n";
               }
               fp << "};" << endl;
               break;

            case 2:
               element_loop(&mesh) {
                  size_t nb_en = The_element.getNbNodes();
                  fp << tt;
                  if (nb_en==3)
                     fp << "T";
                  else if (nb_en==4)
                     fp << "Q";
                  fp << "(";
                  for (k=1; k<nb_en; k++)
                     fp << The_element(k)->getX() << "," << 
                           The_element(k)->getY() << ",0.,";
                  fp << The_element(nb_en)->getX() << ","
                     << The_element(nb_en)->getY() << ",0.) {" << endl;
                  for (k=1; k<=nb_en; k++) {
                     fp << v(The_element(k)->n(),1);
                     if (nb_dof > 1)
                        fp << "," << v(The_element(k)->n(),2) << ",0.0";
                     if (k<nb_en)
                        fp << ",";
                  }
                  fp << "\n};\n";
               }
               fp << "};" << endl;
               break;

            case 3:
               element_loop(&mesh) {
                  size_t nb_en = The_element.getNbNodes();
                  fp << tt;
                  if (nb_en==4)
                     fp << "S";
                  else if (nb_en==8)
                     fp << "H";
                  else if (nb_en==6)
                     fp << "I";
                  else if (nb_en==5)
                     fp << "Y";
                  fp << "(";
                  for (k=1; k<=nb_en-1; k++)
                     fp << The_element(k)->getX() << ","
                        << The_element(k)->getY() << ","
                        << The_element(k)->getZ() << ",";
                  fp << The_element(nb_en)->getX() << ","
                     << The_element(nb_en)->getY() << ","
                     << The_element(nb_en)->getZ() << ") {" << endl;
                  for (k=1; k<=nb_en; k++) {
                     fp << v(The_element(k)->n(),1);
                     if (nb_dof > 1)
                        fp << "," << v(The_element(k)->n(),2) << ","
                                  << v(The_element(k)->n(),3) << endl;
                     if (k<nb_en)
                        fp << ",";
                  }
                  fp << "};" << endl;
               }
               fp << "};" << endl;
               break;
         }
      }
      break;

      case TECPLOT:
      {
         map<int,string> sh = {{LINE,"LINESEG"},{TRIANGLE,"TRIANGLE"},{QUADRILATERAL,"QUADRILATERAL"},
                               {TETRAHEDRON,"TETRAHEDRON"},{HEXAHEDRON,"BRICK"}};
         string shape = sh[mesh.getShape()];
         fp << "TITLE = \" \"\n" << endl;
         fp << "VARIABLES = \"X\", \"Y\"";
         if (mesh.getDim()==3)
            fp << ", \"Z\"";
         if (nb_dof == 1)
            fp << ", \"T\"";
         else if (nb_dof == 2)
            fp << ", \"UX\", \"UY\"";
         else if (nb_dof == 3)
            fp << ", \"UX\", \"UY\", \"UZ\"";
         fp << endl;
         fp << "\nZONE T=\"" << "step-" << 1 << "\", N=" << mesh.getNbNodes() << ", E="
            << mesh.getNbElements() << ", F=FEPOINT, ET=" << shape
            << ", SOLUTIONTIME=" << 0.;
         if (mts) {
            fp << ", D=(1,";
            if (mesh.getDim()>1)
               fp << "2,";
            if (mesh.getDim()==3)
               fp << "3,";
            fp << "FECONNECT)";
         }
         fp << endl;
         node_loop(&mesh) {
            for (size_t i=1; i<=mesh.getDim(); i++)
               fp << "  " << The_node.getCoord(i);
            for (size_t j=0; j<nb_dof; j++)
               fp << "  " << v(node_label,j);
            fp << endl;
         }
         element_loop(&mesh) {
            for (size_t n=1; n<=The_element.getNbNodes(); n++)
               fp << "  " << The_element(n)->n();
            fp << endl;
         }
      }
      break;

      case MATLAB:
         break;
   }
   ts++;
   fp.close();
}


void saveField(Vect<real_t>& v,
               const Grid&   g,
               string        output_file,
               int           opt)
{
   bool scalar = true;
   size_t nb_dof = v.getNbDOF();
   if (nb_dof>1)
      scalar = false;
   if (scalar==false)
      throw OFELIException("In saveField(Vect<real_t>,Grid,string,int): "
                           " This function is not implemented for vector data");
   ofstream fp(output_file.c_str());
   if (fp.fail())
      throw OFELIException("In saveField(Vect<real_t>,Grid,string,int): "
                           " Failed to open file "+output_file);
   fp.setf(ios::right|ios::scientific);
   size_t nx=g.getNx(), ny=g.getNy(), nz=g.getNz();

   switch (opt) {

      case GMSH:
      {
         fp << "View \"" << "#" << "\" {" << endl;
         switch (g.getDim()) {

            case 1:
               for (size_t i=1; i<=nx; i++) {
                  fp << "SL(";
                  fp << g.getX(i) <<  ", 0., 0., " << g.getX(i+1) <<  ", 0., 0. ) {" << endl;
                  fp << v(i) << "," << v(i+1) << "};" << endl;
               }
               fp << "};" << endl;
               break;

            case 2:
               for (size_t i=1; i<=nx; ++i) {
                  for (size_t j=1; j<=ny; ++j) {
                     fp << "ST(";
                     fp << g.getXY(i  ,j  ).x << "," << g.getXY(i  ,j  ).y << ",0.,";
                     fp << g.getXY(i+1,j  ).x << "," << g.getXY(i+1,j  ).y << ",0.,";
                     fp << g.getXY(i+1,j+1).x << "," << g.getXY(i+1,j+1).y << ",0.) {" << endl;
                     fp << v(i,j) << "," << v(i+1,j) << "," << v(i+1,j+1) << "};" << endl;
                     fp << "ST(";
                     fp << g.getXY(i+1,j+1).x << "," << g.getXY(i+1,j+1).y << ",0.,";
                     fp << g.getXY(i  ,j+1).x << "," << g.getXY(i  ,j+1).y << ",0.,";
                     fp << g.getXY(i  ,j  ).x << "," << g.getXY(i  ,j  ).y << ",0.) {" << endl;
                     fp << v(i+1,j+1) << "," << v(i,j+1) << "," << v(i,j) << "};" << endl;
                  }
               }
               fp << "};" << endl;
               break;

            case 3:
               for (size_t i=1; i<=nx; ++i) {
                  for (size_t j=1; j<=ny; ++j) {
                     for (size_t k=1; k<=nz; ++k) {
                        fp << "SH(";
                        fp << g.getXYZ(i  ,j  ,k  ).x << "," << g.getXYZ(i  ,j  ,k  ).y << "," << g.getXYZ(i  ,j  ,k  ).z << ",";
                        fp << g.getXYZ(i+1,j  ,k  ).x << "," << g.getXYZ(i+1,j  ,k  ).y << "," << g.getXYZ(i+1,j  ,k  ).z << ",\n";
                        fp << g.getXYZ(i+1,j+1,k  ).x << "," << g.getXYZ(i+1,j+1,k  ).y << "," << g.getXYZ(i+1,j+1,k  ).z << ",";
                        fp << g.getXYZ(i  ,j  ,k  ).x << "," << g.getXYZ(i  ,j  ,k  ).y << "," << g.getXYZ(i  ,j  ,k  ).z << ",\n";
                        fp << g.getXYZ(i  ,j  ,k+1).x << "," << g.getXYZ(i  ,j  ,k+1).y << "," << g.getXYZ(i  ,j  ,k+1).z << ",";
                        fp << g.getXYZ(i+1,j  ,k+1).x << "," << g.getXYZ(i+1,j  ,k+1).y << "," << g.getXYZ(i+1,j  ,k+1).z << ",\n";
                        fp << g.getXYZ(i+1,j+1,k+1).x << "," << g.getXYZ(i+1,j+1,k+1).y << "," << g.getXYZ(i+1,j+1,k+1).z << ",";
                        fp << g.getXYZ(i  ,j  ,k+1).x << "," << g.getXYZ(i  ,j  ,k+1).y << "," << g.getXYZ(i  ,j  ,k+1).z << ") {";
                        fp << v(i,j,k  ) << "," << v(i+1,j,k  ) << "," << v(i+1,j+1,k  ) << "," << v(i,j,k  ) << ",";
                        fp << v(i,j,k+1) << "," << v(i+1,j,k+1) << "," << v(i+1,j+1,k+1) << "," << v(i,j,k+1) << "};" << endl;
                     }
                  }
               }
               fp << "};" << endl;
               break;
         }
      }
      break;
      
      case VTK:
      {
         if (nz > 1)
            throw OFELIException("In saveField(Vect<real_t>,Grid,string,int): "
                                 "This function is not implemented for 3-D");
         fp << "# vtk DataFile Version 2.0\nImported from OFELI files\nASCII\nDATASET POLYDATA\n";
         fp << "POINTS " << (nx+1)*(ny+1)*(nz+1) << " double\n";
         for (size_t i=1; i<=nx+1; i++) {
            for (size_t j=1; j<=ny+1; j++)
               fp << g.getX(i) << "  " << g.getY(j) << "  0.000" << endl;
         }
         fp << "POLYGONS  " << nx*ny << "  " << 5*nx*ny << endl;
         size_t nn = 1;
         for (size_t i=1; i<=nx; i++) {
            for (size_t j=1; j<=ny; j++)
               fp << "4  " << nn+j-2 << "  " << nn+j+ny-1 << "  " << nn+j+ny << "  " << nn+j-1 << endl;
            nn += ny+1;
         }
         if (v.getDOFType()==NODE_DOF) {
            nx++, ny++;
            fp << "\nPOINT_DATA  " << nx*ny << endl;
         }
         else
            fp << "\nCELL_DATA  " << nx*ny << endl;
         fp << "SCALARS  u  double  1\nLOOKUP_TABLE  default" << endl;
         for (size_t i=1; i<=nx+1; i++) {
            for (size_t j=1; j<=ny+1; j++)
               fp << v(i,j) << "  ";
            fp << endl;
         }
      }
      break;
   }
   fp.close();
}


void saveFormat(Mesh&  mesh,
                string input_file,
                string output_file,
                int    format,
                int    f)
{
   if (f<0)
      return;

   switch (format) {

      case GNUPLOT:
         saveGnuplot(mesh,input_file,output_file,f);
         break;

      case VTK:
         saveVTK(mesh,input_file,output_file,f);
         break;

      case GMSH:
         saveGmsh(mesh,input_file,output_file,f);
         break;

      case TECPLOT:
         saveTecplot(mesh,input_file,output_file,f);
         break;
   }
}


void saveGnuplot(string input_file,
                 string output_file,
                 string mesh_file,
                 int    f)
{
   Mesh mesh(mesh_file);
   saveGnuplot(mesh,input_file,output_file,f);
}


void saveGnuplot(Mesh&  mesh,
                 string input_file,
                 string output_file,
                 int    f)
{
   size_t nb_dof=0;
   cout << "Converting file: " << input_file << " to Gnuplot format." << endl;
   vector<vector<real_t> > field;
   vector<real_t> tm;
   string name;
   getfields(input_file,mesh,nb_dof,tm,field,name);

   Vect<real_t> v(mesh);
   size_t dim = mesh.getDim();
   ofstream fp(output_file.c_str());
   fp.setf(ios::right|ios::scientific);
   element_loop(&mesh) {
      size_t m = 0;
      switch (dim) {

         case 1: if (The_element.getShape()==LINE)
                    m = 2;
                 if (!m)
                    throw OFELIException("In saveGnuplot(string,string,string,int): "
                                          " Unknown element shape " + to_string(element_label));
                 break;

         case 2: if (The_element.getShape()==QUADRILATERAL)
                    m = 4;
                 if (The_element.getShape()==TRIANGLE)
                    m = 3;
                 if (!m)
                    throw OFELIException("In saveGnuplot(string,string,string,int): "
                                         " Unknown element shape " + to_string(element_label));
      }
      for (size_t n=1; n<=m; n++) {
         the_node = The_element(n);
         fp << setprecision(4) << setw(18) << The_node.getX() << " "
            << setprecision(4) << setw(18) << The_node.getY();
      }
      the_node = The_element(1);
      fp << setprecision(4) << setw(18) << The_node.getX() << " "
         << setprecision(4) << setw(18) << The_node.getY() << endl << endl;
   }
   fp.close();
}


void saveTecplot(string input_file,
                 string output_file,
                 string mesh_file,
                 int    f)
{
   Mesh mesh(mesh_file);
   saveTecplot(mesh,input_file,output_file,f);
}


void saveTecplot(Mesh&  mesh,
                 string input_file,
                 string output_file,
                 int    f)
{
   size_t nb_dof=0;
   cout << "Converting file: " << input_file << " to Tecplot format." << endl;
   vector<vector<real_t> > v;
   vector<real_t> tm;
   string name;
   getfields(input_file,mesh,nb_dof,tm,v,name);
   map<int,string> sh = {{LINE,"LINESEG"},{TRIANGLE,"TRIANGLE"},{QUADRILATERAL,"QUADRILATERAL"},
                         {TETRAHEDRON,"TETRAHEDRON"},{HEXAHEDRON,"BRICK"}};
   string shape = sh[mesh.getShape()];
   ofstream fp(output_file.c_str());
   fp.setf(ios::right|ios::scientific);
   size_t count = 0, nb_time = tm.size();

   for (size_t n=0; n<nb_time; n++) {
      if (count==0) {
         fp << "TITLE = \" \"\n" << endl;
         fp << "VARIABLES = \"X\", \"Y\"";
         if (mesh.getDim()==3)
            fp << ", \"Z\"";
         if (nb_dof == 1)
            fp << ", \"T\"";
         else if (nb_dof == 2)
            fp << ", \"UX\", \"UY\"";
         else if (nb_dof == 3)
            fp << ", \"UX\", \"UY\", \"UZ\"";
         fp << endl;
      }

      if ((n+1)%f==0) {
         fp << "\nZONE T=\"" << "step-" << n << "\", N=" << mesh.getNbNodes() << ", E="
            << mesh.getNbElements() << ", F=FEPOINT, ET=" << shape
            << ", SOLUTIONTIME=" << tm[n];
         if (count) {
            fp << ", D=(1,";
            if (mesh.getDim()>1)
               fp << "2,";
            if (mesh.getDim()==3)
               fp << "3,";
            fp << "FECONNECT)";
         }
         fp << endl;

         node_loop(&mesh) {
            if (count==0)
               for (size_t i=1; i<=mesh.getDim(); i++)
                  fp << "  " << The_node.getCoord(i);
            for (size_t j=0; j<nb_dof; j++)
               fp << "  " << v[n][nb_dof*(node_label-1)+j];
            fp << endl;
         }
         if (count==0) {
            element_loop(&mesh) {
               for (size_t i=1; i<=The_element.getNbNodes(); i++)
                  fp << setw(10) << The_element(i)->n();
               fp << endl;
            }
         }
         count++;
      }
   }
   fp.close();
}


void saveVTK(string input_file,
             string output_file,
             string mesh_file,
             int    f)
{
   Mesh mesh(mesh_file);
   saveVTK(mesh,input_file,output_file,f);
}


void saveVTK(Mesh&  mesh,
             string input_file,
             string output_file,
             int    f)
{
   if (f<=0)
      return;

   ofstream *pf=nullptr;
   size_t k=0, nb_dof;
   string of;

   cout << "Converting file: " << input_file << " to VTK format." << endl;
   string proj=output_file.substr(0,output_file.rfind("."));
   vector<vector<real_t> > field;
   vector<real_t> tm;
   string name;
   getfields(input_file,mesh,nb_dof,tm,field,name);

   size_t nb_st=0;
   bool scalar = true;
   if (nb_dof>1)
      scalar = false;
   size_t size=0, nb_time=tm.size();
   map<int,int> sh = {{LINE,2},{TRIANGLE,3},{QUADRILATERAL,4},{TETRAHEDRON,4},
                      {HEXAHEDRON,8},{PENTAHEDRON,6}};
   element_loop(&mesh)
      size += sh[The_element.getShape()] + 1;
   for (size_t n=0; n<nb_time; n+=f) {
      nb_st++;
      string tt = "Time=" + to_string(tm[n]);
      if (nb_time==1)
         of = proj + ".vtk";
      else
         of = proj + "-" + zeros(n) + ".vtk";
      cout << "   Storing time step " << n+1 << " in file " << of << endl;
      pf = new ofstream(of.c_str());
      *pf << setprecision(16) << std::scientific;
      *pf << "# vtk DataFile Version 2.0\nImported from OFELI files\nASCII" << endl;
      *pf << "DATASET UNSTRUCTURED_GRID\nPOINTS " << mesh.getNbNodes() << " double" << endl;
      node_loop(&mesh)
         *pf << The_node.getX() << "  " << The_node.getY() << "  " << The_node.getZ() << endl;
      *pf << "\nCELLS " << mesh.getNbElements() << setw(10) << size << endl;
      element_loop(&mesh) {
         *pf << setw(8) << sh[The_element.getShape()];
         for (int i=0; i<sh[The_element.getShape()]; i++)
           *pf << setw(10) << The_element(i+1)->n()-1;
         *pf << endl;
      }
      *pf << "\nCELL_TYPES  " << mesh.getNbElements() << endl;
      k = 0;
      sh = {{LINE,3},{TRIANGLE,5},{QUADRILATERAL,9},{TETRAHEDRON,10},
            {HEXAHEDRON,12},{PENTAHEDRON,13}};
      element_loop(&mesh) {
         *pf << setw(4) << sh[The_element.getShape()];
         if (++k%30 == 0)
            *pf << endl;
      }
      *pf << "\nPOINT_DATA  " << mesh.getNbNodes() << endl;

      if (scalar)
         *pf << "SCALARS  u  double  1\nLOOKUP_TABLE  default" << endl;
      else
         *pf << "VECTORS  u  double" << endl;

      node_loop(&mesh) {
         *pf << field[n][nb_dof*(node_label-1)] << " ";
         if (!scalar) {
            *pf << field[n][nb_dof*(node_label-1)+1] << " ";
            if (nb_dof > 2)
               *pf << field[n][nb_dof*(node_label-1)+2] << " ";
            else
               *pf << 0. << " ";
         }
         *pf << endl;
      }
   }
   pf->close();
   delete pf;
   cout << "Number of stored time steps: " << nb_st << endl;
}


void saveGmsh(string input_file,
              string output_file,
              string mesh_file,
              int    f)
{
   Mesh mesh(mesh_file);
   saveGmsh(mesh,input_file,output_file,f);
}


void saveGmsh(Mesh&  mesh,
              string input_file,
              string output_file,
              int    f)
{
   size_t nb_dof=0, nb_en=0;
   ofstream pf(output_file.c_str());
   pf << setprecision(16);

   cout << "Converting file: " << input_file << " to Gmsh format." << endl;
   vector<vector<real_t> > v;
   vector<real_t> tm;
   string name;
   getfields(input_file,mesh,nb_dof,tm,v,name);
   char tt = 'S';
   if (nb_dof == mesh.getDim())
      tt = 'V';

   pf << "View \"" << name << "\" {" << endl;
   size_t nb_time=tm.size();

   switch (mesh.getDim()) {

      case 1:
         element_loop(&mesh) {
            pf << "SL(";
            pf << The_element(1)->getX() <<  ", 0., 0., "
               << The_element(2)->getX() <<  ", 0., 0. ) {" << endl;
            for (size_t n=0; n<nb_time; n+=f) {
               pf << v[n][The_element(1)->n()-1] << "," << v[n][The_element(2)->n()-1];
               if (n<nb_time-1)
                  pf << ",";
               pf << endl;
            }
            pf << "};" << endl;
         }
         pf << "};" << endl;
         break;

      case 2:
         element_loop(&mesh) {
            pf << tt;
            if ((nb_en=The_element.getNbNodes())==3)
               pf << "T(";
            else if (nb_en==4)
               pf << "Q(";
            for (size_t k=1; k<nb_en; ++k)
               pf << The_element(k)->getX() << "," << The_element(k)->getY() << ",0.,";
            pf << The_element(nb_en)->getX() << "," << The_element(nb_en)->getY() << ",0.) {" << endl;
            for (size_t n=0; n<nb_time; n+=f) {
               for (size_t k=1; k<nb_en; ++k) {
                  pf << v[n][nb_dof*(The_element(k)->n()-1)];
                  if (nb_dof > 1)
                     pf << "," << v[n][nb_dof*(The_element(k)->n()-1)+1] << ",0.0";
                  pf << ",";
               }
               pf << v[n][nb_dof*(The_element(nb_en)->n()-1)];
               if (nb_dof > 1)
                  pf << "," << v[n][nb_dof*(The_element(nb_en)->n()-1)+1] << ",0.0";
               if (n<nb_time-1 && n+f<nb_time)
                  pf << ",";
               pf << endl;
            }
            pf << "};" << endl;
         }
         pf << "};" << endl;
         break;

      case 3:
         element_loop(&mesh) {
            pf << tt;
            if ((nb_en=The_element.getNbNodes())==4)
               pf << "S";
            else if (nb_en==8)
               pf << "H";
            else if (nb_en==6)
               pf << "I";
            else if (nb_en==5)
               pf << "Y";
            pf << "(";
            for (size_t k=1; k<=nb_en-1; ++k)
               pf << The_element(k)->getX() << ","
                  << The_element(k)->getY() << ","
                  << The_element(k)->getZ() << ",";
            pf << The_element(nb_en)->getX() << ","
               << The_element(nb_en)->getY() << ","
               << The_element(nb_en)->getZ() << ") {" << endl;
            for (size_t n=0; n<nb_time; n+=f) {
               for (size_t k=1; k<nb_en; ++k) {
                  pf << v[n][nb_dof*(The_element(k)->n()-1)];
                  if (nb_dof > 1)
                     pf << "," << v[n][nb_dof*(The_element(k)->n()-1)+1] << ","
                        << v[n][nb_dof*(The_element(k)->n()-1)+2] << endl;
                  pf << ",";
               }
               pf << v[n][nb_dof*(The_element(nb_en)->n()-1)];
               if (nb_dof > 1)
                  pf << "," << v[n][nb_dof*(The_element(nb_en)->n()-1)+1] << ","
                     << v[n][nb_dof*(The_element(nb_en)->n()-1)+2] << endl;
               if (n<nb_time-1 && n+f<nb_time)
                  pf << ",";
            }
            pf << "};" << endl;
         }
         pf << "};" << endl;
         break;
   }
}


void getfields(string                   file,
               Mesh&                    ms,
               size_t&                  nb_dof,
               vector<real_t>&          t,
               vector<vector<real_t> >& u,
               string&                  name
              )
{
   {
      IOField temp(file,ms,IOField::IN);
      temp.scan(t,NODE_DOF);
      nb_dof = temp.getNbDOF();
   }
   size_t n = t.size();
   if (n)
      cout << "Number of found time steps: " << n << ", Running from "
           << "time: " << t[0] << " to time: " << t[n-1] << "." << endl;
   else {
      cout << "No time step was found." << endl;
      return;
   }

   IOField io(file,ms,IOField::IN);
   u.resize(n);
   io.get(ms,u,name);
}


void saveMatrix(const Matrix<real_t>* A,
                string                file)
{
   ofstream of(file.c_str(),ios::out);
   of.setf(ios::right|ios::scientific);
   of << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" ?>" << endl;
   of << "<OFELI_File>" << endl;
   of << "<info>\n   <title></title>" << endl;
   of << "   <date>" << __DATE__ << "</date>" << endl;
   of << "   <author></author>\n</info>" << endl;
   size_t nr = A->getNbRows(), nc=A->getNbColumns();
   of << "<Matrix nb_rows=\"" << nr << "\" nb_cols=\"" << nc << "\">" << endl;
   for (size_t i=1; i<=nr; ++i) {
      for (size_t j=1; j<=nc; ++j)
         of << "  " << (*A)(i,j);
      of << endl;
   }
   of << "</Matrix>\n</OFELI_File>" << endl;
}

} /* namespace OFELI */
