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
using std::setw;
using std::ios;
using std::setprecision;

#include "OFELI_Config.h"
#include "mesh/MeshUtil.h"
#include "io/saveField.h"
#include "io/IOField.h"

namespace OFELI {

void saveField(Vect<real_t>& v,
               string        output_file,
               int           opt)
{
   const Mesh *mesh = &(v.getMesh());
   if (mesh==NULL)
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
   int type = v.getDOFType();
   if (nb_dof>1)
      scalar = false;
   ofstream fp(output_file.c_str());
   fp.setf(ios::right|ios::scientific);
   fp << setprecision(16);

   switch (opt) {

      case GNUPLOT:
         if (dim==1) {
            mesh_nodes(*mesh)
              fp << The_node.getX() << " " << v(node_label) << endl;
            break;
         }
         mesh_elements(*mesh) {
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
         mesh_elements(*mesh) {
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
         mesh_elements(*mesh) {
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
         mesh_elements(*mesh) {
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
         if (type==NODE_FIELD)
            fp << "\nPOINT_DATA  " << mesh->getNbNodes() << endl;
         else if (type==ELEMENT_FIELD)
            fp << "\nCELL_DATA  " << mesh->getNbElements() << endl;
         if (scalar)
            fp << "SCALARS  u  double  1\nLOOKUP_TABLE  default" << endl;
         else
            fp << "VECTORS  u  double" << endl;

         if (type==NODE_FIELD) {
            mesh_nodes(*mesh) {
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
         else if (type==ELEMENT_FIELD) {
            mesh_elements(*mesh) {
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
               mesh_elements(*mesh) {
                  fp << "SL(";
                  fp << The_element(1)->getX() <<  ", 0., 0., "
                     << The_element(2)->getX() <<  ", 0., 0. ) {\n";
                  fp << v(The_element(1)->n(),1) << ","
                     << v(The_element(2)->n(),1) << "};\n";
               }
               fp << "};" << endl;
               break;

            case 2:
               mesh_elements(*mesh) {
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
               mesh_elements(*mesh) {
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
         mesh_nodes(*mesh) {
            for (size_t i=1; i<=dim; i++)
               fp << "  " << The_node.getCoord(i);
            for (size_t j=0; j<nb_dof; j++)
               fp << "  " << v(node_label,j);
            fp << endl;
         }
         mesh_elements(*mesh) {
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
   int type = v.getDOFType();
   if (nb_dof>1)
      scalar = false;
   ofstream fp(output_file.c_str());
   fp.setf(ios::right|ios::scientific);
   fp << setprecision(16);

   switch (opt) {

      case GNUPLOT:
         if (dim==1) {
            mesh_nodes(mesh) {
               fp << setprecision(4) << setw(18) << The_node.getX() << " "
                  << setprecision(4) << setw(18) << v(node_label) << endl;
            }
            break;
         }
         mesh_elements(mesh) {
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
         mesh_elements(mesh) {
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
         mesh_elements(mesh) {
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
         if (type==NODE_FIELD)
            fp << "\nPOINT_DATA  " << mesh.getNbNodes() << endl;
         else if (type==ELEMENT_FIELD)
            fp << "\nCELL_DATA  " << mesh.getNbElements() << endl;
         if (scalar)
            fp << "SCALARS  u  double  1\nLOOKUP_TABLE  default" << endl;
         else
            fp << "VECTORS  u  double" << endl;

         if (type==NODE_FIELD) {
            mesh_nodes(mesh) {
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
         else if (type==ELEMENT_FIELD) {
            mesh_elements(mesh) {
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
               mesh_elements(mesh) {
                  fp << "SL(";
                  fp << The_element(1)->getX() <<  ", 0., 0., "
                     << The_element(2)->getX() <<  ", 0., 0. ) {\n";
                  fp << v(The_element(1)->n(),1) << ","
                     << v(The_element(2)->n(),1) << "};\n";
               }
               fp << "};" << endl;
               break;

            case 2:
               mesh_elements(mesh) {
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
               mesh_elements(mesh) {
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
         mesh_nodes(mesh) {
            for (size_t i=1; i<=dim; i++)
               fp << "  " << The_node.getCoord(i);
            for (size_t j=0; j<nb_dof; j++)
               fp << "  " << v(node_label,j);
            fp << endl;
         }
         mesh_elements(mesh) {
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
   int type = v.getDOFType();
   if (nb_dof>1)
      scalar = false;
   ofstream fp(output_file.c_str());
   fp.setf(ios::right|ios::scientific);

   switch (opt) {

      case GNUPLOT:
         if (mesh.getDim()==1) {
            mesh_nodes(mesh) {
               fp << setprecision(4) << setw(18) << The_node.getX() << " "
                  << setprecision(4) << setw(18) << v(node_label) << endl;
            }
            break;
         }
         mesh_elements(mesh) {
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
         mesh_elements(mesh) {
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
               fp << setw(9) << The_element(i+1)->n()-1;
            fp << endl;
         }
         fp << "\nCELL_TYPES  " << mesh.getNbElements() << endl;
         k = 0;
         mesh_elements(mesh) {
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
         if (type==NODE_FIELD)
            fp << "\nPOINT_DATA  " << mesh.getNbNodes() << endl;
         else if (type==ELEMENT_FIELD)
            fp << "\nCELL_DATA  " << mesh.getNbElements() << endl;
         if (scalar)
            fp << "SCALARS  u  double  1\nLOOKUP_TABLE  default" << endl;
         else
            fp << "VECTORS  u  double" << endl;

         if (type==NODE_FIELD) {
            mesh_nodes(mesh) {
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
         else if (type==ELEMENT_FIELD) {
            mesh_elements(mesh) {
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
               mesh_elements(mesh) {
                  fp << "SL(";
                  fp << The_element(1)->getX() <<  ", 0., 0., "
                     << The_element(2)->getX() <<  ", 0., 0. ) {\n";
                  fp << v(The_element(1)->n(),1) << ","
                     << v(The_element(2)->n(),1) << "};\n";
               }
               fp << "};" << endl;
               break;

            case 2:
               mesh_elements(mesh) {
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
               mesh_elements(mesh) {
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
         mesh_nodes(mesh) {
            for (size_t i=1; i<=mesh.getDim(); i++)
               fp << "  " << The_node.getCoord(i);
            for (size_t j=0; j<nb_dof; j++)
               fp << "  " << v(node_label,j);
            fp << endl;
         }
         mesh_elements(mesh) {
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
   static int ts = 0;
   size_t i, j;
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

      case VTK:
      {
         if (nz > 1)
            throw OFELIException("In saveField(Vect<real_t>,Grid,string,int): "
                                 "This function is not implemented for 3-D");
         fp << "# vtk DataFile Version 2.0\nImported from OFELI files\nASCII\nDATASET POLYDATA\n";
         fp << "POINTS " << (nx+1)*(ny+1)*(nz+1) << " double\n";
         for (i=1; i<=nx+1; i++) {
            for (j=1; j<=ny+1; j++) {
               fp << g.getX(i) << "  " << g.getY(j) << "  0.000"<< endl;
            }
         }
         fp << "POLYGONS  " << nx*ny << "  " << 5*nx*ny << endl;
         size_t nn = 1;
         for (i=1; i<=nx; i++) {
            for (j=1; j<=ny; j++) {
               fp << "4  " << nn+j-2 << "  " << nn+j+ny-1 << "  " << nn+j+ny
                  << "  " << nn+j-1 << endl;
            }
            nn += ny+1;
         }
         if (v.getDOFType()==NODE_FIELD) {
            nx++, ny++;
            fp << "\nPOINT_DATA  " << nx*ny << endl;
         }
         else
            fp << "\nCELL_DATA  " << nx*ny << endl;
         fp << "SCALARS  u  double  1\nLOOKUP_TABLE  default" << endl;
         for (size_t i=1; i<=nx; i++) {
            for (size_t j=1; j<=ny; j++) {
               fp << v(i,j) << "  ";
            }
            fp << endl;
         }
      }
      break;

      case TECPLOT:
      {
      }
      break;
   }
   ts++;
   fp.close();
}


void saveGnuplot(string input_file,
                 string output_file,
                 string mesh_file)
{
   Mesh mesh(mesh_file);
   vector<real_t> tm;
   cout << "Converting file: " << input_file << " to Gnuplot format." << endl;
   {
      IOField temp(mesh_file,input_file,mesh,IOField::IN);
      temp.scan(tm,NODE_FIELD);
   }
   size_t nb_time = tm.size();
   if (nb_time)
      cout << "Number of found time steps: " << nb_time << ", Running from "
           << "time: " << tm[0] << " to time: " << tm[nb_time-1] << "." << endl;
   else
      cout << "No time step was found." << endl;

   Vect<real_t> v(mesh);
   IOField io(mesh_file,input_file,mesh,IOField::IN);
   io.get(v,tm[0]);
   size_t n, m;
   ofstream fp(output_file.c_str());
   fp.setf(ios::right|ios::scientific);
   mesh_elements(mesh) {
      m = 0;
      if (The_element.getShape()==QUADRILATERAL)
         m = 4;
      if (The_element.getShape()==TRIANGLE)
         m = 3;
      if (!m)
         throw OFELIException("In saveGnuplot(string,string,string): "
                              " Unknown element shape " + itos(element_label));
      for (n=1; n<=m; n++) {
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
                 string mesh_file)
{
   cout << "Converting file: " << input_file << " to Tecplot format." << endl;
   Mesh mesh(mesh_file);
   vector<real_t> tm;
   {
      IOField temp(mesh_file,input_file,mesh,IOField::IN);
      temp.scan(tm,NODE_FIELD);
   }
   size_t nb_time = tm.size();
   cout << "Number of found time steps: " << nb_time << ", Running from "
        << "time: " << tm[0] << " to time: " << tm[nb_time-1] << "." << endl;

   vector<vector<real_t> > field(nb_time);
   IOField io(mesh_file,input_file,mesh,IOField::IN);
   string name;
   io.get(mesh,field,name);
   size_t nb_dof = io.getNbDOF();
   string shape;
   if (mesh.getShape()==LINE)
      shape = "LINESEG";
   else if (mesh.getShape()==QUADRILATERAL)
      shape = "QUADRILATERAL";
   else if (mesh.getShape()==TRIANGLE)
      shape = "TRIANGLE";
   else if (mesh.getShape()==TETRAHEDRON)
      shape = "TETRAHEDRON";
   else if (mesh.getShape()==HEXAHEDRON)
      shape = "BRICK";
   else
      ;

   ofstream fp(output_file.c_str());
   fp.setf(ios::right|ios::scientific);
   size_t count = 0;

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

      mesh_nodes(mesh) {
         if (count==0)
            for (size_t i=1; i<=mesh.getDim(); i++)
               fp << "  " << The_node.getCoord(i);
         for (size_t j=0; j<nb_dof; j++)
            fp << "  " << field[n][nb_dof*(node_label-1)+j];
         fp << endl;
      }
      if (count==0) {
         mesh_elements(mesh) {
            for (n=1; n<=The_element.getNbNodes(); n++)
               fp << setw(10) << The_element(n)->n();
            fp << endl;
         }
      }
      count++;
   }
   fp.close();
}


void saveVTK(string input_file,
             string output_file,
             string mesh_file)
{
   Mesh mesh(mesh_file);
   ofstream *pf=NULL;
   size_t m=0, k=0, nb_dof;
   string of;

   cout << "Converting file: " << input_file << " to VTK format." << endl;
   string proj=output_file.substr(0,output_file.rfind("."));
   vector<real_t> tm;
   {
      IOField temp(mesh_file,input_file,mesh,IOField::IN);
      temp.scan(tm,NODE_FIELD);
      nb_dof = temp.getNbDOF();
   }
   size_t nb_time = tm.size();
   if (nb_time)
      cout << "Number of found time steps: " << nb_time << ", Running from "
           << "time: " << tm[0] << " to time: " << tm[nb_time-1] << "." << endl;
   else {
      cout << "No time step was found." << endl;
      return;
   }

   vector<vector<real_t> > field(nb_time);
   for (size_t i=0; i<nb_time; i++)
      field[i].resize(nb_dof*mesh.getNbNodes());
   IOField io(mesh_file,input_file,mesh,IOField::IN);
   string name;
   io.get(mesh,field,name);

   bool scalar = true;
   if (nb_dof>1)
      scalar = false;
   size_t size=0;
   mesh_elements(mesh) {
      int s = The_element.getShape();
      if (s==LINE)
         m = 2;
      else if (s==TRIANGLE)
         m = 3;
      else if (s==QUADRILATERAL || s==TETRAHEDRON)
         m = 4;
      else if (s==HEXAHEDRON)
         m = 8;
      else if (s==PENTAHEDRON)
         m = 6;
      else ;
      size += m + 1;
   }
   for (size_t n=0; n<nb_time; n++) {
      string tt = "Time=" + dtos(tm[n]);
      if (nb_time==1)
         of = proj + ".vtk";
      else
         of = proj + "-" + zeros(n) + ".vtk";
      cout << "   Storing time step " << n+1 << " in file " << of << endl;
      pf = new ofstream(of.c_str());
      *pf << setprecision(16) << scientific;
      *pf << "# vtk DataFile Version 2.0\nImported from OFELI files\nASCII" << endl;
      *pf << "DATASET UNSTRUCTURED_GRID\nPOINTS " << mesh.getNbNodes() << " double" << endl;
      mesh_nodes(mesh)
         *pf << The_node.getX() << "  " << The_node.getY() << "  " << The_node.getZ() << endl;
      *pf << "\nCELLS " << mesh.getNbElements() << setw(10) << size << endl;
      mesh_elements(mesh) {
         int s = The_element.getShape();
         if (s==LINE)
            m = 2;
         else if (s==TRIANGLE)
            m = 3;
         else if (s==QUADRILATERAL || s==TETRAHEDRON)
            m = 4;
         else if (s==HEXAHEDRON)
            m = 8;
         else if (s==PENTAHEDRON)
            m = 6;
         else;
         *pf << setw(8) << m;
         for (size_t i=0; i<m; i++)
           *pf << setw(10) << The_element(i+1)->n()-1;
         *pf << endl;
      }
      *pf << "\nCELL_TYPES  " << mesh.getNbElements() << endl;
      k = 0;
      mesh_elements(mesh) {
         switch (The_element.getShape()) {
            case LINE:            m =  3; break;
            case TRIANGLE:        m =  5; break;
            case QUADRILATERAL:   m =  9; break;
            case TETRAHEDRON:     m = 10; break;
            case HEXAHEDRON:      m = 12; break;
            case PENTAHEDRON:     m = 13; break;
         }
         *pf << setw(4) << m;
         if (++k%30 == 0)
            *pf << endl;
      }
      *pf << "\nPOINT_DATA  " << mesh.getNbNodes() << endl;

      if (scalar)
         *pf << "SCALARS  u  double  1\nLOOKUP_TABLE  default" << endl;
      else
         *pf << "VECTORS  u  double" << endl;

      mesh_nodes(mesh) {
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
}


void saveGmsh(string input_file,
              string output_file,
              string mesh_file)
{
   Mesh mesh(mesh_file);
   char tt;
   size_t i, k, nb_dof=0, nb_en;
   ofstream pf(output_file.c_str());
   pf << setprecision(16);

   cout << "Converting file: " << input_file << " to Gmsh format." << endl;
   vector<real_t> tm;
   {
      IOField temp(input_file,mesh,IOField::IN);
      temp.scan(tm,NODE_FIELD);
      nb_dof = temp.getNbDOF();
   }

   size_t nb_time = tm.size();
   if (nb_time)
      cout << "Number of found time steps: " << nb_time << ", Running from "
           << "time: " << tm[0] << " to time: " << tm[nb_time-1] << "." << endl;
   else {
      cout << "No time step was found." << endl;
      return;
   }

   IOField io(mesh_file,input_file,mesh,IOField::IN);
   vector<vector<real_t> > field(nb_time);
   string name;
   io.get(mesh,field,name);
   tt = 'S';
   if (nb_dof == mesh.getDim())
      tt = 'V';

   pf << "View \"" << name << "\" {" << endl;
   size_t j;
   switch (mesh.getDim()) {

      case 1:
         mesh_elements(mesh) {
            pf << "SL(";
            pf << The_element(1)->getX() <<  ", 0., 0., "
               << The_element(2)->getX() <<  ", 0., 0. ) {" << endl;
            for (i=0; i<nb_time-1; i++)
               pf << field[i][The_element(1)->n()-1] << ","
                  << field[i][The_element(2)->n()-1] << ",";
            pf << field[nb_time-1][The_element(1)->n()-1] << ","
               << field[nb_time-1][The_element(2)->n()-1] << "};" << endl;
         }
         pf << "};" << endl;
         break;

      case 2:
         mesh_elements(mesh) {
            pf << tt;
            if ((nb_en=The_element.getNbNodes())==3)
               pf << "T";
            else if (nb_en==4)
               pf << "Q";
            pf << "(";
            for (k=1; k<nb_en; k++)
               pf << The_element(k)->getX() << "," << The_element(k)->getY() << ",0.,";
            pf << The_element(nb_en)->getX() << "," << The_element(nb_en)->getY() 
               << ",0.) {" << endl;
            for (i=0; i<nb_time; i++) {
               for (k=1; k<=nb_en; k++) {
                  j = nb_dof*(The_element(k)->n()-1);
                  pf << field[i][j];
                  if (nb_dof > 1)
                     pf << "," << field[i][j+1] << ",0.0";
                  if (i<nb_time-1 || k<nb_en)
                     pf << ",";
               }
               pf << endl;
            }
            pf << "};" << endl;
         }
         pf << "};" << endl;
         break;

      case 3:
         mesh_elements(mesh) {
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
            for (k=1; k<=nb_en-1; k++)
               pf << The_element(k)->getX() << ","
                  << The_element(k)->getY() << ","
                  << The_element(k)->getZ() << ",";
            pf << The_element(nb_en)->getX() << ","
               << The_element(nb_en)->getY() << ","
               << The_element(nb_en)->getZ() << ") {" << endl;
            for (i=0; i<nb_time; i++) {
               for (k=1; k<=nb_en; k++) {
                  pf << field[i][nb_dof*(The_element(k)->n()-1)];
                  if (nb_dof > 1)
                     pf << "," << field[i][nb_dof*(The_element(k)->n()-1)+1] << ","
                        << field[i][nb_dof*(The_element(k)->n()-1)+2] << endl;
                  if (i<nb_time-1 || k<nb_en)
                     pf << ",";
               }
            }
            pf << "};" << endl;
         }
         pf << "};" << endl;
         break;
   }
}

} /* namespace OFELI */
