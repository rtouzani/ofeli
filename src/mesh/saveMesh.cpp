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

                 Functions to save mesh in various file formats

  ==============================================================================*/

#include "mesh/saveMesh.h"
#include "OFELIException.h"

namespace OFELI {

void saveMesh(const string&      file,
              const Mesh&        mesh,
              ExternalFileFormat form)
{
   switch(form) {

      case GMSH:
         saveGmsh(file,mesh);
         break;

      case GNUPLOT:
         saveGnuplot(file,mesh);
         break;

      case MATLAB:
         saveMatlab(file,mesh);
         break;

      case TECPLOT:
         saveTecplot(file,mesh);
         break;

      case VTK:
         saveVTK(file,mesh);
         break;

      case OFELI_FF:
         mesh.put(file);
         break;

      case EASYMESH:
         throw OFELIException("saveMesh(string,Mesh,int): No saving to Easymesh format available");
         break;

      case GAMBIT:
         throw OFELIException("saveMesh(string,Mesh,int): No saving to Gambit format available");
         break;

      case BAMG:
         if (mesh.getNbSides()==0) {
            throw OFELIException("saveMesh(string,Mesh,int): Mesh sides must be extracted before "
                                 "converting to Bamg format");
            return;
         }
         saveBamg(file,mesh);
         break;

      case NETGEN:
         throw OFELIException("saveMesh(string,Mesh,int): No saving to Netgen format available");
         break;

      case TETGEN:
         throw OFELIException("saveMesh(string,Mesh,int): No saving to Tetgen format available");
         break;

      case TRIANGLE_FF:
         throw OFELIException("saveMesh(string,Mesh,int): No saving to Triangle format available");
         break;
   }
}


int getType(Element& e)
{
   int t = 0, sh = e.getShape(), n = e.getNbNodes();
   if (sh==LINE) {
      if (n==2)
         t = 1;
      else if (n==3)
         t = 8;
   }
   else if (sh==TRIANGLE) {
      if (n==3)
         t = 2;
      else if (n==6)
         t = 9;
   }
   else if (sh==QUADRILATERAL) {
      if (n==4)
         t = 3;
      else if (n==9)
         t = 10;
   }
   else if (sh==TETRAHEDRON) {
      if (n==4)
         t = 4;
      else if (n==10)
         t = 11;
   }
   else if (sh==HEXAHEDRON) {
      if (n==8)
         t = 5;
      else if (n==27)
         t = 12;
   }
   else if (sh==PENTAHEDRON) {
      if (n==6)
         t = 6;
      else if (n==18)
         t = 13;
   }
   return t;
}


int getType(Side& s)
{
   int t = 0, sh = s.getShape(), n = s.getNbNodes();
   if (sh==LINE) {
      if (n==2)
         t = 1;
      else if (n==3)
         t = 8;
   }
   else if (sh==TRIANGLE) {
      if (n==3)
         t = 2;
      else if (n==6)
         t = 9;
   }
   else if (sh==QUADRILATERAL) {
      if (n==4)
         t = 3;
      else if (n==9)
         t = 10;
   }
   return t;
}


void saveGmsh(const string& file,
              const Mesh&   mesh)
{
   static vector<int> type {0,0,1,2,3,4,5,6};
   size_t i;
   ofstream pf(file.c_str());
   if (pf.fail())
      throw OFELIException("saveGmsh(string,Mesh): cannot open file " + file);

   pf << "$MeshFormat\n2.2  0  8\n$EndMeshFormat" << endl;
   size_t n = 0;
   size_t nd=mesh.getNbNodes(), ne=mesh.getNbElements(), ns=mesh.getNbSides();
   pf << "$Nodes\n" << nd << endl;
   node_loop(&mesh)
      pf << node_label << "  " << the_node->getX() << "  "
         << the_node->getY() << "  " << the_node->getZ() << endl;
   pf << "$EndNodes" << endl;
   if (ne+ns>0)
      pf << "$Elements\n" << ne+ns << endl;
   element_loop(&mesh) {
      pf << ++n << " " << getType(The_element) << " 2 " << The_element.getCode() << " 1";
      for (i=1; i<=The_element.getNbNodes(); i++)
         pf << " " << The_element(i)->n();
      pf << endl;
   }
   side_loop(&mesh) {
      pf << ++n << " " << getType(The_side) << " 2 " << The_side.getCode() << " 1";
      for (i=1; i<=The_side.getNbNodes(); i++)
         pf << " " << The_side(i)->n();
      pf << endl;
   }
   pf << "$EndElements" << endl;
   pf.close();
}


void saveGnuplot(const string& file,
                 const Mesh&   mesh)
{
   if (mesh.getDim() != 2)
      throw OFELIException("saveGnuplot(string,Mesh):  This function is valid for 2-D only.");
   ofstream pf(file.c_str());
   if (pf.fail())
      throw OFELIException("saveGnuplot(string,Mesh): Cannot open file " + file);

   const Node *nd;
   std::map<int,int> sh = {{LINE,2},{QUADRILATERAL,4},{TRIANGLE,3}};
   element_loop(&mesh) {
      size_t m = sh[The_element.getShape()];
      if (m==0)
         throw OFELIException("saveGnuplot(string,Mesh): Illegal element geometry.");
      for (size_t n=1; n<=m; n++) {
         nd = The_element(n);
         pf << std::setprecision(5) << setw(18) << nd->getX()
            << std::setprecision(5) << setw(18) << nd->getY() << endl;
      }
      nd = The_element(1);
      pf << std::setprecision(5) << setw(18) << nd->getX()
         << std::setprecision(5) << setw(18) << nd->getY() << endl << endl;
   }
   pf.close();
}


void saveMatlab(const string& file,
                const Mesh&   mesh)
{
   if (mesh.getDim() != 2)
      throw OFELIException("saveMatlab(string,Mesh): This function is valid for 2-D only.");
   ofstream pf(file.c_str());
   if (pf.fail())
      throw OFELIException("saveMatlab(string,Mesh): Cannot open file " + file);

   element_loop(&mesh) {
      pf << "x = [ ";
      for (size_t n=1; n<=The_element.getNbNodes(); n++)
         pf << std::setprecision(5) << setw(18) << The_element(n)->getX() << " ";
      pf << std::setprecision(5) << setw(18) << The_element(1)->getX() << " ]; ";
      pf << "y = [ ";
      for (size_t m=1; m<=The_element.getNbNodes(); m++)
         pf << std::setprecision(5) << setw(18) << The_element(m)->getY();
      pf << std::setprecision(5) << setw(18) << The_element(1)->getY() << " ]; ";
      pf << "line(x,y);" << endl;
   }
   pf.close();
}


void saveTecplot(const string& file,
                 const Mesh&   mesh)
{
   string shape;
   static vector<string> sh {"","point","line","tria","quad","tetra","hexa","penta"};
   ofstream pf(file.c_str());
   if (pf.fail())
      throw OFELIException("saveTecplot(string,Mesh): Cannot open file " + file);

   size_t n, m;
   pf << "TITLE = \" \"\n";
   pf << "VARIABLES = \"X\", \"Y\"";
   if (mesh.getDim()==3)
      pf << ", \"Z\"";
   pf << endl;
   the_element = mesh.getPtrElement(1);
   if (mesh.getDim()==2) {
      if (the_element->getShape()==QUADRILATERAL) {
         m = 4;
         shape = "QUADRILATERAL";
      }
      if (the_element->getShape()==TRIANGLE) {
         m = 3;
         shape = "TRIANGLE";
      }
      else if (the_element->getShape()==LINE) {
         m = 2;
         shape = "LINESEG";
      }
   }
   else {
      if (the_element->getShape()==TETRAHEDRON) {
         m = 4;
         shape = "TETRAHEDRON";
      }
      else if (the_element->getShape()==HEXAHEDRON) {
         m = 8;
         shape = "BRICK";
      }
      else if (the_element->getShape()==LINE) {
         m = 2;
         shape = "LINESEG";
      }
   }
   pf << "ZONE N=" << mesh.getNbNodes() << ", E=" << mesh.getNbElements() 
      << ", F=FEPOINT, ET=" << shape << endl;
   node_loop(&mesh) {
      for (size_t i=1; i<=mesh.getDim(); i++)
         pf << std::setprecision(5) << setw(18) << The_node.getCoord(i) << " ";
      pf << endl;
   }

   element_loop(&mesh) {
      m = 0;
      if (The_element.getShape()==QUADRILATERAL)
         m = 4;
      else if (The_element.getShape()==TRIANGLE)
         m = 3;
      else if (The_element.getShape()==TETRAHEDRON)
         m = 4;
      else if (The_element.getShape()==HEXAHEDRON)
         m = 8;
      else if (The_element.getShape()==LINE)
         m = 2;
      else if (!m)
         throw OFELIException("saveTecplot(string,Mesh): Illegal element geometry");
      for (n=1; n<=m; n++)
         pf << The_element(n)->n() << "  ";
      pf << endl;
   }
   pf.close();
}


void saveVTK(const string& file,
             const Mesh&   mesh)
{
   size_t m=0, size, i;
   ofstream pf(file.c_str());
   static vector<int> nnd {0,1,2,3,4,4,8,6};
   static vector<int> etype {0,3,5,9,10,12,13};
   if (pf.fail())
      throw OFELIException("saveVTK(string,Mesh): Cannot open file " + file);

   pf << "# vtk DataFile Version 2.0\n";
   pf << "Imported from OFELI file" << endl;
   pf << "ASCII" << endl;
   pf << "DATASET UNSTRUCTURED_GRID\n";
   pf << "POINTS " << mesh.getNbNodes() << " float" << endl;
   node_loop(&mesh) {
      pf << std::setprecision(4) << setw(18) << The_node.getX()
         << std::setprecision(4) << setw(18) << The_node.getY()
         << std::setprecision(4) << setw(18) << The_node.getZ();
      pf << endl;
   }

   size = 0;
   element_loop(&mesh)
      size += nnd[The_element.getShape()] + 1;
   pf << "CELLS " << mesh.getNbElements() << " " << size << endl;
   element_loop(&mesh) {
      pf << nnd[The_element.getShape()];
      for (i=1; i<=m; i++)
         pf << "  " << The_element(i)->n()-1;
      pf << endl;
   }
   pf << "CELL_TYPES  " << mesh.getNbElements() << endl;
   size_t k = 0;
   element_loop(&mesh) {
      pf << "  " << etype[The_element.getShape()];
      if (++k%15 == 0)
         pf << endl;
   }
   pf << endl;
   pf.close();
}


void saveBamg(const string& file,
              const Mesh&   mesh)
{
   ofstream of(file.c_str());
   of << "MeshVersionFormatted  0\n\nDimension\n2\n" << endl;
   of << "Vertices\n" << mesh.getNbNodes() << endl;
   
   node_loop(&mesh)
      of << The_node.getX() << "  " << The_node.getY() << "  " << The_node.getCode(1) << endl;

   size_t nbs = 0;
   side_loop(&mesh)
      if (The_side.isOnBoundary())
         nbs++;
   of << "\nEdges\n" << nbs << endl;
   side_loop(&mesh) {
      if (The_side.isOnBoundary())
      of << The_side(1)->n() << "  " << The_side(2)->n() << "  " << The_side.getCode(1) << endl;
   }
   
   of << "\nTriangles\n" << mesh.getNbElements() << endl;
   element_loop(&mesh)
      of << The_element(1)->n() << "  " << The_element(2)->n() << "  "
         << The_element(3)->n() << "  " << The_element.getCode() << endl;
   of << "\nSubdomainFromMesh\n1" << endl;
   of << "3  1  1" << endl;
   of << "\nEnd" << endl;
}


} /* namespace OFELI */
