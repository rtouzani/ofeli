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
         break;

      case EASYMESH:
         break;

      case GAMBIT:
         break;

      case BAMG:
         break;

      case NETGEN:
         break;

      case TETGEN:
         break;

      case TRIANGLE_FF:
         break;
   }
}


void saveGmsh(const string& file,
              const Mesh&   mesh)
{
  int type=0;
   size_t i;
   ofstream pf(file.c_str());
   if (pf.fail())
      throw OFELIException("saveGmsh(string,Mesh): cannot open file " + file);

   pf << "$NOD\n" << mesh.getNbNodes() << endl;
   mesh_nodes(mesh)
      pf << node_label << "  " << the_node->getX() << "  "
         << the_node->getY() << "  " << the_node->getZ() << endl;
   pf << "$ENDNOD" << endl;
   pf << "$ELM\n" << mesh.getNbElements() << endl;
   mesh_elements(mesh) {
      if (The_element.getShape()==LINE)
         type = 1;
      else if (The_element.getShape()==TRIANGLE)
         type = 2;
      else if (The_element.getShape()==QUADRILATERAL)
         type = 3;
      else if (The_element.getShape()==TETRAHEDRON)
         type = 4;
      else if (The_element.getShape()==HEXAHEDRON)
         type = 5;
      else if (The_element.getShape()==PENTAHEDRON)
         type = 6;
      else
         throw OFELIException("saveGmsh(string,Mesh): Element shape not recognized.");
      pf << element_label << "  " << type << "  " << The_element.getCode() << "  0  "
         << The_element.getNbNodes();
      for (i=1; i<=The_element.getNbNodes(); i++)
         pf << "  " << The_element.getNodeLabel(i);
      pf << endl;
   }
   mesh_sides(mesh) {
      if (The_side.getShape()==LINE)
         type = 1;
      else if (The_side.getShape()==TRIANGLE)
         type = 2;
      else if (The_side.getShape()==QUADRILATERAL)
         type = 3;
      else
         type = 0;
      pf << side_label << "  " << type << "  " << The_side.getCode(1) << "  0  "
         << The_side.getNbNodes();
      for (i=1; i<=The_side.getNbNodes(); i++)
         pf << "  " << The_side.getNodeLabel(i);
      pf << endl;
   }
   pf << "$ENDELM" << endl;
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
   size_t n, m;

   const Node *nd;
   mesh_elements(mesh) {
      m = 0;
      if (The_element.getShape()==LINE)
         m = 2;
      if (The_element.getShape()==QUADRILATERAL)
         m = 4;
      if (The_element.getShape()==TRIANGLE)
         m = 3;
      if (!m)
         throw OFELIException("saveGnuplot(string,Mesh): Illegal element geometry.");

      for (n=1; n<=m; n++) {
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

   mesh_elements(mesh) {
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
   string shape, sh[10];
   sh[LINE] = "line";
   sh[TRIANGLE] = "tria";
   sh[QUADRILATERAL] = "quad";
   sh[TETRAHEDRON] = "tetra";
   sh[HEXAHEDRON] = "hexa";
   sh[PENTAHEDRON] = "penta";
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
   mesh_nodes(mesh) {
      for (size_t i=1; i<=mesh.getDim(); i++)
         pf << std::setprecision(5) << setw(18) << The_node.getCoord(i) << " ";
      pf << endl;
   }

   mesh_elements(mesh) {
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
   if (pf.fail())
      throw OFELIException("saveVTK(string,Mesh): Cannot open file " + file);

   pf << "# vtk DataFile Version 2.0\n";
   pf << "Imported from XML (OFELI) file" << endl;
   pf << "ASCII" << endl;
   pf << "DATASET UNSTRUCTURED_GRID\n";
   pf << "POINTS " << mesh.getNbNodes() << " float" << endl;
   mesh_nodes(mesh) {
      pf << std::setprecision(4) << setw(18) << The_node.getX()
         << std::setprecision(4) << setw(18) << The_node.getY()
         << std::setprecision(4) << setw(18) << The_node.getZ();
      pf << endl;
   }

   size = 0;
   mesh_elements(mesh) {
      if (The_element.getShape()==LINE)
         m = 2;
      else if (The_element.getShape()==TRIANGLE)
         m = 3;
      else if (The_element.getShape()==QUADRILATERAL)
         m = 4;
      else if (The_element.getShape()==TETRAHEDRON)
         m = 4;
      else if (The_element.getShape()==HEXAHEDRON)
         m = 8;
      else if (The_element.getShape()==PENTAHEDRON)
         m = 6;
      else
         throw OFELIException("saveVTK(string,Mesh): Illegal element Geometry.");
      size += m+1;
   }
   pf << "CELLS " << mesh.getNbElements() << " " << size << endl;
   mesh_elements(mesh) {
      if (The_element.getShape()==LINE)
         m = 2;
      else if (The_element.getShape()==TRIANGLE)
         m = 3;
      else if (The_element.getShape()==QUADRILATERAL)
         m = 4;
      else if (The_element.getShape()==TETRAHEDRON)
         m = 4;
      else if (The_element.getShape()==HEXAHEDRON)
         m = 8;
      else if (The_element.getShape()==PENTAHEDRON)
         m = 6;
      else ;
      pf << m;
      for (i=1; i<=m; i++)
         pf << "  " << The_element(i)->n()-1;
      pf << endl;
   }
   pf << "CELL_TYPES  " << mesh.getNbElements() << endl;
   size_t k = 0;
   mesh_elements(mesh) {
      if (The_element.getShape()==LINE)
         m = 3;
      else if (The_element.getShape()==TRIANGLE)
         m = 5;
      else if (The_element.getShape()==QUADRILATERAL)
         m = 9;
      else if (The_element.getShape()==TETRAHEDRON)
         m = 10;
      else if (The_element.getShape()==HEXAHEDRON)
         m = 12;
      else if (The_element.getShape()==PENTAHEDRON)
         m = 13;
      else ;
      pf << "  " << m;
      if (++k%15 == 0)
         pf << endl;
   }
   pf << endl;
   pf.close();
}
   

void saveBamg(const string& file,
              Mesh&         mesh)
{
   mesh.getAllSides();
   ofstream of(file.c_str());
   of << "MeshVersionFormatted  0\n\nDimension\n2\n" << endl;
   of << "Vertices\n" << mesh.getNbNodes() << endl;
   
   mesh_nodes(mesh)
      of << The_node.getX() << "  " << The_node.getY() << "  " << The_node.getCode(1) << endl;

   size_t nbs = 0;
   mesh_sides(mesh)
      if (The_side.isOnBoundary())
         nbs++;
   of << "\nEdges\n" << nbs << endl;
   mesh_sides(mesh) {
      if (The_side.isOnBoundary())
      of << The_side(1)->n() << "  " << The_side(2)->n() << "  " << The_side.getCode(1) << endl;
   }
   
   of << "\nTriangles\n" << mesh.getNbElements() << endl;
   mesh_elements(mesh)
      of << The_element(1)->n() << "  " << The_element(2)->n() << "  "
         << The_element(3)->n() << "  " << The_element.getCode() << endl;
   of << "\nSubdomainFromMesh\n1" << endl;
   of << "3  1  1" << endl;
   of << "\nEnd" << endl;
}


} /* namespace OFELI */
