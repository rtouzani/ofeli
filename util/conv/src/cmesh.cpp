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

                                     c m e s h

                 A Program to convert various formats of mesh files

  ==============================================================================*/

#include <string>
#include "mesh/Mesh.h"
#include "mesh/getMesh.h"
#include "mesh/saveMesh.h"
#include "mesh/Domain.h"
#include "io/tclap/CmdLine.h"

using namespace std;
using namespace TCLAP;
using namespace OFELI;

int main(int argc, char **argv)
{
   void parse(int, char **, string &, string &, string &, string &, int &);
   int nb_dof;
   string project, input_format, output_format, input_file, output_file;

   cout << "\n\n";
   cout << "cmesh, A Program to convert various formats of mesh files\n";
   cout << "cmesh, version 2.0, Copyright (c) 1998 - 2018  Rachid Touzani\n\n";
   cout << "This program is free software: you can redistribute it and/or modify\n";
   cout << "it under the terms of the GNU Lesser General Public License as published by\n";
   cout << "the Free Software Foundation, either version 3 of the License, or\n";
   cout << "(at your option) any later version.\n\n";
   cout << "This program is distributed in the hope that it will be useful,\n";
   cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of\n";
   cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n";
   cout << "GNU Lesser General Public License for more details.\n\n";
   cout << "You should have received a copy of the GNU Lesser General Public License\n";
   cout << "along with this program.  If not, see <http://www.gnu.org/licenses/>.\n";

   parse(argc,argv,input_format,output_format,input_file,output_file,nb_dof);
   Mesh mesh;

//----------
// I N P U T
//----------

// OFELI File
   if (input_format == "xml" || input_format == "ofeli") {
      cout << "Reading mesh file " << input_file << " ...\n";
      mesh.get(input_file);
   }

// EasyMesh Files
   else if (input_format == "em") {
      cout << "Reading mesh files " << argv[1] << " ...\n";
      getEasymesh(input_file,mesh,nb_dof);
   }

// BAMG Files
   else if (input_format == "bamg") {
      cout << "Reading mesh file " << input_file << " ...\n";
      getBamg(input_file,mesh,nb_dof);
   }

// Gambit File
   else if (input_format == "gambit") {
      cout << "Reading mesh file " << input_file << " ...\n";
      getGambit(input_file,mesh,nb_dof);
   }

// Gmsh File
   else if (input_format == "gmsh") {
      cout << "Reading mesh file " << input_file << " ...\n";
      getGmsh(input_file,mesh,nb_dof);
   }

// Netgen File
   else if (input_format == "netgen") {
      cout << "Reading mesh files " << input_file << " ...\n";
      getNetgen(input_file,mesh,nb_dof);
   }

// Tetgen File
   else if (input_format == "tetgen") {
      cout << "Reading mesh files " << input_file << " ...\n";
      getTetgen(input_file,mesh,nb_dof);
   }

// Matlab File
   else if (input_format == "matlab") {
      cout << "Reading matlab file " << input_file << " ...\n";
      getMatlab(input_file,mesh,nb_dof);
   }

// Triangle Files
   else if (input_format == "triangle") {
      cout << "Reading triangle generic files " << input_file << " ...\n";
      getTriangle(input_file,mesh,nb_dof);
   }


//------------
// O U T P U T
//------------

// OFELI File
   if (output_format == "xml" || output_format == "ofeli") {
      cout << "Saving mesh file " << output_file << " ...\n";
      mesh.put(output_file);
   }

// Gnuplot File
   else if (output_format == "gpl" || output_format == "gnuplot") {
      cout << "Saving mesh file " << output_file << " ...\n";
      saveGnuplot(output_file,mesh);
   }

// Tecplot File
   else if (output_format == "tec") {
      cout << "Saving mesh file " << output_file << " ...\n";
      saveTecplot(output_file,mesh);
   }

// Matlab File
   else if (output_format == "matlab") {
      cout << "Saving mesh file " << output_file << " ...\n";
      saveMatlab(output_file,mesh);
   }

// Gmsh File
   else if (output_format == "gmsh") {
      cout << "Saving mesh file " << output_file << " ...\n";
      saveGmsh(output_file,mesh);
   }

// vtk File
   else if (output_format == "vtk") {
      cout << "Saving mesh file " << output_file << " ...\n";
      saveVTK(output_file,mesh);
   }

// Binary vtk File
   else if (output_format == "bvtk") {
      cout << "Saving mesh file " << output_file << " ...\n";
//      MDF2BVTK(output_file,mesh,1);
   }

   cout << "done." << endl;
   return 0;
}


void parse(int argc, char **argv, string &input_format, string &output_format, 
           string &input_file, string &output_file, int &nb_dof)
{
   const char help_in[]=
       "\nAvailable input formats:"
       "\n   ofeli     OFELI XML file (*.m): Same as xml"
       "\n   em        EasyMesh files (*.s *.e *.n)"
       "\n   bamg      BAMG file (*.bamg)"
       "\n   emc2      EMC2 Formatted file (*.msh)"
       "\n   gambit    Gambit Neutral file (*.neu)"
       "\n   gmsh      Gmsh file (*.msh)"
       "\n   netgen    Netgen file (*.vol)"
       "\n   tetgen    Tetgen files (*.node and *.ele)"
       "\n   matlab    Matlab file (*_matlab.m)"
       "\n   triangle  Triangle files (*.node and *.ele)";

   const char help_out[]=
       "\nAvailable output formats:"
       "\n   ofeli    XML file (*.m): Same as xml"
       "\n   gpl      Gnuplot file (*.gpl)"
       "\n   gmsh     Gmsh file (*.msh)"
       "\n   tec      TecPlot file"
       "\n   vtk      vtk file"
       "\n   matlab   Matlab file";

   try {
      CmdLine cmd(" ",' ',"2.3");

//    Nb. of dof
      ValueArg<string> nb("n","nb_dof","Nb. of degrees of freedom per node",false,"1","string");
      cmd.add(nb);

//    Input and output mesh files
      ValueArg<string> output("o","output","Mesh Output File",false,"","string");
      ValueArg<string> input("i","input","Mesh Input File",true,"","string");
      cmd.add(output);
      cmd.add(input);

//    Output format
      vector<string> allowed_out;
      allowed_out.push_back("xml");
      allowed_out.push_back("ofeli");
      allowed_out.push_back("gpl");
      allowed_out.push_back("gnuplot");
      allowed_out.push_back("gmsh");
      allowed_out.push_back("tec");
      allowed_out.push_back("tecplot");
      allowed_out.push_back("vtk");
      allowed_out.push_back("matlab");
      ValuesConstraint<string> allowedOutput(allowed_out);
      ValueArg<string> to("","to",help_out,true,"ofeli",&allowedOutput);
      cmd.add(to);

//    Input format
      vector<string> allowed_in;
      allowed_in.push_back("xml");
      allowed_in.push_back("ofeli");
      allowed_in.push_back("em");
      allowed_in.push_back("bamg");
      allowed_in.push_back("emc2");
      allowed_in.push_back("gambit");
      allowed_in.push_back("gmsh");
      allowed_in.push_back("netgen");
      allowed_in.push_back("tetgen");
      allowed_in.push_back("matlab");
      allowed_in.push_back("triangle");
      ValuesConstraint<string> allowedInput(allowed_in);
      ValueArg<string> from("","from",help_in,true,"ofeli",&allowedInput);
      cmd.add(from);

//    Parse the command line
      cmd.parse(argc,argv);

//    Set variables
      input_file = input.getValue();
      input_format = from.getValue();
      output_format = to.getValue();
      nb_dof = stringTo<int>(nb.getValue());
      if (output.isSet())
         output_file = output.getValue();
      else {
         string project = input_file.substr(0,input_file.rfind("."));
         output_file = project;
         if (output_format == "xml" || output_format == "ofeli")
            output_file += ".m";
         else if (output_format == "gpl")
            output_file += "-gnuplot.dat";
         else if (output_format == "gmsh")
            output_file += ".geo";
         else if (output_format == "tec")
            output_file = project + "-tecplot.dat";
         else if (output_format == "vtk")
            output_file += ".vtk";
        else if (output_format == "matlab")
            output_file += "-matlab.m";
        else
            throw("Error: Unknown output file format");
      }
   } catch ( ArgException& e )
   { cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
}
