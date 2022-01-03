/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani

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

                                    c f i e l d
                 A Program to convert various formats of field files

  ==============================================================================*/

#include "io/IOField.h"
#include "mesh/Mesh.h"
#include "io/saveField.h"
#include "io/tclap/CmdLine.h"

using namespace TCLAP;

using namespace OFELI;

void parse(int     argc,
           char**  argv,
           string& mesh_file,
           string& input,
           string& output,
           string& output_format,
           int&    freq);


int main(int argc, char **argv)
{
   string mesh_file, input_file, output_file, output_format, inp;
   int freq;
   cout << "\n\n";
   cout << "cfield, A Program to convert an XML OFELI file to other field file formats.\n";
   cout << "cfield, version 2.0, Copyright (c) 1998 - 2022 Rachid Touzani\n\n";
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

   parse(argc, argv, mesh_file, input_file, output_file, output_format, freq);
   ifstream inf(input_file.c_str());
   string cc;
   inf >> cc;
   inf.close();
   Mesh mesh(mesh_file);

// Gnuplot File
   if (output_format=="gpl") {
      cout << "Saving gnuplot file: " << output_file << " ...\n";
      saveGnuplot(input_file,output_file,mesh_file,freq);
   }

// Tecplot File
   else if (output_format=="tec" || output_format=="tecplot") {
      cout << "Saving Tecplot file: " << output_file << " ...\n";
      saveTecplot(input_file,output_file,mesh_file,freq);
   }

// vtk File
   else if (output_format=="vtk" || output_format=="paraview") {
      cout << "Saving vtk file: " << output_file << " ..." << endl;
      saveVTK(input_file,output_file,mesh_file,freq);
   }

// Gmsh File
   else if (output_format=="gmsh") {
      cout << "Saving Gmsh file: " << output_file << " ..." << endl;
      saveGmsh(input_file,output_file,mesh_file,freq);
   }
   cout << "done." << endl;
   return EXIT_SUCCESS;
}


void parse(int     argc,
           char**  argv,
           string& mesh_file,
           string& input_file, 
           string& output_file,
           string& output_format,
           int&    freq)
{
   const char help[]=
        "Available output formats:"
        "\n   gmsh: Gmsh Postprocessing File (*.pos)"
        "\n   gpl : Gnuplot File (*_gnuplot.dat)"
        "\n   tec : Tecplot file (*_tecplot.dat)"
        "\n   vtk : vtk file (*.vtk)";

   try {
      CmdLine cmd("cfield",' ',"2.2");
      ValueArg<string> output("o","output","Output field file name",false,"","string");
      ValueArg<string> input("i","input","Input field file name",true,"","string");
      ValueArg<string> mesh("m","mesh","Mesh file name",true,"","string");
      cmd.add(output);
      cmd.add(input);
      cmd.add(mesh);

//    Output format
      vector<string> allowed_out;
      allowed_out.push_back("gmsh");
      allowed_out.push_back("gpl");
      allowed_out.push_back("gnuplot");
      allowed_out.push_back("tec");
      allowed_out.push_back("tecplot");
      allowed_out.push_back("vtk");
      allowed_out.push_back("paraview");
      ValuesConstraint<string> allowedOutVals(allowed_out);
      ValueArg<string> format("f","format",help,true,"gmsh",&allowedOutVals);
      cmd.add(format);
      ValueArg<string> fr("s","saving_frequency","Saving frequency",false,"1","string");
      cmd.add(fr);

//    Parse the command line.
      cmd.parse(argc,argv);

//    Set variables
      output_format = format.getValue();
      mesh_file = mesh.getValue();
      input_file = input.getValue();
      freq = stringTo<int>(fr.getValue());
      string project = mesh_file.substr(0,mesh_file.rfind("."));
      if (mesh_file.substr(mesh_file.rfind(".")+1)=="xml")
         project = mesh_file.substr(0,mesh_file.rfind(".m"));
      if (output.isSet()==true)
         output_file = output.getValue();
      else  {
         if (output_format=="gmsh")
            output_file = project + ".pos";
         else if (output_format=="gpl" || output_format=="gnuplot")
            output_file = project + "-gnuplot.dat";
         else if (output_format=="tec" || output_format=="tecplot")
            output_file = project + "-tecplot.dat";
         else if (output_format=="vtk" || output_format=="paraview")
            output_file = project + ".vtk";
         else
            throw("Error: Output file name must be given.");
      }
   } catch ( ArgException& e )
   { cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
}
