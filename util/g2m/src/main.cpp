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
                                      g 2 m

            Program to generate 2-D meshes and store in OFELI mesh file
  ==============================================================================*/

#include "OFELI.h"
#include "io/tclap/CmdLine.h"

using namespace std;
using namespace TCLAP;
using namespace OFELI;


int main(int argc, char *argv[])
{
   cout << "\n";
   cout << "===========================================================================\n";
   cout << "g2m, A Program to generate 2-D meshes\n";
   cout << "g2m, version 2.0, Copyright (c) 1998 - 2018  Rachid Touzani\n";
   cout << "---------------------------------------------------------------------------\n";
   cout << "This program is free software: you can redistribute it and/or modify\n";
   cout << "it under the terms of the GNU Lesser General Public License as published by\n";
   cout << "the Free Software Foundation, either version 3 of the License, or\n";
   cout << "(at your option) any later version.\n\n";
   cout << "This program is distributed in the hope that it will be useful,\n";
   cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of\n";
   cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n";
   cout << "GNU Lesser General Public License for more details.\n\n";
   cout << "You should have received a copy of the GNU Lesser General Public License\n";
   cout << "along with this program. If not, see <http://www.gnu.org/licenses/>.\n";
   cout << "---------------------------------------------------------------------------\n" << endl;

   bool parse(int argc, char **argv, string &dom_file, string &geo_file, 
              string &output_file, bool& vtk_opt);
   string *file, dom_file=" ", geo_file=" ", output_file;
   bool vtk_opt=false;
   Domain *d = NULL;
   bool inter = parse(argc,argv,dom_file,geo_file,output_file,vtk_opt);
   if (inter) {
      Domain d;
      int ret = d.get();
      if (ret) {
         d.setOutputFile(output_file);
         if (ret>0)
            d.generateMesh();
         d.getMesh().put(output_file);
      }
   }
   else {
      if (geo_file==" ") {
         cout << "Processing with domain file: " << dom_file << endl;
         d = new Domain(dom_file);
         d->genMesh(output_file);
         cout << "Mesh generated and stored in file: " << output_file << endl;
         file = new string(dom_file.substr(0,dom_file.rfind(".")));
         geo_file = *file + ".geo";
         remove(geo_file.c_str());
         cout << "File " << geo_file << " deleted." << endl;
         string bamg_file = *file + ".bamg";
         remove(bamg_file.c_str());
         cout << "File " << bamg_file << " deleted." << endl;
         if (vtk_opt) {
            saveVTK(*file+".vtk",d->getMesh());
            cout << "Mesh stored in vtk file " << *file+".vtk" << endl;
         }
      }
      else {
         file = new string(geo_file.substr(0,geo_file.rfind(".")));
         string bamg_file = *file + ".bamg";
         cout << "Processing with geometry file: " << geo_file << endl;
         main_bamg(geo_file,bamg_file);
         Mesh ms;
         getBamg(bamg_file,ms,1);
         string ofeli_file = *file + ".m";
         ms.put(ofeli_file);
         cout << "Mesh generated and stored in file: " << ofeli_file << endl;
         remove(bamg_file.c_str());
         cout << "File " << bamg_file << " deleted." << endl;
         saveGmsh(*file+".msh",ms);
         cout << "Mesh stored in gmsh file " << *file+".msh" << endl;
         if (vtk_opt) {
            saveVTK(*file+".vtk",ms);
            cout << "Mesh stored in vtk file " << *file+".vtk" << endl;
         }
      }
      if (d)
         delete d;
      delete file;
   }
   return 0;
}


bool parse(int     argc,
           char**  argv,
           string& dom_file,
           string& geo_file,
           string& output_file,
           bool&   vtk_opt)
{
   bool inter = true;
   try {
      CmdLine cmd(" ",' ',"2.0");
      SwitchArg vtk("v","vtk","Save mesh in vtk file");
      ValueArg<string> output("o","output","Mesh Output File",false,"out.m","string");
      ValueArg<string> dom("d","domain","Domain Input File",false,"","string");
      ValueArg<string> geo("g","geo","Geometry Input File (bamg format)",false,"","string");

      cmd.add(vtk);
      cmd.add(dom);
      cmd.add(geo);
      cmd.add(output);
      cmd.parse(argc,argv);
      vtk_opt = vtk.getValue();

      if (dom.isSet() && geo.isSet())
         throw OFELIException("You cannot give domain and geometry input files.");

      if (dom.isSet() || geo.isSet())
         inter = false;

      string project;
      if (dom.isSet()) {
         dom_file = dom.getValue();
         project = dom_file.substr(0,dom_file.rfind("."));
      }
      else if (geo.isSet()) {
         geo_file = geo.getValue();
         project = geo_file.substr(0,geo_file.rfind("."));
      }
      else
         ;

      output_file = project + ".m";
      if (output.isSet())
         output_file = output.getValue();
      else
         output_file = "out.m";

   } catch (ArgException &e)
   { cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
   return inter;
}
