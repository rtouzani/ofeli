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
   cout << "g2m, version 2.0, Copyright (c) 1998 - 2016  Rachid Touzani\n";
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

   void parse(int argc, char **argv, bool &inter, string &dom_file, string &geo_file, 
              string &output_file, string &w);
   string *file, dom_file=" ", geo_file=" ", output_file, w;
   Domain *d = NULL;
   bool inter = false;
   parse(argc,argv,inter,dom_file,geo_file,output_file,w);
   if (inter) {
      Domain d;
      if (d.get())
         d.generateMesh(output_file);
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
         if (w=="gmsh") {
            saveGmsh(*file+".msh",d->getMesh());
            cout << "Mesh stored in gmsh file " << *file+".msh" << endl;
         }
         if (w=="vtk") {
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
         if (w=="gmsh") {
            saveGmsh(*file+".msh",ms);
            cout << "Mesh stored in gmsh file " << *file+".msh" << endl;
         }
         if (w=="vtk") {
            saveVTK(*file+".vtk",ms);
            cout << "Mesh stored in vtk file " << *file+".vtk" << endl;
         }
      }
      cout << "===========================================================================\n";

      if (d)
         delete d;
      delete file;
   }
   return 0;
}


void parse(int     argc,
           char**  argv,
           bool&   inter,
           string& dom_file,
           string& geo_file,
           string& output_file,
           string& w)
{
   const char help_w[]="\nAvailable output formats:"
                       "\n   gmsh  Gmsh file (*.msh)"
                       "\n   vtk   vtk file (*.vtk)";

   try {
      CmdLine cmd(" ",' ',"1.0");
      vector<string> allowed_with;
      allowed_with.push_back("gmsh");
      allowed_with.push_back("vtk");
      allowed_with.push_back("none");
      ValuesConstraint<string> allowedWith(allowed_with);
      ValueArg<string> with("w","with",help_w,false,"none",&allowedWith);
      ValueArg<string> output("o","output","Mesh Output File",false,"out.m","string");

      SwitchArg interactive("i","interactive","Interactive mode",true);
      ValueArg<string> dom("d","domain","Domain Input File",true,"","string");
      ValueArg<string> geo("g","geo","Geometry Input File (bamg format)",true,"","string");
      vector<TCLAP::Arg*> xorlist;
      xorlist.push_back(&interactive);
      xorlist.push_back(&dom);
      xorlist.push_back(&geo);
      cmd.xorAdd(xorlist);

      cmd.add(with);
      cmd.add(output);
      cmd.parse(argc,argv);
      w = with.getValue();
      string project;
      if (dom.isSet()) {
         dom_file = dom.getValue();
         project = dom_file.substr(0,dom_file.rfind("."));
      }
      else if (geo.isSet()) {
         geo_file = geo.getValue();
         project = geo_file.substr(0,geo_file.rfind("."));
      }
      else if (interactive.isSet())
         inter = true;

      if (output.isSet())
         output_file = output.getValue();
      else
         output_file = project + ".m";

   } catch (ArgException &e)
   { cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
}
