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

                                     v m e s h

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
   cout << "\n";
   cout << "===========================================================================\n";
   cout << "vmesh, A Program to visualize an OFELI mesh file by gmsh\n";
   cout << "vmesh, version 1.0, Copyright (c) 1998 - 2022  Rachid Touzani\n";
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
   if (argc<2) {
      cout << "Usage: vmesh <mesh_file>\n" << endl;
      return EXIT_FAILURE;
   }
   string file = string(argv[1]).substr(0,string(argv[1]).rfind(".")) + ".msh";
   Mesh ms(argv[1]);
   ms.save(file);
   string com = "gmsh " + file;
   int ret = system(com.c_str());
   remove(file.c_str());
   cout << "done." << endl;
   return ret;
}
