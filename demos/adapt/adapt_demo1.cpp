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

                       A demo program for mesh adaptation

          This program illustrate a 2-D mesh adaptation example where
              the mesh is adapted to a given solution vector

  ==============================================================================*/

#include "OFELI.h"
using namespace OFELI;

int main(int argc, char **argv)
{
   if (argc < 2) {
      cout << "Usage: ad1 <domain_file>" << endl;
      exit(1);
   }

// Generate background mesh
   try {
      Domain dom(argv[1]);
      MeshAdapt ma(dom);
      Mesh &ms = ma.getMesh();
//    Background mesh stored in file 'mesh1.m'
      ms.put("mesh1.m");

//    Set solution to adapt
      Vect<double> u(ms), v;
      u.set("exp(-100*x)+exp(-100*y)");
//    Solution on the background mesh stored for visualization with GMSH
      saveField(u,"u.pos",GMSH);

//    Generate adapted mesh
      ma.setHMax(0.1);
      ma.setError(0.04);
      ma.run(u,v);
      ms = ma.getMesh();

//    Adapted mesh stored in file 'mesh2.m'
//    Interpolated solution on the adapted mesh is stored for visualization with GMSH
      ms.put("mesh2.m");
      saveField(v,"v.pos",GMSH);
      cout << ma;
   } CATCH_EXCEPTION
   return 0;
}
