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

                       A demo program for mesh adaptation

           This program illustrate a 2-D mesh adaptation example where
          the mesh is adapted for the solution of a 2-D elliptic problem

  ==============================================================================*/

#include "OFELI.h"
#include "Therm.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
   void BC_RHS(int           opt,
               const Mesh&   ms,
               Vect<double>& bc,
               Vect<double>& f);

   void solve(Mesh&         ms,
              Vect<double>& f,
              Vect<double>& bc,
              Vect<double>& u);

// Read and output mesh data
   if (argc <= 1) {
      cout << " Usage: ad2 <opt>" << endl;
      cout << " opt = 0: Rectangle test (file: rect.dom)" << endl;
      cout << " opt = 1: L-shaped domain test (file: L.dom)" << endl;
      exit(1);
   }

   try {
      int opt = atoi(argv[1]);
      string domain_file = "rect.dom";
      if (opt>0)
         domain_file = "L.dom";
      Domain dom(domain_file);
      MeshAdapt ma(dom);
      Mesh ms0 = ma.getMesh();
      ms0.put("test.m");
      Mesh ms("test.m");

      Vect<double> bc(ms), f(ms);
      BC_RHS(opt,ms,bc,f);

//    Output solution
      Vect<double> u(ms);
      solve(ms,f,bc,u);
      saveField(u,"u.pos",GMSH);
      ms.put("mesh1.m");

//    Generate adapted mesh
      if (opt==0)
         ma.setError(0.005);
      else
         ma.setError(0.0002);
      Vect<double> v;
      ma.run(u,v);
      cout << ma;
      Mesh &msa = ma.getMesh();
      msa.put("mesh2.m");

//    New solution
      saveField(v,"v.pos",GMSH);
   } CATCH_EXCEPTION
   return 0;
}


void solve(Mesh&         ms,
           Vect<double>& f,
           Vect<double>& bc,
           Vect<double>& u)
{
   DC2DT3 eq(ms,u);
   eq.setInput(BOUNDARY_CONDITION,bc);
   eq.setInput(SOURCE,f);
   eq.run();
}


void BC_RHS(int           opt,
            const Mesh&   ms,
            Vect<double>& bc,
            Vect<double>& f)
{
   if (opt==0) {
      bc = 0;
      f.set("exp(-100*y^2)");
   }
   else {
      NodeLoop(ms) {
         double x = TheNode.getX(), y = TheNode.getY();
         double rr = pow(x*x+y*y,OFELI_THIRD);
         double t = 2*OFELI_THIRD*atan2(y,x);
         bc(theNodeLabel) = rr*sin(t);
         f(theNodeLabel) = 0;
      }
   }
}
