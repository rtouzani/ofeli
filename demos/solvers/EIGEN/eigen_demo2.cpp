/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

   This program is free software; you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free
   Software Foundation; Version 2 of the License.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
   details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the :

   Free Software Foqundation
   Inc., 59 Temple Place - Suite 330
   Boston, MA  02111-1307, USA

  ==============================================================================

                 An example of a Finite Element Code using OFELI

            Solution of an eigenvalue problem for the Laplace equation
                         in 2-D using P1 finite elements

  ==============================================================================*/

#include "OFELI.h"
#include "Laplace.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
// Expand arguments
   if (argc < 2) {
      cout << "Usage: " << argv[0] << " <project_file>" << endl;
      return 0;
   }

   try {
      banner("eigen_demo2");
      IPF data(argv[1]);
      int nb = data.getInteger("nb");
      Mesh ms(data.getMeshFile());

//    Solve the eigenvalue problem
      Laplace2DT3 eq(ms);
      EigenProblemSolver e(eq);
      e.run(nb);

//    Output eigenvalues, save eigenvectors for Gmsh post-processing
      Vect<double> v(ms);
      for (int i=1; i<=nb; i++) {
         cout << "Eigenvalue #" << i << ": " << e.getEigenValue(i) << endl;
         e.getEigenVector(i,v);
         saveField(v,ms,data.getPlotFile(i),GMSH);
      }
   } CATCH_EXCEPTION
   return 0;
}
