/*==============================================================================

                                ********************
                                *     E C 2 D 1    *
                                ********************

             A Finite Element Code for Quasi-Static Eddy Current Problems
                Model for 2-D Geometries with scalar magnetic field

  ------------------------------------------------------------------------------

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

  ==============================================================================*/


#include <string.h>
#include "OFELI.h"
#include "User.h"
#include "UpdateMF.h"
using namespace OFELI;

void ScaledMF(Mesh &ms, SkMatrix<complex_t> &a, Vect<complex_t> &b,
              complex_t &current, double omega, int flag);

int main(int argc, char *argv[])
{
   Mesh *ms0=NULL, *ms1=NULL;
   SkMatrix<complex_t> *A0=NULL, *A1=NULL;
   Vect<complex_t>     *b0=NULL, *b1=NULL;

   if (argc < 2) {
      cout << "\nUsage:  ec2d1 <parameter_file>\n";
      return 0;
   }

   IPF data("ec2d1 - 1.0",argv[1]);
   int output_flag = data.getOutput();
   int flag_volt = data.getInteger("flag_volt");
   double omega = data.getDouble("omega");
   complex_t volt = data.getComplexPar(1);
   complex_t current(1,0);

   if (!flag_volt)
      current = volt;

   if (output_flag) {
      cout << endl << endl;
      cout << "=====================================================================\n\n";
      cout << "                          E  C  2  D  1\n\n";
      cout << "            A Finite Element Code for Quasi-Static Analysis\n";
      cout << "      of Eddy Currents in 2-D Geometries with Scalar Magnetic Field\n\n";
      cout << "            EC2D1 uses OFELI Library of Finite Element Classes\n\n";
      cout << "                           V E R S I O N   1.0\n\n";
      cout << "                     Copyright R. Touzani, 1999\n\n";
      cout << "=====================================================================\n\n";
   }

//---------------------------------
// Read data
//---------------------------------

// Read Mesh data
   if (output_flag > 1)
      cout << "Reading mesh data ...\n";
   if (flag_volt)
      ms0 = new Mesh(data.getMeshFile(1));
   ms1 = new Mesh(data.getMeshFile(2));
//   IOField pl_file(data.getMesh(2),data.getPlotFile(),ms1,OUT);

   if (output_flag > 1 && flag_volt)
      cout << "Mesh Data of Inductor\n\n" << *ms0;

   if (output_flag > 1)
      cout << "Mesh Data of Conductor\n\n" << *ms1;

// Declare problem data (matrix, rhs, boundary conditions, body forces)
   if (output_flag > 1)
      cout << "Allocating memory for matrices and R.H.S. ...\n";
   A1 = new SkMatrix<complex_t>(*ms1);
   b1 = new Vect<complex_t>(ms1->getNbDOF());
   if (flag_volt) {
      A0 = new SkMatrix<complex_t>(*ms0);
      b0 = new Vect<complex_t>(ms0->getNbDOF());
   }

// Calculate scaled magnetic field
// -------------------------------

   if (output_flag > 1)
      cout << "Calculating scaled magnetic field ...\n";
   if (flag_volt)
      ScaledMF(*ms0,*A0,*b0,current,omega,flag_volt);
   ScaledMF(*ms1,*A1,*b1,current,omega,flag_volt);

// Calculate magnetic field in the free space and update it in the conductors
// --------------------------------------------------------------------------

   if (output_flag > 1)
      cout << "Updating magnetic field ...\n";
   if (flag_volt)
      UpdateMF(*ms0, *ms1, omega, volt, *b0, *b1);

// Output and save fields
// ----------------------

   if (flag_volt && output_flag > 0)
      cout << *b0;
   if (output_flag > 0)
      cout << *b1;
//   pl_file.put(u0);

   delete ms1;
   delete A1;
   delete b1;
   if (flag_volt) {
      delete ms0;
      delete b0;
      delete A0;
   }
#ifdef WITH_PAUSE
   system("PAUSE");
#endif
   return 0;
}
