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

                        Implementation of class 'FMMSolver'

  ==============================================================================*/

#include "equations/interface/FMMSolver.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

FMMSolver::FMMSolver(const Grid&   g,
                     Vect<real_t>& phi,
                     bool          HA)
{
   _phi = &phi;
   _dim = CheckDimension(g);
   if (_dim==3)
      _theFM = new FMM3D(g,*_phi,HA);
   else
      _theFM = new FMM2D(g,*_phi,HA);
}


FMMSolver::~FMMSolver()
{
   if (_theFM)
      delete _theFM;
}


int FMMSolver::CheckDimension(const Grid& g)
{
   return 0;
/*   size_t ng=
   if ((g.getNz()!=0) && (g.getNz()+1==_phi->getNz()))
      return 3;
   else if ((g.getNz()==0) && (_phi->getNz()==0))
      return 2;
   else
      cerr << "The grid dimension doesn't match point dimension." << endl;
   exit(-1);*/
}


void FMMSolver::ExtendSpeed(Vect<real_t>& F)
{
   if (_dim==3)
      cerr << "Error: 3-D is not implemented." << endl;
   else
      _theFM->ExtendSpeed(F);
}

} /* namespace OFELI */
