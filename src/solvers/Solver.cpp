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

  ==============================================================================*/

#include "util/macros.h"

namespace OFELI {

int     theStep=1, theIteration=1, NbTimeSteps, MaxNbIterations=1000, Verbosity=1;
real_t  theTimeStep, theTime=0, theFinalTime;
real_t  theTolerance=1.e-8, theDiscrepancy=1.0;
bool     Converged=false;

} /* namespace OFELI */
