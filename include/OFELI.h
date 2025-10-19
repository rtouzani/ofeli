/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2025 Rachid Touzani

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

                 Header file that includes headers of kernel classes

  ==============================================================================*/

#ifndef __OFELI_H
#define __OFELI_H

/*! \file OFELI.h
 *  \ingroup Util
 *  \brief Header file that includes all kernel classes of the library.
 *
 *  To be included in conjunction with problem dependent header files.
 */

#include "util/banner.h"
#include "util/Timer.h"

#include "linear_algebra/Matrix_impl.h"
#include "linear_algebra/DMatrix_impl.h"
#include "linear_algebra/DSMatrix_impl.h"
#include "linear_algebra/SkMatrix_impl.h"
#include "linear_algebra/SkSMatrix_impl.h"
#include "linear_algebra/TrMatrix_impl.h"
#include "linear_algebra/SpMatrix_impl.h"
#include "linear_algebra/BMatrix_impl.h"
#include "linear_algebra/LocalMatrix_impl.h"
#include "linear_algebra/LocalVect_impl.h"

#include "mesh/Domain.h"
#include "mesh/MeshAdapt.h"
#include "mesh/MeshExtract.h"
#include "mesh/MeshUtil.h"
#include "mesh/Partition.h"
#include "mesh/getMesh.h"
#include "mesh/saveMesh.h"

#include "post/Reconstruction.h"
#include "post/Estimator.h"

#include "io/saveField.h"
#include "io/Prescription.h"
#include "io/IOField.h"
#include "io/IPF.h"

#include "solvers/GMRes.h"
#include "solvers/CG.h"
#include "solvers/CGS.h"
#include "solvers/BiCG.h"
#include "solvers/BiCGStab.h"
#include "solvers/Jacobi.h"
#include "solvers/SSOR.h"

#include "solvers/ODESolver.h"
#include "solvers/OptSolver.h"
#include "solvers/EigenProblemSolver.h"
#include "solvers/NLASSolver.h"
#include "solvers/LPSolver.h"
#include "solvers/TimeStepping.h"
#include "solvers/FuncApprox.h"
#include "solvers/LeastSquare.h"

#include "equations/interface/FastMarching.h"
#include "equations/interface/FastMarching1DG.h"
#include "equations/interface/FastMarching2DG.h"
#include "equations/interface/FastMarching3DG.h"

//#include "io/Tabulation.h"
//#include "io/Funct.h"

#ifdef USE_PETSC
#include "linear_algebra/petsc/PETScWrapper.h"
#include "linear_algebra/petsc/PETScMatrix.h"
#include "linear_algebra/petsc/PETScVect.h"
#endif

#endif