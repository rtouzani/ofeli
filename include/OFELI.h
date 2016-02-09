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

#include "util/macros.h"
#include "util/banner.h"
#include "util/Timer.h"

#include "mesh/Domain.h"
#include "mesh/Mesh.h"
#include "mesh/MeshAdapt.h"
#include "mesh/MeshExtract.h"
#include "mesh/MeshUtil.h"
#include "mesh/Partition.h"
#include "mesh/getMesh.h"
#include "mesh/saveMesh.h"

#include "post/Reconstruction.h"
#include "post/Estimator.h"

#include "io/output.h"
#include "io/UserData.h"
#include "io/IPF.h"
#include "io/saveField.h"
#include "io/Prescription.h"
#include "io/IOField.h"

#include "util/Gauss.h"
#include "shape_functions/Line2.h"
#include "shape_functions/Line2H.h"
#include "shape_functions/Line3.h"
#include "shape_functions/Hexa8.h"
#include "shape_functions/Quad4.h"
#include "shape_functions/Tetra4.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Triang6S.h"
#include "shape_functions/Penta6.h"

#include "linear_algebra/Vect.h"
#include "linear_algebra/DMatrix.h"
#include "linear_algebra/DSMatrix.h"
#include "linear_algebra/SkMatrix.h"
#include "linear_algebra/SkSMatrix.h"
#include "linear_algebra/TrMatrix.h"
#include "linear_algebra/SpMatrix.h"
#include "linear_algebra/BMatrix.h"
#include "linear_algebra/LocalMatrix.h"
#include "linear_algebra/Assembly.h"

#ifndef USE_EIGEN
#include "solvers/GMRes.h"
#include "solvers/CG.h"
#include "solvers/CGS.h"
#include "solvers/BiCG.h"
#include "solvers/BiCGStab.h"
#include "solvers/Jacobi.h"
#include "solvers/SSOR.h"
#include "solvers/TimeStepping.h"
#endif

#include "solvers/ODESolver.h"
#include "solvers/OptimTN.h"
#include "solvers/OptimAux.h"
#include "solvers/EigenProblemSolver.h"
#include "equations/interface/FastMarching2D.h"

#include "io/Tabulation.h"
#include "io/Funct.h"

#ifdef USE_PETSC
#include "linear_algebra/petsc/PETScWrapper.h"
#include "linear_algebra/petsc/PETScMatrix.h"
#include "linear_algebra/petsc/PETScVect.h"
#endif

#endif
