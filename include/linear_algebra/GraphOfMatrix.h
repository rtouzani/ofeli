/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani

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

                  Various functions to determine matrix graph

  ==============================================================================*/


#ifndef __GRAPH_OF_MATRIX_H
#define __GRAPH_OF_MATRIX_H

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>

#include "mesh/Mesh.h"
#include "util/util.h"

using std::vector;
using std::pair;

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
typedef std::pair<size_t,size_t> RC;

using std::vector;
using std::sort;
using std::unique;


/// Set simple skyline structure of matrix (for node based unknowns)
    size_t SimpleSkyline(const Mesh&     mesh,
                         vector<size_t>& ch);

/// Set skyline structure of matrix (for node based unknowns)
    size_t NodeSkyline(const Mesh&     mesh,
                       vector<size_t>& ch);

/// Set skyline structure of matrix (for node based unknowns and selected DOF)
    size_t NodeSkyline(const Mesh&     mesh,
                       vector<size_t>& ch,
                       size_t          dof1,
                       size_t          dof2);

/// Set skyline structure of matrix (for node based unknowns)
    size_t NodeSkyline(const Mesh&     mesh,
                       vector<size_t>& ch,
                       size_t          dof);

/// Set skyline structure of matrix (for side based unknowns)
    size_t SideSkyline(const Mesh&     mesh,
                       vector<size_t>& ch);

/// Set skyline structure of matrix (for side based unknowns)
    size_t SideSkyline(const Mesh&     m,
                       vector<size_t>& ch,
                       size_t          dof1,
                       size_t          dof2);

/// Set skyline structure of matrix (for side based unknowns)
    size_t SideSkyline(const Mesh&     mesh,
                       vector<size_t>& ch,
                       size_t          dof);

/// Set skyline structure of matrix (for element based unknowns)
    size_t ElementSkyline(const Mesh&     mesh,
                          vector<size_t>& ch);

/// Set skyline structure of matrix (for element based unknowns)
    size_t ElementSkyline(const Mesh&     mesh,
                          vector<size_t>& ch,
                          size_t          dof);

/// Set a simple graph of matrix (for node based unknowns)
    size_t SimpleGraph(const Mesh&     mesh,
                       vector<long>&   xadj,
                       vector<size_t>& adjncy);

/// Set graph of matrix (for node based unknowns)
    size_t NodeGraph(const Mesh&     mesh,
                     vector<size_t>& row_ptr,
                     vector<size_t>& col_ind,
                     vector<RC>&     IJ,
                     vector<size_t>& nbc);

/// Set graph of matrix (for node based unknowns)
    size_t NodeGraph(const Mesh&     mesh,
                     size_t          dof1,
                     size_t          dof2,
                     vector<size_t>& row_ptr,
                     vector<size_t>& col_ind,
                     vector<RC>&     IJ,
                     vector<size_t>& nbc);

/// Set graph of matrix (for side based unknowns)
    size_t SideGraph(const Mesh&           mesh,
                           vector<size_t>& row_ptr,
                           vector<size_t>& col_ind,
                           vector<RC>&     IJ,
                           vector<size_t>& nbc);

/// Set graph of matrix (for side based unknowns)
    size_t SideGraph(const Mesh&     mesh,
                     size_t          dof1,
                     size_t          dof2,
                     vector<size_t>& row_ptr,
                     vector<size_t>& col_ind,
                     vector<RC>&     IJ,
                     vector<size_t>& nbc);

/// Set graph of matrix (for element based unknowns)
    size_t ElementGraph(const Mesh&     mesh,
                        vector<size_t>& row_ptr,
                        vector<size_t>& col_ind,
                        vector<RC>&     IJ,
                        vector<size_t>& nbc);

/// Set graph of matrix (for side/node coupling)
    size_t SideNodeGraph(const Mesh&     m,
                         vector<size_t>& row_ptr,
                         vector<size_t>& col_ind,
                         vector<RC>&     IJ,
                         vector<size_t>& nbc);

/// Set graph of matrix (for side/node coupling)
    size_t NodeSideGraph(const Mesh&     mesh,
                         vector<size_t>& row_ptr,
                         vector<size_t>& col_ind,
                         vector<RC>&     IJ,
                         vector<size_t>& nbc);

    size_t NodeGraphScal(const Mesh&     mesh,
                         vector<size_t>& row_ptr,
                         vector<size_t>& col_ind,
                         vector<RC>&     IJ,
                         vector<size_t>& nbc);

    size_t NodeGraphScal(const Mesh&     mesh,
                         size_t          dof,
                         size_t          nb_eq,
                         vector<size_t>& row_ptr,
                         vector<size_t>& col_ind,
                         vector<RC>&     IJ,
                         vector<size_t>& nbc);

    size_t XGraph(const Mesh&     mesh,
                  vector<size_t>& row_ptr,
                  vector<size_t>& col_ind,
                  vector<RC>&     IJ,
                  vector<size_t>& nbc);

/// Set graph of matrix (for DG P0/P0 matrix)
    size_t DG0Graph(const Mesh&     m,
                    vector<RC>&     I,
                    vector<size_t>& nbc);

/// Set graph of matrix for DG matrix
    size_t DGGraph(const Mesh& m,
    	             vector<RC>& I);

    void StoreGraph(const vector<RC>& IJ,
                    vector<size_t>&   row_ptr,
                    vector<size_t>&   col_ind);

    size_t getMatrixLength(const vector<size_t>& ch);

    size_t FinalizeGraph(vector<size_t>& row_ptr,
                         vector<size_t>& col_ind,
                         vector<RC>&     IJ,
                         vector<size_t>& nbc);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
