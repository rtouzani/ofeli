/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

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

                            Definition of Class Muscl2DT
     Class for Muscl finite volumes for 2-D problems using triangular meshes

       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/

#ifndef __MUSCL2DT_H
#define __MUSCL2DT_H

#include <iostream>
using std::endl;

#include <valarray>
using std::valarray;

#include <vector>
using std::vector;

#include <string>
using std::string;

#include "equations/cl/Muscl.h"
#include "mesh/MeshUtil.h"
#include "linear_algebra/Vect.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */


/*! \file Muscl2DT.h
 *  \brief Definition file for class Muscl2DT.
 */

template<class T_,size_t N_> class LocalVect;
template<class T_> class Vect;

/*! \class Muscl2DT
 *  \ingroup ConservationLaws
 *  \brief Class for 2-D hyperbolic solvers with Muscl scheme.
 *
 * \author S. Clain, V. Clauzon
 * \copyright GNU Lesser Public License
 */

class Muscl2DT : public Muscl {

 public:

   using Muscl::_nb_sides;
   using Muscl::_nb_elements;
   using Muscl::_nb_nodes;

/// \brief Constructor using mesh
    Muscl2DT(Mesh &m);

/// \brief Destructor
    ~Muscl2DT();

/** \brief Function to reconstruct by the Muscl method
 *  @param [in] U Field to reconstruct
 *  @param [out] LU Left gradient vector
 *  @param [out] RU Right gradient vector
 *  @param [in] dof Label of dof to reconstruct
 */
    bool setReconstruction(const Vect<real_t>& U,
                           Vect<real_t>&       LU,
                           Vect<real_t>&       RU,
                           size_t              dof);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    real_t getMinimumFaceArea() const { return _MinimumFaceArea; }
    real_t getMinimumElementVolume() const { return _MinimumElementVolume; }
    real_t getMaximumFaceArea() const { return _MaximumFaceArea; }
    real_t getMaximumElementVolume() const { return _MaximumElementVolume; }
    real_t getMeanFaceArea() const { return _MeanFaceArea; }
    real_t getMeanElementVolume() const { return _MeanElementVolume; }
    real_t getMinimumEdgeLength() const { return _MinimumEdgeLength; }
    real_t getMinimumVolumebyArea() const { return _MinimumVolumebyArea; }
    real_t getMaximumEdgeLength() const { return _MaximumEdgeLength; }
    real_t getTauLim() const { return _taulim; }
    real_t getComega() const { return _Comega; }
    void setbetalim(real_t bl) { _betalim = bl; }
    void print_mesh_stat();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   Vect<Point<real_t> > _n;       // Normal vector
   Vect<real_t> _Lrate, _Rrate;   // Ki/Sij
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Construction of normals to sides.
 *  \details Convention: for a given side, getPtrElement(1) is the left element and getPtrElement(2)
 *  is the right element. The normal goes from left to right.
 *  For boundary sides, the normal points outward.
 */
    void Initialize();

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   void FirstOrder(const Vect<real_t> &U, Vect<real_t> &LU, Vect<real_t> &RU, size_t dof);
   void MultiSlopeM(const Vect<real_t> &U, Vect<real_t> &LU, Vect<real_t> &RU, size_t dof);
   void MultiSlopeQ(const Vect<real_t> &U, Vect<real_t> &LU, Vect<real_t> &RU, size_t dof);
   void GradientQ(const Vect<real_t> &U, Vect<real_t> &LU, Vect<real_t> &RU, size_t dof);
   void GradientM(const Vect<real_t> &U, Vect<real_t> &LU, Vect<real_t> &RU, size_t dof);
   void MultiSlopeMStable();

/** \brief Predictor for monoslope methods
 *  @param [in] U contains U1, U2, U3 for the three adjacent triangles
 *  @param [out] slope contains the gradient                      
 *  @param [in] n Label of the treated central triangle
 */
   void Grad(const LocalVect<real_t,3> &U, Point<real_t> &slope, size_t n);

/** \brief Calculate slope using least squares 
 *  @param [in] U contains U1, U2, U3 for the three adjacent triangles
 *  @param [out] slope contains the gradient
 *  @param [in] n Label of the treated central triangle
 */
    void LeastSquare(const LocalVect<real_t,3>& U,
                           Point<real_t>&       slope,
                           size_t               n);

// Algebra
   real_t _determ(Point<real_t> P, Point<real_t> Q) const { return (P.x*Q.y-P.y*Q.x); }

// B_iB_j, B_iQ_j, B_iM_j
   Vect<Point<real_t> > _BBv, _BQv, _BMv;

// norm of BiMij, B_iB_j and B_iQ_j 
   Vect<real_t> _BMl, _BBl, _BQl;
   
// projection coeff for multislope methods Q t1 versus t_2 and t_3
   Vect<real_t> _beta12, _beta13;

// projection coeff for multislope methods M s1 versus s_2 and s_3
   Vect<real_t> _betabis12, _betabis13;

// projection coeff for multislope M method s_i versus t_i
   Vect<real_t> _mu1, _mu2, _mu3;

// If a cell is not multislope convenient, this variable is set to 0
   Vect<int> _order;

// mesh characteristics
   real_t _MinimumFaceArea, _MinimumElementVolume, _MaximumFaceArea, _MaximumElementVolume;
   real_t _MeanFaceArea, _MeanElementVolume, _MinimumEdgeLength, _MinimumVolumebyArea;
   real_t _MaximumEdgeLength, _taulim, _Comega, _betalim;
   real_t getMinSize32();
   real_t gettaulim();
   real_t getComega();
   void getgraphtau();
   void getgraphComega();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
