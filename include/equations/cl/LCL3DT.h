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

                            Definition of Class LCL3DT
              Class to solve the linear conservation law equation in 3-D
                  using MUSCL finite volumes on tetrahedral meshes

       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/

#ifndef  __LCL3DT_H
#define  __LCL3DT_H

#include "equations/cl/Muscl3DT.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file LCL3DT.h
 *  \brief Definition file for class LCL3DT.
 */

template<class T_,size_t N_> class LocalVect;
class Mesh;
class Side;

/*! \class LCL3DT
 *  \ingroup ConservationLaws
 *  \brief Class to solve the linear conservation law equation in 3-D by a MUSCL Finite Volume scheme
 *  on tetrahedra
 *
 * \author S. Clain, V. Clauzon
 * \copyright GNU Lesser Public License
 */
class LCL3DT : public Muscl3DT {

 public:

   using Muscl3DT::_nb_sides;
   using Muscl3DT::_nb_elements;

/// \brief Constructor using mesh
    LCL3DT(Mesh& m);

/** \brief Constructor using mesh and initial field
 *  @param [in] m Reference to Mesh instance
 *  @param [in] U Vector containing initial (elementwise) solution
 */
    LCL3DT(Mesh&         m,
           Vect<real_t>& U);

/// \brief Destructor
    ~LCL3DT();

/// \brief Set elementwise initial condition
/// @param [in] u Vect instance containing initial condition values
    void setInitialCondition(Vect<real_t>& u);

/// \brief Set	a constant initial condition
/// @param [in] u Value of initial condition to assign to all elements
    void setInitialCondition(real_t u);

/// \brief Reconstruct flux using Muscl scheme
    void setReconstruction();

/// \brief Run one time step
    real_t runOneTimeStep();

/// \brief Set Dirichlet boundary condition.
/// Assign a constant value <tt>u</tt> to all boundary sides
    void setBC(real_t u);

/** \brief Set Dirichlet boundary condition.
 *  \details Assign a constant value to a side
 *  @param [in] sd Side to which value is prescibed
 *  @param [in] u Value to prescribe
 */
    void setBC(const Side&  sd,
                     real_t u);

/** \brief Set Dirichlet boundary condition.
 *  \details Assign a constant value sides with a given code
 *  @param [in] code Code of sides to which value is prescibed
 *  @param [in] u Value to prescribe
 */
    void setBC(int    code,
               real_t u);

/// \brief Set convection velocity
/// @param [in] v Vect instance containing velocity
    void setVelocity(const Vect<real_t>& v);

/// \brief Set (constant) convection velocity
/// @param [in] v Vector containing constant velocity to prescribe
    void setVelocity(const LocalVect<real_t,3>& v);

/// \brief Assign reference length value
    void setReferenceLength(real_t dx) { _ReferenceLength=dx; }

/// \brief Return reference length
    real_t getReferenceLength() const { return _ReferenceLength; }

/** \brief Computation of the primal variable n->n+1.
 *  \details Vector <tt>Flux</tt> contains elementwise fluxes issued from the Riemann problem, 
 *  calculated with, as left element, \b getNeighborElement(1) and right element
 *  \b getNeighborElement(2) if \b getNeighborElement(2) doesn't exist, we are on a boundary 
 *  and we prescribe a symmetry condition
 */
    void Forward(const Vect<real_t>& Flux,
                       Vect<real_t>& Field);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
//  add flux to value
    void forward();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    Vect<real_t> *_U, _LU, _RU, _FU, _u, _v;
    void init();
    real_t getFluxUpwind();
    real_t getRiemannUpwind(size_t s);
    bool _init_alloc;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
