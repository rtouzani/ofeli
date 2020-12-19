/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

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

                            Definition of Class LCL2DT
                Class to solve the linear conservation law in 2-D
                 using MUSCL finite volumes on triangular meshes
                
       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/

#ifndef  __LCL2DT_H
#define  __LCL2DT_H

#include "equations/cl/Muscl2DT.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */


/*! \file LCL2DT.h
 *  \brief Definition file for class LCL2DT.
 */

template<class T_,size_t N_> class LocalVect;
class Mesh;
class Side;

/*! \class LCL2DT
 *  \ingroup ConservationLaws
 *  \brief Class to solve the linear hyperbolic equation in 2-D by a MUSCL Finite Volume
 *  scheme on triangles
 *
 * \author S. Clain, V. Clauzon
 * \copyright GNU Lesser Public License
 */
class LCL2DT : public Muscl2DT {

public:

   using Muscl2DT::_nb_sides;
   using Muscl2DT::_nb_elements;

/// \brief Constructor using Mesh instance
    LCL2DT(Mesh& m);

/** \brief Constructor using mesh and initial data
    @param [in] m Reference to Mesh instance
    @param [in] U Vector containing initial (elementwise) solution 
 */
    LCL2DT(Mesh&         m,
           Vect<real_t>& U);

/// \brief Destructor
    ~LCL2DT();

/// \brief Return sidewise flux vector
    Vect<real_t> &getFlux() { return _FU; }

/// \brief Set elementwise initial condition
/// @param [in] u Vect instance containing initial condition values
    void setInitialCondition(Vect<real_t>& u);

/// \brief Set	a constant initial condition
/// @param [in] u Value of initial condition to assign to all elements
    void setInitialCondition(real_t u);

/// \brief Reconstruct flux using Muscl scheme
    void setReconstruction();

/// \brief Run one time step of the linear conservation law
/// @return Value of the time step
    real_t runOneTimeStep();

/// \brief Set Dirichlet boundary condition.
/// \details Assign a constant value <tt>u</tt> to all boundary sides
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
    void setVelocity(const LocalVect<real_t,2>& v);

/** \brief Computation of the primal variable n->n+1.
 *  \details Vector \a Flux contains elementwise fluxes issued from the Riemann problem, calculated
 *  with, as left element, \b getNeighborElement(1) and right element \b getNeighborElement(2)
 *  if \b getNeighborElement(2) doesn't exist, we are on a boundary and we prescribe a symmetry condition
 */
    void Forward(const Vect<real_t>& Flux,
                       Vect<real_t>& Field);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
//  add flux to value
    void forward();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void init();
    real_t getFluxUpwind();
    real_t getRiemannUpwind(size_t s);
    Vect<real_t> *_U, _LU, _RU, _FU, _v;
    bool _init_alloc;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
