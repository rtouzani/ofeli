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

                            Definition of Class LCL1D
                Class to solve the linear conservation law equation
                         in 1-D using MUSCL finite volumes
                
       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/

#ifndef  __LCL1D_H
#define  __LCL1D_H

#include "equations/cl/Muscl1D.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */


/*! \file LCL1D.h
 *  \brief Definition file for class LCL1D.
 */


class Mesh;
class Side;

/*! \class LCL1D
 *  \ingroup ConservationLaws
 *  \brief Class to solve the linear conservation law (Hyperbolic equation)
 *  in 1-D by a MUSCL Finite Volume scheme
 *
 * \author S. Clain, V. Clauzon
 * \copyright GNU Lesser Public License
 */
class LCL1D : public Muscl1D {

public:

   using Muscl1D::_nb_sides;
   using Muscl1D::_nb_elements;
	
/// \brief Constructor using mesh instance
    LCL1D(Mesh& m);

/// \brief Constructor
    LCL1D(Mesh&         m,
          Vect<real_t>& U);

/// \brief Destructor
    ~LCL1D();

/// \brief Return sidewise fluxes
    Vect<real_t> &getFlux() { return _FU; }

/// \brief Assign initial condition by a vector
/// @param [in] u Vector containing initial condition
    void setInitialCondition(Vect<real_t>& u);

/// \brief Assign a constant initial condition
/// @param [in] u Constant value for the initial condition 
    void setInitialCondition(real_t u);

/// \brief Run MUSCL reconstruction
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
    void setVelocity(Vect<real_t>& v);

/// \brief Set (constant) convection velocity
    void setVelocity(real_t v);

/// \brief Assign reference length value
    void setReferenceLength(real_t dx) { _ReferenceLength = dx; }

/// \brief Return reference length
    real_t getReferenceLength() const { return _ReferenceLength; }

/** \brief Computation of the primal variable n->n+1.
 *  \details Vector \b Flux contains elementwise fluxes issued from the Riemann problem, calculated
 *  with, as left element, \b getNeighborElement(1) and right element \b getNeighborElement(2)
 *  if \b getNeighborElement(2) doesn't exist, we are on a boundary and we prescribe a
 *  symmetry condition
 */
    void Forward(const Vect<real_t>& Flux,
                       Vect<real_t>& Field);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
//  add flux to value
    void forward();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 protected:
 
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    bool _init_alloc;
    void _Init();
    real_t getFluxUpwind();
    real_t computeRiemannUpwind(size_t);
    Vect<real_t> *_U, _LU, _RU, _FU, _v;
    real_t _ReferenceLength;
    bool _alloc_init;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
