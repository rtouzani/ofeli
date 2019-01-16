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

                          Definition of class NS2DT3BT3
      for 2-D Navier-Stokes equations using P1-Bubble/P1 (Mini) finite element

  ==============================================================================*/


#ifndef __NS2DT3BT3_H
#define __NS2DT3BT3_H

#include "equations/fluid/Equa_Fluid.h"
#include "io/UserData.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/*! \file NS2DT3BT3.h
 *  \brief Definition file for class NS2DT3BT3.
 */

/*! \class NS2DT3BT3
 *  \ingroup Fluid
 *  \brief Builds finite element arrays for incompressible Navier-Stokes equations in 2-D
 *  domains using P<sub>1</sub>+Bubble/P<sub>1</sub> element.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class NS2DT3BT3 : virtual public Equa_Fluid<real_t,3,9,2,6> {

 public :

/// \brief Default Constructor
/// \details Builds an empty equation
    NS2DT3BT3() { }

/** \brief Constructor for a given mesh and initial solution
 *  @param [in] mesh Mesh instance
 *  @param [out] u Vector that contains velocity at nodes this one is obtained
 *  @param [in] Re Reynolds number. This is used if no viscosity is provided
 */
    NS2DT3BT3(Mesh&         mesh,
              Vect<real_t>& u,
              real_t        Re=1.0);

/// \brief Constructor using element data
/// @param [in] el Pointer to Element instance
    NS2DT3BT3(Element* el);

/// \brief Constructor using side data
/// @param [in] sd Pointer to Side instance
    NS2DT3BT3(Side* sd);

/** \brief Constructor using element and previous time data
 *  @param [in] el Pointer to Element instance
 *  @param [in,out] u Velocity vector. Will be updated with solution once this one is obtained
 *  @param [in,out] time Time value. Updated by time step value after this call
 */
    NS2DT3BT3(      Element*      el,
              const Vect<real_t>& u,
              const real_t&       time=0.);

/// \brief Constructor using side and previous time data
    NS2DT3BT3(      Side*         sd,
              const Vect<real_t>& u,
              const real_t&       time=0.);

/// \brief Destructor
    ~NS2DT3BT3();

/// \brief Add element lumped mass contribution to matrix after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    void LMass(real_t coef=1.);

/// \brief Add element viscous contribution to matrix after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    void Viscous(real_t coef=1.);

/// \brief Add element pressure gradient contribution to matrix after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    void PressureGradient(real_t coef=1.);

/// \brief Add convection contribution to left hand side after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
/// \details First term, explicit velocity, implicit velocity derivatives
    void LHS1_Convection(real_t coef=1.);

/// \brief Add convection contribution to left hand side after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
/// \details Second term, implicit velocity, explicit velocity derivatives
    void LHS2_Convection(real_t coef=1.);

/// \brief Add convection contribution to right hand side after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    void RHS_Convection(real_t coef=1.);

/// \brief Add body right-hand side term to right hand-side
/// @param [in] ud UserData instance that defines data 
    void BodyRHS(UserData<real_t>& ud);

/// \brief Add boundary right-hand side term to right-hand side
/// @param [in] ud UserData instance that defines data 
    void BoundaryRHS(UserData<real_t>& ud);

/// \brief Solve the problem
    int run();

protected :

   void set(const Element *el);
   void set(const Side *sd);

private :

   real_t _Re;
   LocalMatrix<real_t,4,4> _bb, _cc, _aa;
   void Misc();
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
