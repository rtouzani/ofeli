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

                          Definition of class NS2DT3BT3
      for 2-D Navier-Stokes equations using P1-Bubble/P1 (Mini) finite element

  ==============================================================================*/


#ifndef __NS2DT3BT3_H
#define __NS2DT3BT3_H

#include "equations/fluid/Equa_Fluid.h"

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

class NS2DT3BT3 : virtual public Equa_Fluid<3,9,2,6> {

 public :

/// \brief Default Constructor
/// \details Builds an empty equation
    NS2DT3BT3() { }

/** \brief Constructor for a given mesh
 *  @param [in] ms Mesh instance
 */
    NS2DT3BT3(Mesh& ms);

/** \brief Constructor for a given mesh and initial solution
 *  @param [in] ms Mesh instance
 *  @param [in,out] u Vector that contains velocity at nodes this one is obtained
 */
    NS2DT3BT3(Mesh&         ms,
              Vect<real_t>& u);

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

/// \brief Add body force term to right hand-side
/// @param [in] f Vector containing body force at nodes
    void BodyRHS(Vect<real_t>& f);

/// \brief Add boundary traction  term to right-hand side
/// @param [in] f Vector containing body force at nodes
    void BoundaryRHS(Vect<real_t>& f);

 protected:

   void set(const Element *el);
   void set(const Side *sd);

 private:

   LocalMatrix<real_t,4,4> _bb, _cc, _aa;
   void Misc();
   void build();
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
