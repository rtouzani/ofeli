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

                            Definition of class Pres1DL2
                         for 1-D wave propagation of pressure
       
  ==============================================================================*/


#ifndef __PRES1DL2_H
#define __PRES1DL2_H


#include "equations/acoustics/Equa_Acoustics.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Pres1DL2.h
 *  \brief Definition file for class PRES1DL2.
 */

/*! \class Pres1DL2
 *  \ingroup Acoustics
 *  \brief Builds finite element arrays for acoustic propagation in 1-D 
 *  using 2-Node elements.
 *
 *  \details Note that members calculating element arrays have as an argument a real
 *  <tt>coef</tt> that will be multiplied by the contribution of the current element. 
 *  This makes possible testing different algorithms.
 */

class Pres1DL2 : public Equa_Acoustics<2,2,1,1>
{

 public:

/// \brief Default Constructor.
/// \details Constructs an empty equation.
    Pres1DL2();

/// Constructor using mesh instance
/// @param [in] ms Mesh instance
    Pres1DL2(Mesh& ms);

/** Constructor using mesh instance and solution vector
 *  @param [in] ms Mesh instance
 *  @param [in,out] u Vect instance containing solution vector
 */ 
    Pres1DL2(Mesh&         ms,
             Vect<real_t>& u);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// Constructors for a mesh
    Pres1DL2(Mesh&         mesh,
             Vect<real_t>& b,
             real_t&       init_time,
             real_t&       final_time,
             real_t&       time_step)
    : Equation<2,2,1,1>(mesh,b,init_time,final_time,time_step) { }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Destructor
    ~Pres1DL2();

/** \brief Add lumped mass matrix to element matrix after multiplying it by coefficient
 *  <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void LMass(real_t coef=1);

/** \brief Add Consistent mass matrix to element matrix after multiplying it by 
 *  coefficient \c coef.
 *  @param [in] coef Coefficient to multiply by added term [default: 1]
 */
    void Mass(real_t coef=1);

/** \brief Add diffusion matrix to element matrix after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void Diffusion(real_t coef=1);

/** \brief Add body right-hand side term to right hand side.
 *  @param [in] f Vector containing source at nodes.
 */
    void BodyRHS(const Vect<real_t>& f);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Add boundary right-hand side flux to right hand side.
/// @param [in] flux Vector containing source at side nodes.
    void BoundaryRHS(real_t flux) { }

/** \brief Add boundary right-hand side term to right hand side after multiplying it by 
 *  coefficient <tt>coef</tt>
 *  @param [in] f Vector containing source at nodes.
 */
    void BoundaryRHS(const Vect<real_t>& f) { }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Return (constant) heat flux in element.
    real_t Flux() const;

/** \brief Set equation input data
 *  @param [in] opt Parameter that selects data type for input. This parameter
 *  is to be chosen in the enumerated variable EqDataType
 *  <ul>
 *     <li><tt>INITIAL_FIELD</tt>: Initial temperature
 *     <li><tt>BOUNDARY_CONDITION_DATA</tt>: Boundary condition (Dirichlet)
 *     <li><tt>SOURCE_DATA</tt>: Heat source
 *     <li><tt>FLUX_DATA</tt>: Heat flux (Neumann boundary condition)
 *     <li><tt>VELOCITY</tt>: Velocity vector (for the convection term)
 *  </ul>
 *  @param [in] u Vector containing input data
 */
    void setInput(EType         opt,
                  Vect<real_t>& u);

 private:
    void set(const Element* el);
    void set(const Side* sd) { }
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
