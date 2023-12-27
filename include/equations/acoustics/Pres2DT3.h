/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2024 Rachid Touzani

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

                            Definition of class Pres2DT3
                for 2-D wave propagation of pressure using P1-triangles

  ==============================================================================*/


#ifndef __PRES2DT3_H
#define __PRES2DT3_H

#include "equations/acoustics/Equa_Acoustics.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Pres2DT3.h
 *  \brief Definition file for class Pres2DT3.
 */

/*! \class Pres2DT3
 *  \ingroup Acoustics
 *  \brief Builds finite element arrays for wave propagation in 2-D 
 *  using 3-Node elements.
 *
 *  \details Note that members calculating element arrays have as an argument a real
 *  <tt>coef</tt> that will be multiplied by the contribution of the current element. 
 *  This makes possible testing different algorithms.
 *
 */

class Pres2DT3 : public Equa_Acoustics<3,3,2,2>
{

 public:

    using Equa_Acoustics<3,3,2,2>::run;

/// \brief Default Constructor.
/// Constructs an empty equation.
    Pres2DT3();

/// \brief Constructor using Mesh data
/// @param [in] ms Mesh instance
    Pres2DT3(Mesh& ms);

/** \brief Constructor using Mesh and initial condition
 *  @param [in] ms Mesh instance
 *  @param [in] u Vect instance containing initial solution
 */
    Pres2DT3(Mesh&         ms,
             Vect<real_t>& u);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Constructors for a mesh
    Pres2DT3(Mesh&         mesh,
             Vect<real_t>& b,
             real_t        init_time,
             real_t        final_time,
             real_t        time_step)
   : Equation<3,3,2,2>(mesh,b,init_time,final_time,time_step) { }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Destructor
    ~Pres2DT3();

/** \brief Add lumped mass matrix to element matrix after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void LMass(real_t coef=1);

/** \brief Add Consistent mass matrix to element matrix after multiplying it 
 *  by coefficient <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
    void Mass(real_t coef=1);

/** \brief Add diffusion matrix to element matrix after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
    void Diffusion(real_t coef=1);

/** \brief Add body right-hand side term to right hand side.
 *  @param [in] f Vector containing source at nodes.
 */
    void BodyRHS(const Vect<real_t>& f);

/** \brief Add body right-hand side term to right hand side.
 *  \details Case where the body right-hand side is piecewise constant.
 *  @param [in] f Value of thermal source (Constant in element).
 */
    void BodyRHS(real_t f);

/// \brief Add boundary right-hand side flux to right hand side.
/// @param [in] flux Vector containing source at side nodes.
    void BoundaryRHS(real_t flux);

/** \brief Add boundary right-hand side term to right hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  @param [in] f Vector containing source at nodes
 */
    void BoundaryRHS(const Vect<real_t>& f);

/// \brief Return (constant) heat flux in element.
    Point<real_t>& Flux() const;

/// \brief Compute gradient of solution.
/// @param [in] g Elementwise vector containing gradient of solution.
    void Grad(Vect<Point<real_t> >& g);

/** \brief Return gradient of a vector in element.
 *  @param [in] u Global vector for which gradient is computed.
 *  Vector <tt>u</tt> has as size the total number of nodes
 */
    Point<real_t>& Grad(const Vect<real_t>& u) const;

/** \brief Set equation input data
 *  @param [in] opt Parameter to select type of input (enumerated values)
 *  <ul>
 *     <li><tt>INITIAL_FIELD</tt>: Initial temperature
 *     <li><tt>BOUNDARY_CONDITION_DATA</tt>: Boundary condition (Dirichlet)
 *     <li><tt>SOURCE_DATA</tt>: Heat source
 *     <li><tt>FLUX_DATA</tt>: Heat flux (Neumann boundary condition)
 *     <li><tt>SPEED</tt>: Speed vector (for the convection term)
 *  </ul>
 *  @param [in] u Vector containing input data
 */
    void setInput(EqDataType    opt,
                  Vect<real_t>& u);

 private:
    mutable Point<real_t> _f;
    void set(const Element* el);
    void set(const Side* sd);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
