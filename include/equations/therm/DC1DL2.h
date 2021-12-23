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

                            Definition of class DC1DL2
       for 1-D diffusion-convection scalar equation using 2-node line element

  ==============================================================================*/


#ifndef __DC1DL2_H
#define __DC1DL2_H


#include "equations/therm/Equa_Therm.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file DC1DL2.h
 *  \brief Definition file for class DC1DL2.
 */

/*! \class DC1DL2
 *  \ingroup Therm
 *  \brief Builds finite element arrays for thermal diffusion and convection in 1-D 
 *  using 2-Node elements.
 *
 *  \details Note that members calculating element arrays have as an argument a real
 *  <tt>coef</tt> that will be multiplied by the contribution of the current element. 
 *  This makes possible testing different algorithms.
 */

class DC1DL2 : public Equa_Therm<2,2,1,1>
{

 public:

/// \brief Default Constructor.
/// \details Constructs an empty equation.
    DC1DL2();

/// Constructor using mesh instance
/// @param [in] ms Mesh instance
    DC1DL2(Mesh& ms);

/** Constructor using mesh instance and solution vector
 *  @param [in] ms Mesh instance
 *  @param [in,out] u Vect instance containing solution vector
 */ 
    DC1DL2(Mesh&         ms,
           Vect<real_t>& u);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// Constructors for a mesh
    DC1DL2(Mesh&         mesh,
           Vect<real_t>& b,
           real_t&       init_time,
           real_t&       final_time,
           real_t&       time_step)
    : Equation<2,2,1,1>(mesh,b,init_time,final_time,time_step) { }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Destructor
    ~DC1DL2();

/** \brief Add lumped capacity matrix to element matrix after multiplying it by coefficient
 *  <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void LCapacity(real_t coef=1);

/** \brief Add Consistent capacity matrix to element matrix after multiplying it by 
 *  coefficient \c coef.
 *  @param [in] coef Coefficient to multiply by added term [default: 1]
 */
    void Capacity(real_t coef=1);

/** \brief Add diffusion matrix to element matrix after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void Diffusion(real_t coef=1);

/** \brief Add convection matrix to element matrix after multiplying it by coefficient
 *  <tt>coef</tt>
 *  @param [in] v Constant velocity vector
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void Convection(const real_t& v,
                    real_t        coef=1);

/** \brief Add convection matrix to element matrix after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  \details Case where velocity field is given by a vector <tt>v</tt>
 *  @param [in] v Velocity vector
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void Convection(const Vect<real_t>& v,
                    real_t              coef=1);

/** \brief Add convection matrix to element matrix after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  \details Case where velocity field has been previouly defined
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void Convection(real_t coef=1);

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
    void setInput(EqDataType    opt,
                  Vect<real_t>& u);

 private:
    void set(const Element* el);
    void set(const Side* sd) { }
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
