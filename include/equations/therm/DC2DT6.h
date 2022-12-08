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

                            Definition of class DC2DT6
    for 2-D diffusion-convection scalar equation using 6-node triangular element

  ==============================================================================*/

#ifndef __DC2DT6_H
#define __DC2DT6_H

#include "equations/therm/Equa_Therm.h"
#include "shape_functions/Triang6S.h"
#include "shape_functions/Line3.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file DC2DT6.h
 *  \brief Definition file for class DC2DT6.
 */

/*! \class DC2DT6
 *  \ingroup Therm
 *  \brief Builds finite element arrays for thermal diffusion and convection in 2-D domains
 *  using 6-Node triangles.
 *
 *  \details Note that members calculating element arrays have as an argument a real
 *  <tt>coef</tt> that will be multiplied by the contribution of the current element.
 *  This makes possible testing different algorithms.
 *
 */

class DC2DT6 : virtual public Equa_Therm<6,6,3,3>
{

 public :

/// \brief Default Constructor.
/// \details Constructs an empty equation.
    DC2DT6();

/// \brief Constructor using Mesh data
/// @param [in] ms Mesh instance
    DC2DT6(Mesh& ms);

/** \brief Constructor using Mesh data and solution vector
 *  @param [in] ms Mesh instance
 *  @param [in,out] u Vect instance containing solution vector
 */
    DC2DT6(Mesh&         ms,
           Vect<real_t>& u);

/// \brief Destructor
    ~DC2DT6() { }

/// \brief Add lumped capacity matrix to element matrix after multiplying it by coefficient \c coef.
/// @param [in] coef Coefficient to multiply by added term (default value = 1).
    void LCapacity(real_t coef=1);

/// \brief Add Consistent capacity matrix to element matrix after multiplying it by coefficient \c coef.
/// @param [in] coef Coefficient to multiply by added term (default value = 1).
    void Capacity(real_t coef=1);

/// \brief Add diffusion matrix to element matrix after multiplying it by coefficient <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void Diffusion(real_t coef=1);

/** \brief Add convection matrix to left-hand side after multiplying it by coefficient <tt>coef</tt>
 *  \details Case where velocity field has been previouly defined
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Convection(real_t coef=1) { coef = 1; }

/** \brief Add convection matrix to left hand side after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  @param [in] v Constant velocity vector.
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Convection(Point<real_t>& v,
                    real_t         coef=1);

/** \brief Add convection matrix to left-hand side after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  \details Case where velocity field is given by a vector <tt>v</tt>
 *  @param [in] v Velocity vector.
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Convection(const Vect<real_t>& v,
                    real_t              coef=1);

/** \brief Add body right-hand side term to right hand side
 *  @param [in] f Local vector (of size <tt>6</tt>) containing source at nodes
 */
    void BodyRHS(const Vect<real_t>& f);

/** \brief Add boundary right-hand side term to right hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  @param [in] f Vector containing source at nodes
 */
    void BoundaryRHS(const Vect<real_t>& f);

 private:
   void set(const Element* el);
   void set(const Side* sd);
   Point<real_t> _x[6], _s[3];
   real_t        _a3;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
