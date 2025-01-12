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

                           Definition of class DC3DAT3
       for 3-D diffusion-convection scalar equation in axisymmetric geometry
                         using 3-node triangular element

  ==============================================================================*/


#ifndef __DC3DAT3_H
#define __DC3DAT3_H


#include "equations/therm/Equa_Therm.h"
#include "linear_algebra/LocalMatrix.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file DC3DAT3.h
 *  \brief Definition file for class DC3DAT3.
 */

/*! \class DC3DAT3
 *  \ingroup Therm
 *  \brief Builds finite element arrays for thermal diffusion and convection in 3-D domains 
 *  with axisymmetry using 3-Node triangles.
 *
 *  \details Note that members calculating element arrays have as an argument a real
 *  \c coef that will be multiplied by the contribution of the current element. 
 *  This makes possible testing different algorithms.
 *
 */

class DC3DAT3 : virtual public Equa_Therm<3,3,2,2>
{

 public:

/// \brief Default Constructor.
/// \details Constructs an empty equation.
    DC3DAT3();

/// \brief Constructor using Mesh data
/// @param [in] ms Mesh instance
    DC3DAT3(Mesh& ms);

/** \brief Constructor using Mesh data and solution vector
 *  @param [in] ms Mesh instance
 *  @param [in,out] u Vect instance containing solution vector
 */
    DC3DAT3(Mesh&         ms,
            Vect<real_t>& u);

/// \brief Destructor
    ~DC3DAT3();

/// \brief Add lumped capacity matrix to element matrix after multiplying it by coefficient \c coef.
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void LCapacity(real_t coef=1);

/** \brief Add Consistent capacity matrix to element matrix after multiplying it by
 *  coefficient<tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Capacity(real_t coef=1);

/** \brief Add diffusion matrix to left-hand side after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Diffusion(real_t coef=1);

/** \brief Add diffusion matrix to left-hand side after multiplying it by coefficient <tt>coef</tt>
 *  \details Case where the diffusivity matrix is given as an argument
 *  @param [in] diff Instance of class DMatrix containing diffusivity matrix
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
    void Diffusion(const LocalMatrix<real_t,2,2>& diff,
                   real_t                         coef=1);

/** \brief Add body right-hand side term to right hand side.
 *  @param [in] f Local vector (of size <tt>3</tt>) containing source at odes.
 */
    void BodyRHS(const Vect<real_t>& f);

/// \brief Add boundary right-hand side term to right hand side.
/// @param [in] flux Value of flux to impose on the side
    void BoundaryRHS(real_t flux);

/** \brief Add boundary right-hand side term to right hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  @param [in] f Vector containing source at nodes
 */
    void BoundaryRHS(const Vect<real_t>& f);

/// \brief Return gradient of a vector in element.
/// @param [in] u Vector for which gradient is computed.
    Point<real_t> & Grad(const Vect<real_t>& u);

 private:
    void set(const Element* el);
    void set(const Side* sd);
    real_t        _r[3];
    Point<real_t> _grad;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif

