/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2018 Rachid Touzani

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

                           Prototypes for class 'Triang6S'
                    (Shape functions for 6-node straight triangle)

  ==============================================================================*/


#ifndef __TRIANG6S_H
#define __TRIANG6S_H

#include "shape_functions/FEShape.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Triang6S.h
 *  \brief Definition file for class Triang6S.
 */

class Element;
class Side;

/*! \class Triang6S
 *  \ingroup Shape
 *  \brief Defines a 6-Node straight triangular finite element using <tt>P<sub>2</sub></tt>
 *  interpolation
 *
 *  \details The reference element is the rectangle triangle with two unit edges.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class Triang6S : public triangle
{

   using FEShape::_el;
   using FEShape::_sd;

 public:

/// \brief Default Constructor
    Triang6S();

/** \brief Constructor for an element
 *  \details The constructed triangle is an element in a 2-D mesh.
 *  @param [in] el Pointer to Element instance
 */
    Triang6S(const Element* el);

/// \brief Destructor
    ~Triang6S() { }

/** \brief Calculate shape function of a node
 *  @param [in] i Local label of the node <tt>1&le;i&le;6</tt>
 *  @param [in] s Local coordinates of the point where the shape function is evaluated
 */
    real_t Sh(      size_t         i,
              const Point<real_t>& s) const;

/** \brief Calculate derivatives of shape function of a node
 *  @param [in] i Local label of node
 *  @param [in] s Local coordinates of the point where the gradient of the shape function 
 *                is evaluated
 */
    Point<real_t> DSh(      size_t         i,
                      const Point<real_t>& s) const;

/// \brief Return coordinates of center of element.
    Point<real_t> getCenter() const { return _c; }

/** \brief Return gradient vector in triangle at a given point
 *  @param [in] s Local coordinates of the point where the gradient of the shape function 
 *                is evaluated
 *  @param [in] u Local vector for which the gradient is evaluated
 */
    Point<real_t> Grad(const LocalVect<real_t,6>& u,
                       const Point<real_t>&       s) const;

/// \brief Return maximal edge length of triangle
    real_t getMaxEdgeLength() const;

/// \brief Return minimal edge length of triangle
    real_t getMinEdgeLength() const;

 private:
   Point<real_t> _x21, _x31, _x32;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
