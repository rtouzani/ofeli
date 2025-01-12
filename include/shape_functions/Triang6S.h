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

                           Prototypes for class 'Triang6S'
                    (Shape functions for 6-node straight triangle)

  ==============================================================================*/


#ifndef __TRIANG6S_H
#define __TRIANG6S_H

#include "shape_functions/FEShape.h"
#include "linear_algebra/LocalVect.h"

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

/** \brief Calculate shape functions
 *  @param [in] s Local first coordinate of the point where the gradient of the shape functions are evaluated
 *  @param [in] t Local second coordinate of the point where the gradient of the shape functions are evaluated
 *  @param [out] sh Array of of shape functions at \c (s,t)
 */
    void Sh(real_t  s,
            real_t  t,
            real_t* sh) const;

/// \brief Return coordinates of center of element.
    Point<real_t> getCenter() const { return _c; }

/// \brief Return maximal edge length of triangle
    real_t getMaxEdgeLength() const;

/// \brief Return minimal edge length of triangle
    real_t getMinEdgeLength() const;

/** \brief Initialize local point coordinates in element.
 *  @param [in] s Local first coordinate of the point where the gradient of the shape functions are evaluated
 *  @param [in] t Local second coordinate of the point where the gradient of the shape functions are evaluated
 */
    void setLocal(real_t s,
                  real_t t);

/** \brief Compute partial derivatives of shape functions at mid edges of triangles
 *  \details This member function can be called for integrations using partial derivatives
 *  of shape functions and approximated by midedge integration formula
 *  @param [out] dsh Vector containing partial derivatives of shape functions 
 *  @param [out] w Vector containing weights for the integration formula
 */
    void atMidEdges(std::vector<Point<real_t> >& dsh,
                    std::vector<real_t>&         w);

/** \brief Return partial derivatives of shape functions of element nodes
 *  @return LocalVect instance of partial derivatives of shape functions
 *          <i>e.g.</i> \c dsh(i).x, \c dsh(i).y, are partial derivatives of the <i>i</i>-th
 *          shape function. 
 *  @note The local point at which the derivatives are computed must be chosen before by using the
 *  member function setLocal
 */
    std::vector<Point<real_t> > DSh() const;

 private:
    Point<real_t> _x21, _x31, _x32;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
