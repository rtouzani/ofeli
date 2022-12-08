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

     Definition of class 'Quad4' for shape functions of 4-node quadrilateral

  ==============================================================================*/


#ifndef __QUAD4_H
#define __QUAD4_H

#include <vector>

#include "shape_functions/FEShape.h"
#include "linear_algebra/LocalVect.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Quad4.h
 *  \brief Definition file for class Quad4.
 */

/*! \class Quad4
 *  \ingroup Shape
 *  \brief Defines a 4-node quadrilateral finite element using <tt>Q<sub>1</sub></tt>
 *  isoparametric interpolation.
 *
 * \details The reference element is the square <tt>[-1,1]x[-1,1]</tt>. The user must 
 * take care to the
 * fact that determinant of jacobian and other quantities depend on the point in the
 * reference element where they are calculated. For this, before any utilization of
 * shape functions or jacobian, function \b setLocal() must be invoked.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class Quad4 : public FEShape
{

 public:

/// \brief Default Constructor
    Quad4();

/// \brief Constructor when data of Element <tt>el</tt> are given
/// @param [in] element Pointer to Element
    Quad4(const Element* element);

/// \brief Constructor when data of Side <tt>sd</tt> are given
/// @param [in] side Pointer to Side
    Quad4(const Side *side);

/// \brief Destructor
    ~Quad4() { }

/// \brief Choose element by giving its pointer
    void set(const Element* el);

/// \brief Choose side by giving its pointer
    void set(const Side* sd);

/** \brief Initialize local point coordinates in element.
 *  @param [in] s Point in the reference element
 *  This function computes jacobian, shape functions and their partial
 *  derivatives at <tt>s</tt>. Other member functions only return these values.
 */
    void setLocal(const Point<real_t>& s);

/** \brief Calculate shape functions and their partial derivatives and integration weights
 *  @param [in] n Number of Gauss-Legendre integration points in each direction
 *  @param [in] sh Vector of shape functions at Gauss points
 *  @param [in] dsh Vector of shape function derivatives at Gauss points
 *  @param [in] w Weights of integration formula at Gauss points
 */
    void atGauss(int                          n,
                 std::vector<real_t>&         sh,
                 std::vector<Point<real_t> >& dsh,
                 std::vector<real_t>&         w);

/** \brief Return gradient of a function defined at element nodes
 *  @param [in] u Vector of values at nodes
 *  @param [in] s Local coordinates (in <tt>[-1,1]*[-1,1]</tt>) of point where
 *  the gradient is evaluated
 *  @return Value of gradient
 *  @note If the derivatives of shape functions were not computed before calling
 *  this function (by calling setLocal), this function will compute them
 */
    Point<real_t> Grad(const LocalVect<real_t,4>& u,
                       const Point<real_t>&       s);

/// \brief Return maximal edge length of quadrilateral
    real_t getMaxEdgeLength() const;

/// \brief Return minimal edge length of quadrilateral
    real_t getMinEdgeLength() const;

 private:
    bool _localized;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
