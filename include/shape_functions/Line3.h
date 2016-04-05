/*==============================================================================

                                    O  F  E  L  I

                          Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzani

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

     Definition of class 'Line3' for shape functions of 3-Node Line Element

  ==============================================================================*/


#ifndef __LINE3_H
#define __LINE3_H

#include "shape_functions/FEShape.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "linear_algebra/Point.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Line3.h
 *  \brief Definition file for class Line3.
 */

class Element;
class Side;

/*! \class Line3
 *  \ingroup Shape
 * \brief To describe a 3-Node quadratic planar line finite element.
 *
 * Defines geometric quantities associated to 3-node quadratic
 * element <tt>P<sub>2</sub></tt> in the space. The reference element is the segment
 * <tt>[-1,1]</tt>. The user must take care to the fact that determinant
 * of jacobian and other quantities depend on the point in the
 * reference element where they are calculated. For this, before
 * any utilization of shape functions or jacobian, function \b setLocal()
 * must be invoked.\n
 * Element nodes are ordered as the following: the left one, the central one
 * and the right one.
 */

class Line3 : public FEShape
{

 public :

/// \brief Default Constructor
    Line3();

/// \brief Constructor for an element
    Line3(const Element* el);

/// \brief Constructor for a side
    Line3(const Side* sd);

/// \brief Destructor
    ~Line3() { }

/// \brief Initialize local point coordinates in element.
    void setLocal(real_t s);

/// \brief Return derivatives of shape function of node \a i at a given point
    real_t DSh(size_t i) const { return _dsh[i-1].x; }

/// \brief Return actual coordinates of localized point
    Point<real_t> getLocalPoint() const { return (_sh[0]*_x[0] + _sh[1]*_x[1]); }

};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
