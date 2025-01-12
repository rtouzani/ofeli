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

     Definition of class 'Line2' for shape functions of 2-Node Line Element

  ==============================================================================*/


#ifndef __LINE2_H
#define __LINE2_H

#include "shape_functions/FEShape.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "mesh/Node.h"
#include "linear_algebra/LocalVect.h"
#include "linear_algebra/Point.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/** \defgroup Shape Shape Function
 *  \brief Shape function classes
 */

/*! \file Line2.h
 *  \brief Definition file for class Line2.
 */

/*! \class Line2 Line2.h "Line2.h"
 *  \brief To describe a 2-Node planar line finite element.
 *  \ingroup Shape
 *
 * \details Defines geometric quantities associated to 2-node linear segment element
 * <tt>P<sub>1</sub></tt> in the space. The reference element is the segment <tt>[-1,1]</tt>.
 * Note that the line length is not checked unless the function check is called.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class Line2 : public FEShape
{

 public:

/// \brief Default Constructor
    Line2();

/// \brief Constructor for an element
/// @param [in] el Pointer to element
    Line2(const Element* el);

/// \brief Constructor for a side
/// @param [in] side Pointer to side
    Line2(const Side* side);
   
/// \brief Constructor for an edge
/// @param [in] edge Pointer to edge
    Line2(const Edge* edge);
   
/// \brief Destructor
    ~Line2();

/// \brief Return element length
    real_t getLength() const { return _length; }

/// \brief Return unit normal vector to line
    Point<real_t> getNormal() const;

/// \brief Return unit tangent vector to line
    Point<real_t> getTangent() const;

/** \brief Calculate shape function of a given node at a given point.
 *  @param [in] i Node number (1 or 2).
 *  @param [in] s Localization of point in natural coordinates (must be between <tt>-1</tt>
 *  and <tt>1</tt>).
 */
    real_t Sh(size_t i,
              real_t s) const;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Calculate shape function of node <tt>i</tt> at a given point <tt>s</tt>
/// \details <tt>s</tt> must be between <tt>-1</tt> and <tt>1</tt>.
    real_t Sh(size_t i, Point<real_t> s) const { return Sh(i,s.x); }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Return partial derivatives of shape functions of element nodes
 *  @return LocalVect instance of partial derivatives of shape functions
 *          <i>e.g.</i> \c dsh(i).x, \c dsh(i).y, are partial derivatives of the <i>i</i>-th
 *          shape function. 
 */
    std::vector<Point<real_t> > DSh() const;

/// \brief Return reference coordinates of a point <tt>x</tt> in element
/// \details Only the x-coordinate of the returned value has a meaning
    Point<real_t> getRefCoord(const Point<real_t>& x);

/// \brief Check whether point <tt>x</tt> is in current line element or not
    bool isIn(const Point<real_t>& x);

/** \brief Return interpolated value at a given point
 *  @param [in] x Point where interpolation is evaluated (in the reference element).
 *  @param [out] v Computed value.
 */
    real_t getInterpolate(const Point<real_t>&       x, 
                          const LocalVect<real_t,2>& v);

 private:
    real_t _length;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
