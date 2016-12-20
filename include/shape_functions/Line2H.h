/*==============================================================================

                                    O  F  E  L  I

                          Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2017 Rachid Touzani

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

             Definition of class 'Line2H' for shape functions of 2-Node
                          Hermite Line Element in Space

  ==============================================================================*/


#ifndef __LINE2H_H
#define __LINE2H_H

#include "shape_functions/FEShape.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Line2H.h
 *  \brief Definition file for class Line2H.
 */

class Element;
class Side;

/*! \class Line2H
 *  \ingroup Shape
 *  \brief To describe a 2-Node Hermite planar line finite element.
 *
 * \details Defines geometric quantities associated to 2-node segment element in
 * the space using Hermite (C<sup>1</sup>) interpolation. The interpolation functions
 * are polynomials of degree <tt>3</tt>. The reference element is the segment <tt>[-1,1]</tt>.
 * The unknowns are supported by extremities of the interval: each
 * extremity supports two unknowns, the function and its line derivative.
 *
 */

class Line2H : public FEShape
{

 public :

/// \brief Default Constructor
    Line2H();

/// \brief Constructor for an element
    Line2H(const Element* el);

/// \brief Constructor for a side
    Line2H(const Side* side);

/// \brief Destructor
    ~Line2H() { }

/** \brief Localize a point in the element.
 *  \details For a point <tt>s</tt> in the reference element, return coordinates
 *  in the real element.
 */
    Point<real_t> getLocalPoint(real_t s) const;

/// \brief Return shape function value of node <tt>i</tt> at given point <tt>s</tt>
    real_t Sh(size_t i,
              real_t s) const;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// Calculate shape function of node <tt>i</tt> at a given point <tt>s</tt>.
/// <tt>s</tt> must be between <tt>-1</tt> and <tt>1</tt>.
    real_t Sh(size_t i, Point<real_t> s) const { return Sh(i,s.x); }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Return first derivative (along the abscissa) of shape function of node
/// <tt>i</tt> at a given point
    real_t DSh(size_t i,
               real_t s) const;

/// \brief Return second derivatives (along the abscissa) of shape function of node <tt>i</tt>
    real_t D2Sh(size_t i,
                real_t s) const;

/// \brief Return determinant of jacobian
    real_t getDet() const { return _det; }

/// \brief Return element length
    real_t getLength() { return (_length=2*_det); }

/** \brief Check line length and number of line nodes
 *  @return
 *  <ul>
 *    <li><tt> > 0</tt>: <tt>m</tt> is the length
 *    <li><tt> = 0</tt>: zero length (=> Error)
 *  </ul>
 */
    real_t check() const;

 private:
    real_t _length;

};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
