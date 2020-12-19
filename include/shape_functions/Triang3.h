/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

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

                            Definition of class 'Triang3'
                        (Shape functions for 3-node triangle)

  ==============================================================================*/


#ifndef __TRIANG3_H
#define __TRIANG3_H

#include <vector>
#include "shape_functions/FEShape.h"
#include "linear_algebra/LocalVect.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Triang3.h
 *  \brief Definition file for class Triang3.
 */

template<class T_,size_t N_> class LocalVect;
class Element;
class Side;


/*! \class Triang3
 *  \ingroup Shape
 *  \brief Defines a 3-Node (<tt>P<sub>1</sub></tt>) triangle.
 *  \details The reference element is the rectangle triangle with two unit edges.
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */

class Triang3 : public triangle
{

   using FEShape::_el;
   using FEShape::_sd;
 
public:

/// \brief Default Constructor
    Triang3();

/// \brief Constructor for an element.
/// \details The constructed triangle is an element in a 2-D mesh.
    Triang3(const Element* el);

/// \brief Constructor for a side.
/// \details The constructed triangle is a side in a 3-D mesh.
    Triang3(const Side* sd);

/// \brief Destructor
    ~Triang3();

/// \brief Choose element by giving its pointer
    void set(const Element* el);

/// \brief Choose side by giving its pointer
    void set(const Side* sd);

/** \brief Calculate shape function of node at a given point.
 *  @param [in] i Label (local) of node
 *  @param [in] s Natural coordinates of node where to evaluate
 */
    real_t Sh(size_t        i,
              Point<real_t> s) const
    {
       switch (i) {
          case 1 : return (1.-s.x-s.y);
          case 2 : return s.x;
          case 3 : return s.y;
       }
       return 0;
    }

/** \brief Return partial derivatives of shape functions of element nodes
 *  @return Vector of partial derivatives of shape functions
 *  <i>e.g.</i> \c dsh[i-1].x, \c dsh[i-1].y, are partial derivatives of the <i>i</i>-th
 *  shape function. 
 */
    std::vector<Point<real_t> > DSh() const;

/// \brief Return interpolated value at point of coordinate <tt>x</tt>
    real_t getInterpolate(const Point<real_t>&       x,
                          const LocalVect<real_t,3>& v);

/** \brief Check element area and number of nodes
 *  @return
 *  <ul>
 *    <li><tt> > 0</tt>: \a m is the length
 *    <li><tt> = 0</tt>: zero length (=> Error)
 *  </ul>
 */
    real_t check() const;

/** \brief Return constant gradient vector in triangle
 *  @param [in] u Local vector for which the gradient is evaluated
 */
    Point<real_t> Grad(const LocalVect<real_t,3>& u) const;

/// \brief Return maximal edge length of triangle
    real_t getMaxEdgeLength() const;

/// \brief Return minimal edge length of triangle
    real_t getMinEdgeLength() const;

};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
