/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2019 Rachid Touzani

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

     Definition of class 'Penta6' for shape functions of 6-node pentahedron

  ==============================================================================*/


#ifndef __PENTA6_H
#define __PENTA6_H

#include "shape_functions/FEShape.h"
#include "mesh/Element.h"
#include "mesh/Node.h"
#include "linear_algebra/Point.h"
#include "linear_algebra/LocalVect.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Penta6.h
 *  \brief Definition file for class Penta6.
 */

/*! \class Penta6
 *  \ingroup Shape
 *  \brief Defines a 6-node pentahedral finite element using <tt>P<sub>1</sub></tt>
 *  interpolation in local coordinates <tt>(s.x,s.y)</tt> and <tt>Q<sub>1</sub></tt>
 *  isoparametric interpolation in local coordinates <tt>(s.x,s.z)</tt> and <tt>(s.y,s.z)</tt>.
 *
 * \details The reference element is the cartesian product of the standard reference triangle 
 * with the line <tt>[-1,1]</tt>.
 * The nodes are ordered as follows:
 *    Node 1 in reference element is at s=(1,0,0) 
 *    Node 2 in reference element is at s=(0,1,0) 
 *    Node 3 in reference element is at s=(0,0,0) 
 *    Node 4 in reference element is at s=(1,0,1) 
 *    Node 5 in reference element is at s=(0,1,1) 
 *    Node 6 in reference element is at s=(0,0,1) 
 *
 * The user must take care to the
 * fact that determinant of jacobian and other quantities depend on the point in the
 * reference element where they are calculated. For this, before any utilization of
 * shape functions or jacobian, function \b setLocal() must be invoked.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class Penta6 : public FEShape
{

 public :

/// \brief Default Constructor
    Penta6();

/// \brief Constructor when data of Element <tt>el</tt> are given
/// @param [in] element Pointer to Element
    Penta6(const Element* element);

/// \brief Destructor
    ~Penta6() { }

/// \brief Choose element by giving its pointer
    void set(const Element* el);

/** \brief Initialize local point coordinates in element.
 *  @param [in] s Point in the reference element
 *  This function computes jacobian, shape functions and their partial
 *  derivatives at <tt>s</tt>. Other member functions only return these values.
 */
    void setLocal(const Point<real_t>& s);

/** \brief Return partial derivatives of shape functions of element nodes
 *  @return LocalVect instance of partial derivatives of shape functions
 *          <i>e.g.</i> \c dsh(i).x, \c dsh(i).y, are partial derivatives of the <i>i</i>-th
 *          shape function. 
 *  @note The local point at which the derivatives are computed must be chosen before by using the
 *  member function setLocal
 */
    vector<Point<real_t> > DSh() const;

/// \brief Return Maximum length of pentahedron edges
    real_t getMaxEdgeLength() const;

/// \brief Return Mimimum length of pentahedron edges
    real_t getMinEdgeLength() const;

};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
