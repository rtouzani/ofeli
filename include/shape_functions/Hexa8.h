/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani

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

        Definition of class Hexa8 for shape functions of 8-Node hexahedral
                            element in three dimensions

  ==============================================================================*/


#ifndef __HEXA8_H
#define __HEXA8_H

#include <vector>

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

/*! \file Hexa8.h
 *  \brief Definition file for class Hexa8.
 */

/*! \class Hexa8
 *  \ingroup Shape
 *  \brief Defines a three-dimensional 8-node hexahedral finite element using Q1-isoparametric
 *  interpolation.
 *
 *  \details The reference element is the cube <tt>[-1,1]x[-1,1]x[-1,1]</tt>.
 *  The user must take care to the fact
 *  that determinant of jacobian and other quantities depend on the point in the
 *  reference element where they are calculated. For this, before any utilization of
 *  shape functions or jacobian, function \b getLocal(s) must be invoked.
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */

class Hexa8 : public FEShape
{

 public :

/// \brief Default Constructor
    Hexa8();

/// \brief Constructor when data of Element \a el are given
    Hexa8(const Element* el);

/// \brief Destructor
    ~Hexa8() { }

/** \brief Initialize local point coordinates in element.
 *  @param [in] s Point in the reference element
 *  This function computes jacobian, shape functions and their partial
 *  derivatives at <tt>s</tt>. Other member functions only return these values.
 */
    void setLocal(const Point<real_t>& s);

/** \brief Calculate shape function derivatives and integration weights
 *  @param [in] n Number of Gauss-Legendre integration points in each direction
 *  @param [in] dsh Vector of shape function derivatives at the Gauss points
 *  @param [in] w Weights of integration formula at Gauss points
 */
    void atGauss(int                          n,
                 std::vector<Point<real_t> >& dsh,
                 std::vector<real_t>&         w);

/** \brief Calculate shape functions and integration weights
 *  @param [in] n Number of Gauss-Legendre integration points in each direction
 *  @param [in] sh Vector of shape functions at the Gauss points
 *  @param [in] w Weights of integration formula at Gauss points
 */
    void atGauss(int                  n,
                 std::vector<real_t>& sh,
                 std::vector<real_t>& w);

/// \brief Return maximal edge length
    real_t getMaxEdgeLength() const;

/// \brief Return minimal edge length
    real_t getMinEdgeLength() const;

/** \brief Return gradient of a function defined at element nodes
 *  @param [in] u Vector of values at nodes
 *  @param [in] s Local coordinates (in <tt>[-1,1]*[-1,1]*[-1,1]</tt>) of point where
 *  the gradient is evaluated
 *  @return Value of gradient
 *  @note If the derivatives of shape functions were not computed before calling
 *  this function (by calling setLocal), this function will compute them
 */
    Point<real_t> Grad(const LocalVect<real_t,8>& u,
                       const Point<real_t>&       s);

 private:
    std::vector<Point<real_t> > _dsh;

};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
