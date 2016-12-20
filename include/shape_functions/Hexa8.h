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

      Definition of class Hexa8 for shape functions of 8-Node hexahedral
                            element in three dimensions

  ==============================================================================*/


#ifndef __HEXA8_H
#define __HEXA8_H

#include "shape_functions/FEShape.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Hexa8.h
 *  \brief Definition file for class Hexa8.
 */

template<class T_,size_t N_> class LocalVect;
template<class T_,size_t NR_,size_t NC_> class LocalMatrix;
class Element;
class Side;


/*! \class Hexa8
 *  \ingroup Shape
 *  \brief Defines a three-dimensional 8-node hexahedral finite element using Q1-isoparametric interpolation.
 *
 *  \details The reference element is the cube <tt>[-1,1]*[-1,1]*[-1,1]</tt>.
 *  The user must take care to the fact
 *  that determinant of jacobian and other quantities depend on the point in the
 *  reference element where they are calculated. For this, before any utilization of
 *  shape functions or jacobian, function \b getLocal(s) must be invoked.
 *
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

/** \brief Return <tt>x</tt>, <tt>y</tt> and <tt>z</tt> partial derivatives of shape 
 *  function of node <tt>i</tt> at a given point.
 *  \details Member function \a setLocal must have been called before in order to 
 *  calculate relevant quantities.
 */
    Point<real_t> DSh(size_t i) { return _dsh[i-1]; }

/// \brief Calculate shape function derivatives and integration weights
/// for 1-point Gauss rule
    void atGauss1(LocalVect<Point<real_t>,8>& dsh, real_t& w);

/// \brief Calculate shape function derivatives and integration weights
/// for 2x2x2-point Gauss rule
    void atGauss2(LocalMatrix<Point<real_t>,8,8>& dsh,
                  LocalVect<real_t,8>&            w);

/// \brief Return maximal edge length
    real_t getMaxEdgeLength() const;

/// \brief Return minimal edge length
    real_t getMinEdgeLength() const;

};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
