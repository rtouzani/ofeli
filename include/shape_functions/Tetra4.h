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

      Definition of class 'Tetra4' for shape functions of 4-node tetrahedron

  ==============================================================================*/


#ifndef __TETRA4_H
#define __TETRA4_H

#include "shape_functions/FEShape.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Tetra4.h
 *  \brief Definition file for class Tetra4.
 */

template<class T_,size_t N_> class LocalVect;
class Element;
class Side;

/*! \class Tetra4
 *  \ingroup Shape
 *  \brief  Defines a three-dimensional 4-node tetrahedral finite element
 *  using <tt>P<sub>1</sub></tt> interpolation.
 *
 * The reference element is the right tetrahedron with four unit edges interpolation.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class Tetra4 : public FEShape
{

public :

/// \brief Default Constructor
    Tetra4();

/// \brief Constructor when data of Element <tt>el</tt> are given
    Tetra4(const Element* el);

/// \brief Destructor
    ~Tetra4() { }

/// \brief Choose element by giving its pointer
    void set(const Element* el);

/// \brief Calculate shape function of node <tt>i</tt> at a given point <tt>s</tt>.
/// \details <tt>s</tt> is a point in the reference tetrahedron.
    real_t Sh(size_t        i,
              Point<real_t> s) const;

/// \brief Return <tt>x</tt>, <tt>y</tt> and <tt>z</tt> partial derivatives of shape 
/// function associated to node <tt>i</tt>.
/// \details Note that these are constant in element.
    Point<real_t> DSh(size_t i) const { return _dsh[i-1]; }

/// \brief Return volume of element
    real_t getVolume() const { return OFELI_SIXTH*_det; }

/// \brief Return reference coordinates of a point <tt>x</tt> in element.
    Point<real_t> getRefCoord(const Point<real_t>& x) const;

/// \brief Check whether point <tt>x</tt> is in current tetrahedron or not
    bool isIn(const Point<real_t>& x);

/// \brief Return interpolated value at point of coordinate <tt>x</tt>
    real_t getInterpolate(const Point<real_t>&       x,
                          const LocalVect<real_t,4>& v);

/** \brief Return edge shape function
 *  @param [in] k Local edge number for which the edge shape function is computed 
 *  @param [in] s Local coordinates in element
 *  @remark Element edges are ordered as follows: Edge <tt>k</tt> has end vertices
 *  <tt>k</tt> and <tt>k+1</tt>
 */
    Point<real_t> EdgeSh(size_t        k,
                         Point<real_t> s);

/** \brief Return curl of edge shape function
 *  @param [in] k Local edge number for which the curl of the edge shape function is computed 
 *  @remark Element edges are ordered as follows: Edge <tt>k</tt> has end vertices <tt>k</tt>
 *  and  <tt>k+1</tt>
 */
    Point<real_t> CurlEdgeSh(size_t k);

/// \brief Return maximal edge length of tetrahedron
    real_t getMaxEdgeLength() const;

/// \brief Return minimal edge length of tetrahedron
    real_t getMinEdgeLength() const;

 private:
   void CalculateShape();

};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
