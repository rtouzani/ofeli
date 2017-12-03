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

               Definition of class 'FMM' for Fast Marching Methods

  ==============================================================================*/

#ifndef __FMM_H
#define __FMM_H

/*!
 * \file FMM.h
 * \brief Fast Marching Method
 * \author Boris Meden, Mohamed Sylla
 * \date 14 March 2009
 */

#include "equations/interface/Heap.h"
#include <limits>

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS

class Grid;

/*!
 * \class FMM
 * \ingroup Interface
 * \brief Abstract class for the fast marching algorithm
 *
 * This abstract class manages the fast marching method
 *
 * \author M. Sylla, B. Meden
 * \copyright GNU Lesser Public License
 */
class FMM
{

 protected:
   Grid *_gd;
   real_t _inf;
   Vect<real_t> _AlivePt, _TAlive, *_phi;
   bool _high_accuracy;
   int _nx, _ny, _nz;
   int MaxQuadratic(real_t a, real_t b, real_t c, real_t& x);

 public:

/*!
 *  \brief Constructor.
 *  \details Constructor using grid and and initial interface
 *  @param [in] g Instance of class Grid
 *  @param [in] phi Vector containing the level set function at grid nodes. The values are 
 *  0 on the interface (from which the distance is computed), positive on one side and 
 *  negative on the other side. They must contain the signed distance on the nodes 
 *  surrounding the interface but can take any value on other nodes, provided they have 
 *  the right sign.
 *  @param [in] HA <tt>true</tt> if the program must be executed with high accuracy,
 *  \a false otherwise
 */
    FMM(const Grid&   g,
        Vect<real_t>& phi,
        bool          HA=false);

/// \brief Destructor
/// \details FMM class destructor
    virtual ~FMM();

/// \brief Initialize the heap
/// @param [in,out] NarrowPt Heap containing Narrow points
    virtual void InitHeap(Heap& NarrowPt) = 0;

/// \brief Execute Fast Marching Procedure
    virtual void solve() = 0;

/*!
 *  \brief compute the distance from node to interface
 *  @param [in] pt node to treat
 *  @param [in] sign node sign
 *  @return distance from node <tt>pt</tt> to interface
 */
    virtual void Evaluate(IPoint& pt,
                          int     sign) = 0;

/** \brief Extend speed by Sethian's method.
 *  \details The method consists in calculating a speed <tt>F</tt> such that
 *  its gradient is orthogonal to the gradient of the level set function
 */
    virtual void ExtendSpeed(Vect<real_t>& F) = 0;

/** \brief Check error by comparing with the gradient norm.
 *  \details This function prints discrete L<sup>2</sup> and Max. errors
 */
    virtual real_t check_error() = 0;
};

/** \brief Solve the quadratic equation and compute the maximal solution
 *  \details The quadratic equation is <tt>a*x<sup>2</sup> + b*x + c = 0</tt>
 *  @param [in] a coefficient <tt>a</tt>
 *  @param [in] b coefficient <tt>b</tt>
 *  @param [in] c coefficient <tt>c</tt>
 *  @param [out] x largest solution (if it exists)
 *  @return <tt>0</tt> if solution doesn't exist, <tt>1</tt> if it does
*/
int MaxQuadratic(real_t  a,
                 real_t  b,
                 real_t  c,
                 real_t& x);

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
