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

             Definition of class 'FMM2D' for 2-D Fast Marching Methods

  ==============================================================================*/

#ifndef __FMM2D_H
#define __FMM2D_H

/*!
 * \file FMM2D.h
 * \brief Fast Marching Method
 * \author Boris Meden, Mohamed Sylla
 */

#include "equations/interface/FMM.h"

namespace OFELI {

/*!
 * \class FMM2D
 * \brief class for the fast marching 2-D algorithm
 *
 * \details This class manages the 2-D Fast Marching method 
 */

class FMM2D : public FMM
{
 public:
/*!
 * \brief Constructor
 * \param [in] g Instance of class Grid
 * \param [in] phi Vector containing the level set function at grid nodes.
 * The values are <tt>0</tt> on the interface (from which the distance is computed),
 * positive on one side and negative on the other side. They must contain the
 * signed distance on the nodes surrounding the interface but can take any value on
 * other nodes, provided they have the right sign.
 * \param [in] HA true if the program must be executed with high accuracy,
 * false otherwise
 */
    FMM2D(const Grid&         g,
                Vect<real_t>* phi,
                bool          HA);

/// Initialize the heap
/// \param [in,out] NarrowPt Heap containing Narrow points
    void InitHeap(Heap& NarrowPt);

/// Execute Fast Marching Procedure
    void solve();

/** \brief compute the distance from node to interface
 * \param [in] pt node to treat
 * \param [in] sign Node sign
 * \return  distance from node <tt>pt</tt> to interface
 */
   void Evaluate(IPoint& pt,
                 int     sign);

/** \brief Extend the speed function to the whole grid
 *  \param [in,out] F Vector containing the speed at interface nodes on input
 *  and extended speed to all grid nodes
 */
   void ExtendSpeed(Vect<real_t>& F);

/// \brief Check error by comparing with the gradient norm.
/// \details This function returns discrete L<sup>2</sup> and Max. errors
    real_t check_error();

 private:
    
    void UpdateExt(size_t        i,
                   size_t        j,
                   real_t        dx,
                   real_t        dy,
                   Vect<real_t>& F);
};

} /* namespace OFELI */

#endif
