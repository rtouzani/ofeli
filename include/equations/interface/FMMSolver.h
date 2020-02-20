/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

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

           Definition of class 'FMMSolver' for Fast Marching Methods

  ==============================================================================*/

#ifndef __FMMSOLVER_H
#define __FMMSOLVER_H


#include "equations/interface/FMM2D.h"
#include "equations/interface/FMM3D.h"
#include "mesh/Grid.h"

/*!
 * \file FMM.h
 * \brief Fast Marching Method
 * \author Boris Meden, Mohamed Sylla
 */

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/** \defgroup Interface Interface Problems
 *  \brief Interface problems, including image processing
 */

/** \class FMMSolver
 *  \ingroup Interface
 *  \brief The Fast Marching Method solver
 *
 *  \details This class enables computing the signed distance function with respect to
 *  an interface. It works in 2-D and 3-D on a structured grid.
 *  The class is an interface for client. It points to FMM
 *
 * \author M. Sylla, B. Meden
 * \copyright GNU Lesser Public License
 */
class FMMSolver
{
 private:
   FMM *_theFM;
   Vect<real_t> *_phi;
   int _dim;

 public:

/** \brief Constructor
 * @param [in] g Instance of class Grid defining the grid on which the distance is computed.
 * @param [in] phi Vector containing the level set function at grid nodes.
 * The vector entries are 0 on the interface (from which the distance is computed),
 * positive on one side and negative on the other side.
 * They must contain the signed distance on the nodes surrounding the interface. These values
 * identify by linear interpolation the interface position.
 * The vector entries can take any value on other grid nodes, provided they have the right sign.
 * @param [in] ha true if high accuracy FMM is active. The high accuracy version is more accurate
 * but requires more accurate values on the nodes neighbouring the interface.
 */
    FMMSolver(const Grid&         g,
                    Vect<real_t>& phi,
                    bool          ha=false);

/// Destructor
    ~FMMSolver();

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** \brief check the compatibility between the grid and the vector
 *  \details The function throws an exception if the grid and vector dimension don't mach,
 *  and stops the program.
 *  @param [in] g Instance of class Grid
 */
    int CheckDimension(const Grid& g);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Execute the fast marching program
    void run() { this->_theFM->run(); }

/** \brief Extend speed by Sethian's method.
 *  \details The method consists in calculating a speed <tt>F</tt> such that
 *  its gradient is orthogonal to the gradient of the level set function
 *  @param [in,out] F Speed function where on input the value of the function
 *  is meaningful on the interface. On output <tt>F</tt> contains the extended speed
 */
    void ExtendSpeed(Vect<real_t>& F);

/** \brief Return the consistency error of the method
 *  \details Consistency is measured by computing the discrete value of the norm of the
 *  gradient of the signed distance and subtracting the obtained norm from 1.
 *  The absolute value of the result is returned.
 */
    real_t check_error() { return this->_theFM->check_error(); }
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
