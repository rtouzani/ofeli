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

                        Definition of function BSpline
                   implementing the B-Spline curve algorithm

               Copied and adapted from the program by Keith Vertanen

  ==============================================================================*/

#ifndef __BSPLINE_H
#define __BSPLINE_H


#include "OFELI_Config.h"
#include "linear_algebra/Point.h"
#include "linear_algebra/Vect.h"


namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

void compute_intervals(Vect<size_t>& u,
                       size_t        n,
                       size_t        t);
real_t blend(size_t              k,
             size_t              t,
             const Vect<size_t>& u,
             real_t              v);
void compute_point(const Vect<size_t>&         u,
                   size_t                      n,
                   size_t                      t,
                   real_t                      v,
                   const Vect<Point<real_t> >& control,
                   Point<real_t>&              output);


/*! \file BSpline.h
 *  \brief Function to perform a B-spline interpolation
 */

/** \fn BSpline(size_t n, size_t t, Vect<Point<real_t> >& control, 
 *              Vect<Point<real_t> >& output, size_t num_output)
 *  \ingroup Util
 *  \brief Function to perform a B-spline interpolation.
 *  \details This program is adapted from a free program ditributed by
 *  Keith Vertanen (vertankd@cda.mrs.umn.edu) in 1994.
 *  @param [in] n Number of control points minus 1.
 *  @param [in] t Degree of the polynomial plus 1.
 *  @param [in] control Control point array made up of Point stucture.
 *  @param [out] output Vector in which the calculated spline points are to be put.
 *  @param [in] num_output How many points on the spline are to be calculated.
 *
 *  \note Condition: <tt>n+2>t</tt> (No curve results if <tt>n+2<=t</tt>)
 *  Control vector contains the number of points specified by <tt>n</tt>
 *  Output array is the proper size to hold <tt>num_output</tt> point structures
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void BSpline(size_t                n,
             size_t                t,
             Vect<Point<real_t> >& control, 
             Vect<Point<real_t> >& output,
             size_t                num_output);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
