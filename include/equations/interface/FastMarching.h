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

          Definition of class 'FastMarching' for Fast Marching Methods

  ==============================================================================*/

#ifndef __FAST_MARCHING_H
#define __FAST_MARCHING_H

#include "equations/Equa.h"
#include "equations/interface/FastMarching1DG.h"
#include "equations/interface/FastMarching2DG.h"
#include "equations/interface/FastMarching3DG.h"
#include "linear_algebra/Vect.h"

/*!
 * \file FastMarching.h
 * \brief Fast Marching Method
 */

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*!
 * \class FastMarching
 * \brief class for the fast marching algorithm on uniform grids
 *
 * \details This class implements the Fast Marching method to solve
 * the eikonal equation in a uniform grid (1-D, 2-D or 3-D).
 * In other words, the class solves the partial differential equation
 *    |&nabla;u|F = 1
 * with u = 0 on the interface, where \c F is the velocity 
 */

class FastMarching : virtual public Equa
{

 public:

/*!
 * \brief Default Constructor
 * \details Initializes to default value grid data
 */
    FastMarching();

/*!
 * \brief Constructor using grid data
 * \details Constructor using Grid instance
 * \param [in] g Instance of class Grid
 * \param [in] T Vector containing the on input an initialization of the distance function and once
 * the function \c run is invoked the distance at grid nodes.
 * The initialization vector must use the following rules:
 * <ul>
 *    <li>The solution must be supplied at all grid points in the vicinity of the interface(s).
 *    <li>All other grid nodes must have the value \c INFINITY wth positive value if the node is in an outer
 *        domain and negative if it is in an inner domain
 * </ul>
 */
    FastMarching(const Grid&   g,
                 Vect<real_t>& T);

/*!
 * \brief Constructor
 * \details Constructor using Grid instance and propagation speed
 * \param [in] g Instance of class Grid
 * \param [in] T Vector containing the on input an initialization of the distance function and once
 * the function \c run is invoked the distance at grid nodes.
 * The initialization vector must use the following rules:
 * <ul>
 *    <li>The solution must be supplied at all grid points in the vicinity of the interface(s).
 *    <li>All other grid nodes must have the value \c INFINITY wth positive value if the node is in an outer
 *        domain and negative if it is in an inner domain
 * </ul>
 * \param [in] F Vector containing propagation speed at grid nodes
 */
    FastMarching(const Grid&   g,
                 Vect<real_t>& T,
                 Vect<real_t>& F);

/// \brief Destructor
    ~FastMarching();

/**
 * @brief Define grid and solution vector
 * @details This function is to be used if the default constructor has been used
 * @param [in] g Instance of class Grid
 * @param [in] T Vector containing the on input an initialization of the distance function and once
 * the function \c run is invoked the distance at grid nodes.
 * The initialization vector must use the following rules:
 * <ul>
 *    <li>The solution must be supplied at all grid points in the vicinity of the interface(s).
 *    <li>All other grid nodes must have the value \c INFINITY wth positive value if the node is in an outer
 *        domain and negative if it is in an inner domain
 * </ul>
 */
    void set(const Grid&   g,
             Vect<real_t>& T);

/**
 * @brief Define grid, solution vector and prppagation speed
 * @details This function is to be used if the default constructor has been used
 * @param [in] g Instance of class Grid
 * @param [in] T Vector containing the on input an initialization of the distance function and once
 * the function \c run is invoked the distance at grid nodes.
 * The initialization vector must use the following rules:
 * <ul>
 *    <li>The solution must be supplied at all grid points in the vicinity of the interface(s).
 *    <li>All other grid nodes must have the value \c INFINITY wth positive value if the node is in an outer
 *        domain and negative if it is in an inner domain
 * </ul>
 * @param [in] F Vector containing propagation speed at grid nodes
 */
    void set(const Grid&   g,
             Vect<real_t>& T,
             Vect<real_t>& F);

/** \brief Execute Fast Marching Procedure
 *  \details Once this function is invoked, the vector \c T in the constructor or in the member function \c set
 *  contains the solution.
 *  @return Return value:
 *  <ul>
 *    <li> = 0 if solution has been normally computed
 *    <li> != 0 An error has occurred
 *  </ul>
 */
    int run();

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void build() { }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Check consistency by computing the discrete residual.
/// \details This function returns residual error (||&nabla; u|^2|F|-1|)
    real_t getResidual();
   
 private:

   int _dim;
   FastMarching1DG *_theFMM1D;
   FastMarching2DG *_theFMM2D;
   FastMarching3DG *_theFMM3D;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
