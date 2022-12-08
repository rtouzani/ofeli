/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani

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
                           for 3-D uniform grids

  ==============================================================================*/

#ifndef __FAST_MARCHING_3DG_H
#define __FAST_MARCHING_3DG_H

#include "equations/interface/FMHeap.h"
#include "mesh/Grid.h"
#include "linear_algebra/Vect.h"
#include "util/util.h"
#include "equations/Equa.h"

/*!
 * \file FastMarching3DG.h
 * \brief Fast Marching Method for 3-D uniform grids
 */

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! 
 * \class FastMarching3DG
 * \brief class for the fast marching algorithm on 3-D uniform grids
 *
 * \details This class implements the Fast Marching method to solve
 * the eikonal equation in a 3-D uniform grid.
 * In other words, the class solves the partial differential equation
 *    |&nabla;T|F = 1
 * with T = 0 on the interface, where \c F is the velocity 
 */

class FastMarching3DG : virtual public Equa
{

 public:

   using Equa::_u;
   using Equa::_v;

/*!
 * \brief Default Constructor
 * \details Initializes to default value grid data
 */
    FastMarching3DG();

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
    FastMarching3DG(const Grid&   g,
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
    FastMarching3DG(const Grid&   g,
                    Vect<real_t>& T,
                    Vect<real_t>& F);

/// \brief Destructor
    ~FastMarching3DG();

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
 *  contains the solution
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
/// \details This function returns residual error (||&nabla;u|<sup>2</sup>|F|<sup>2</sup>-1|)
    real_t getResidual();

 private:

   FMHeap _Narrow;
   Pt *_p, *_np;
   Vect<real_t> _b;
   vector<Pt> _U;
   vector<Pt *> _neigs;
   const Grid *_theGrid;
   size_t _nx, _ny, _nz;
   real_t _hx, _hy, _hz;
   void init();
   real_t eval();
   size_t Neigs();
   size_t IJK(int i, int j, int k) const { return (_ny+1)*(_nz+1)*(i-1)+(_nz+1)*(j-1)+k-1; }
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
