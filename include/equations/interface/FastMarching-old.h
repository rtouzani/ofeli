/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

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

#include "equations/interface/FMHeap.h"
#include "mesh/Grid.h"
#include "linear_algebra/Vect.h"
#include "util/util.h"

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
 *    |\nabla u|F = 1
 * with u = 0 on the interface, where \c F is the velocity 
 */

class FastMarching
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
    FastMarching(const Grid&         g,
                 Vect<real_t>&       T,
                 const Vect<real_t>& F);

/// \brief Destructor
    ~FastMarching() { }

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
    void set(const Grid&         g,
             Vect<real_t>&       T,
             const Vect<real_t>& F);


/** Execute Fast Marching Procedure
 *  Once this function is invoked, the vector \c phi in the constructor or in the member function \c set
 *  contains the solution
 */
    void run();

/// \brief Check consistency by computing the discrete residual.
/// \details This function returns residual error (||\nabla u|^2|F|-1|)
    real_t getResidual();

/**
 * @brief Extend speed vector to whole domain
 * 
 * @param[out] v Vector containing speed at each grid node
 */
   void ExtendSpeed(Vect<real_t>& v);
   
 private:

   FMHeap _Narrow;
   Pt *_p, *_np;
   vector<Pt> _u;
   vector<Pt *> _neigs;
   const Grid *_theGrid;
   Vect<real_t> *_T, _F;
   size_t _dim, _nx, _ny, _nz;
   real_t _hx, _hy, _hz;
   void run1D();
   void run2D();
   void run3D();
   void init1D();
   void init2D();
   void init3D();
   real_t eval1D();
   real_t eval2D();
   real_t eval3D();
   int Neigs1D();
   int Neigs2D();
   int Neigs3D();
   void ExtendSpeed1D(Vect<real_t>& v);
   void ExtendSpeed2D(Vect<real_t>& v);
   void ExtendSpeed3D(Vect<real_t>& v);
   void UpdateExt1D(int i, Vect<real_t>& v);
   void UpdateExt2D(int i, int j, Vect<real_t>& v);
   void UpdateExt3D(int i, int j, int k, Vect<real_t>& v);

   int IJ(int i, int j) const { return (_ny+1)*(i-1)+j-1; }
   int IJK(int i, int j, int k) const { return (_ny+1)*(_nz+1)*(i-1)+(_nz+1)*(j-1)+k-1; }

   int MaxQuad(const real_t& a, const real_t& b, const real_t& c, real_t& x)
   {
      real_t d = b*b - a*c;
      if (d<0.)
         return -1;
      d = sqrt(d);
      x = fmax((-b-d)/a,(-b+d)/a);
      return 0;
   }
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
