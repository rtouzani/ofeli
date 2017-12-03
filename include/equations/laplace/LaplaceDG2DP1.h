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

                        Definition of class LaplaceDG2DP1
              for 2-D Laplace equation using P1 triangular DG element

  ==============================================================================*/


#ifndef __LAPLACE_DG_2DP1_H
#define __LAPLACE_DG_2DP1_H

#include "equations/DG.h"
#include "linear_algebra/LocalMatrix.h"

namespace OFELI {

/*! \file LaplaceDG2DP1.h
 *  \brief Definition file for class LaplaceDG2DP1.
 */

/*! \class LaplaceDG2DP1
 *  \ingroup Laplace
 *  \brief To build and solve the linear system for the Poisson problem
 *  using the DG P<sub>1</sub> 2-D triangle element.
 *  \details This class build the linear system of equations for a standard elliptic
 *  equation using the Discontinuous Galerkin P<sub>1</sub> finite element method.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class LaplaceDG2DP1 : public DG {

 public:


/** \brief Constructor with mesh and vector data
 *  @param [in] ms Mesh instance
 *  @param [in] f Vector containing the right-hand side of the elliptic equation at triangle vertices
 *  @param [in] Dbc Vector containing prescribed values of the solution (Dirichlet boundary condition)
 *  at nodes having a positive code 
 *  @param [in] Nbc Vector containing prescribed values of the flux (Neumann boundary condition) at
 *  each side having a positive code
 *  @param [in] u Vector where the solution is stored once the linear system is solved
 */
    LaplaceDG2DP1(Mesh&         ms,
                  Vect<real_t>& f,
                  Vect<real_t>& Dbc,
                  Vect<real_t>& Nbc,
                  Vect<real_t>& u);

/// \brief Destructor
    ~LaplaceDG2DP1();


/** \brief Set parameters for the DG method
 * @param [in] sigma Penalty parameters to enforce continuity at nodes (Must be positive) [Default: <tt>100</tt>]
 * @param [in] eps Epsilon value of the DG method to choose among the values:
 * <ul>
 *    <li> 0 Incomplete Interior Penalty Galerkin method (IIPG)
 *    <li> -1 Symmetric Interior Penalty Galerkin method (SIPG)
 *    <li> 1 Non symmetric interior penalty Galerkin method (NIPG)
 * </ul>
 * For a user not familiar with the method, please choose the value of <tt>eps=-1</tt> and <tt>sigma>100</tt>
 * which leads to a symmetric positive definite matrix [Default: <tt>-1</tt>]
 */
    void set(real_t sigma,
             real_t eps);

/** \brief Set diffusivity matrix
 *  \details This function provides the diffusivity matrix as instance of class
 *  LocalMatrix. The default diffusivity matrix is the identity matrix
 *  @param [in] K Diffusivity matrix
 */
    void set(const LocalMatrix<real_t,2,2>& K);

/** \brief Build global matrix and right-hand side.
 *  \details The problem matrix and right-hand side are the ones used in the constructor.
 *  They are updated in this member function.
 */
    void build();

/** \brief Perform post calculations
 *  \details This function gives an averaged solution given at mesh nodes (triangle vertices) by
 *  a standard L<sub>2</sub>-projection method.
 *  @param [in] u Solution at nodes
 */
    void Smooth(Vect<real_t>& u);

/** \brief Build and solve the linear system of equations using an iterative method.
 *  \details The matrix is preconditioned by the diagonal ILU method.
 *  The linear system is solved either by the Conjugate Gradient method if the matrix is symmetric
 *  positive definite (<tt>eps=-1</tt>) or the GMRES method if not. The solution is stored in the vector
 *  <tt>u</tt> given in the constructor.
 *  @return Number of performed iterations. Note that the maximal number
 *  is 1000 and the tolerance is 1.e-8
 */
    int run();

 private:

   LocalMatrix<real_t,2,2> _K;
   LocalVect<real_t,3> _F1, _F2, _z;
   LocalVect<Point<real_t>,3> _dSh;
   LocalVect<size_t,2> _is;
   LocalVect<size_t,3> _ls1, _ls2;
   real_t _sigma, _eps;
   void set(Element& el);
   void setSide(size_t k);
   int solve();
};

} /* namespace OFELI */

#endif
