/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2019 Rachid Touzani

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

                         Definition of class Laplace2DT3
              for 2-D Laplace equation using 3-node triangular element

  ==============================================================================*/


#ifndef __LAPLACE_2DT3_H
#define __LAPLACE_2DT3_H

#include "equations/laplace/Equa_Laplace.h"
#include "linear_algebra/Assembly.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Laplace2DT3.h
 *  \brief Definition file for class Laplace2DT3.
 */

/*! \class Laplace2DT3
 *  \ingroup Laplace
 *  \brief To build element equation for the Laplace equation
 *  using the 2-D triangle element (<tt>P<sub>1</sub></tt>).
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class Laplace2DT3 : virtual public Equa_Laplace<real_t,3,3,2,2> {

 public:

   using Equation<real_t,3,3,2,2>::_dSh;

/// \brief Default constructor.
    Laplace2DT3();

/** \brief Constructor with mesh.
 *  @param [in] ms Mesh instance
 */
    Laplace2DT3(Mesh& ms);

/** \brief Constructor using mesh and solution vector
 *  @param [in] ms Mesh instance
 *  @param [in] u Problem right-hand side
 */
    Laplace2DT3(Mesh&         ms,
                Vect<real_t>& u);

/** \brief Constructor that initializes a standard Poisson equation
 *  \details This constructor sets data for the Poisson equation with mixed (Dirichlet and Neumann) boundary conditions.
 *  @param [in] ms Mesh instance
 *  @param [in] b Vector containing the source term (right-hand side of the equation) at
 *  mesh nodes
 *  @param [in] Dbc Vector containing prescribed values of the solution (Dirichlet
 *  boundary condition) at nodes with positive code. Its size is the total number of nodes
 *  @param [in] Nbc Vector containing prescribed fluxes (Neumann boundary conditions)
 *  at sides, its size is the total number of sides
 *  @param [in] u Vector to contain the finite element solution at nodes once the member function run() is called.
 */
    Laplace2DT3(Mesh&         ms,
                Vect<real_t>& b,
                Vect<real_t>& Dbc,
                Vect<real_t>& Nbc,
                Vect<real_t>& u);

/// \brief Destructor
    ~Laplace2DT3() { }

/// \brief Add finite element matrix to left-hand side
    void LHS();

/// \brief Add body source term to right-hand side
/// @param [in] f Vector containing the source given function at mesh nodes
    void BodyRHS(const Vect<real_t>& f);

/// \brief Add boundary source term to right-hand side
/// @param [in] h Vector containing the source given function at mesh nodes
    void BoundaryRHS(const Vect<real_t>& h);

/** \brief Build global stiffness and mass matrices for the eigen system
 *  @param [in] opt Flag to choose a lumped mass matrix (0) or consistent (1) [Default: <tt>0</tt>]
 */
    void buildEigen(int opt=0);

/** \brief Perform post calculations
 *  @param [in] u Solution at nodes
 *  @param [out] p Vector containing gradient at elements
 */
    void Post(const Vect<real_t>&   u,
              Vect<Point<real_t> >& p);

 private:

   void set(const Element* el);
   void set(const Side* sd);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
