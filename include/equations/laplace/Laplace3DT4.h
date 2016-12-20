/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2017 Rachid Touzani

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

                         Definition of class Laplace3DT4
              for 3-D Laplace equation using 4-node tetrahedral element

  ==============================================================================*/


#ifndef __LAPLACE_3DT4_H
#define __LAPLACE_3DT4_H

#include "equations/laplace/Equa_Laplace.h"
#include "linear_algebra/Assembly.h"
#include "io/UserData.h"
#include "shape_functions/Tetra4.h"
#include "shape_functions/Triang3.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Laplace3DT4.h
 *  \brief Definition file for class Laplace3DT4.
 */

/*! \class Laplace2DT3
 *  \ingroup Laplace
 *  \brief To build element equation for the Laplace equation
 *  using the 3-D tetrahedral element (<tt>P<sub>1</sub></tt>).
 */

class Laplace3DT4 : virtual public Equa_Laplace<real_t,4,4,3,3> {

 public:

/** \brief Constructor with mesh.
 *  @param [in] ms Mesh instance
 */
    Laplace3DT4(Mesh& ms);

/** \brief Constructor with problem data
 *  @param [in] ms Mesh instance
 *  @param [in] A Problem matrix in Sparse format.
 *  This matrix must be zeroed before calling the constructor
 *  @param [in] b Problem right-hand side
 */
    Laplace3DT4(Mesh&             ms,
                SpMatrix<real_t>& A,
                Vect<real_t>&     b);

/** \brief Constructor using mesh and solution vector
 *  @param [in] ms Mesh instance
 *  @param [in] b Problem right-hand side
 */
    Laplace3DT4(Mesh&         ms,
                Vect<real_t>& b);

/// \brief Constructor for an element
    Laplace3DT4(Element* el);

/// \brief Constructor for a side
    Laplace3DT4(Side* sd);

/// \brief Destructor
    ~Laplace3DT4() { }

/// \brief Add finite element matrix to left-hand side
/// @param [in] coef Value to multiply by the added matrix
    void LHS(real_t coef=1.);

/// \brief Add body source term to right-hand side
/// @param [in] f Vector containing the source given function at mesh nodes
    void BodyRHS(const Vect<real_t>& f);

/// \brief Add boundary source term to right-hand side
/// @param [in] h Vector containing the source given function at mesh nodes
    void BoundaryRHS(const Vect<real_t>& h);

/// \brief Define Source right-hand side of the equation
/// @param f [in] Vector containing source values at nodes
    void setSource(const Vect<real_t>& f);

/** \brief Build global matrix and right-hand side.
 *  \details The problem matrix and right-hand side are the ones used in the constructor.
 *  They are updated in this member function.
 */
    void build();

/** \brief Build global stiffness and mass matrices for the eigen system
 *  @param [in] opt Flag to choose a lumed mass matrix (0) or consistent (1) [Default: <tt>0</tt>]
 */
    void buildEigen(int opt=0);

/** \brief Perform post calculations
 *  @param [in] u Solution at nodes
 *  @param [out] p Vector containing gradient at elements
 */
    void Post(const Vect<real_t>&   u,
              Vect<Point<real_t> >& p);

/** \brief Solve the linear system of equations using the Conjugate Gradient iterative method.
 *  \details The matrix is preconditioned by an ILU method.
 *  @param [in] u Vector containing the solution at all sides (Sides where boundary conditions 
 *  are prescribed are included).
 *  @return Number of performed iterations in the CG method. Note that the maximal number
 *  is <tt>1000</tt> and the tolerance is <tt>1.e-8</tt>
 */
    int solve(Vect<real_t>& u);

/** \brief Compute the product of the stiffness matrix by a given vector
 *  @param [in] x Vector by which the matrix is multiplied
 *  @param [in] b Product vector
 */
    void Axb(const Vect<real_t>& x,
             Vect<real_t>&       b);

 private:

   SpMatrix<real_t> _A;
   void set(const Element* el);
   void set(const Side* sd);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
