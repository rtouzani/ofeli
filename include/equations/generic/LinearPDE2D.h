/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2024 Rachid Touzani

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

	   Class LinearPDE2D : A generic class projected to contain all PDE's in 2D
                   solved using 3-Node Triangular Finite element

  ==============================================================================*/


#ifndef __LINEAR_PDE_2D_H
#define __LINEAR_PDE_2D_H

#include "equations/generic/Equa_LinearPDE.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file LinearPDE2D.h
 *  \brief Definition file for class LinearPDE2D.
 */

/*! \class LinearPDE2D
 *  \ingroup Generic
 *  \brief Solves a generic linear PDE in 2-D using 3-Node triangular finite elements.
 *
 *  \details This class handles the scalar linear equation:
 *       a d^2u/dt^2 + b du/dt - div(c Grad u) + d.Grad u + eu = f
 *   in a given interval, with Dirichlet and/or Neumann boundary conditions.
 *   The coefficients a, b, c, d, e can be constants or given functions of space and time
 *   
 *   Numerical solution uses the P1 finite element method for space discretization.
 *   For transient (time-dependent) problems, the class TimeStepping must be used.
 *   For eigenproblems, the class EigenProblemSolver must be used.
 */
class LinearPDE2D : public Equa_LinearPDE<3,2>
{

 public:

    using Equa::getPDECoef;

/// \brief Default Constructor.
/// Constructs an empty equation.
    LinearPDE2D();

/// \brief Constructor using Mesh data
/// @param [in] ms Mesh instance
    LinearPDE2D(Mesh& ms);

/** \brief Constructor using Mesh and initial condition
 *  @param [in] ms Mesh instance
 *  @param [in] u Vect instance containing initial solution
 */
    LinearPDE2D(Mesh&         ms,
                Vect<real_t>& u);

/// \brief Destructor
    ~LinearPDE2D() { }

/** \brief Add 0th order term, in time and space, to left-hand side
 *  <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void Mat_00(real_t coef=1.0);

/** \brief Add 1st order term in time, 0th in space to left-hand side
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void Mat_10(real_t coef=1.0);

/** \brief Add Oth order term in time, 1st in space to left-hand side
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void Mat_01(real_t coef=1.0);

/** \brief Add 2nd order term in time, 0th in space to left-hand side
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void Mat_20(real_t coef=1.0);

/** \brief Add 0th order term in time, 2nd in space to left-hand side
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void Mat_02(real_t coef=1.0);

/** \brief Add body right-hand side term to right hand side.
 *  @param [in] f Vector containing source at nodes.
 */
    void BodyRHS(const Vect<real_t>& f);

/** \brief Add body right-hand side term to right hand side.
 *  \details Case where the body right-hand side is piecewise constant.
 *  @param [in] f Value of thermal source (Constant in element).
 */
    void BodyRHS(real_t f);

/// \brief Add boundary right-hand side flux to right hand side.
/// @param [in] flux Vector containing source at side nodes.
    void BoundaryRHS(real_t flux);

/** \brief Add boundary right-hand side term to right hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  @param [in] f Vector containing source at nodes
 */
    void BoundaryRHS(const Vect<real_t>& f);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// Return element contribution to energy.
    real_t Energy(Vect<real_t>& u);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// Return energy gradient in element.
    void EnergyGrad(Vect<real_t>& u,
                    Vect<real_t>& g);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Return (constant) heat flux in element.
    Point<real_t>& Flux() const;

/// \brief Compute gradient of solution.
/// @param [in] g Elementwise vector containing gradient of solution.
    void Grad(Vect<Point<real_t> >& g);

/** \brief Return gradient of a vector in element.
 *  @param [in] u Global vector for which gradient is computed.
 *  Vector <tt>u</tt> has as size the total number of nodes
 */
    Point<real_t>& Grad(const Vect<real_t>& u) const;

/** \brief Set equation input data
 *  @param [in] opt Parameter to select type of input (enumerated values)
 *  <ul>
 *     <li><tt>INITIAL_FIELD</tt>: Initial temperature
 *     <li><tt>BOUNDARY_CONDITION_DATA</tt>: Boundary condition (Dirichlet)
 *     <li><tt>SOURCE_DATA</tt>: Heat source
 *     <li><tt>FLUX_DATA</tt>: Heat flux (Neumann boundary condition)
 *  </ul>
 *  @param [in] u Vector containing input data
 */
    void setInput(EqDataType    opt,
                  Vect<real_t>& u);

 private:
    mutable Point<real_t> _f;
    void set(const Element* el);
    void set(const Side* sd);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif