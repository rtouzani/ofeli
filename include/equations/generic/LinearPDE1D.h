/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2026 Rachid Touzani

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

    Class LinearPDE1D : A generic class projected to contain all PDE's in 1D
                 solved using 2-Node Line (P1) Finite element

  ==============================================================================*/

#ifndef __LINEAR_PDE_1D_H
#define __LINEAR_PDE_1D_H

#include "equations/generic/Equa_LinearPDE.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file LinearPDE1D.h
 *  \brief Definition file for class LinearPDE1D.
 */

/*! \class LinearPDE1D
 *  \ingroup GenericPDE
 *  \brief Solves a generic linear PDE in 1-D using 2-Node line finite elements.
 *
 *  \details This class handles the scalar linear equation:
 *       a d^2u/dt^2 + b du/dt - d/dx(c du/dx) + d du/dx + eu = f
 *   in a given interval, with Dirichlet and/or Neumann boundary conditions.
 *   The coefficients a, b, c, d, e can be constants or given functions of space and time
 *   
 *   Numerical solution uses the P1 finite element method for space discretization.
 *   For transient (time-dependent) problems, the class TimeStepping must be used.
 *   For eigenproblems, the class EigenProblemSolver must be used.
 */
class LinearPDE1D : public Equa_LinearPDE<2,1>
{

 public:

    using Equa::getPDECoef;

/// \brief Default Constructor.
/// \details Constructs an empty equation.
    LinearPDE1D();

/// Constructor using mesh instance
/// @param [in] ms Mesh instance
    LinearPDE1D(Mesh& ms);

/** Constructor using mesh instance and solution vector
 *  @param [in] ms Mesh instance
 *  @param [in,out] u Vect instance containing solution vector
 */ 
    LinearPDE1D(Mesh&         ms,
                Vect<real_t>& u);

/// \brief Destructor
    ~LinearPDE1D() { }

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

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Add boundary right-hand side flux to right hand side.
/// @param [in] flux Vector containing source at side nodes.
    void BoundaryRHS(real_t flux) { }

/** \brief Add boundary right-hand side term to right hand side after multiplying it by 
 *  coefficient <tt>coef</tt>
 *  @param [in] f Vector containing source at nodes.
 */
    void BoundaryRHS(const Vect<real_t>& f) { }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Return (constant) heat flux in element.
    real_t Flux() const;

/** \brief Set equation input data
 *  @param [in] opt Parameter that selects data type for input. This parameter
 *  is to be chosen in the enumerated variable EqDataType
 *  <ul>
 *     <li><tt>INITIAL_FIELD</tt>: Initial temperature
 *     <li><tt>BOUNDARY_CONDITION_DATA</tt>: Boundary condition (Dirichlet)
 *     <li><tt>SOURCE_DATA</tt>: Heat source
 *     <li><tt>FLUX_DATA</tt>: Heat flux (Neumann boundary condition)
 *  </ul>
 *  @param [in] u Vector containing input data
 */
//    void setInput(EqDataType    opt,
//                  Vect<real_t>& u);

 private:
    void set(const Element* el);
    void set(const Side* sd) { }
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
