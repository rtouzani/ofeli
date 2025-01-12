/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2025 Rachid Touzani

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

                         Definition of class Bar2DL2
                for Planar Bar element with 2 d.o.f. per node

  ==============================================================================*/


#ifndef __BAR2DL2_H
#define __BAR2DL2_H

#include "equations/solid/Equa_Solid.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Bar2DL2.h
 *  \brief Definition file for class Bar2DL2.
 */

/*! \class Bar2DL2
 *  \ingroup Solid
 *  \brief To build element equations for Planar Elastic Bar element with 2 DOF (Degrees of Freedom) per node.
 * 
 *  \details This class implements a planar (two-dimensional) elastic bar using
 *  2-node lines.
 *  Note that members calculating element arrays have as an argument a real
 *  <tt>coef</tt> that is multiplied by the contribution of the current element.
 *  This makes possible testing different algorithms.
 */

class Bar2DL2 : public Equa_Solid<2,4,1,2>
{

 public :

/// \brief Default Constructor.
/// \details Constructs an empty equation.
    Bar2DL2() { }

/** \brief Constructor using a Mesh instance.
 *  @param [in] ms Reference Mesh instance
 */
    Bar2DL2(Mesh& ms);

/** \brief Constructor using a Mesh instance and a solution vector instance.
 *  @param [in] ms Reference Mesh instance
 *  @param [in,out] u Reference to solution vector
 */
    Bar2DL2(Mesh&         ms,
            Vect<real_t>& u);

/// \brief Destructor
    ~Bar2DL2() { }

/// \brief Define bar section
    void setSection(real_t A);

/// \brief Add lumped mass matrix to element matrix after multiplying it by coefficient <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void LMass(real_t coef=1);

/// \brief Add consistent mass matrix to element matrix after multiplying it by coefficient <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void Mass(real_t coef=1);

/// \brief Add element stiffness to left hand side.
/// @param [in] coef Coefficient to multuply by added term [Default: <tt>1</tt>].
    void Stiffness(real_t coef=1.);

/// \brief Return stresses in bar.
    real_t Stress() const;

/** \brief Return stresses in the truss structure (elementwise)
 *  @param [in] s Vect instance containing axial stresses in elements
 */
    void getStresses(Vect<real_t>& s);

/** \brief Build the linear system of equations
 *  \details Before using this function, one must have properly selected 
 *  appropriate options for:
 *  <ul>
 *     <li>The choice of a steady state or transient analysis
 *     <li>In the case of transient analysis, the choice of a time integration scheme
 *         and a lumped or consistent mass matrix
 *     <li>The choice of desired linear system solver
 *  </ul>
 */
    void build();

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void BodyRHS(const Vect<real_t> &f);
    void BoundaryRHS(const Vect<real_t> &f) { }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 private:
   real_t _section, _cc, _ss, _sc;
   void Load();
   void set(const Element *el);
   void set(const Side *sd) { }
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
