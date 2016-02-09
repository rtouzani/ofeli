/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzani

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
#include "io/UserData.h"

namespace OFELI {

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

class Bar2DL2 : public Equa_Solid<real_t,2,4,1,2>
{

 public :

/// \brief Default Constructor.
/// \details Constructs an empty equation.
    Bar2DL2() { }

/** \brief Constructor using element data.
 *  @param [in] el Pointer to Element
 *  @param [in] section Section of bar at present element
 */
    Bar2DL2(Element* el,
            real_t   section);

/// \brief Destructor
    ~Bar2DL2() { }

/// \brief Add element consistent mass contribution to matrix and right-hand side after multiplication by <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void Mass(real_t coef=1.) { MassToLHS(coef); MassToRHS(coef); }

/// \brief Add element lumped mass contribution to matrix ans right-hand side after multiplication by <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void LMass(real_t coef=1.) { LMassToLHS(coef); LMassToRHS(coef); }

/// \brief Add lumped mass matrix to left-hand side after multiplying it by coefficient <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void LMassToLHS(real_t coef=1);

/// \brief Add lumped mass contribution to right-hand side after multiplying it by coefficient <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void LMassToRHS(real_t coef=1); 

/// \brief Add consistent mass matrix to left-hand side after multiplying it by coefficient <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void MassToLHS(real_t coef=1);

/// \brief Add consistent mass contribution to right-hand side after multiplying it by coefficient <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void MassToRHS(real_t coef=1); 

/// \brief Add element stiffness to left hand side.
/// @param [in] coef Coefficient to multuply by added term [Default: <tt>1</tt>].
    void Stiffness(real_t coef=1.);

/// \brief Add body right-hand side term to right hand side.
/// @param [in] ud instance containing user data with prescribes loads
    void BodyRHS(UserData<real_t>& ud);

/// \brief Return stresses in bar.
    real_t Stress() const;

/** \brief Return stresses in the truss structure (elementwise)
 *  @param [in] u Vect instance containing displacements at nodes
 *  @param [in] s Vect instance containing axial stresses in elements
 */
    void getStresses(const Vect<real_t>& u,
                           Vect<real_t>& s);

/** \brief Run one time step
 *  \details This function performs one time step in equation solving.
 *  It is to be used only if a <tt>TRANSIENT</tt> analysis is required.
 *  @return Return error from the linear system solver
 */
    int runOneTimeStep();

/** \brief Solve the equation
 *  \details If the analysis (see function \b setAnalysis) is <tt>STEADY_STATE</tt>, then
 *  the function solves the stationary equation.\n
 *  If the analysis is <tt>TRANSIENT</tt>, then the function performs time stepping
 *  until the final time is reached.
 */
    int run();

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

/** \brief Build global stiffness and mass matrices for the eigen system
 *  \details Case where the mass matrix is consistent
 *  @param [in] K Stiffness matrix
 *  @param [in] M Consistent mass matrix
 */
    void buildEigen(SkSMatrix<real_t>& K,
                    SkSMatrix<real_t>& M);

/** \brief Build global stiffness and mass matrices for the eigen system
 *  \details Case where the mass matrix is lumped
 *  @param [in] K Stiffness matrix
 *  @param [in] M Vector containing diagonal mass matrix
 */
    void buildEigen(SkSMatrix<real_t>& K,
                    Vect<real_t>&      M);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void BodyRHS(const Vect<real_t> &f, int opt=GLOBAL_ARRAY) { }
    void BoundaryRHS(const Vect<real_t> &f) { }
    void BoundaryRHS(UserData<real_t> &ud) { }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 protected:
   void set(const Element *el);

 private:
   real_t  _section, _cc, _ss, _sc;
   SkSMatrix<real_t> _A;
};

} /* namespace OFELI */

#endif
