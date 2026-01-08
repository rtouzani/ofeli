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

                         Definition of class Elas2DT3
             for Two-Dimensional Linearized Elasticity equations using
                           3-node triangular finite element

  ==============================================================================*/


#ifndef __ELAS2DT3_H
#define __ELAS2DT3_H

#include "equations/solid/Equa_Solid.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Elas2DT3.h
 *  \brief Definition file for class Elas2DT3.
 */

/*! \class Elas2DT3
 *  \ingroup Solid
 *  \brief To build element equations for 2-D linearized elasticity using 3-node triangles.
 *
 *  \details This class enables building finite element arrays for linearized isotropic
 *  elasticity problem in 2-D domains using 3-Node triangles.\n
 *  Unilateral contact is handled using a penalty function.
 *  Note that members calculating element arrays have as an argument a real
 *  <tt>coef</tt> that is multiplied by the contribution of the current element.
 *  This makes possible testing different algorithms.
 *
 */

class Elas2DT3 : public Equa_Solid<3,6,2,4>
{

 public:

    using Equa_Solid<3,6,2,4>::run;

/// \brief Default Constructor.
/// \details Constructs an empty equation.
    Elas2DT3();

/// \brief Constructor using Mesh data
/// @param [in] ms Mesh instance
    Elas2DT3(Mesh& ms);

/** \brief Constructor using Mesh data and solution vector
 *  @param [in] ms Mesh instance
 *  @param [in,out] u Reference to solution vector
 */ 
    Elas2DT3(Mesh&         ms,
             Vect<real_t>& u);

/// Destructor
    ~Elas2DT3();

/// \brief Set media properties.
/// \details Useful to override material properties deduced from mesh file.
    void Media(real_t E,
               real_t nu,
               real_t rho);

/// \brief Set plane strain hypothesis
    void PlaneStrain();

/// \brief Set plane strain hypothesis by giving values of Young's modulus <tt>E</tt>
/// and Poisson ratio <tt>nu</tt>
    void PlaneStrain(real_t E,
                     real_t nu);

/// \brief Set plane stress hypothesis.
    void PlaneStress();

/// \brief Set plane stress hypothesis by giving values of Young's modulus <tt>E</tt>
/// and Poisson ratio <tt>nu</tt>
    void PlaneStress(real_t E,
                     real_t nu);

/// \brief Add element lumped mass contribution to element matrix after multiplication by <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void LMass(real_t coef=1.);

/// \brief Add element consistent mass contribution to element matrix after multiplication by <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void Mass(real_t coef=1.);

/// \brief Add element deviatoric matrix to element matrix after multiplication by <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void Deviator(real_t coef=1.);

/// \brief Add element dilatational contribution to element matrix after multiplication by 
/// <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void Dilatation(real_t coef=1.);

/** \brief Add body right-hand side term to right hand side
 *  @param [in] f Vector containing source at nodes (DOF by DOF)
 */
    void BodyRHS(const Vect<real_t>& f);

/// \brief Add body right-hand side term to right hand side
    void BodyRHS();

/// \brief Add boundary right-hand side term to right hand side.
/// @param [in] f Vect instance that contains constant traction to impose to side.
    void BoundaryRHS(const Vect<real_t>& f);

/// \brief Add boundary right-hand side term to right hand side.
    void BoundaryRHS();

/** \brief Penalty Signorini contact side contribution to matrix and right-hand side.
 *  @param [in] coef Penalty value by which the added term is multiplied
 *              [Default: <tt>1.e07</tt>]
 *  @return = <tt>0</tt> if no contact is achieved on this side, <tt>1</tt> otherwise
 */
    int Contact(real_t coef=1.e07);

/** \brief Calculate reactions
 *  \details This function can be invoked in postprocessing
 *  @param [in] r Reaction on the side
 */
    void Reaction(Vect<real_t>& r);

/** \brief Calculate contact pressure
 *  \details This function can be invoked in postprocessing
 *  @param [in] f 
 *  @param [in] penal Penalty parameter that was used to impose contact condition
 *  @param [out] p Contact pressure
 */
    void ContactPressure(const Vect<real_t>& f,
                         real_t              penal,
                         Point<real_t>&      p);

/** \brief Calculate strains in element.
 *  \details This function can be invoked in postprocessing.
 *  @param [out] eps vector of strains in elements
 */
    void Strain(Vect<real_t>& eps);

/** \brief Calculate principal stresses and Von-Mises stress in element.
 *  @param [out] s vector of principal stresses in elements
 *  @param [out] vm Von-Mises stresses in elements
 *  This function can be invoked in postprocessing.
 */
    void Stress(Vect<real_t>& s,
                Vect<real_t>& vm);

/** \brief Add contribution of periodic boundary condition (by a penalty technique).
 *  \details Boundary nodes where periodic boundary conditions are to be imposed must
 *  have codes equal to <tt>PERIODIC_A</tt> on one side and <tt>PERIODIC_B</tt> on the opposite side.
 *  @param [in] coef Value of penalty parameter [Default: <tt>1.e20</tt>]
 */
    void Periodic(real_t coef=1.e20);

 private:
   void set(const Element *el);
   void set(const Side *sd);
   real_t _E1, _E2, _E3, _E6;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
