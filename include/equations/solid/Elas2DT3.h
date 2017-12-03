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

class Elas2DT3 : public Equa_Solid<real_t,3,6,2,4>
{

 public :

/// \brief Default Constructor.
/// \details Constructs an empty equation.
    Elas2DT3();

/// \brief Constructor using Mesh data
/// @param [in] ms Mesh instance
    Elas2DT3(Mesh& ms);

/// \brief Constructor using element data
/// @param el Pointer to Element instance
    Elas2DT3(const Element* el);

/// \brief Constructor using side data
/// @param sd Pointer to Side instance
    Elas2DT3(const Side* sd);

/** \brief Constructor using element, previous time solution <tt>u</tt> and time value
 *  @param [in] el Pointer to element.
 *  @param [in] u Vect instance that contains solution at previous time step.
 *  @param [in] time Current time value [Default: <tt>0</tt>].
 */
    Elas2DT3(const Element*      el,
             const Vect<real_t>& u,
                   real_t        time=0.);

/** \brief Constructor for an element (transient case) with specification of time integration scheme.
 *  @param [in] el Pointer to element.
 *  @param [in] u Vect instance that contains solution at previous time step.
 *  @param [in] time Current time value.
 *  @param [in] deltat Time step.
 *  @param [in] scheme Time Integration Scheme: To be chosen among the enumerated values:
 *  <ul>
 *     <li><tt>FORWARD_EULER</tt>: Forward Euler scheme
 *     <li><tt>BACKWARD_EULER</tt>: Backward Euler scheme,
 *     <li><tt>CRANK_NICOLSON</tt>: Crank-Nicolson Euler scheme.
 *  </ul>
 */
    Elas2DT3(const Element*      el,
             const Vect<real_t>& u,
                   real_t        time,
                   real_t        deltat,
                   int           scheme);

/** \brief Constructor using side, previous time solution <tt>u</tt> and time value.
 *  @param [in] sd Pointer to side.
 *  @param [in] u Vect instance that contains solution at previous time step.
 *  @param [in] time Current time value [Default: <tt>0</tt>].
 */
    Elas2DT3(const Side*         sd,
             const Vect<real_t>& u,
                   real_t        time=0.);

/** \brief Constructor for a side (transient case) with specification of time integration scheme.
 *  @param [in] sd Pointer to side.
 *  @param [in] u Vect instance that contains solution at previous time step.
 *  @param [in] time Current time value [Default: <tt>0</tt>].
 *  @param [in] deltat Time step.
 *  @param [in] scheme Time Integration Scheme
 */
    Elas2DT3(const Side*         sd,
             const Vect<real_t>& u,
                   real_t        time,
                   real_t        deltat,
                   int           scheme);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Constructors for a mesh, matrix and right-hand side
/*    Elas2DT3(Mesh &mesh, SkSMatrix<real_t> &A, Vect<real_t> &b) : Equation<real_t,3,6,2,4>(mesh, A, b) { }
    Elas2DT3(Mesh &mesh, SkMatrix<real_t>  &A, Vect<real_t> &b) : Equation<real_t,3,6,2,4>(mesh, A, b) { }
    Elas2DT3(Mesh &mesh, SpMatrix<real_t>  &A, Vect<real_t> &b) : Equation<real_t,3,6,2,4>(mesh, A, b) { }*/

// Constructors for a mesh, matrix, right-hand side and previous solution
/*   Elas2DT3(Mesh &mesh, SkSMatrix<real_t> &A, Vect<real_t> &b, Vect<real_t> &u) :
            Equation<real_t,3,6,2,4>(mesh, A, b, u) { }
            Elas2DT3(Mesh &mesh, SkMatrix<real_t>  &A, Vect<real_t> &b, Vect<real_t> &u) :
            Equation<real_t,3,6,2,4>(mesh, A, b, u) { }
   Elas2DT3(Mesh &mesh, SpMatrix<real_t>  &A, Vect<real_t> &b, Vect<real_t> &u) :
            Equation<real_t,3,6,2,4>(mesh, A, b, u) { }*/
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

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

/// \brief Add element lumped mass contribution to matrix after multiplication by <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void LMassToLHS(real_t coef=1.);

/// \brief Add element lumped mass contribution to right-hand side after multiplication by
/// <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void LMassToRHS(real_t coef=1.);

/// \brief Add element lumped mass contribution to matrix and right-hand side after
/// multiplication by <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void LMass(real_t coef=1.) { LMassToLHS(coef); LMassToRHS(coef); };

/// \brief Add element consistent mass contribution to matrix after multiplication by <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void MassToLHS(real_t coef=1.);

/// \brief Add element consistent mass contribution to right-hand side after multiplication
/// by <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void MassToRHS(real_t coef=1.);

/// \brief Add element consistent mass contribution to matrix and right-hand side after
/// multiplication by <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void Mass(real_t coef=1.) { MassToLHS(coef); MassToRHS(coef); };

/// \brief Add element deviatoric matrix to left-hand side after multiplication by <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void Deviator(real_t coef=1.);

/// \brief Add element deviatoric contribution to right-hand side after multiplication by 
/// <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void DeviatorToRHS(real_t coef=1.);

/// \brief Add element dilatational contribution to left-hand side after multiplication by 
/// <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void Dilatation(real_t coef=1.);

/** \brief Add element dilatational contribution to right-hand side after multiplication by 
 *  <tt>coef</tt>
 *  \details To use for explicit formulations
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
    void DilatationToRHS(real_t coef=1.);

/// \brief Add body right-hand side term to right-hand side after multiplication by <tt>coef</tt>
/// \details Body forces are deduced from UserData instance <tt>ud</tt>
    void BodyRHS(UserData<real_t>& ud);

/** \brief Add body right-hand side term to right hand side
 *  @param [in] f Vector containing source at element nodes (DOF by DOF)
 *  @param [in] opt Vector is local (<tt>LOCAL_ARRAY</tt>) with size <tt>6</tt> or global
 *  (<tt>GLOBAL_ARRAY</tt>) with size = Number of element DOF [Default: <tt>GLOBAL_ARRAY</tt>].
 */
    void BodyRHS(const Vect<real_t>& f,
                       int           opt=GLOBAL_ARRAY);

/// \brief Add boundary right-hand side term to right hand side after multiplication by 
/// <tt>coef</tt>
/// @param [in] ud UserData instance defining boundary forces
    void BoundaryRHS(UserData<real_t>& ud);

/// \brief Add boundary right-hand side term to right hand side.
/// @param [in] f Vect instance that contains constant traction to impose to side.
    void BoundaryRHS(const Vect<real_t>& f);

/** \brief Penalty Signorini contact side contribution to matrix and right-hand side.
 *  @param [in] ud UserData instance defining contact data
 *  @param [in] coef Penalty value by which the added term is multiplied [Default: <tt>1.e07</tt>]
 *  @return = <tt>0</tt> if no contact is achieved on this side, <tt>1</tt> otherwise
 */
    int SignoriniContact(UserData<real_t>& ud,
                         real_t            coef=1.e07);

/** \brief Penalty Signorini contact side contribution to matrix and right-hand side.
 *  @param [in] f Vect instance that contains contact data
 *  @param [in] coef Penalty value by which the added term is multiplied [Default: <tt>1.e07</tt>]
 *  @return = <tt>0</tt> if no contact is achieved on this side, <tt>1</tt> otherwise
 */
    int SignoriniContact(Vect<real_t>& f,
                         real_t        coef=1.e07);

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
    void ContactPressure(const Vect<real_t>&  f,
                               real_t         penal,
                               Point<real_t>& p);

/// \brief Calculate strains in element.
/// \details This function can be invoked in postprocessing.
    void Strain(Vect<real_t>& eps);

/** \brief Calculate principal stresses and Von-Mises stress in element.
 *  @param [in] s vector of principal stresses
 *  @param [in] vm Von-Mises stress.
 *  This function can be invoked in postprocessing.
 */
    void Stress(Vect<real_t>& s,
                real_t&       vm);

/** \brief Add contribution of periodic boundary condition (by a penalty technique).
 *  \details Boundary nodes where periodic boundary conditions are to be imposed must
 *  have codes equal to <tt>PERIODIC_A</tt> on one side and <tt>PERIODIC_B</tt> on the opposite side.
 *  @param [in] coef Value of penalty parameter [Default: <tt>1.e20</tt>]
 */
    void Periodic(real_t coef=1.e20);

 protected:
   void set(const Element *el);
   void set(const Side *sd);

 private:
   real_t _E1, _E2, _E3, _E6;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
