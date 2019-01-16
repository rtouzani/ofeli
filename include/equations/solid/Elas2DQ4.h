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

                         Definition of class Elas2DQ4
             for Two-Dimensional Linearized Elasticity equations using
                       4-node quadrilateral finite element

  ==============================================================================*/


#ifndef __ELAS2DQ4_H
#define __ELAS2DQ4_H


#include "equations/solid/Equa_Solid.h"
#include "shape_functions/Quad4.h"
#include "shape_functions/Line2.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Elas2DQ4.h
 *  \brief Definition file for class Elas2DQ4.
 */

/*! \class Elas2DQ4
 *  \ingroup Solid
 *  \brief To build element equations for 2-D linearized elasticity using 4-node quadrilaterals.
 *
 *  \details This class enables building finite element arrays for linearized isotropic
 *  elasticity problem in 2-D domains using 4-Node quadrilaterals.\n
 *  Unilateral contact is handled using a penalty function.
 *  Note that members calculating element arrays have as an argument a real
 *  <tt>coef</tt> that is multiplied by the contribution of the current element.
 *  This makes possible testing different algorithms.
 */

class Elas2DQ4 : virtual public Equa_Solid<real_t,4,8,2,4>
{

 public:

/// \brief Default Constructor
/// \details Constructs an empty equation.
    Elas2DQ4()
    {
       _quad = NULL;
       _ln = NULL;
    }

/// \brief Constructor using element data
    Elas2DQ4(const Element* el);

/// \brief Constructor using side data
    Elas2DQ4(const Side* sd);

/** \brief Constructor using element and previous time data
 *  @param [in] element Pointer to element
 *  @param [in] u Vect instance containing solution at previous time step
 *  @param [in] time Current time value [Default: <tt>0</tt>]
 */
    Elas2DQ4(const Element*      element,
             const Vect<real_t>& u,
             const real_t&       time=0.);

/** \brief Constructor using side and previous time data
 *  @param [in] side Pointer to side
 *  @param [in] u Vect instance containing solution at previous time step
 *  @param [in] time Current time value [Default: <tt>0</tt>]
 */
    Elas2DQ4(const Side*         side,
             const Vect<real_t>& u,
             const real_t&       time=0.);

/// \brief Destructor
    ~Elas2DQ4();

/// \brief Set plane strain hypothesis
    void PlaneStrain();

/** \brief Set plane strain hypothesis by giving values of Young's modulus and Poisson ratio
 *  @param [in] E Young's modulus
 *  @param [in] nu Poisson ratio
 */
    void PlaneStrain(real_t E,
                     real_t nu);

/// \brief Set plane stress hypothesis
    void PlaneStress();

/** \brief Set plane stress hypothesis by giving values of Young's modulus and Poisson ratio
 *  @param [in] E Young's modulus
 *  @param [in] nu Poisson ratio
 */
    void PlaneStress(real_t E,
                     real_t nu);

/// \brief Add element lumped mass contribution to matrix after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    void LMassToLHS(real_t coef=1.);

/// \brief Add element lumped mass contribution to right-hand side after multiplication by 
/// <tt>coef</tt> [Default: <tt>1</tt>]
    void LMassToRHS(real_t coef=1.);

/// \brief Add element lumped mass contribution to matrix and right-hand side after
/// multiplication by <tt>coef</tt> [Default: <tt>1</tt>]
    void LMass(real_t coef=1.) { LMassToLHS(coef); LMassToRHS(coef); }

/// \brief Add element consistent mass contribution to matrix and right-hand side after 
/// multiplication by <tt>coef</tt> [Default: <tt>1</tt>]
    void Mass(real_t coef=1.) { coef=1; std::cerr << "Sorry, consistent mass matrix is not implemented !\n"; }

/// \brief Add element deviatoric matrix to left-hand side after multiplication by 
/// <tt>coef</tt> [Default: <tt>1</tt>]
    void Deviator(real_t coef=1.);

/// \brief Add element deviatoric contribution to right-hand side after multiplication by 
/// <tt>coef</tt> [Default: <tt>1</tt>]
//  \details To use for explicit formulations
    void DeviatorToRHS(real_t coef=1.);

/// \brief Add element dilatational contribution to left-hand side after multiplication by 
/// <tt>coef</tt> [Default: <tt>1</tt>]
    void Dilatation(real_t coef=1.);

/// \brief Add element dilatational contribution to right hand side after multiplication by
/// <tt>coef</tt> [Default: <tt>1</tt>]
/// \details To use for explicit formulations
    void DilatationToRHS(real_t coef=1.);

/// \brief Add body right-hand side term to right hand side after multiplication by <tt>coef</tt>
/// \details Body forces are deduced from UserData instance <tt>ud</tt>
    void BodyRHS(UserData<real_t>& ud);

/** \brief Add body right-hand side term to right hand side.
 *  @param [in] bf Vector containing source at element nodes (DOF by DOF).
 *  @param [in] opt Vector is local (<tt>LOCAL_ARRAY</tt>) with size 8 or global
 *  (<tt>GLOBAL_ARRAY</tt>) with size = Total number of DOF [Default: <tt>GLOBAL_ARRAY</tt>].
 */
    void BodyRHS(const Vect<real_t>& bf,
                       int           opt=GLOBAL_ARRAY);

/// \brief Add boundary right-hand side term to right hand side after multiplication by
/// <tt>coef</tt>
/// \details Boundary forces are deduced from UserData instance <tt>ud</tt>
    void BoundaryRHS(UserData<real_t>& ud);

/** \brief Add boundary right-hand side term to right hand side.
 *  @param [in] sf Vector containing source at element nodes (DOF by DOF).
 *  @warning The vector <tt>sf</tt> is sidewise constant, <i>i.e.</i> its size
 *  is twice the number of sides.
 */
    void BoundaryRHS(const Vect<real_t>& sf);

/** \brief Penalty Signorini contact side contribution to matrix and right-hand side.
 *  @param [in] ud UserData instance defining contact data
 *  @param [in] coef Penalty value by which the added term is multiplied [Default: <tt>1.e07</tt>]
 *  @return <tt>0</tt> if no contact was realized on this side, <tt>1</tt> otherwise
 */
    int SignoriniContact(UserData<real_t>& ud,
                         real_t            coef=1.e07);

/// \brief Calculate strains at element barycenter
/// @param [out] eps Vector containing strains in element
    void Strain(LocalVect<real_t,3>& eps);

/** \brief Calculate principal stresses and Von-Mises stress at element barycenter.
 *  @param [out] s LocalVect containing principal stresses in element
 *  @param [out] vm Value of Von-Mises stress in element
 */
    void Stress(LocalVect<real_t,3>& s,
                real_t&              vm);

/** \brief Calculate principal stresses and Von-Mises stress at element barycenter.
 *  @param [out] sigma Vector containing principal stresses in element
 *  @param [out] s Vector containing principal stresses in element
 *  @param [out] vm Value of Von-Mises stress in element
 */
    void Stress(LocalVect<real_t,3>& sigma,
                LocalVect<real_t,3>& s,
                real_t&              vm);

 private:

   real_t        _E1, _E2, _E3, _E6;
   Quad4         *_quad;
   Line2         *_ln;
   real_t        _g[2], _w[2], _xl[4], _yl[4];
   Point<real_t> _cg;
   void set(const Element *el);
   void set(const Side *sd);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
