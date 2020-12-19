/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

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
    Elas2DQ4() : _quad(nullptr), _ln(nullptr)
    { }

/// \brief Constructor using Mesh instance
/// @param [in] ms Reference to Mesh instance
    Elas2DQ4(Mesh& ms);

/** \brief Constructor using Mesh instance and solution vector
 *  @param [in] ms Reference to Mesh instance
 *  @param [in,out] u Solution vector
 */
    Elas2DQ4(Mesh&         ms,
             Vect<real_t>& u);

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

/// \brief Add element lumped mass contribution to element matrix after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    void LMass(real_t coef=1.);

/// \brief Add element consistent mass contribution to matrix and right-hand side after 
/// multiplication by <tt>coef</tt> [Default: <tt>1</tt>]
    void Mass(real_t coef=1.) { coef=1; std::cerr << "Sorry, consistent mass matrix is not implemented !\n"; }

/// \brief Add element deviatoric matrix to element matrix after multiplication by 
/// <tt>coef</tt> [Default: <tt>1</tt>]
    void Deviator(real_t coef=1.);

/// \brief Add element dilatational contribution to element matrix after multiplication by 
/// <tt>coef</tt> [Default: <tt>1</tt>]
    void Dilatation(real_t coef=1.);

/** \brief Add body right-hand side term to right hand side.
 *  @param [in] f Vector containing source at nodes (DOF by DOF).
 */
    void BodyRHS(const Vect<real_t>& f);

/// \brief Add body right-hand side term to right hand side.
    void BodyRHS();

/** \brief Add boundary right-hand side term to right hand side.
 *  @param [in] f Vector containing source at nodes (DOF by DOF).
 */
    void BoundaryRHS(const Vect<real_t>& f);

/// \brief Add boundary right-hand side term to right hand side.
    void BoundaryRHS();

/** \brief Calculate strains at element barycenters
 *  @param [out] eps Vector containing strains in elements
 *  @remark The instance of Elas2DQ4 must have been constructed using the constructor
 *  with Mesh instance and solution vector
 */
    void Strain(Vect<real_t>& eps);

/** \brief Calculate principal stresses and Von-Mises stress at element barycenter.
 *  @param [out] st Vector containing principal stresses in elements
 *  @param [out] vm Vector containing Von-Mises stresses in elements
 *  @remark The instance of Elas2DQ4 must have been constructed using the constructor
 *  with Mesh instance and solution vector
 */
    void Stress(Vect<real_t>& st,
                Vect<real_t>& vm);

/** \brief Calculate principal stresses and Von-Mises stress at element barycenter.
 *  @param [out] sigma Vector containing principal stresses in elements
 *  @param [out] s Vector containing principal stresses in elements
 *  @param [out] st Value of Von-Mises stress in elements
 *  @remark The instance of Elas2DQ4 must have been constructed using the constructor
 *  with Mesh instance and solution vector
 */
    void Stress(Vect<real_t>& sigma,
                Vect<real_t>& s,
                Vect<real_t>& st);

 private:

   real_t _E1, _E2, _E3, _E6;
   Quad4  *_quad;
   Line2  *_ln;
   real_t _xl[4], _yl[4], _g[2], _ww[2];
   void set(const Element *el);
   void set(const Side *sd);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
