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

                         Definition of class Elas3DT4
             for Three-Dimensional Linearized Elasticity equations using
                       4-node tetrahedral finite element

  ==============================================================================*/


#ifndef __ELAS3DT4_H
#define __ELAS3DT4_H


#include "equations/solid/Equa_Solid.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Elas3DT4.h
 *  \brief Definition file for class Elas3DT4.
 */

/*! \class Elas3DT4
 *  \ingroup Solid
 *  \brief To build element equations for 3-D linearized elasticity using 4-node tetrahedra.
 *
 * \details This class enables building finite element arrays for linearized isotropic
 *  elasticity problem in 3-D domains using 4-Node tetrahedra.\n
 *
 */

class Elas3DT4 : virtual public Equa_Solid<real_t,4,12,3,9>
{

 public :

/// \brief Default Constructor
    Elas3DT4() { }

/// \brief Constructor using element data
    Elas3DT4(const Element *el);

/// \brief Constructor using side data
    Elas3DT4(const Side *sd);

/// \brief Constructor using element and previous time data
    Elas3DT4(const Element*      element,
             const Vect<real_t>& u,
	     const real_t&       time=0.);

/// \brief Constructor using side and previous time data
    Elas3DT4(const Side*         side,
	     const Vect<real_t>& u,
	     const real_t&       time=0.);

/// \brief Destructor
    ~Elas3DT4();

/** \brief Set Media properties.
 *  @param [in] E Young's modulus
 *  @param [in] nu Poisson ratio
 *  @param [in] rho Density
 */
    void Media(real_t E,
	       real_t nu,
	       real_t rho);

/// \brief Add element lumped mass contribution to matrix after multiplication by \a coef
    void LMassToLHS(real_t coef=1);

/// \brief Add element lumped mass contribution to right-hand side after multiplication by \a coef
    void LMassToRHS(real_t coef=1);

/// \brief Add element lumped mass contribution to matrix and right-hand side after multiplication by \a coef
    void LMass(real_t coef) { LMassToLHS(coef); LMassToRHS(coef); }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Add element consistent Mass contribution to matrix after multiplication by coef
    void Mass(real_t coef=1.) { coef=1; std::cerr << "Sorry, consistent mass matrix is not implemented !\n"; }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Add element deviatoric matrix to left hand-side after multiplication by \a coef.
    void Deviator(real_t coef=1.);

/// \brief Add element deviatoric matrix to right hand-side after multiplication by \a coef.
    void DeviatorToRHS(real_t coef=1.);

/// \brief Add element dilatational contribution to left-hand side after multiplication by \a coef.
    void Dilatation(real_t coef=1.);

/// \brief Add element dilatational contribution to right-hand side after multiplication by \a coef.
    void DilatationToRHS(real_t coef=1.);

/// \brief Add body right-hand side term to right hand side after multiplication by \a coef.
/// \details Body forces are deduced from UserData instance \a ud.
    void BodyRHS(UserData<real_t> &ud);

/** \brief Add body right-hand side term to right hand side.
 *  @param [in] f Vect instance containing source at element nodes (DOF by DOF).
 *  @param [in] opt Vector is local (LOCAL_ARRAY) with size 12 or global
 * (GLOBAL_ARRAY) with size = Number of element DOF.
 */
    void BodyRHS(const Vect<real_t>& f,
		       int           opt=LOCAL_ARRAY);

/// \brief Add boundary right-hand side term to right hand side.
/// @param [in] f Vect instance that contains constant traction to impose to side.
    void BoundaryRHS(const Vect<real_t>& f);

/** \brief Build global stiffness and mass matrices for the eigen system
 *  \details Case where the mass matrix is lumped
 *  @param [in] K Stiffness matrix
 *  @param [in] M Vector containing diagonal mass matrix
 */
    void buildEigen(SkSMatrix<real_t>& K,
		    Vect<real_t>&      M);

 private:
   void set(const Element *el);
   void set(const Side *sd);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
