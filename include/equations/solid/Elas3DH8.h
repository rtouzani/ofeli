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

                         Definition of class Elas3DH8
             for Three-Dimensional Linearized Elasticity equations using
                       8-node hexahedral finite element

  ==============================================================================*/


#ifndef __ELAS3DH8_H
#define __ELAS3DH8_H


#include "equations/solid/Equa_Solid.h"
#include "shape_functions/Hexa8.h"
#include "shape_functions/Quad4.h"

namespace OFELI {

/*! \file Elas3DH8.h
 *  \brief Definition file for class Elas3DH8.
 */

/*! \class Elas3DH8
 *  \ingroup Solid
 *  \brief To build element equations for 3-D linearized elasticity using 8-node hexahedra.
 *
 *  \details This class enables building finite element arrays for linearized isotropic
 *  elasticity problem in 3-D domains using 8-Node hexahedra.\n
 *
 *  Note that members calculating element arrays have as an argument a double
 *  \c coef that is multiplied by the contribution of the current element.
 *  This makes possible testing different algorithms.
 *
 */

class Elas3DH8 : virtual public Equa_Solid<double,8,24,4,12>
{

public :

/// \brief Default Constructor.
/// \details Constructs an empty equation.
    Elas3DH8()
    {
       _hexa = NULL;
       _quad = NULL;
    }

/// \brief Constructor using element data.
    Elas3DH8(const Element *el);

/// \brief Constructor using side data
    Elas3DH8(const Side *sd);

/// \brief Constructor using element, previous time solution \a u and time value.
    Elas3DH8(const Element*      element,
             const Vect<real_t>& u,
             const real_t&       time=0.);

/// \brief Constructor using side, previous time solution \a u and time value.
    Elas3DH8(const Side*         side,
             const Vect<real_t>& u,
             const real_t&       time=0.);

/// \brief Destructor
    ~Elas3DH8();

/// \brief Add element lumped mass contribution to matrix after multiplication by \a coef.
    void LMassToLHS(real_t coef=1.);

/// \brief Add element lumped mass contribution to right-hand side after multiplication by \a coef.
    void LMassToRHS(real_t coef=1.);

/// \brief Add element lumped mass contribution to right-hand and left-hand sides after multiplication by \a coef.
    void LMass(real_t coef) { LMassToLHS(coef); LMassToRHS(coef); }

/// \brief Add element lumped mass contribution to matrix and right-hand side after multiplication by \a coef.
    void Mass(real_t coef=1.) { coef=1; std::cerr << "Sorry, consistent mass matrix is not implemented !\n"; }

/// \brief Add element deviatoric matrix to left-hand side after multiplication by \a coef
    void Deviator(real_t coef=1.);

/// \brief Add element deviatoric matrix to right-hand side after multiplication by \a coef
    void DeviatorToRHS(real_t coef=1.);

/// \brief Add element dilatational contribution to left hand-side after multiplication by \a coef.
    void Dilatation(real_t coef=1.);

/// \brief Add element dilatational contribution to right hand-side after multiplication by \a coef.
    void DilatationToRHS(real_t coef=1.);

/** \brief Add body right-hand side term to right hand side.
 *  \details Body forces are deduced from UserData instance \a ud.
 */
    void BodyRHS(UserData<real_t> &ud);

/** \brief Add boundary right-hand side term to right hand side
 *  @param [in] f Vector containing traction (boundary force) at sides
 */
    void BoundaryRHS(const Vect<real_t> &f);

/** \brief Add body right-hand side term to right hand side.
    @param [in] bf Vector containing source at element nodes (DOF by DOF).
    @param [in] opt Vector is local (\a LOCAL_ARRAY) with size 24 or global
    (\a GLOBAL_ARRAY) with size = Number of element DOF.
 */
    void BodyRHS(const Vect<real_t>& bf,
                       int           opt=LOCAL_ARRAY);

private :

   Hexa8                           *_hexa;
   Quad4                           *_quad;
   LocalVect<real_t,8>             _w;
   real_t                          _xg[2], _wg[2];
   Point<real_t>                   _cg, _xl[8];
   LocalMatrix<Point<real_t> ,8,8> _dsh;
   void set(const Element *el);
   void set(const Side *sd);
};

} /* namespace OFELI */

#endif
