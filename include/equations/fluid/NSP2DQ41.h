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

                        Definition of class NSP2DQ41
            for 2-D Navier-Stokes equations using penalty formulation and
                        Quadrilateral finite element
                    (4-nodes velocity, 1-node pressure)

  ==============================================================================*/


#ifndef __NSP2DQ41_H
#define __NSP2DQ41_H

#include "equations/fluid/Equa_Fluid.h"
#include "shape_functions/Quad4.h"
#include "shape_functions/Line2.h"
#include "io/UserData.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file NSP2DQ41.h
 *  \brief Definition file for class NSP2DQ41.
 */

/*! \class NSP2DQ41
 *  \ingroup Fluid
 *  \brief Builds finite element arrays for incompressible Navier-Stokes equations in 2-D
 *  domains using Q<sub>1</sub>/P<sub>0</sub> element and a penaly formulation for the incompressibility condition.
 */

class NSP2DQ41 : virtual public Equa_Fluid<real_t,4,8,2,4>
{

 public:

/// \brief Default Constructor
/// \details Builds an empty equation
    NSP2DQ41()
    {
       _quad = NULL;
       _ln = NULL;
    }

/// \brief Constructor using Element data
/// @param [in] el Pointer to Element instance
    NSP2DQ41(const Element* el);

/// \brief Constructor using Side data
/// @param [in] sd Pointer to Side instance
    NSP2DQ41(const Side* sd);

/** \brief Constructor using element and previous time data
 *  @param [in] el Pointer to Element instance
 *  @param [in] u Vector that contains velocity at previous time step
 *  @param [in] time Time value [Default: <tt>0.</tt>]
 */
    NSP2DQ41(const Element*      el,
             const Vect<real_t>& u,
             const real_t&       time=0.);

/** \brief Constructor using side and previous time data
 *  @param [in] sd Pointer to Side instance
 *  @param [in] u Vector that contains velocity at previous time step
 *  @param [in] time Time value [Default: <tt>0.</tt>]
 */
    NSP2DQ41(const Side*         sd,
             const Vect<real_t>& u,
             const real_t&       time=0.);

/// \brief Destructor
    ~NSP2DQ41();

/// \brief Define constant viscosity
    void Viscosity(real_t visc) { _visc = visc; }

/// \brief Define constant density
    void Density(real_t dens) { _dens = dens; }

/// \brief Add element lumped mass contribution to matrix after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    void LMass(real_t coef=1.);

/// \brief Add element consistent mass contribution to matrix after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    void Mass(real_t coef=1.);

/// \brief Add element viscous contribution to matrix after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    void Viscous(real_t coef=1.);

/// \brief Add element viscous contribution to right-hand side after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    void RHS_Viscous(real_t coef=1.);

/// \brief Add element penalty contribution to matrix after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    void Penal(real_t coef=1.);

/// \brief Add convection contribution to left-hand side after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
/// \details First term, explicit velocity, implicit velocity derivatives
    void LHS1_Convection(real_t coef=1.);

/// \brief Add convection contribution to left-hand side after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
/// \details Second term, implicit velocity, explicit velocity derivatives
    void LHS2_Convection(real_t coef=1.);

/// \brief Add convection contribution to right-hand side after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    void RHS_Convection(real_t coef=1.);

/// \brief Add body right-hand side term to right-hand side
/// @param [in] ud UserData instance that defines data 
    void BodyRHS(UserData<real_t>& ud);

/// \brief Add boundary right-hand side term to right-hand side
/// @param [in] ud UserData instance that defines data 
    void BoundaryRHS(UserData<real_t>& ud);

/** \brief Add contribution of periodic boundary condition (by a penalty technique).
 *  \details Boundary nodes where periodic boundary conditions are to be imposed must
 *  have codes equal to <tt>PERIODIC_A</tt> on one side and <tt>PERIODIC_B</tt> on the opposite side.
 *  @param [in] coef Value of penalty parameter [Default: <tt>1.e20</tt>].
 */
    void Periodic(real_t coef=1.e20);

/// \brief Calculate element pressure by penalization after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    real_t Pressure(real_t coef=1.);

 private:

   Quad4               *_quad;
   Line2               *_ln;
   real_t              _gauss[2], _xl[4], _yl[4];
   real_t              _cgx, _cgy;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
