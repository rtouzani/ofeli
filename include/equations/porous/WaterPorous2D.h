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

                          Definition of class WaterPorous2D

  ==============================================================================*/


#ifndef __WATER_POROUS_2D_H
#define __WATER_POROUS_2D_H

#include "equations/porous/Equa_Porous.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file WaterPorous2D.h
 *  \brief Definition file for class WaterPorous2D.
 */

/*! \class WaterPorous2D
 *  \brief To solve water flow equations in porous media (2-D)
 *  \details Class WaterPorous2D solves the fluid flow equations of water or any
 *  incompressible or slightly compressible fluid in a porous medium
 *  in two-dimensional configurations.
 *
 *  Porous media flows are modelled here by the Darcy law. The water, or any
 *  other fluid is considered as slightly compressible, i.e., its compressibility
 *  coefficient is constant.
 *
 *  Space discretization uses the P<sub>1</sub> (3-Node triangle) finite element 
 *  method. Time integration uses class TimeStepping that provides various
 *  well known time integration schemes.
 */


class WaterPorous2D : public Equa_Porous<real_t,3,3,2,2>
{

 public:

/// \brief Default Constructor.
/// \details Constructs an empty equation.
    WaterPorous2D();

/** \brief Constructor
 *  \details This constructor uses mesh and reservoir information
 *  @param [in] ms Mesh instance
 */
    WaterPorous2D(Mesh& ms);

/// \brief Destructor
    ~WaterPorous2D();

/** \brief Set constant coefficients
 * @param [in] cw Compressibility coefficient
 * @param [in] phi Porosity
 * @param [in] rho Density
 * @param [in] Kx x-Absolute permeability
 * @param [in] Ky y-Absolute permeability
 * @param [in] mu Viscosity
 */
    void setCoef(real_t cw,
                 real_t phi,
                 real_t rho,
                 real_t Kx,
                 real_t Ky,
                 real_t mu);

/// \brief Add mass term contribution the element matrix
    void Mass();

/// \brief Add mobility term contribution the element matrix
    void Mobility();

/** \brief Add source right-hand side term to right-hand side.
 *  @param [in] bf Vector containing source at nodes.
 */
    void BodyRHS(const Vect<real_t>& bf);

/** \brief Add boundary right-hand side term to right-hand side.
 *  @param [in] sf Vector containing source at nodes.
 */
    void BoundaryRHS(const Vect<real_t>& sf);

 private:
    void set(const Element *el);
    void set(const Side *sd);
    real_t _Kxe, _Kye;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
