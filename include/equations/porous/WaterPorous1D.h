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

                          Definition of class WaterPorous1D

  ==============================================================================*/


#ifndef __WATER_POROUS_1D_H
#define __WATER_POROUS_1D_H

#include "equations/porous/Equa_Porous.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file WaterPorous1D.h
 *  \brief Definition file for class WaterPorous1D.
 */

/*! \class WaterPorous2D
 *  \brief To solve water flow equations in porous media (1-D)
 *  \details Class WaterPorous2D solves the fluid flow equations of water or any
 *  incompressible or slightly compressible fluid in a porous medium
 *  in two-dimensional configurations.
 *
 *  Porous media flows are modelled here by the Darcy law. The water, or any
 *  other fluid is considered as slightly compressible, i.e., its compressibility
 *  coefficient is constant.
 *
 *  Space discretization uses the P<sub>1</sub> (2-Node line) finite element 
 *  method. Time integration uses class TimeStepping that provides various
 *  well known time integration schemes.
 */


class WaterPorous1D : public Equa_Porous<real_t,2,2,1,1>
{

 public:

/// \brief Default Constructor.
/// \details Constructs an empty equation.
    WaterPorous1D();

/** \brief Constructor
 *  \details This constructor uses mesh and reservoir information
 *  @param [in] ms Mesh instance
 *  @param [in] verb Verbosity parameter
 */
    WaterPorous1D(Mesh&  ms,
                  size_t verb=1);

/// \brief Destructor
    ~WaterPorous1D();

/** \brief Set constant coefficients
 * @param [in] cw Compressibility coefficient
 * @param [in] phi Porosity
 * @param [in] rho Density
 * @param [in] K Absolute permeability
 * @param [in] mu Viscosity
 */
    void setCoef(real_t cw,
                 real_t phi,
                 real_t rho,
                 real_t K,
                 real_t mu);

/// \brief Add mass term contribution the element matrix
    void Mass();

/// \brief Add mobility term contribution the element matrix
    void Mobility();

/** \brief Add source right-hand side term to right-hand side.
 *  @param [in] bf Vector containing source at element nodes.
 *  @param [in] opt Vector is local (<tt>LOCAL_ARRAY</tt>) with size <tt>3</tt> or global
 *  (<tt>GLOBAL_ARRAY</tt>) with size = Number of nodes [Default: <tt>GLOBAL_ARRAY</tt>].
 */
    void BodyRHS(const Vect<real_t>& bf,
                 int                 opt=GLOBAL_ARRAY);

 private:
    void set(const Element *el);
    size_t _nb_nodes, _nb_elements, _verb;
    real_t _Ke, _length;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
