/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani

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
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class NSP2DQ41 : virtual public Equa_Fluid<4,8,2,4>
{

 public:

/** \brief Constructor using mesh data
 *  @param [in] ms Mesh instance
 */
    NSP2DQ41(Mesh& ms);

/** \brief Constructor using mesh data and velocity vector
 *  @param [in] ms Mesh instance
 *  @param [in,out] u Velocity vector
 */
    NSP2DQ41(Mesh&         ms,
             Vect<real_t>& u);

/// \brief Destructor
    ~NSP2DQ41();

/** \brief Define penalty parameter
 *  \details Penalty parameter is used to enforce the incompressibility constraint
 *  @param [in] lambda Penaly parameter: Large value [Default: <tt>1.e07</tt>]
 */  
    void setPenalty(real_t lambda) { _penal = lambda; }

/** \brief Set equation input data
 *  @param [in] opt Parameter that selects data type for input. This parameter
 *  is to be chosen in the enumerated variable EqDataType
 *  @param [in] u Vect instance that contains input vector data
 *  List of data types contains <tt>INITIAL_FIELD</tt>, <tt>BOUNDARY_CONDITION_DATA</tt>, 
 *  <tt>SOURCE_DATA</tt> or <tt>FLUX</tt> with obvious meaning
 */
    void setInput(EqDataType    opt,
                  Vect<real_t>& u);

/** \brief Add contribution of periodic boundary condition (by a penalty technique).
 *  \details Boundary nodes where periodic boundary conditions are to be imposed must
 *  have codes equal to <tt>PERIODIC_A</tt> on one side and <tt>PERIODIC_B</tt> on the opposite side.
 *  @param [in] coef Value of penalty parameter [Default: <tt>1.e20</tt>].
 */
    void Periodic(real_t coef=1.e20);

/** \brief Build the linear system of equations
 *  \details Before using this function, one must have properly selected 
 *  appropriate options for:
 *  <ul>
 *     <li>The choice of a steady state or transient analysis. By default, the analysis is stationary
 *     <li>In the case of transient analysis, the choice of a time integration scheme
 *         and a lumped or consistent capacity matrix. If transient analysis is chosen, the lumped
 *         capacity matrix option is chosen by default, and the implicit Euler scheme is used
 *         by default for time integration.
 *  </ul>
 */
    void build();

/** \brief Run one time step 
 *  \details This function performs one time step, once a time integration scheme has been 
 *  selected.
 */
    int runOneTimeStep();


 private:

   Vect<real_t> *_p;
   Quad4 *_quad;
   Line2 *_ln;
   real_t _xl[4], _yl[4];
   real_t _penal;
   void set(const Element *el);
   void set(const Side *sd);

   void getPressure();

/// \brief Add element lumped mass contribution to matrix after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    void LMass(real_t coef=1.);

/// \brief Add element consistent mass contribution to matrix after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    void Mass(real_t coef=1.);

/// \brief Add element viscous contribution to matrix after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    void Viscous(real_t coef=1.);

/// \brief Add element penalty term after multiplication by <tt>coef</tt>
/// [Default: <tt>1.e07</tt>]
    void Penal(real_t coef=1.e07);

/// \brief Add element viscous contribution to right-hand side after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    void RHS_Viscous(real_t coef=1.);

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
/// @param [in] f Vector containing body forces at nodes 
    void BodyRHS(Vect<real_t>& f);

/// \brief Add boundary right-hand side term to right-hand side
/// @param [in] f Vector containing boundary forces (tractions) at nodes 
    void BoundaryRHS(Vect<real_t>& f);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
