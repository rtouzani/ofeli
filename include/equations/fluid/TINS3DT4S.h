/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani

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

                        Definition of class TINS3DT4S
             for 3-D Transient Incompressible Navier-Stokes Equations
             using P1 Stabilized elements for velocity, P1 for pressure
         Time integration is done using a second-order fractional step method

  ==============================================================================*/


#ifndef __TINS3DT4S_H
#define __TINS3DT4S_H


#include "equations/fluid/Equa_Fluid.h"
#include "linear_algebra/Vect.h"
#include "solvers/Prec.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file TINS3DT4S.h
 *  \brief Definition file for class TINS3DT4S.
 */

/*! \class TINS3DT4S
 *  \ingroup Fluid
 *  \brief Builds finite element arrays for transient incompressible fluid flow using Navier-Stokes
 *  equations in 3-D domains. 
 *  Numerical approximation uses stabilized 4-node tatrahedral finite elements for velocity and pressure.
 *  2nd-order projection scheme is used for time integration.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */


class TINS3DT4S : virtual public Equa_Fluid<4,12,3,9> {

 public:

/// Default Constructor
    TINS3DT4S();

/** \brief Constructor using mesh
 *  @param [in] ms Mesh instance
 */
    TINS3DT4S(Mesh& ms);

/** \brief Constructor using mesh and velocity
 *  @param [in] ms Mesh instance
 *  @param [in,out] u Vect instance containing initial velocity. This vector is updated
 *  during computations and will therefore contain velocity at each time step
 */
    TINS3DT4S(Mesh&         ms,
              Vect<real_t>& u);

/// \brief Destructor
    ~TINS3DT4S();

/** \brief Set equation input data
 *  @param [in] opt Parameter to select type of input (enumerated values)
 *  <ul>
 *     <li><tt>INITIAL_FIELD</tt>: Initial temperature
 *     <li><tt>BOUNDARY_CONDITION_DATA</tt>: Boundary condition (Dirichlet)
 *     <li><tt>SOURCE_DATA</tt>: Heat source
 *     <li><tt>FLUX_DATA</tt>: Heat flux (Neumann boundary condition). NOT IMPLEMENTED
 *     <li><tt>VELOCITY_FIELD</tt>: Velocity vector (for the convection term)
 *  </ul>
 *  @param [in] u Vector containing input data (Vect instance)
 */
    void setInput(EqDataType    opt,
                  Vect<real_t>& u);

/// \brief Run one time step
    int runOneTimeStep();

/// \brief Run (in the case of one step run)
    int run() { return runOneTimeStep(); }

private:

   bool             _constant_matrix;
   real_t           _cr, _c24, _vol, _cfl;
   size_t           _ne, _en[4];
   Vect<real_t>     *_p, _MM, _c, _q, _bp;
   vector<size_t>   _col_ind, _row_ptr;
   SpMatrix<real_t> _PM;
#ifndef USE_EIGEN
   Prec<real_t>     _PP, _PV;
#endif

   void build();
   void set(Element *el);
   void set(Side *sd);
   void init();
   void PressureMatrix();
   int getPressure();
   void getMomentum();
   void updateVelocity();
   void setBuyoancy();
   void ElementVelocityMatrix();
   void SideVelocityMatrix();
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
