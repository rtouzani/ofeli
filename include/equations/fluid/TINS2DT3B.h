/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2017 Rachid Touzani

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

                            Definition of class TINS2DT3B
                for 2-D Transient Incompressible Navier-Stokes Equations
                using P1+Bubble elements for velocity, P1 for pressure
       Time integration is done using a second-order fractional step method

  ==============================================================================*/


#ifndef __TINS2DT3B_H
#define __TINS2DT3B_H


#include "equations/fluid/Equa_Fluid.h"
#include "linear_algebra/Assembly.h"
#include "linear_algebra/Vect.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file TINS2DT3B.h
 *  \brief Definition file for class TINS2DT3B.
 */

/*! \class TINS2DT3B
 *  \ingroup Fluid
 *  \brief Builds finite element arrays for thermal diffusion and convection in 2-D domains 
 *  using 3-Node triangles.
 *
 * Note that members calculating element arrays have as an argument a double
 * \c coef that will be multiplied by the contribution of the current element. 
 * This makes possible testing different algorithms.
 *
 */


class TINS2DT3B : virtual public Equa_Fluid<real_t,3,6,2,4> {

 public:

/// Default Constructor
    TINS2DT3B();

/** \brief Constructor using mesh
 *  @param [in] mesh Mesh instance
 *  @param [in,out] u Vect instance containing initial velocity. This vector is updated
 *  during computations and will therefore contain velocity at each time step
 *  @param [out] p Vect instance that will contain pressure at nodes. This vector is updated
 *  during computations and will therefore contain pressure at each time step
 *  @param [in] ts Time step
 *  @param [in] Re Reynolds number. The default value (<tt>0</tt>) means that no Reynolds number
 *  is given and problem data are supplied by material properties. If Re has any other
 *  value, then nondimensional form of the equations is assumed and material properties
 *  are ignored.
 */
    TINS2DT3B(Mesh&         mesh,
              Vect<real_t>& u,
              Vect<real_t>& p,
              real_t&       ts,
              real_t        Re=0.);

/// \brief Destructor
    ~TINS2DT3B();

/** \brief Set equation input data
 *  @param [in] opt Parameter to select type of input (enumerated values)
 *  <ul>
 *     <li><tt>INITIAL_FIELD</tt>: Initial temperature
 *     <li><tt>BOUNDARY_CONDITION_DATA</tt>: Boundary condition (Dirichlet)
 *     <li><tt>SOURCE_DATA</tt>: Heat source
 *     <li><tt>FLUX_DATA</tt>: Heat flux (Neumann boundary condition)
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

/// \brief Add surface tension interface force
    int setSurfaceTension(real_t sigma);

private:

   real_t           _cfl, _Re;
   Vect<real_t>     _ub, *_p, _MM, _c, _q;
   vector<size_t>   _col_ind, _row_ptr;
   SpMatrix<real_t> _VM, _PM;
#ifndef USE_EIGEN
   Prec<real_t>     _PP, _PV;
#endif
   Point<real_t>    _g;

   void build();
   void set(Element *el);
   void set(Side *sd);
   void initEquation(Mesh& mesh, real_t ts);
   void PressureMatrix();
   void VelocityMatrix();
   int getPressure();
   void getMomentum();
   void updateVelocity();
   void setBuyoancy();
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
