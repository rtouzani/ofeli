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

                            Definition of Class ICPG1D
        Class to solve the Inviscid compressible fluid flows (Euler equations)
                    for perfect gas in 1-D using MUSCL finite volumes
                
       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/

#ifndef __ICPG1D_H
#define __ICPG1D_H

#include "equations/cl/Muscl1D.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file ICPG1D.h
 *  \brief Definition file for class ICPG1D.
 */

template<class T_,size_t N_> class LocalVect;
class Mesh;
class Side;

/*! \class ICPG1D
 *  \ingroup ConservationLaws
 *  \brief Class to solve the Inviscid compressible fluid flows (Euler equations) for 
 *  perfect gas in 1-D
 *
 *  Solution method is a second-order MUSCL Finite Volume scheme
 *
 */
class ICPG1D : public Muscl1D {

 public:

   using Muscl1D::_nb_sides;
   using Muscl1D::_nb_elements;

/// \brief Constructor using Mesh instance
    ICPG1D(Mesh& ms);

/** \brief Constructor using mesh and initial data
 *  @param [in] ms Reference to Mesh instance
 *  @param [in] r Vector containing initial (elementwise) density
 *  @param [in] v Vector containing initial (elementwise) velocity 
 *  @param [in] p Vector containing initial (elementwise) pressure 
 */
    ICPG1D(Mesh&         ms,
           Vect<real_t>& r,
           Vect<real_t>& v,
           Vect<real_t>& p);

/// \brief Destructor
    ~ICPG1D();

/// \brief Set reconstruction from class Muscl
    void setReconstruction();

/// \brief Advance one time step
    real_t runOneTimeStep();

/** \brief Add flux to field.
    \details If this function is used, the user must call getFlux himself
    @param [in] flux Vector containing fluxes at sides (points)
    @param [out] field Vector containing solution vector
 */
    void Forward(const Vect<real_t>& flux,
                       Vect<real_t>& field);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef real_t (ICPG1D::* e_fptr)(int) ;  // type member function pointer
    static e_fptr fsolver[7];                 // array of member function pointer. MUST be initialized in .cpp only!!!
    e_fptr mySolver;                          // function pointer on solver is used to avoid "if" condition
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Choose solver type
    void setSolver(SolverType solver)
    {
       if (solver>6)
          mySolver = NULL;
       else
          mySolver = fsolver[int(solver)];
    }

/// \brief Set value of constant Gamma for gases
    void setGamma(real_t gamma) { _Gamma = gamma; }

/// \brief Set value of Cv (specific heat at constant volume)
    void setCv(real_t Cv) { _Cv = Cv; }

/// \brief Set value of C<sub>p</sub> (specific heat at constant pressure)
    void setCp(real_t Cp) { _Cp = Cp; }

/// \brief Set value of constant Kappa
    void setKappa(real_t Kappa) { _Kappa = Kappa; }

/// \brief Return value of constant Gamma
    real_t getGamma() const { return _Gamma; }

/// \brief Return value of C<sub>v</sub> (specific heat at constant volume)
    real_t getCv() const { return _Cv; }

/// \brief Return value of C<sub>p</sub> (specific heat at constant pressure)
    real_t getCp() const { return _Cp; }

/// \brief Return value of constant Kappa
    real_t getKappa() const { return _Kappa; }

/// \brief Get vector of momentum at elements
/// @param [in,out] m Vect instance that contains on output element momentum
    void getMomentum(Vect<real_t>& m) const;

/// \brief Get vector of internal energy at elements
/// @param [in,out] ie Vect instance that contains on output element internal energy
    void getInternalEnergy(Vect<real_t>& ie) const;

/// \brief Get vector of total energy at elements
/// @param [in,out] te Vect instance that contains on output element total energy
    void getTotalEnergy(Vect<real_t>& te) const;

/// \brief Get vector of sound speed at elements
/// @param [in,out] s Vect instance that contains on output element sound speed
    void getSoundSpeed(Vect<real_t>& s) const;

/// \brief Get vector of elementwise Mach number
/// @param [in,out] m Vect instance that contains on output element Mach number
    void getMach(Vect<real_t>& m) const;

/// \brief Initial condition corresponding to the shock tube
    void setInitialCondition_shock_tube(const LocalVect<real_t,3>& BcG,
                                        const LocalVect<real_t,3>& BcD,
                                              real_t               x0);

/// \brief A constant initial condition
/// @param [in] u LocalVect instance containing density, velocity and pressure 
    void setInitialCondition(const LocalVect<real_t,3>& u);

/** \brief Assign a boundary condition as a constant to a given side
 *  @param [in] sd Side to which the value is assigned
 *  @param [in] u Value to assign
 */
    void setBC(const Side&  sd,
                     real_t u);

/** \brief Assign a boundary condition value
 *  @param [in] code Code value to which boundary condition is assigned
 *  @param [in] a Value to assign to sides that have code <tt>code</tt>
 */
    void setBC(int    code,
               real_t a);

/// \brief Assign a boundary condition value
/// @param [in] a Value to assign to all boundary sides
    void setBC(real_t a);

/** \brief Assign a Dirichlet boundary condition vector
    @param [in] sd Side instance to which the values are assigned
    @param [in] u LocalVect instance that contains values to assign to the side
 */
    void setBC(const Side&                sd,
               const LocalVect<real_t,3>& u);

/** \brief Assign a Dirichlet boundary condition vector
    @param [in] code Side code for which the values are assigned
    @param [in] U LocalVect instance that contains values to assign to sides with code \a code
 */
    void setBC(      int                  code,
               const LocalVect<real_t,3>& U);

/// \brief Assign a Dirichlet boundary condition vector
/// @param [in] u LocalVect instance that contains values to assign to all boundary sides
    void setBC(const LocalVect<real_t,3>& u);

/** \brief Impose a constant inflow or outflow boundary condition on a given side
 *  @param [in] sd Instance of Side on which the condition is prescribed
 *  @param [in] u LocalVect instance that contains values to assign to the side
 */
    void setInOutflowBC(const Side&                sd,
                        const LocalVect<real_t,3>& u);

/** \brief Impose a constant inflow or outflow boundary condition on sides with a given code
 *  @param [in] code Value of code for which the condition is prescribed
 *  @param [in] u LocalVect instance that contains values to assign to the sides
 */
    void setInOutflowBC(      int                  code,
                        const LocalVect<real_t,3>& u);

/// \brief Impose a constant inflow or outflow boundary condition on boundary sides
/// @param [in] u LocalVect instance that contains values to assign to the sides
    void setInOutflowBC(const LocalVect<real_t,3>& u);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
//  add flux to value
    void forward();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS

//  primary variables are used here for positivity
    Vect<real_t> *_r, *_v, *_p;

//  memory allocation and default paramaters
    void init();

//  solvers
    real_t RiemannSolver_ROE(int s);
    real_t RiemannSolver_VFROE(int s);
    real_t RiemannSolver_LF(int s);
    real_t RiemannSolver_RUSANOV(int s);
    real_t RiemannSolver_HLL(int s);
    real_t RiemannSolver_HLLC(int s);

/// compute flux
    real_t getFlux();

//  conversion function
    void fromPrimalToConservative();
    void fromConservativeToPrimal();

//  temp variables used by all solvers
    real_t _Gc, _Dc;

//  primary variables, left-hand side
    Vect<real_t> _Gr, _Gv, _Gp;

//  primary variables, right-hand side
    Vect<real_t> _Dr, _Dv, _Dp;

//  !!!conservative!!! numerical flux
    Vect<real_t> _Fr, _FrU, _FE;

    real_t _Gamma, _Cv, _Cp, _Kappa, _min_density;
    bool _init_alloc;

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
