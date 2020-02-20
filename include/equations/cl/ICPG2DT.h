/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

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

                            Implementation of Class ICPG2DT
                Class to solve compressible Euler equations for oerfect gas
                         in 2-D using MUSCL finite volumes
                
       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/

#ifndef  __ICPG2DT_H
#define  __ICPG2DT_H

#include "equations/cl/Muscl2DT.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file ICPG2DT.h
 *  \brief Definition file for class ICPG2DT.
 */

template<class T_,size_t N_> class LocalVect;
class Mesh;
class Side;

/*! \class ICPG2DT
 *  \ingroup ConservationLaws
 *  \brief Class to solve the Inviscid compressible fluid flows (Euler equations) for perfect gas in 2-D
 *
 *  \details Solution method is a second-order MUSCL Finite Volume scheme on triangles
 *
 * \author S. Clain, V. Clauzon
 * \copyright GNU Lesser Public License
 */
class ICPG2DT : public Muscl2DT {

 public:

   using Muscl2DT::_nb_sides;
   using Muscl2DT::_nb_elements;

/// \brief Constructor using mesh instance
    ICPG2DT(Mesh &ms);

/** \brief Constructor using mesh and initial data
    @param [in] ms Mesh instance
    @param [in] r Initial density vector (as instance of Vect)
    @param [in] v Initial velocity vector (as instance of Vect)
    @param [in] p Initial pressure vector (as instance of Vect)
 */
    ICPG2DT(Mesh&         ms,
            Vect<real_t>& r,
            Vect<real_t>& v,
            Vect<real_t>& p);

/// \brief Destructor
    ~ICPG2DT();

/// \brief Reconstruct.
/// \details exit(3) if reconstruction fails
    void setReconstruction();

/// \brief Advance one time step
    real_t runOneTimeStep();

/// \brief Add Flux to Field
/// \details If this function is used, the function getFlux must be called
    void Forward(const Vect<real_t>& Flux,
                       Vect<real_t>& Field);

/// \brief Get flux
    real_t getFlux();

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef real_t (ICPG2DT::* e_fptr)(int) ;    // type member function pointer
    static e_fptr fsolver[7];	                 // Array of member function pointer. MUST be initialized in .cpp only!!!
    e_fptr mySolver;                             // function pointer on solver is used to avoid a heavy "if"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Choose solver
 *  @param [in] s Index of solver in the enumerated variable SolverType
 *  Available values are: <tt>ROE_SOLVER</tt>, <tt>VFROE_SOLVER</tt>, <tt>LF_SOLVER</tt>,
 *  <tt>RUSANOV_SOLVER</tt>, <tt>HLL_SOLVER</tt>, <tt>HLLC_SOLVER</tt>, <tt>MAX_SOLVER</tt>
 */
    void setSolver(SolverType s);

/// \brief Set Gamma value
    void setGamma(real_t gamma) { _Gamma = gamma; }

/// \brief Set value of heat capacity at constant volume
    void setCv(real_t Cv) { _Cv = Cv; }

/// \brief Set value of heat capacity at constant pressure
    void setCp(real_t Cp) { _Cp = Cp; }

/// \brief Set Kappa value
    void setKappa(real_t Kappa) { _Kappa = Kappa; }
	
/// \brief Return value of Gamma
    real_t getGamma() const { return _Gamma; }

/// \brief Return value of heat capacity at constant volume
    real_t getCv() const { return _Cv; }

/// \brief Return value of heat capacity at constant pressure
    real_t getCp() const { return _Cp; }

/// \brief Return value of Kappa
    real_t getKappa() const { return _Kappa; }

/// \brief Return reference to mesh instance
    Mesh &getMesh() { return* _theMesh; }

/// \brief Calculate elementwise momentum
    void getMomentum(Vect<real_t>& m) const;

/// \brief Calculate elementwise internal energy
    void getInternalEnergy(Vect<real_t>& e) const;

/// \brief Return elementwise total energy
    void getTotalEnergy(Vect<real_t>& e) const;

/// \brief Return elementwise sound speed
    void getSoundSpeed(Vect<real_t>& s) const;

/// \brief Return elementwise Mach number
    void getMach(Vect<real_t>& m) const;

/// \brief Set initial condition for the schock tube problem
    void setInitialConditionShockTube(const LocalVect<real_t,4>& BcL,
                                      const LocalVect<real_t,4>& BcR,
                                      real_t                     x0);

/// \brief Set initial condition
    void setInitialCondition(const LocalVect<real_t,4>& u);

/** \brief Prescribe a constant boundary condition at given side
 *  @param [in] sd Reference to Side instance
 *  @param [in] a Value to prescribe
 */
    void setBC(const Side&  sd,
               real_t      a);

/** \brief Prescribe a constant boundary condition for a given code
 *  @param [in] code Code for which value is imposed
 *  @param [in] a Value to prescribe
 */
    void setBC(int    code,
               real_t a);

/// \brief Prescribe a constant boundary condition on all boundary sides
/// @param [in] u Value to prescribe
    void setBC(real_t u);

/** \brief Prescribe a constant boundary condition at a given side
 *  @param [in] sd Reference to Side instance
 *  @param [in] u Vector (instance of class LocalVect) with as components the 
 *  constant values to prescribe for the four fields (<tt>r</tt>, <tt>vx</tt>,
 *  <tt>vy</tt>, <tt>p</tt>)
 */
    void setBC(const Side&                sd,
               const LocalVect<real_t,4>& u);

/** \brief Prescribe a constant boundary condition for a given code
 *  @param [in] code Code for which value is imposed
 *  @param [in] u Vector (instance of class LocalVect) with as components the 
 *  constant values to prescribe for the four fields (<tt>r</tt>, <tt>vx</tt>,
 *  <tt>vy</tt>, <tt>p</tt>)
 */
    void setBC(int                        code,
               const LocalVect<real_t,4>& u);

/** \brief Prescribe a constant boundary condition at all boundary sides
 *  @param [in] u Vector (instance of class LocalVect) with as components the 
 *  constant values to prescribe for the four fields (<tt>r</tt>, <tt>vx</tt>, 
 *  <tt>vy</tt>, <tt>p</tt>)
 */
    void setBC(const LocalVect<real_t,4>& u);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setInOutFlowBC(const Side&                sd,
                        const LocalVect<real_t,4>& u);
    void setInOutFlowBC(int                        code,
                        const LocalVect<real_t,4>& u);
    void setInOutFlowBC(const LocalVect<real_t,4>& u);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Return density at given element label
/// @param [in] i Element label
    real_t getR(size_t i) const { return (*_r)(i); }

/// Return velocity at given element label
/// @param [in] i Element label
/// @param [in] j component index (<tt>1</tt> or <tt>2</tt>)
    real_t getV(size_t i, size_t j) const { return (*_v)(i,j); }

/// \brief Return pressure at given element label
/// @param [in] i Element label
    real_t getP(size_t i) const { return (*_p)(i); }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
//  add flux to value
    void forward();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS

//  primary variables are used here for positivity reason
    Vect<real_t> *_r, *_v, *_p;

//  memory allocation and default paramaters
    void init();

/// \brief set Roe Riemann solver
    real_t RiemannSolver_ROE(int s);

/// \brief set Roe Riemann solver
    real_t RiemannSolver_VFROE(int s);

/// \brief set Roe Riemann solver
    real_t RiemannSolver_LF(int s);

/// \brief set Roe Riemann solver
    real_t RiemannSolver_RUSANOV(int s);

/// \brief set Roe Riemann solver
    real_t RiemannSolver_HLL(int s);

/// \brief set Roe Riemann solver
    real_t RiemannSolver_HLLC(int s);

/// \brief Set passage for primal to conservative variables
    void fromPrimalToConservative();

/// \brief Set passage for conservative to primal variables
    void fromConservativeToPrimal();

//  temp variables used by all solvers
    real_t _tLr, _tLp, _Lc, _LH, _Le;
    real_t _tRr, _tRp, _Rc, _RH, _Re;
    real_t _tLv[2], _tRv[2], _Lrv[2], _Rrv[2];

//  primary variables, left and right sides
    Vect<real_t> _Lr, _Lp, _Lv, _Rr, _Rp, _Rv;

// !!!conservative!!! numerical flux
    Vect<real_t> _Fr, _FE, _Frv;
    real_t _Gamma, _Cv, _Cp, _Kappa, _min_density;
    bool _init_alloc;	

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
