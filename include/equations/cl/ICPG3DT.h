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

                            Definition of Class ICPG3DT
        Class to solve the Inviscid compressible fluid flows (Euler equations)
           for perfect gas in 3-D using MUSCL triangular finite volumes
                
       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/
  
#ifndef __ICPG3DT_H
#define __ICPG3DT_H

#include "equations/cl/Muscl3DT.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file ICPG3DT.h
 *  \brief Definition file for class ICPG3DT.
 */

template<class T_,size_t N_> class LocalVect;
class Mesh;
class Side;

/*! \class ICPG3DT
 *  \ingroup ConservationLaws
 *  \brief Class to solve the Inviscid compressible fluid flows (Euler equations) for perfect gas in 3-D
 *
 *  \details Solution method is a second-order MUSCL Finite Volume scheme with tetrahedra
 *
 * \author S. Clain, V. Clauzon
 * \copyright GNU Lesser Public License
 */
class ICPG3DT : public Muscl3DT {

public:

   using Muscl3DT::_nb_sides;
   using Muscl3DT::_nb_elements;

/// \brief Constructor using mesh data
/// @param [in] ms Mesh instance
    ICPG3DT(Mesh& ms);

/** \brief Constructor using mesh and initial data
    @param [in] ms Mesh instance
    @param [in] r Elementwise initial density vector (as instance of Element Vect)
    @param [in] v Elementwise initial velocity vector (as instance of Element Vect)
    @param [in] p Elementwise initial pressure vector (as instance of Element Vect)
 */
    ICPG3DT(Mesh&         ms,
            Vect<real_t>& r,
            Vect<real_t>& v,
            Vect<real_t>& p);

/// \brief Destructor
    ~ICPG3DT();

/// \brief Reconstruct.
/// \details exit(3) if reconstruction failed
    void setReconstruction();

/// \brief Advance one time step
    real_t runOneTimeStep();

/// \brief Add flux to field
//  \details If this function is used, user must call getFlux himself
    void Forward(const Vect<real_t>& flux,
                       Vect<real_t>& field);

/// \brief Return flux
    real_t getFlux();

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef real_t (ICPG3DT::*e_fptr)(int); // type member function pointer
    static e_fptr fsolver[7];               // Array of member function pointer. MUST be initialized in .cpp only!!!
    e_fptr mySolver; // function pointer on solver is used to get rid of an embarrassing "if"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Choose solver
    void setSolver(SolverType solver)
    {
       if (solver>=6)
          mySolver = nullptr;
       else 
          mySolver = fsolver[int(solver)];
	}

/// \brief Assign a reference length
    void setReferenceLength(real_t dx) { _ReferenceLength=dx; }

/// \brief Assign a time step
    void setTimeStep(real_t dt) { _TimeStep = dt; }

/// \brief Assign CFL value
    void setCFL(real_t CFL) { _CFL = CFL; };

/// \brief Return reference length
    real_t getReferenceLength() const { return _ReferenceLength; }

/// \brief Return time step
    real_t getTimeStep() const { return _TimeStep; }

/// \brief Return CFL
    real_t getCFL() const { return _CFL; }

/// \brief Set &gamma; value
    void setGamma(real_t gamma) { _Gamma = gamma; }

/// \brief Set value of C<sub>v</sub> (Heat capacity at constant volume)
    void setCv(real_t Cv) { _Cv = Cv; }

/// \brief Set value of C<sub>p</sub> (Heat capacity at constant pressure)
    void setCp(real_t Cp) { _Cp = Cp; }

/// \brief Set Kappa value
    void setKappa(real_t Kappa) { _Kappa = Kappa; }

/// \brief Return value of &gamma;
    real_t getGamma() const { return _Gamma; }

/// \brief Return value of C<sub>v</sub> (Heat capacity at constant volume)
    real_t getCv() const { return _Cv; }

/// \brief Return value of C<sub>p</sub> (Heat capacity at constant pressure)
    real_t getCp() const { return _Cp; }

/// \brief Return value of \f$\kappa\f$
    real_t getKappa() const { return _Kappa; }

/// \brief Return reference to mesh instance
    Mesh &getMesh() { return *_theMesh; }

/// \brief Return pointer to mesh
    Mesh *getPtrMesh() { return _theMesh; }

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
    void setInitialConditionShockTube(const LocalVect<real_t,5>& BcG,
                                      const LocalVect<real_t,5>& BcD,
                                      real_t                     x0);

/// \brief Set initial condition
    void setInitialCondition(const LocalVect<real_t,5>& u);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setBC(const Side &sd, real_t u);
    void setBC(int code, real_t u);
    void setBC(real_t u);
    void setBC(const Side &sd, const LocalVect<real_t,5> &u);
    void setBC(int code, const LocalVect<real_t,5> &u);
    void setBC(const LocalVect<real_t,5> &u);
    void setInOutFlowBC(const Side &sd, const LocalVect<real_t,5> &u);
    void setInOutFlowBC(int code, const LocalVect<real_t,5> &u);
    void setInOutFlowBC(const LocalVect<real_t,5> &u);
    void forward();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS

// primary variables are used here for positivity reason
   Vect<real_t> *_r, *_v, *_p;

// memory allocation and default paramaters
   void Init();

// solvers
   real_t RiemannSolver_ROE(int s);
   real_t RiemannSolver_VFROE(int s);
   real_t RiemannSolver_LF(int s);
   real_t RiemannSolver_RUSANOV(int s);
   real_t RiemannSolver_HLL(int s);
   real_t RiemannSolver_HLLC(int s);

// temp variables used by all solvers
   real_t _tGr, _Grv[3], _tGv[3], _tGp, _GT, _Gc, _GH, _Ge, _GE;
   real_t _tDr, _Drv[3], _tDv[3], _tDp, _DT, _Dc, _DH, _De, _DE;

// primary variables, left side
   Vect<real_t> _Gr, _Gv, _Gp;

// primary variables, right side
   Vect<real_t> _Dr, _Dv, _Dp;

// !!!conservative!!! numerical flux
   Vect<real_t> _Fr, _Frv, _FE;

// conversion function
   void fromPrimalToConservative();
   void fromConservativeToPrimal();

/// physical constants
    real_t _Gamma, _Cv, _Cp, _Kappa;

//  dt, dx, CFL
    real_t _TimeStep, _ReferenceLength, _CFL;
    bool _init_alloc;
    real_t _min_density;

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
