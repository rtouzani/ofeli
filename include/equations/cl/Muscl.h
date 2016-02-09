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

                            Definition of Class Muscl
                      Mother class for Muscl finite volumes
                
       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/


#ifndef __MUSCL_H
#define __MUSCL_H

#include "OFELI_Config.h"
#include "linear_algebra/Vect.h"

namespace OFELI {

/*! \file Muscl.h
 *  \brief Definition file for class Muscl.
 */

/** \defgroup ConservationLaws Conservation Law Equations
 *  \brief For solvers of conservation law equations
 */


class Mesh;

/*! \class Muscl
 *  \ingroup ConservationLaws
 *  \brief Parent class for hyperbolic solvers with %Muscl scheme.
 *
 * Everything here is common for both 2D and 3D muscl methods !
 * Virtual functions are implemented in Muscl2D and Muscl3D classes
 *
 */
 
class Muscl {

 public:

/// \brief Constructor using mesh instance
    Muscl(Mesh& m);

/// \brief Destructor
    virtual ~Muscl() { }

/// \enum Method
/// \brief Enumeration for flux choice
    enum Method {
       FIRST_ORDER_METHOD   =  0,        /*!< First Order upwind method */
       MULTI_SLOPE_Q_METHOD =  1,        /*!< Multislope Q method       */
       MULTI_SLOPE_M_METHOD =  2,        /*!< Multislope M method       */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
       M_FORCE_WORD         = 0xFFFFFFFF /*!< for compatibility reasons, gcc may return a warning anyway */
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
    };

/// \enum Limiter
/// \brief Enumeration of flux limiting methods
    enum Limiter {
       MINMOD_LIMITER    = 0,          /*!< MinMod limiter              */
       VANLEER_LIMITER   = 1,          /*!< Van Leer limiter            */
       SUPERBEE_LIMITER  = 2,          /*!< Superbee limiter            */
       VANALBADA_LIMITER = 3,          /*!< Van Albada limiter          */
       MAX_LIMITER       = 4,          /*!< Max limiter                 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
       L_FORCE_WORD      = 0xFFFFFFFF  /*!< for compatibility reason, gcc may return a warning anyway */
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
    };

/// \enum SolverType
/// \brief Enumeration of various solvers for the Riemann problem
    enum SolverType {
       ROE_SOLVER     = 0,              /*!< Roe solver                  */
       VFROE_SOLVER   = 1,              /*!< Finite Volume Roe solver    */
       LF_SOLVER      = 2,              /*!< LF solver                   */
       RUSANOV_SOLVER = 3,              /*!< Rusanov solver              */
       HLL_SOLVER     = 4,              /*!< HLL solver                  */
       HLLC_SOLVER    = 5,              /*!< HLLC solver                 */
       MAX_SOLVER     = 6,              /*!< Max solver                  */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
       FORCE_WORD     = 0xFFFFFFFF      /*!< for compatibility reason, gcc may return a warning anyway */
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
    };

/// \brief Assign time step value
/// @param [in] dt Time step value
    void setTimeStep(real_t dt) { _TimeStep = dt; }

/// \brief Return time step value
    real_t getTimeStep() const { return _TimeStep; }

/// \brief Assign CFL value
/// @param [in] CFL Value of CFL
    void setCFL(real_t CFL) { _CFL = CFL; }

/// \brief Return CFL value
    real_t getCFL() const { return _CFL; }

/// \brief Assign reference length value
/// @param [in] dx Value of reference length
    void setReferenceLength(real_t dx) { _ReferenceLength = dx; }

/// \brief Return reference length
    real_t getReferenceLength() const { return _ReferenceLength; }

/// \brief Return reference to Mesh instance
    Mesh &getMesh() const { return *_theMesh; }

/// \brief Set verbosity parameter
/// @param [in] v Value of verbosity parameter
    void setVerbose(int v) { _verbose = v; }

/** \brief Function to reconstruct by the Muscl method
 *  @param [in] U Field to reconstruct
 *  @param [out] LU Left gradient vector
 *  @param [out] RU Right gradient vector
 *  @param [in] dof Label of dof to reconstruct
 */
    bool setReconstruction(const Vect<real_t>& U,
                                 Vect<real_t>& LU,
                                 Vect<real_t>& RU,
                                 size_t        dof);

/// \brief Choose a flux solver
/// @param [in] s Solver to choose
    void setMethod(const Method &s) { _method = s; }

/// \brief Choose a code for solid zone
    void setSolidZoneCode(int c) { _solid_zone_code = c; _solid_zone = true; }

/// \brief Return flag for presence of solid zones
    bool getSolidZone() const { return _solid_zone; }

/// \brief Return code of solid zone, 0 if this one is not present
    int getSolidZoneCode() const { if (_solid_zone) return _solid_zone_code; else return 0; }

/// \brief Choose a flux limiter
/// @param [in] l Limiter to choose
    void setLimiter(Limiter l) {
         if (l>=2)
            myLimiter = NULL;
         else
            myLimiter = flimiter[int(l)];
    }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
//  type member function pointer
    typedef real_t (Muscl::* m_fptr) (real_t,real_t);
//  Array of member function pointers. Must be initialized in .cpp only!!! 
    static m_fptr flimiter[5];
//  function pointer on solver is used to avoid "if" statements
    m_fptr myLimiter;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    Mesh    *_theMesh;
    real_t  _TimeStep, _CFL, _ReferenceLength;
    int     _verbose, _solid_zone_code;
    size_t  _nb_sides, _nb_elements;

//  this is the main function here, compute all geometric stuff. Could take a while
    virtual void Initialize() = 0;
    Method _method;
    Limiter _limiter;
    bool _solid_zone;

// limiters // !!!!! ++ tau_lim !!!!

    typedef unsigned long long ull;

//  strongly optimized !
//  Take care for portability !!!!
    real_t minmod(real_t val_plus, real_t val_minus);
    real_t superbee(real_t val_plus, real_t val_minus);
    real_t vanleer(real_t val_plus, real_t val_minus);
    real_t vanalbada(real_t val_plus, real_t val_minus);

//  contains for multislope M method
    real_t m_limiter(real_t a, real_t lim1, real_t lim2);
    real_t m_limiter2(real_t a, real_t s, real_t p, real_t d);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

} /* namespace OFELI */

#endif
