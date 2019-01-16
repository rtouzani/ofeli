/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2019 Rachid Touzani

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

    Definition of abstract class 'Equa_Therm' to solve heat transfer problems

  ==============================================================================*/


#ifndef __EQUA_THERM_H
#define __EQUA_THERM_H

#include "equations/Equation.h"
#include "equations/therm/PhaseChange.h"
#include "solvers/TimeStepping.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \defgroup Therm Heat Transfer
 *  \brief Heat Transfer equations
 */

/*! \file Equa_Therm.h
 *  \brief Definition file for class Equa_Therm.
 */


class Element;
class Side;

/*! \class Equa_Therm
 *  \ingroup Therm
 * \brief Abstract class for Heat transfer Finite %Element classes.
 *
 * \tparam <T_> data type (real_t, float, ...)
 * \tparam <NEN> Number of element nodes
 * \tparam <NEE_> Number of element equations
 * \tparam <NSN_> Number of side nodes
 * \tparam <NSE_> Number of side equations
 */

template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
class Equa_Therm : virtual public Equation<T_,NEN_,NEE_,NSN_,NSE_>
{

 public:

   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theta;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::setMaterialProperty;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theMesh;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theElement;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theSide;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_terms;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_analysis;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_time_scheme;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_time_step;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_time;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_final_time;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_A;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_b;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_uu;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_u;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_bc;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_bf;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_sf;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eA0;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eA1;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eRHS;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::sRHS;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_bf_given;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_bc_given;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_sf_given;


/// \brief Default constructor.
/// \details Constructs an empty equation.
    Equa_Therm()
    {
       _terms = DIFFUSION;
       _analysis = STEADY_STATE;
       _time_scheme = STATIONARY;
       _rhocp_is_set = _conductivity_is_set = false;
       _terms = LUMPED_CAPACITY|DIFFUSION;
    }

/// \brief Destructor
    virtual ~Equa_Therm() { }

/** \brief Set stabilized formulation
 *  \details Stabilized variational formulations are to be used when the Péclet number
 *  is large.\n
 *  By default, no stabilization is used.
 */
    virtual void setStab() { _stab = true; }

/// \brief Add lumped capacity contribution to left-hand side
/// @param [in] coef coefficient to multiply by the matrix before adding [Default: <tt>1</tt>]
    virtual void LCapacityToLHS(real_t coef=1) { _coef = coef; }

/// \brief Add lumped capacity contribution to right-hand side
/// @param [in] coef coefficient to multiply by the vector before adding [Default: <tt>1</tt>]
    virtual void LCapacityToRHS(real_t coef=1) { _coef = coef; }

/// \brief Add consistent capacity contribution to left-hand side
/// @param [in] coef coefficient to multiply by the matrix before adding [Default: <tt>1</tt>]
    virtual void CapacityToLHS(real_t coef=1) { _coef = coef; }

/// \brief Add consistent capacity contribution to right-hand side
/// @param [in] coef coefficient to multiply by the vector before adding [Default: <tt>1</tt>]
    virtual void CapacityToRHS(real_t coef=1) { _coef = coef; }

/// \brief Add lumped capacity contribution to left and right-hand sides taking into
/// account time integration scheme.
    void setLumpedCapacity() 
    {
       LCapacityToLHS(1./AbsEqua<T_>::_time_step);
       LCapacityToRHS(1./AbsEqua<T_>::_time_step);
    }

/// \brief Add consistent capacity contribution to left and right-hand sides taking
/// into account time integration scheme.
    void setCapacity() 
    {
       CapacityToLHS(1./AbsEqua<T_>::_time_step);
       CapacityToRHS(1./AbsEqua<T_>::_time_step);
    }

/// \brief Add diffusion term to left-hand side
    virtual void Diffusion(real_t coef=1.) { _coef = coef; }

/// \brief Add diffusion term to right-hand side
    virtual void DiffusionToRHS(real_t coef=1.) { _coef = coef; }

/// \brief Add diffusion contribution to left and/or right-hand side taking into
/// account time integration scheme.
    void setDiffusion()
    {
       if (AbsEqua<T_>::_analysis==STEADY_STATE || AbsEqua<T_>::_time_scheme==BACKWARD_EULER)
          Diffusion(1.);
       else if (AbsEqua<T_>::_time_scheme==CRANK_NICOLSON)
          Diffusion(0.5);
       if (AbsEqua<T_>::_time_scheme==CRANK_NICOLSON)
          DiffusionToRHS(0.5);
       else if (AbsEqua<T_>::_time_scheme==FORWARD_EULER)
          DiffusionToRHS();
       else
          ;
    }

/// \brief Add convection term to left-hand side
    virtual void Convection(real_t coef=1.) { _coef = coef; }

/// \brief Add convection term to right-hand side
    virtual void ConvectionToRHS(real_t coef=1.) { _coef = coef; }

/// \brief Add convection contribution to left and/or right-hand side taking into
/// account time integration scheme.
    void setConvection()
    {
       if (AbsEqua<T_>::_time_scheme==STEADY_STATE || 
           AbsEqua<T_>::_time_scheme==BACKWARD_EULER)
          Convection(1.);
       else if (AbsEqua<T_>::_time_scheme==CRANK_NICOLSON)
          Convection(0.5);
       if (AbsEqua<T_>::_time_scheme==CRANK_NICOLSON)
          ConvectionToRHS(0.5);
       else if (AbsEqua<T_>::_time_scheme==FORWARD_EULER)
          ConvectionToRHS();
       else;
    }

/** \brief Add body right-hand side term to right-hand side.
 *  @param [in] bf Vector containing source at element nodes.
 *  @param [in] opt Vector is local (<tt>LOCAL_ARRAY</tt>) with size <tt>3</tt> or global
 *  (<tt>GLOBAL_ARRAY</tt>) with size = Number of nodes [Default: <tt>GLOBAL_ARRAY</tt>].
 */
    virtual void BodyRHS(const Vect<real_t>& bf,
                         int                 opt=GLOBAL_ARRAY) { }

/** \brief Add boundary right-hand side term to right-hand side.
 *  @param [in] sf Vector containing source at side nodes.
 *  @param [in] opt Vector is local (<tt>LOCAL_ARRAY</tt>) with size <tt>3</tt> or global
 *  (<tt>GLOBAL_ARRAY</tt>) with size = Number of nodes [Default: <tt>GLOBAL_ARRAY</tt>].
 */
    virtual void BoundaryRHS(const Vect<real_t>& sf,
                             int                 opt=GLOBAL_ARRAY) { }

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
    void build()
    {
      *_A = 0;
       if (_time_scheme==FORWARD_EULER)
          _theta = 0;
       else if (_time_scheme==BACKWARD_EULER)
          _theta = 1;
       else if (_time_scheme==CRANK_NICOLSON)
          _theta = 0.5;
       MESH_EL {
          set(theElement);
          this->ElementVector(_uu);
          if (_terms&CAPACITY)
             setCapacity();
          if (_terms&LUMPED_CAPACITY)
             setLumpedCapacity();
          if (_terms&DIFFUSION)
             setDiffusion();
          if (_terms&CONVECTION)
             setConvection();
          this->ElementAssembly(_A);
          if (_terms&SOURCE)
             BodyRHS(*_bf,GLOBAL_ARRAY);
          this->updateBC(*_bc);
          this->ElementAssembly(*_b);
       }
       if (_sf_given) {
          MESH_SD {
             set(theSide);
             BoundaryRHS(*_sf);
             this->SideAssembly(*_b);
          }
       }
    }


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
 *  @param [in] s Reference to used TimeStepping instance
 */
    void build(TimeStepping& s)
    {
       MESH_EL {
          set(theElement);
          this->ElementVector(*_u);
          if (_terms&CAPACITY)
             CapacityToLHS(1.);
          if (_terms&LUMPED_CAPACITY)
             LCapacityToLHS(1.);
          if (_terms&DIFFUSION)
             Diffusion();
          if (_terms&CONVECTION)
             Convection();
          if (_terms&SOURCE && _bf)
             BodyRHS(*_bf,GLOBAL_ARRAY);
          s.Assembly(TheElement,eRHS.get(),eA0.get(),eA1.get());
       }
       if (_sf_given) {
          MESH_SD {
             if (TheSide.isReferenced()) {
                set(theSide);
                this->SideVector(*_u);
                if (_terms&FLUX && _sf)
                   BoundaryRHS(*_sf,GLOBAL_ARRAY);
                s.SAssembly(TheSide,sRHS.get());
             }
          }
       }
    }


/** \brief Build the linear system for an eigenvalue problem
 *  @param [in] e Reference to used EigenProblemSolver instance
 */
    void build(EigenProblemSolver& e)
    {
       MESH_EL {
          set(theElement);
          this->ElementVector(*_u);
          if (_terms&CAPACITY)
             CapacityToLHS(1.);
          if (_terms&LUMPED_CAPACITY)
             LCapacityToLHS(1.);
          Diffusion();
          if (_terms&CONVECTION)
             Convection();
          e.Assembly(TheElement,eA0.get(),eA1.get());
       }
    }

/** \brief Run one time step
 *  \details This function performs one time step in equation solving.
 *  It is to be used only if a \a TRANSIENT analysis is required.
 *  @return Return error from the linear system solver
 */
    int runTransient()
    {
       *_b = 0;
       build();
       int ret=AbsEqua<T_>::solveLinearSystem(*_b,_uu);
       _u->insertBC(*_theMesh,_uu,*_bc);
       return ret;
    }

/** \brief Run one time step
 *  \details This function performs one time step in equation solving.
 *  It is identical to the function runTransient.
 *  @return Return error from the linear system solver
 */
    int runOneTimeStep() { return runTransient(); }

/** \brief Run the equation
 *  \details If the analysis (see function setAnalysis) is <tt>STEADY_STATE</tt>, then
 *  the function solves the stationary equation.\n
 *  If the analysis is <tt>TRANSIENT</tt>, then the function performs time stepping
 *  until the final time is reached.
 */
    int run()
    {
       int ret=0;
       if (_b==NULL)
          _b = new Vect<T_>(_theMesh->getNbEq());
       _uu.setSize(_theMesh->getNbEq());
       if (_analysis==STEADY_STATE) {
          build();
          ret = AbsEqua<T_>::solveLinearSystem(*_b,_uu);
          _u->insertBC(*_theMesh,_uu,*_bc);
       }
       else {
          for (; _time<=_final_time; _time+=_time_step)
             runOneTimeStep();
       }
       return ret;
    }

/// \brief Set product of Density by Specific heat (constants)
    void setRhoCp(const real_t& rhocp) { _rhocp = rhocp; }

/// \brief Set (constant) thermal conductivity
    void setConductivity(const real_t& diff) { _diff = diff; }

/// \brief Set product of Density by Specific heat given by an algebraic expression
    void RhoCp(const string& exp) { _rhocp = setMaterialProperty(exp.c_str(),"RhoCp"); }

/// \brief Set thermal conductivity given by an algebraic expression
    void Conduc(const string& exp) { _diff = setMaterialProperty(exp.c_str(),"Conductivity"); }

 protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   PhaseChange *_phase;
   virtual void set(const Element *el)=0;
   virtual void set(const Side *sd)=0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Set material properties
    void setMaterial()
    {
       theMaterial.setCode(_theElement->getCode());
       _rhocp = theMaterial.Density() * theMaterial.SpecificHeat();
       _diff  = theMaterial.ThermalConductivity();
    }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   bool     _stab, _rhocp_is_set, _conductivity_is_set;
   real_t   _diff, _rhocp, _coef;
   Vect<T_> *_vel;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
