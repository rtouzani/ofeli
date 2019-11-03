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
#include "mesh/Material.h"

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
extern Material theMaterial;

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

   using AbsEqua<T_>::setInput;
   using AbsEqua<T_>::setTerms;
   using AbsEqua<T_>::_u;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::setMaterialProperty;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theMesh;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theElement;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theSide;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_terms;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_analysis;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_TimeInt;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_eu;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eA0;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eA1;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::sA0;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eMat;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eRHS;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::sRHS;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_nb_nodes;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_nb_sides;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_nb_el;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_nb_eq;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_nb_dof_total;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_nb_dof;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_el_geo;


/// \brief Default constructor.
/// \details Constructs an empty equation.
    Equa_Therm()
    {
       _terms = DIFFUSION;
       _analysis = STEADY_STATE;
       _TimeInt.scheme = NONE;
       _rhocp_is_set = _conductivity_is_set = false;
       _terms = DIFFUSION;
    }

/// \brief Destructor
    virtual ~Equa_Therm() { }

/** \brief Set stabilized formulation
 *  \details Stabilized variational formulations are to be used when the Péclet number
 *  is large.\n
 *  By default, no stabilization is used.
 */
    virtual void setStab() { _stab = true; }

/// \brief Add lumped capacity contribution to element matrix
/// @param [in] coef coefficient to multiply by the matrix before adding [Default: <tt>1</tt>]
    virtual void LCapacity(real_t coef=1) { _coef = coef; }

/// \brief Add consistent capacity contribution to left-hand side
/// @param [in] coef coefficient to multiply by the matrix before adding [Default: <tt>1</tt>]
    virtual void Capacity(real_t coef=1) { _coef = coef; }

/// \brief Add diffusion term to element matrix
    virtual void Diffusion(real_t coef=1.) { _coef = coef; }

/// \brief Add convection term to element matrix
    virtual void Convection(real_t coef=1.) { _coef = coef; }

/** \brief Add body right-hand side term to right-hand side.
 *  @param [in] f Vector containing source at nodes.
 */
    virtual void BodyRHS(const Vect<real_t>& f) { }

/** \brief Add boundary right-hand side term to right-hand side.
 *  @param [in] f Vector containing source at nodes.
 */
    virtual void BoundaryRHS(const Vect<real_t>& f) { }

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
       static bool matrix_set = false;
       real_t dti = 1./AbsEqua<T_>::_TimeInt.delta;
       if (AbsEqua<T_>::_A==nullptr && !matrix_set) {
          AbsEqua<T_>::setMatrixType(SPARSE);
          if (_terms&CONVECTION)
             AbsEqua<real_t>::setSolver(BICG_STAB_SOLVER,DILU_PREC);
          else
             AbsEqua<T_>::setSolver(CG_SOLVER,DILU_PREC);
          matrix_set = true;
       }
       AbsEqua<T_>::_A->clear();
       _TimeInt.theta = 1;
       if (_TimeInt.scheme==FORWARD_EULER)
          _TimeInt.theta = 0;
       else if (_TimeInt.scheme==CRANK_NICOLSON)
          _TimeInt.theta = 0.5;

       mesh_elements(*_theMesh) {
          set(the_element);
          if (_terms&CAPACITY)
             Capacity();
          if (_terms&LUMPED_CAPACITY)
             LCapacity();
          if (_terms&DIFFUSION)
             Diffusion();
          if (_terms&CONVECTION)
             Convection();
          if (_TimeInt.scheme==BACKWARD_EULER)
             eMat = dti*eA1 + eA0;
          else if (_TimeInt.scheme==FORWARD_EULER)
             eMat = dti*eA1;
          else if (_TimeInt.scheme==CRANK_NICOLSON)
             eMat = dti*eA1 + 0.5*eA0;
          AbsEqua<T_>::_A->Assembly(The_element,eMat.get());
          if (AbsEqua<T_>::_bf!=nullptr)
             BodyRHS(*AbsEqua<T_>::_bf);
          if (AbsEqua<T_>::_bc!=nullptr)
             this->updateBC(The_element,*AbsEqua<T_>::_bc);
          if (_TimeInt.scheme==BACKWARD_EULER)
             eRHS += dti*eA1*_eu;
          else if (_TimeInt.scheme==FORWARD_EULER)
             eRHS += (dti*eA1-eA0)*_eu;
          else if (_TimeInt.scheme==CRANK_NICOLSON)
             eRHS += (dti*eA1-0.5*eA0)*_eu;
          AbsEqua<T_>::_b->Assembly(The_element,eRHS.get());
       }
       if (AbsEqua<T_>::_sf!=nullptr) {
          mesh_boundary_sides(*_theMesh) {
             set(the_side);
             BoundaryRHS(*AbsEqua<T_>::_sf);
             AbsEqua<T_>::_A->Assembly(The_side,sA0.get());
             AbsEqua<T_>::_b->Assembly(The_side,sRHS.get());
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
      if (_u==nullptr)
         throw OFELIException("In Equa_Therm::build: No solution vector given.");
       mesh_elements(*_theMesh) {
          set(the_element);
          if (_terms&CAPACITY)
             Capacity(1.);
          if (_terms&LUMPED_CAPACITY)
             LCapacity();
          if (_terms&DIFFUSION)
             Diffusion();
          if (_terms&CONVECTION)
             Convection();
          if (_terms&SOURCE && AbsEqua<T_>::_bf!=nullptr)
             BodyRHS(*AbsEqua<T_>::_bf);
          s.Assembly(*_theElement,eRHS.get(),eA0.get(),eA1.get());
       }
       if (AbsEqua<T_>::_sf!=nullptr) {
          mesh_sides(*_theMesh) {
             if (The_side.isReferenced()) {
                set(the_side);
                this->SideVector(*AbsEqua<T_>::_u);
                if (_terms&FLUX && AbsEqua<T_>::_sf!=nullptr)
                   BoundaryRHS(*AbsEqua<T_>::_sf);
                s.SAssembly(*_theSide,sRHS.get());
             }
          }
       }
    }

/** \brief Build the linear system for an eigenvalue problem
 *  @param [in] e Reference to used EigenProblemSolver instance
 */
    void build(EigenProblemSolver& e)
    {
       mesh_elements(*_theMesh) {
          set(the_element);
          this->ElementVector(*AbsEqua<T_>::_u);
          if (_terms&CAPACITY)
             Capacity(1.);
          if (_terms&LUMPED_CAPACITY)
             LCapacity(1.);
          Diffusion();
          if (_terms&CONVECTION)
             Convection();
          e.Assembly(*_theElement,eA0.get(),eA1.get());
       }
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
