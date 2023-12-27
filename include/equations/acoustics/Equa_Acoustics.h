/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2024 Rachid Touzani

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

     Definition of abstract class 'Equa_Wave' to solve acoustic wave equations

  ==============================================================================*/


#ifndef __EQUA_ACOUSTICS_H
#define __EQUA_ACOUSTICS_H

#include "equations/Equation.h"
#include "solvers/TimeStepping.h"
#include "mesh/Material.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \defgroup Acoustics Acoustic wave propagation problems
 *  \brief Acoustic wave propagation problems
 */

/*! \file Equa_Acoustics.h
 *  \brief Definition file for class Equa_Acoustics.
 */

class Element;
class Side;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
extern Material theMaterial;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! \class Equa_Wave
 *  \ingroup Wave
 * \brief Abstract class for Wave Finite %Element classes.
 *
 * \tparam <NEN> Number of element nodes
 * \tparam <NEE_> Number of element equations
 * \tparam <NSN_> Number of side nodes
 * \tparam <NSE_> Number of side equations
 */

template<size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
class Equa_Acoustics : virtual public Equation<NEN_,NEE_,NSN_,NSE_>
{

 public:

   using Equa::setInput;
   using Equa::setTerms;
   using Equa::_u;
   using Equation<NEN_,NEE_,NSN_,NSE_>::setMaterialProperty;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_theMesh;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_theElement;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_theSide;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_terms;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_analysis;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_TimeInt;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_eu;
   using Equation<NEN_,NEE_,NSN_,NSE_>::eA0;
   using Equation<NEN_,NEE_,NSN_,NSE_>::eA1;
   using Equation<NEN_,NEE_,NSN_,NSE_>::eA2;
   using Equation<NEN_,NEE_,NSN_,NSE_>::sA0;
   using Equation<NEN_,NEE_,NSN_,NSE_>::eMat;
   using Equation<NEN_,NEE_,NSN_,NSE_>::eRHS;
   using Equation<NEN_,NEE_,NSN_,NSE_>::sRHS;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_nb_nodes;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_nb_sides;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_nb_el;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_nb_eq;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_nb_dof_total;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_nb_dof;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_el_geo;

/// \brief Default constructor.
/// \details Constructs an empty equation.
    Equa_Acoustics()
    {
       _terms = DIFFUSION;
       _analysis = TRANSIENT;
       _TimeInt.scheme = NONE;
    }

/// \brief Destructor
    virtual ~Equa_Acoustics() { }

/// \brief Add lumped mass contribution to element matrix
/// @param [in] coef coefficient to multiply by the matrix before adding [Default: <tt>1</tt>]
    virtual void LMass(real_t coef=1) { _coef = coef; }

/// \brief Add consistent mass contribution to left-hand side
/// @param [in] coef coefficient to multiply by the matrix before adding [Default: <tt>1</tt>]
    virtual void Mass(real_t coef=1) { _coef = coef; }

/// \brief Add diffusion term to element matrix
    virtual void Diffusion(real_t coef=1.) { _coef = coef; }

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
 *         and a lumped or consistent mass matrix. If transient analysis is chosen, the lumped
 *         mass matrix option is chosen by default, and the implicit Euler scheme is used
 *         by default for time integration.
 *  </ul>
 */
    void build()
    {
       static bool matrix_set = false;
       real_t dti = 1./Equa::_TimeInt.delta;
       if (Equa::_A==nullptr && !matrix_set) {
          Equa::setMatrixType(SPARSE);
          Equa::setSolver(CG_SOLVER,DILU_PREC);
          matrix_set = true;
       }
       Equa::_A->clear();
       _TimeInt.theta = 1;
       if (_TimeInt.scheme==FORWARD_EULER)
          _TimeInt.theta = 0;
       else if (_TimeInt.scheme==CRANK_NICOLSON)
          _TimeInt.theta = 0.5;

       MESH_EL {
          set(the_element);
          if (_terms&MASS)
             Mass();
          if (_terms&LUMPED_MASS)
             LMass();
          if (_terms&DIFFUSION)
             Diffusion();
          if (_TimeInt.scheme==BACKWARD_EULER)
             eMat = dti*eA1 + eA0;
          else if (_TimeInt.scheme==FORWARD_EULER)
             eMat = dti*eA1;
          else if (_TimeInt.scheme==CRANK_NICOLSON)
             eMat = dti*eA1 + 0.5*eA0;
          Equa::_A->Assembly(The_element,eMat.get());
          if (Equa::_bf!=nullptr)
             BodyRHS(*Equa::_bf);
          if (Equa::_bc!=nullptr)
             this->updateBC(The_element,*Equa::_bc);
          if (_TimeInt.scheme==BACKWARD_EULER)
             eRHS += dti*eA1*_eu;
          else if (_TimeInt.scheme==FORWARD_EULER)
             eRHS += (dti*eA1-eA0)*_eu;
          else if (_TimeInt.scheme==CRANK_NICOLSON)
             eRHS += (dti*eA1-0.5*eA0)*_eu;
          Equa::_b->Assembly(The_element,eRHS.get());
       }
       if (Equa::_sf!=nullptr) {
          MESH_BD_SD {
             set(the_side);
             BoundaryRHS(*Equa::_sf);
             Equa::_A->Assembly(The_side,sA0.get());
             Equa::_b->Assembly(The_side,sRHS.get());
          }
       }
    }

/** \brief Build the linear system of equations
 *  \details Before using this function, one must have properly selected 
 *  appropriate options for:
 *  <ul>
 *     <li>The choice of a steady state or transient analysis. By default, the analysis is stationary
 *     <li>In the case of transient analysis, the choice of a time integration scheme
 *         and a lumped or consistent mass matrix. If transient analysis is chosen, the lumped
 *         mass matrix option is chosen by default, and the implicit Euler scheme is used
 *         by default for time integration.
 *  </ul>
 *  @param [in] s Reference to used TimeStepping instance
 */
    void build(TimeStepping& s)
    {
       MESH_EL {
          set(the_element);
          if (_terms&MASS)
             Mass();
          if (_terms&LUMPED_MASS)
             LMass();
          if (_terms&DIFFUSION)
             Diffusion();
          if (_terms&SOURCE && Equa::_bf!=nullptr)
             BodyRHS(*Equa::_bf);
          s.Assembly(The_element,eRHS.get(),eA0.get(),eA1.get());
       }
       if (Equa::_sf!=nullptr) {
          MESH_SD {
             if (The_side.isReferenced()) {
                set(the_side);
                if (_terms&FLUX && Equa::_sf!=nullptr)
                   BoundaryRHS(*Equa::_sf);
                s.SAssembly(The_side,sRHS.get());
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
          set(the_element);
          this->ElementVector(*Equa::_u);
          if (_terms&MASS)
             Mass(1.);
          if (_terms&LUMPED_MASS)
             Mass(1.);
          Diffusion();
          e.Assembly(*_theElement,eA0.get(),eA1.get());
       }
    }

/// \brief Set speed of sound (constant)
    void setSpeed(const real_t& c) { _speed = c; }

 protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   virtual void set(const Element *el)=0;
   virtual void set(const Side *sd)=0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    real_t  _speed, _coef;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
