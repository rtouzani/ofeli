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

     Definition of abstract class 'Equa_Porous' to solve porous media problems

  ==============================================================================*/


#ifndef __EQUA_POROUS_H
#define __EQUA_POROUS_H

#include "mesh/Material.h"
#include "equations/Equation.h"
#include "solvers/TimeStepping.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \defgroup Porous Porous Media problems
 *  \brief Porous Media equation classes
 */

/*! \file Equa_Porous.h
 *  \brief Definition file for class Equa_Porous.
 */


#ifndef DOXYGEN_SHOULD_SKIP_THIS
extern Material theMaterial;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
 
class Element;
class Side;

/*! \class Equa_Porous
 *  \ingroup Porous
 * \brief Abstract class for Porous Media Finite %Element classes.
 *
 * \tparam <T_> data type (real_t, float, ...)
 * \tparam <NEN> Number of element nodes
 * \tparam <NEE_> Number of element equations
 * \tparam <NSN_> Number of side nodes
 * \tparam <NSE_> Number of side equations
 */

template<size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
class Equa_Porous : virtual public Equation<NEN_,NEE_,NSN_,NSE_>
{

 public:

   using Equation<NEN_,NEE_,NSN_,NSE_>::setMaterialProperty;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_theMesh;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_theElement;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_theSide;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_terms;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_analysis;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_TimeInt;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_A;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_b;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_uu;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_u;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_bc;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_bf;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_sf;
   using Equation<NEN_,NEE_,NSN_,NSE_>::eA0;
   using Equation<NEN_,NEE_,NSN_,NSE_>::eA1;
   using Equation<NEN_,NEE_,NSN_,NSE_>::eRHS;
   using Equation<NEN_,NEE_,NSN_,NSE_>::sRHS;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_nb_nodes;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_nb_el;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_nb_eq;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_nb_dof_total;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_nb_dof;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_el_geo;

/// \brief Default constructor.
/// \details Constructs an empty equation.
    Equa_Porous()
    {
       _terms = int(PDE_Terms::MASS)|int(PDE_Terms::MOBILITY);
       _analysis = STEADY_STATE;
       _TimeInt.scheme = NONE;
       _phi_is_set = _permeability_is_set = false;
    }

/// \brief Destructor
    virtual ~Equa_Porous() { }

/// \brief Add mobility term to the 0-th order element matrix
    virtual void Mobility() { }

/// \brief Add porosity term to the 1-st order element matrix
    virtual void Mass() { }

/** \brief Add source right-hand side term to right-hand side.
 *  @param [in] bf Vector containing source at nodes.
 */
    virtual void BodyRHS(const Vect<real_t>& bf) { }

/** \brief Add boundary right-hand side term to right-hand side.
 *  @param [in] sf Vector containing source at nodes.
 */
    virtual void BoundaryRHS(const Vect<real_t>& sf) { }

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
    }

/** \brief Build the linear system of equations
 *  \details Before using this function, one must have properly selected 
 *  appropriate options for:
 *  <ul>
 *     <li>The choice of a steady state or transient analysis. By default, the analysis is stationary
 *     <li>In the case of transient analysis, the choice of a time integration scheme.
 *         If transient analysis is chosen, the implicit Euler scheme is used
 *         by default for time integration.
 *  </ul>
 *  @param [in] s Reference to used TimeStepping instance
 */
    void build(TimeStepping& s)
    {
       MESH_EL {
          set(the_element);
          if (_terms&int(PDE_Terms::MASS))
             Mass();
          if (_terms&int(PDE_Terms::MOBILITY))
             Mobility();
          if (_terms&int(PDE_Terms::SOURCE) && _bf!=nullptr)
             BodyRHS(*_bf);
          s.Assembly(The_element,eRHS.get(),eA0.get(),eA1.get());
       }
       MESH_SD {
          if (The_side.isReferenced()) {
             set(the_side);
             if (_terms&int(PDE_Terms::FLUX) && _bf!=nullptr)
                BoundaryRHS(*_bf);
             s.SAssembly(The_side,sRHS.get());
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
          if (_terms&int(PDE_Terms::MASS))
             Mass();
          Mobility();
          e.Assembly(The_element,eA0.get(),eA1.get());
       }
    }

/** \brief Run the equation
 *  \details If the analysis (see function setAnalysis) is <tt>STEADY_STATE</tt>, then
 *  the function solves the stationary equation.\n
 *  If the analysis is <tt>TRANSIENT</tt>, then the function performs time stepping
 *  until the final time is reached.
 */
    int run()
    {
       return 0;
    }

/// \brief Set viscosity given by an algebraic expression
    void Mu(const string& exp) { _mu = setMaterialProperty(exp.c_str(),"Viscosity"); }

 protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   virtual void set(const Element *el)=0;
   virtual void set(const Side *sd)=0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Set material properties
    void setMaterial()
    {
       theMaterial.setCode(_theElement->getCode());
       _mu = theMaterial.Viscosity();
    }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   bool         _density_is_set, _phi_is_set,  _mu_is_set, _permeability_is_set;
   real_t       _density, _mu, _cw, _Mw;
   Vect<real_t> _phi, _Kx, _Ky, _Kz;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
