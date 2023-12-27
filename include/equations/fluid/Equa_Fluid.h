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

      Definition of abstract class 'Equa_Fluid' for Finite Element Equations
                         to solve fluid flow problems

  ==============================================================================*/


#ifndef __EQUA_FLUID_H
#define __EQUA_FLUID_H

#include "solvers/TimeStepping.h"
#include "mesh/Material.h"
#include "equations/Equation_impl.h"


namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \defgroup Fluid Fluid Dynamics
 *  \brief Fluid Dynamics equations
 */

/*! \file Equa_Fluid.h
 *  \brief Definition file for class Equa_Fluid.
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
extern Material theMaterial;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! \class Equa_Fluid
 *  \ingroup Fluid
 * \brief Abstract class for Fluid Dynamics Equation classes.
 *
 * \tparam <T_> data type (double, float, ...)
 * \tparam <NEN> Number of element nodes
 * \tparam <NEE_> Number of element equations
 * \tparam <NSN_> Number of side nodes
 * \tparam <NSE_> Number of side equations
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

template<size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
class Equa_Fluid : virtual public Equation<NEN_,NEE_,NSN_,NSE_> {

 public:

   using Equation<NEN_,NEE_,NSN_,NSE_>::setMaterialProperty;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_theMesh;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_theElement;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_theSide;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_terms;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_analysis;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_TimeInt;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_nb_nodes;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_nb_sides;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_nb_el;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_nb_eq;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_nb_dof_total;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_nb_dof;
   using Equation<NEN_,NEE_,NSN_,NSE_>::_el_geo;
   using Equation<NEN_,NEE_,NSN_,NSE_>::updateBC;
   using Equa::_rho_set;
   using Equa::_mu_set;
   using Equa::_beta_set;


/// \brief Default constructor.
/// \details Constructs an empty equation.
    Equa_Fluid() : Equation<NEN_,NEE_,NSN_,NSE_>()
    {
       _terms = 0;
    }

/// \brief Destructor
    virtual ~Equa_Fluid() { }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void build() { }
    int runOneTimeStep() { return 0; }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Set Reynolds number
    void Reynolds(const real_t& Re) { _Re = Re; }

/// \brief Set (constant) Viscosity
    void Viscosity(const real_t& visc) { _mu = visc; }

/// \brief Set viscosity given by an algebraic expression
    void Viscosity(const string& exp) { setMaterialProperty(exp.c_str(),"Viscosity"); }

/// \brief Set (constant) Viscosity
    void Density(const real_t &dens) { _rho = dens; }

/// \brief Set Density given by an algebraic expression
    void Density(const string& exp) { setMaterialProperty(exp.c_str(),"Density"); }

/// \brief Set (constant) thermal expansion coefficient
    void ThermalExpansion(const real_t *e) { _beta = e; }

/// \brief Set thermal expansion coefficient given by an algebraic expression
    void ThermalExpansion(const string& exp) { setMaterialProperty(exp.c_str(),"ThermalExpansion"); }

/// \brief Set material properties
    void setMaterial()
    {
       theMaterial.setCode(_theElement->getCode());
       _mu = theMaterial.Viscosity();
       _rho = theMaterial.Density();
    }
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    real_t _mu, _rho, _beta, _Re, _cfl;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
