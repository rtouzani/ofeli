/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

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

      Definition of abstract class 'FE_Electromagnetics' for Finite Element
                  Equations to solve electromagnetic problems

  ==============================================================================*/


#ifndef __FE_ELECTROMAGNETICS_H
#define __FE_ELECTROMAGNETICS_H

#include "solvers/TimeStepping.h"
#include "mesh/Material.h"

#include "equations/Equa_impl.h"
#include "equations/Equation_impl.h"
using std::complex;
using std::abs;

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \defgroup Electromagnetics Electromagnetics
 *  \brief Electromagnetic equations
 */

/*! \file Equa_Electromagnetics.h
 *  \brief Definition file for class FE_Electromagnetics.
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
extern Material theMaterial;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! \class Equa_Electromagnetics
 *  \ingroup Electromagnetics
 * \brief Abstract class for Electromagnetics Equation classes.
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

template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_> class Equa_Electromagnetics;

template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
class Equa_Electromagnetics : virtual public Equation<T_,NEN_,NEE_,NSN_,NSE_>
{

 public:

   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::setMaterialProperty;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theMesh;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theElement;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theSide;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_terms;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_analysis;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_TimeInt;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_nb_nodes;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_nb_sides;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_nb_el;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_nb_eq;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_nb_dof_total;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_nb_dof;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_el_geo;
   using Equa<T_>::_Mu_set;
   using Equa<T_>::_sigma_set;
   using Equa<T_>::_omega_set;


#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void build() { }
    int runOneTimeStep() { return 0; }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    int runTransient()
    {
       Equa<T_>::_b->clear();
       build();
       int ret=Equa<T_>::solveLinearSystem(*Equa<T_>::_b,Equa<T_>::_uu);
       Equa<T_>::_u->insertBC(*_theMesh,Equa<T_>::_uu,*Equa<T_>::_bc);
       return ret;
    }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 protected:

/// \brief Set (constant) magnetic permeability
    void MagneticPermeability(const real_t& mu) { _Mu = mu; }

/// \brief Set magnetic permeability given by an algebraic expression
    void MagneticPermeability(const string& exp) { setMaterialProperty(exp.c_str(),"Magnetic Permeability"); }

/// \brief Set (constant) electric conductivity
    void ElectricConductivity(const real_t& sigma) { _sigma = sigma; }

/// \brief set electric conductivity given by an algebraic expression
    void ElectricConductivity(const string& exp) { setMaterialProperty(exp.c_str(),"Electric Conductivity"); }

/// \brief Set material properties
    void setMaterial()
    {
       theMaterial.setCode(_theElement->getCode());
       _sigma = theMaterial.ElectricConductivity();
       _Mu    = theMaterial.MagneticPermeability();
    }
#ifndef DOXYGEN_SHOULD_SKIP_THIS
   Point<real_t> _x[NEN_];
   real_t        _Mu, _sigma, _omega, _body_source, _bound_source;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
