/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2018 Rachid Touzani

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

#include "equations/Equation.h"
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
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theElement;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theSide;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_analysis;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_time_scheme;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_time_step;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_final_time;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_A;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_b;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_uu;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_u;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_bc;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_dSh;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Run one time step
/*    int runOneTimeStep()
    {
       build();
       int ret = SolveLinearSystem(_A,*_b,_uu);
       _u->insertBC(*_theMesh,_uu,*_bc);
       return ret;
       }*/

/** \brief Solve the equation
 *  \details If the analysis (see function setAnalysis) is \a STEADY_STATE, then
 *  the function solves the stationary equation.\n
 *  If the analysis is \a TRANSIENT, then the function performs time stepping
 *  until the final time is reached.
 */
/*    int run()
    {
       int ret=0;
       _b = new Vect<T_>(_theMesh->getNbEq());
      if (_analysis==STEADY_STATE) {
          build();
          ret = SolveLinearSystem(_A,*_b,_uu);
          _u->insertBC(*_theMesh,_uu,*_bc);
       }
       else {
          for (real_t t=0.; t<=_final_time; t+=_time_step)
             _time_step = runOneTimeStep();
       }
       delete _b;
       return ret;
    }
*/
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 protected:

/// \brief Set (constant) magnetic permeability
    void MagneticPermeability(const real_t& mu) { _mu = mu; }

/// \brief Set magnetic permeability given by an algebraic expression
    void MagneticPermeability(const string& exp) { setMaterialProperty(exp.c_str(),"Magnetic Permeability"); }

/// \brief Set (constant) electric conductivity
    void ElectricConductivity(const real_t& sigma) { _sigma = sigma; }

/// \brief set electric conductivity given by an algebraic expression
    void ElectricConductivity(const string& exp) { setMaterialProperty(exp.c_str(),"Electric Conductivity"); }

/// \brief Set (constant) electric resistivity
    void ElectricResistivity(const real_t& rho) { _rho = rho; }

/// \brief Set electric resistivity given by an algebraic expression
    void ElectricResistivity(const string& exp) { setMaterialProperty(exp.c_str(),"Electric Resistivity"); }

/// \brief Set material properties
    void setMaterial()
    {
       theMaterial.setCode(_theElement->getCode());
       _sigma = theMaterial.ElectricConductivity();
       _rho   = theMaterial.ElectricResistivity();
       _mu    = theMaterial.MagneticPermeability();
    }
#ifndef DOXYGEN_SHOULD_SKIP_THIS
   Point<real_t> _x[NEN_];
   real_t        _mu, _sigma, _rho, _body_source, _bound_source;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
