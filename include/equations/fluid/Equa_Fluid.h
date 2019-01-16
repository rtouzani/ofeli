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

      Definition of abstract class 'Equa_Fluid' for Finite Element Equations
                         to solve fluid flow problems

  ==============================================================================*/


#ifndef __EQUA_FLUID_H
#define __EQUA_FLUID_H

#include "equations/Equation.h"

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

template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_> class Equa_Fluid;

template<class T_=real_t, size_t NEN_=3, size_t NEE_=3, size_t NSN_=2, size_t NSE_=2>
class Equa_Fluid : virtual public Equation<T_,NEN_,NEE_,NSN_,NSE_> {

public:

   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::setMaterialProperty;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theMesh;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theElement;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theSide;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_terms;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_analysis;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_time_scheme;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_verbose;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_time;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_final_time;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_time_step;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_A;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_b;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_bc;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_bf;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_uu;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_u;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_sf;


/// \brief Default constructor.
/// \details Constructs an empty equation.
    Equa_Fluid()
    {
       _terms = 0;
    }

/// \brief Destructor
    virtual ~Equa_Fluid() { }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/*
/// \brief Run one time step
    int runOneTimeStep()
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
    }*/
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

protected:

/// \brief Set (constant) Viscosity
    void Viscosity(const real_t& visc) { _visc = visc; }

/// \brief Set viscosity given by an algebraic expression
    void Viscosity(const string& exp) { setMaterialProperty(exp.c_str(),"Viscosity"); }

/// \brief Set (constant) Viscosity
    void Density(const real_t &dens) { _dens = dens; }

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
       _visc = theMaterial.Viscosity();
       _dens = theMaterial.Density();
    }
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    real_t _visc, _dens, _beta;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
