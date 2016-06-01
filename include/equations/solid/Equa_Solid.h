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

                Definition of abstract class 'Equa_Solid' to
                to solve Solid and Structural Mechanicsr problems

  ==============================================================================*/


#ifndef __EQUA_SOLID_H
#define __EQUA_SOLID_H

#include "equations/Equation.h"
#include "solvers/TimeStepping.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \defgroup Solid Solid Mechanics
 *  \brief Solid Mechanics finite element equations
 */

/*! \file Equa_Solid.h
 *  \brief Definition file for class Equa_Solid.
 */

/*! \class Equa_Solid
 *  \ingroup Solid
 * \brief Abstract class for Solid Mechanics Finite %Element classes.
 *
 * \tparam <T_> data type (double, float, ...)
 * \tparam <NEN> Number of element nodes
 * \tparam <NEE_> Number of element equations
 * \tparam <NSN_> Number of side nodes
 * \tparam <NSE_> Number of side equations
 */


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_> class Equa_Solid;

template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
class Equa_Solid : virtual public Equation<T_,NEN_,NEE_,NSN_,NSE_>
{

 public:

   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::setMaterialProperty;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theMesh;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theElement;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theSide;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_terms;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_analysis;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_time_scheme;
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
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eA0;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eA1;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eA2;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eRHS;

/// \brief Default constructor.
/// \details Constructs an empty equation.
    Equa_Solid() { }

/// \brief Destructor
    virtual ~Equa_Solid() { }

/// \brief Add lumped mass contribution to left-hand side
/// @param [in] coef coefficient to multiply by the matrix before adding [Default: <tt>1</tt>]
    virtual void LMassToLHS(real_t coef=1) { _coef = coef; }

/// \brief Add lumped mass contribution to right-hand side
/// @param [in] coef coefficient to multiply by the vector before adding [Default: <tt>1</tt>]
    virtual void LMassToRHS(real_t coef=1) { _coef = coef; }

/// \brief Add consistent mass contribution to left-hand side
/// @param [in] coef coefficient to multiply by the matrix before adding [Default: <tt>1</tt>]
    virtual void MassToLHS(real_t coef=1) { _coef = coef; }

/// \brief Add consistent mass contribution to right-hand side
/// @param [in] coef coefficient to multiply by the vector before adding [Default: <tt>1</tt>]
    virtual void MassToRHS(real_t coef=1) { _coef = coef; }

/// \brief Add lumped mass contribution to left and right-hand sides taking into account 
/// time integration scheme.
    void setLumpedMass() 
    {
       LMassToLHS(1./_time_step);
       LMassToRHS(1./_time_step);
    }

/// \brief Add consistent mass contribution to left and right-hand sides taking into account 
/// time integration scheme.
    void setMass() 
    {
       MassToLHS(1./_time_step);
       MassToRHS(1./_time_step);
    }

/// \brief Add consistent mass matrix to left-hand side after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    virtual void Mass(real_t coef=1) { coef=1; }

/// \brief Add lumped mass matrix to left-hand side after multiplication by <tt>coef</tt>
/// [Default: <tt>1</tt>]
    virtual void LMass(real_t coef=1) { coef=1; }

/// \brief Add deviator matrix to left-hand side taking into account time integration scheme,
/// after multiplication by <tt>coef</tt> [Default: <tt>1</tt>]
    virtual void Deviator(real_t coef=1) { coef=1; }

/// \brief Add dilatation matrix to left-hand side taking into account time integration scheme,
/// after multiplication by <tt>coef</tt> [Default: <tt>1</tt>]
    virtual void Dilatation(real_t coef=1) { coef=1; }

/// \brief Add dilatation vector to right-hand side taking into account time integration scheme,
///  after multiplication by <tt>coef</tt> [Default: <tt>1</tt>]
    virtual void DilatationToRHS(real_t coef=1) { coef=1; }

/// \brief Add deviator vector to right-hand side taking into account time integration scheme,
/// after multiplication by <tt>coef</tt> [Default: <tt>1</tt>]
    virtual void DeviatorToRHS(real_t coef=1) { coef=1; }

/// \brief Add stiffness matrix to left-hand side taking into account time integration scheme,
/// after multiplication by <tt>coef</tt> [Default: <tt>1</tt>]
    virtual void Stiffness(real_t coef=1) { coef=1; }

/// \brief Add stiffness matrix to right-hand side taking into account time integration scheme,
/// after multiplication by <tt>coef</tt> [Default: <tt>1</tt>]
    virtual void StiffnessToRHS(real_t coef=1) { coef=1; }

/// \brief Add dilatation matrix to left and/or right-hand side taking into account time 
//  integration scheme
    void setDilatation()
    {
       if (_time_scheme==STEADY_STATE || 
           _time_scheme==FORWARD_EULER || 
           _time_scheme==BACKWARD_EULER)
          Dilatation(1.);
       else if (_time_scheme==CRANK_NICOLSON)
          Dilatation(0.5);
       if (_time_scheme==CRANK_NICOLSON)
          DilatationToRHS(0.5);
    }

/// \brief Add deviator matrix to left and/or right-hand side taking into account time
/// integration scheme.
    void setDeviator()
    {
       if (_time_scheme==STEADY_STATE || 
           _time_scheme==FORWARD_EULER || 
           _time_scheme==BACKWARD_EULER)
          Deviator(1.);
       else if (_time_scheme==CRANK_NICOLSON)
          Deviator(0.5);
       if (_time_scheme==CRANK_NICOLSON)
          DeviatorToRHS(0.5);
    }

/// \brief Add convection contribution to left and/or right-hand side taking into account 
/// time integration scheme.
    void setStiffness()
    {
       if (_time_scheme==STEADY_STATE || _time_scheme==BACKWARD_EULER)
          Stiffness(1.);
       else if (_time_scheme==FORWARD_EULER)
          StiffnessToRHS();
       else
          ;
    }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    virtual void BodyRHS(UserData<T_> &ud)  { }
    virtual void BodyRHS(const Vect<T_>& f, int opt=GLOBAL_ARRAY)  { }
    virtual void BoundaryRHS(const Vect<T_>& f) { }
    virtual void BoundaryRHS(UserData<T_>& ud) { }
    virtual void Periodic(real_t coef=1) { }
    int SignoriniContact(UserData<T_>& ud, real_t coef=1.e07) { return 0; }
    int SignoriniContact(Vect<T_>& f, real_t coef=1.e07) { return 0; }
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
       MESH_EL {
          set(theElement);
          this->ElementVector(_uu);
          if (_terms&MASS)
             setMass();
          if (_terms&LUMPED_MASS)
             setLumpedMass();
          if (_terms&DEVIATORIC)
             setDeviator();
          if (_terms&DILATATION)
             setDilatation();
          this->ElementAssembly(_A);
          if (_terms&SOURCE)
             BodyRHS(*_bf,GLOBAL_ARRAY);
          this->updateBC(*_bc);
          this->ElementAssembly(*_b);
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
          if (_terms&MASS)
             MassToLHS(1.);
          if (_terms&LUMPED_MASS)
             LMassToLHS(1.);
          if (_terms&DEVIATORIC)
             Deviator(1.);
          if (_terms&DILATATION)
             Dilatation(1.);
          if (_terms&LOAD && _bf)
             BodyRHS(*_bf,GLOBAL_ARRAY);
          s.Assembly(TheElement,eRHS.get(),eA0.get(),eA1.get(),eA2.get());
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
          if (_terms&MASS)
             MassToLHS();
          if (_terms&LUMPED_MASS)
             LMassToLHS();
          Deviator();
          Dilatation();
          e.Assembly(TheElement,eA0.get(),eA2.get());
       }
    }

/** \brief Build global stiffness and mass matrices for the eigen system
 *  \details Case where the mass matrix is lumped
 *  @param [in] K Stiffness matrix
 *  @param [in] M Vector containing diagonal mass matrix
 */
/*    void buildEigen(SkSMatrix<T_>& K,
                      Vect<T_>&      M)
    {
       MESH_EL {
          set(theElement);
          Deviator();
          Dilatation();
          this->ElementAssembly(K);
       }
       MESH_EL {
          set(theElement);
          LMassToLHS();
          M.Assembly(theElement,eRHS.get());
       }
       }*/
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 protected:

/// \brief Set (constant) Young modulus
    void Young(const real_t& E) { _E = E; }

/// \brief Set (constant) Poisson ratio
    void Poisson(const real_t& nu) { _nu = nu; }

/// \brief Set (constant) density
    void Density(const real_t& rho) { _rho = rho; }

/// \brief Set Young modulus given by an algebraic expression
    void Young(const string& exp) { _E = setMaterialProperty(exp.c_str(),"Young Modulus"); }

/// \brief Set Poisson ratio given by an algebraic expression
    void Poisson(const string& exp) { _nu = setMaterialProperty(exp.c_str(),"Poisson Coefficient"); }

/// \brief Set density given by an algebraic expression
    void Density(const string& exp) { _rho = setMaterialProperty(exp.c_str(),"Density"); }

/// \brief Set material properties
    void setMaterial()
    {
       theMaterial.setCode(_theElement->getCode());
       _rho = theMaterial.Density();
       _nu  = theMaterial.PoissonRatio();
       _E   = theMaterial.YoungModulus();
    }
    virtual void set(const Element* el) { }
    virtual void set(const Side* sd) { }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   real_t  _E, _nu, _lambda, _G, _rho, _coef;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
