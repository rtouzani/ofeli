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

                Definition of abstract class 'Equa_Solid' to
                to solve Solid and Structural Mechanics problems

  ==============================================================================*/


#ifndef __EQUA_SOLID_H
#define __EQUA_SOLID_H

#include "solvers/TimeStepping.h"
#include "mesh/Material.h"
#include "equations/AbsEqua_impl.h"
#include "equations/Equation_impl.h"


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

class Element;
class Side;
extern Material theMaterial;


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
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_TimeInt;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eA0;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eA1;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eA2;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::sA0;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eMat;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::sMat;
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
    Equa_Solid() : _cd(nullptr) { }

/// \brief Destructor
    virtual ~Equa_Solid() { }

/// \brief Add lumped mass contribution to left-hand side
/// @param [in] coef coefficient to multiply by the matrix before adding [Default: <tt>1</tt>]
    virtual void LMass(real_t coef=1) { _coef = coef; }

/// \brief Add consistent mass contribution to left-hand side
/// @param [in] coef coefficient to multiply by the matrix before adding [Default: <tt>1</tt>]
    virtual void Mass(real_t coef=1) { _coef = coef; }

/// \brief Add deviator matrix to left-hand side taking into account time integration scheme,
/// after multiplication by <tt>coef</tt> [Default: <tt>1</tt>]
    virtual void Deviator(real_t coef=1) { coef=1; }

/// \brief Add dilatation matrix to left-hand side taking into account time integration scheme,
/// after multiplication by <tt>coef</tt> [Default: <tt>1</tt>]
    virtual void Dilatation(real_t coef=1) { coef=1; }

/// \brief Add stiffness matrix to left-hand side taking into account time integration scheme,
/// after multiplication by <tt>coef</tt> [Default: <tt>1</tt>]
    virtual void Stiffness(real_t coef=1) { coef=1; }

/// \brief Set specific input data to solid mechanics
    void setInput(EqDataType    opt,
                  Vect<real_t>& u)
    {
       AbsEqua<real_t>::setInput(opt,u);
       if (opt==CONTACT_DISTANCE)
          _cd = &u;
    }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    virtual void BodyRHS(const Vect<T_>& f) { }
    virtual void BoundaryRHS(const Vect<T_>& f) { }
    virtual void Periodic(real_t coef=1) { }
    virtual int Contact(real_t coef=1.e07) { return 0; }

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
       AbsEqua<T_>::_A->clear();
       mesh_elements(*_theMesh) {
          set(the_element);
          if (_terms&MASS)
             Mass();
          if (_terms&LUMPED_MASS)
             LMass();
          if (_terms&DEVIATORIC)
             Deviator();
          if (_terms&DILATATION)
             Dilatation();
          eMat = eA0;
          AbsEqua<T_>::_A->Assembly(*_theElement,eMat.get());
          if (AbsEqua<T_>::_bf!=nullptr)
             BodyRHS(*AbsEqua<T_>::_bf);
          if (AbsEqua<T_>::_bc!=nullptr)
             this->updateBC(*_theElement,*AbsEqua<T_>::_bc);
          AbsEqua<T_>::_b->Assembly(*_theElement,eRHS.get());
       }
       if (AbsEqua<T_>::_sf!=nullptr) {
          mesh_sides(*_theMesh) {
             set(the_side);
             if (_terms&CONTACT)
                Contact(1.e07);
             BoundaryRHS();
             AbsEqua<T_>::_A->Assembly(The_side,sMat.get());
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
       mesh_elements(*_theMesh) {
          set(the_element);
          this->ElementVector(*AbsEqua<T_>::_u);
          if (_terms&MASS)
             Mass();
          if (_terms&LUMPED_MASS)
             LMass();
          if (_terms&DEVIATORIC)
             Deviator();
          if (_terms&DILATATION)
             Dilatation();
          if (_terms&CONTACT)
             Contact(1.e07);
          if ((_terms&LOAD) && (AbsEqua<T_>::_bf!=nullptr))
             BodyRHS();
          s.Assembly(The_element,eRHS.get(),eA0.get(),eA1.get(),eA2.get());
       }
       mesh_sides(*_theMesh) {
          set(the_side);
          this->SideVector(*AbsEqua<T_>::_u);
          if (_terms&CONTACT)
             Contact(1.e07);
          s.SAssembly(The_side,sRHS.get(),sA0.get());
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
          if (_terms&MASS)
             Mass();
          if (_terms&LUMPED_MASS)
             LMass();
          Deviator();
          Dilatation();
          e.Assembly(The_element,eA0.get(),eA2.get());
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

/** \brief Run one time step
 *  \details This function performs one time step in equation solving.
 *  It is to be used only if a \a TRANSIENT analysis is required.
 *  @return Return error from the linear system solver
 */
    int runTransient()
    {
       *AbsEqua<T_>::_b = 0;
       build();
       int ret=AbsEqua<T_>::solveLinearSystem(*AbsEqua<T_>::_b,AbsEqua<T_>::_uu);
       AbsEqua<T_>::_u->insertBC(*_theMesh,AbsEqua<T_>::_uu,*AbsEqua<T_>::_bc);
       return ret;
    }
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
    virtual void set(const Element* el) = 0;
    virtual void set(const Side* sd) = 0;
    void BodyRHS() { BodyRHS(*AbsEqua<T_>::_bf); }
    void BoundaryRHS() { BoundaryRHS(*AbsEqua<T_>::_sf); }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    real_t _E, _nu, _lambda, _G, _rho, _coef;
    Vect<real_t> *_cd;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
