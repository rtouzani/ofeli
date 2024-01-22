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

                 Definition of abstract class 'Equa_LinearPDE' 

  ==============================================================================*/


#ifndef __EQUA_LINEAR_PDE_H
#define __EQUA_LINEAR_PDE_H

#include "equations/Equation.h"
#include "solvers/TimeStepping.h"
#include "io/Fct.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \defgroup Generic
 *  \brief Generic partial differential equation
 */

/*! \file Equa_LinearPDE.h
 *  \brief Definition file for class Equa_LinearPDE.
 */


class Element;
class Side;

/*! \class Equa_LinearPDE
 *  \ingroup Generic
 * \brief Abstract class for Finite %Element classes for lienar PDEs'.
 *
 * \tparam <NEN> Number of element nodes
 * \tparam <NSN_> Number of side nodes
 */

template<size_t NEN_, size_t NSN_>
class Equa_LinearPDE : virtual public Equation<NEN_,NEN_,NSN_,NSN_>
{

 public:

   using Equa::setInput;
   using Equa::setTerms;
   using Equa::_u;
   using Equa::_terms;
   using Equation<NEN_,NEN_,NSN_,NSN_>::_theMesh;
   using Equation<NEN_,NEN_,NSN_,NSN_>::_theElement;
   using Equation<NEN_,NEN_,NSN_,NSN_>::_theSide;
   using Equation<NEN_,NEN_,NSN_,NSN_>::_analysis;
   using Equation<NEN_,NEN_,NSN_,NSN_>::_TimeInt;
   using Equation<NEN_,NEN_,NSN_,NSN_>::_eu;
   using Equation<NEN_,NEN_,NSN_,NSN_>::eA0;
   using Equation<NEN_,NEN_,NSN_,NSN_>::eA1;
   using Equation<NEN_,NEN_,NSN_,NSN_>::eA2;
   using Equation<NEN_,NEN_,NSN_,NSN_>::sA0;
   using Equation<NEN_,NEN_,NSN_,NSN_>::eRHS;
   using Equation<NEN_,NEN_,NSN_,NSN_>::sRHS;
   using Equation<NEN_,NEN_,NSN_,NSN_>::_nb_nodes;
   using Equation<NEN_,NEN_,NSN_,NSN_>::_nb_sides;
   using Equation<NEN_,NEN_,NSN_,NSN_>::_nb_el;
   using Equation<NEN_,NEN_,NSN_,NSN_>::_nb_eq;
   using Equation<NEN_,NEN_,NSN_,NSN_>::_nb_dof_total;
   using Equation<NEN_,NEN_,NSN_,NSN_>::_nb_dof;
   using Equation<NEN_,NEN_,NSN_,NSN_>::_el_geo;
   using Equation<NEN_,NEN_,NSN_,NSN_>::setPDECoef;

/// \brief Default constructor.
/// \details Constructs an empty equation.
    Equa_LinearPDE()
    {
       _analysis = STEADY_STATE;
       _TimeInt.scheme = NONE;
       _type00 = _type10 = _type01 = _type20 = _type02 = 1;
       _c00 = _c10 = _c01 = _c20 = _c02 = 1.0;
       _terms = int(PDE_Terms::NOTERM);
    }

/// \brief Destructor
    virtual ~Equa_LinearPDE() { }

/// \brief Set no lumping
/// @remark Default is lumping
    void setNoLumping() { _lump = false; }

/// \brief Set stabilization for convection term
/// @remark Default is no stabilization
    void setStab() { _stab = true; }

/// \brief Set coefficient for term (0,0): 0th order in time and space
/// @param [in] a Constant coefficient to multiply by 0-th order term [Default: <tt>1.</tt>]
    void set_00(real_t a=1.0) { _c00 = a; _type00 = 1; _terms = (_terms|int(PDE_Terms::L00)); }

/// \brief Set coefficient for term (0,0): 0th order in time and space
/// @param [in] f Function to multiply by 0-th order term (Function of \c x and \c t )
    void set_00(Fct& f) { _f00 = f; _type00 = 3; _terms = (_terms|int(PDE_Terms::L00)); }

/// \brief Set coefficient for term (0,0): 0th order in time and space
/// @param [in] f Function to multiply by 0-th order term (Function of \c x and \c t)
    void set_00(const string& f) { _f00.set(f,_var); _type00 = 3; _terms = (_terms|int(PDE_Terms::L00)); }

/// \brief Set coefficient for term (1,0): 1st order in time, 0th order in space
/// @param [in] a Constant coefficient to multiply by (1,0)-order term [Default: <tt>1.</tt>]
    void set_10(real_t a=1.0) { _c10 = a; _type10 = 1; _terms = (_terms|int(PDE_Terms::L10)); }

/// \brief Set coefficient for term (1,0): 1st order in time, 0th order in space
/// @param [in] f Function to multiply by (1,0)-order term (Function of \c x and \c t )
    void set_10(Fct& f) { _f10 = f; _type10 = 3; _terms = (_terms|int(PDE_Terms::L10)); }

/// \brief Set coefficient for term (2,0): 2nd order in time, 0th order in space
/// @param [in] a Constant coefficient to multiply by (2,0)-order term [Default: <tt>1.</tt>]
    void set_20(real_t a=1.0) { _c20 = a; _type20 = 1; _terms = (_terms|int(PDE_Terms::L20)); }

/// \brief Set coefficient for term (2,0): 2nd order in time, 0th order in space
/// @param [in] f Function to multiply by (2,0)-order term (Function of \c x and \c t)
    void set_20(Fct& f) { _f20 = f; _type20 = 3; _terms = (_terms|int(PDE_Terms::L20)); }

/// \brief Set coefficient for term (0,1): 0th order in time, 1st order in space
/// @param [in] a Constant coefficient to multiply by (0,1)-order term [Default: <tt>1.</tt>]
    void set_01(real_t a=1.0) { _c01 = a; _type01 = 1; _terms = (_terms|int(PDE_Terms::L01)); }

/// \brief Set coefficient for term (0,1): 0th order in time, 1st order in space
/// @param [in] a Constant coefficient to multiply by (0,1)-order term [Default: <tt>1.</tt>]
    void set_01(Point<real_t> &a) { _d01 = a; _type01 = 2; _terms = (_terms|int(PDE_Terms::L01)); }

/// \brief Set coefficient for term (0,1): 0th order in time, 1st order in time and space
/// @param [in] f Function to multiply by (0,1)-order term (Function of \c x and \c t)
    void set_01(Fct& f) { _f01 = f; _type01 = 3; _terms = (_terms|int(PDE_Terms::L01)); }

/// \brief Set coefficient for term (0,1): 0th order in time, 1st order in space
/// @param [in] a Constant coefficient to multiply by (0,1)-order term [Default: <tt>1.</tt>]
    void set_01(Vect<real_t> &a) { _e01 = a; _type01 = 4; _terms = (_terms|int(PDE_Terms::L01)); }

/// \brief Set coefficient for term (0,2): 0th order in time, 2nd order in space
/// @param [in] a Constant coefficient to multiply by (0,2)-order term [Default: <tt>1.</tt>]
    void set_02(real_t a=1.0) { _c02 = a; _type02 = 1; _terms = (_terms|int(PDE_Terms::L02)); }

/// \brief Set coefficient for term (0,2): 0th order in time, 2nd order in time and space
/// @param [in] f Function to multiply by (0,2)-order term (Function of \c x and \c t)
    void set_02(Fct& f) { _f02 = f; _type02 = 3; _terms = (_terms|int(PDE_Terms::L02)); }

/// \brief Add 0th order term, in time and space, to left-hand side
/// @param [in] coef coefficient to multiply by the matrix before adding [Default: <tt>1</tt>]
    virtual void Mat_00(real_t coef=1.0) { _coef = coef; }

/// \brief Add 1st order term in time, 0th in space to left-hand side
/// @param [in] coef coefficient to multiply by the matrix before adding [Default: <tt>1</tt>]
    virtual void Mat_10(real_t coef=1.0) { _coef = coef; }

/// \brief Add 2nd order term in time, 0th in space to left-hand side
/// @param [in] coef coefficient to multiply by the matrix before adding [Default: <tt>1</tt>]
    virtual void Mat_20(real_t coef=1.0) { _coef = coef; }

/// \brief Add 0th order term in time, 1st in space to left-hand side
/// @param [in] coef coefficient to multiply by the matrix before adding [Default: <tt>1</tt>]
    virtual void Mat_01(real_t coef=1.0) { _coef = coef; }

/// \brief Add 0th order term in time, 2nd in space to left-hand side
/// @param [in] coef coefficient to multiply by the matrix before adding [Default: <tt>1</tt>]
    virtual void Mat_02(real_t coef=1.0) { _coef = coef; }

/// \brief Add body right-hand side term to right-hand side.
/// @param [in] f Vector containing source at nodes.
    virtual void BodyRHS(const Vect<real_t>& f) { }

/// \brief Add boundary right-hand side term to right-hand side.
/// @param [in] f Vector containing source at nodes.
    virtual void BoundaryRHS(const Vect<real_t>& f) { }

/// \brief Build the linear system of equations for the steady state case
    void build()
    {
       static bool matrix_set = false;
       if (Equa::_A==nullptr && !matrix_set) {
          Equa::setMatrixType(SPARSE);
          if (_terms&int(PDE_Terms::L01))
             Equa::setSolver(BICG_STAB_SOLVER,DILU_PREC);
          else
             Equa::setSolver(CG_SOLVER,DILU_PREC);
          matrix_set = true;
       }
       Equa::_A->clear();
       MESH_EL {
          set(the_element);
          if (_terms&int(PDE_Terms::L00))
             Mat_00();
          if (_terms&int(PDE_Terms::L01))
             Mat_01();
          if (_terms&int(PDE_Terms::L02))
             Mat_02();
          Equa::_A->Assembly(The_element,eA0.get());
          if (Equa::_bf!=nullptr)
             BodyRHS(*Equa::_bf);
          if (Equa::_bc!=nullptr)
             this->updateBC(The_element,eA0.get(),*Equa::_bc);
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
 *         and a lumped or consistent capacity matrix. If transient analysis is chosen, the lumped
 *         capacity matrix option is chosen by default, and the implicit Euler scheme is used
 *         by default for time integration.
 *  </ul>
 *  @param [in] s Reference to used TimeStepping instance
 */
    void build(TimeStepping& s)
    {
       MESH_EL {
          set(the_element);
          if (_terms&int(PDE_Terms::L10))
             Mat_10();
          if (_terms&int(PDE_Terms::L20))
             Mat_20();
          if (_terms&int(PDE_Terms::L00))
             Mat_00();
          if (_terms&int(PDE_Terms::L01))
             Mat_01();
          if (_terms&int(PDE_Terms::L02))
             Mat_02();
          if (_terms&int(PDE_Terms::BODY_RHS))
             BodyRHS(*Equa::_bf);
          s.Assembly(The_element,eRHS.get(),eA0.get(),eA1.get(),eA2.get());
       }
       if (Equa::_sf!=nullptr) {
          MESH_SD {
             if (The_side.isReferenced()) {
                set(the_side);
                if (_terms&int(PDE_Terms::BOUNDARY_RHS))
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
          if (_terms&int(PDE_Terms::L10))
             Mat_10();
          if (_terms&int(PDE_Terms::L20))
             Mat_20();
          if (_terms&int(PDE_Terms::L02))
             Mat_02();
          if (_terms&int(PDE_Terms::L01))
             Mat_01();
          e.Assembly(*_theElement,eA0.get(),eA1.get());
       }
    }

 protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   virtual void set(const Element *el)=0;
   virtual void set(const Side *sd)=0;
   bool _lump, _stab;
   real_t _coef, _c00, _c10, _c01, _c20, _c02;
   Fct _f00, _f10, _f01, _f20, _f02;
   Point<real_t> _d01, _d02;
   Vect<real_t> _e01, _e02;
   int _type00, _type10, _type01, _type20, _type02;
   const vector<string> _var {"x","y","z","t"};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif