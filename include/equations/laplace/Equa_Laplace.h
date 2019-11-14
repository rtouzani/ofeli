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

      Definition of abstract class 'Equa_Laplace' for the Laplace equation

  ==============================================================================*/


#ifndef __EQUA_LAPLACE_H
#define __EQUA_LAPLACE_H

#include "equations/Equation.h"
#include "util/macros.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \defgroup Laplace Laplace equation
 *  \brief Laplace and Poisson equations
 */

/*! \file Equa_Laplace.h
 *  \brief Definition file for class Equa_Laplace.
 */

/*! \class Equa_Laplace
 *  \ingroup Laplace
 * \brief Abstract class for classes about the Laplace equation.
 *
 * \tparam T_ Data type (real_t, float, complex<real_t>, ...)
 * \tparam NEN_ Number of element nodes
 * \tparam NEE_ Number of element equations
 * \tparam NSN_ Number of side nodes
 * \tparam NSE_ Number of side equations
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
class Equa_Laplace : virtual public Equation<T_,NEN_,NEE_,NSN_,NSE_> {

 public:
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theMesh;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theElement;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theSide;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_A;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_b;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_bf;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_sf;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_u;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_bc;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eMat;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eA0;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eA1;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eRHS;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::sRHS;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_nb_nodes;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_nb_el;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_nb_eq;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_el_geo;

/// \brief Default constructor.
/// \details Constructs an empty equation.
    Equa_Laplace() { }

/// \brief Destructor
    virtual ~Equa_Laplace() { }

/// \brief Add finite element matrix to left-hand side
    virtual void LHS() { }

/// \brief Add body source term to right-hand side
/// @param [in] f Vector containing the source given function at mesh nodes
    virtual void BodyRHS(const Vect<real_t>& f) { }

/// \brief Add boundary source term to right-hand side
/// @param [in] h Vector containing the source given function at mesh nodes
    virtual void BoundaryRHS(const Vect<real_t>& h) { }

/** \brief Build global matrix and right-hand side.
 *  \details The problem matrix and right-hand side are the ones used in the constructor.
 *  They are updated in this member function.
 */
    void build()
    {
      if (_u==nullptr)
         throw OFELIException("In Equa_Therm::build(): No solution vector given.");
       _A->clear();
       mesh_elements(*_theMesh) {
          set(the_element);
          LHS();
          if (_bf!=nullptr)
             BodyRHS(*_bf);
          if (_bc!=nullptr)
             this->updateBC(*_theElement,*_bc);
          AbsEqua<T_>::_A->Assembly(*_theElement,eMat.get());
          AbsEqua<T_>::_b->Assembly(*_theElement,eRHS.get());
       }
       if (_sf!=nullptr) {
          mesh_boundary_sides(*_theMesh) {
             set(the_side);
             if (_sf!=nullptr)
                BoundaryRHS(*_sf);
             AbsEqua<T_>::_b->Assembly(*_theSide,sRHS.get());
          }
       }
    }

/// \brief Build matrices for an eigenvalue problem
    virtual void buildEigen(int opt=0) { }

/** \brief Build the linear system for an eigenvalue problem
 *  @param [in] e Reference to used EigenProblemSolver instance
 */
    void build(EigenProblemSolver& e)
    {
       mesh_elements(*_theMesh) {
          set(the_element);
          buildEigen();
          e.Assembly(The_element,eA0.get(),eA1.get());
       }
    }

 protected:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
   virtual void set(const Element *el)=0;
   virtual void set(const Side *sd)=0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
