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

      Definition of abstract class 'Equa_Laplace' for the Laplace equation

  ==============================================================================*/


#ifndef __EQUA_LAPLACE_H
#define __EQUA_LAPLACE_H

#include "equations/Equation.h"

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

template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_> class Equa_Laplace;

template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
class Equa_Laplace : virtual public Equation<T_,NEN_,NEE_,NSN_,NSE_> {

 public:
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theMesh;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theElement;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_theSide;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_terms;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_A;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_b;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_bf;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_sf;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_uu;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_u;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_bc;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eA0;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eA1;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eRHS;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_bf_given;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_bc_given;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_sf_given;

/// \brief Default constructor.
/// \details Constructs an empty equation.
    Equa_Laplace()
    {
       _terms = 0;
    }

/// \brief Destructor
    virtual ~Equa_Laplace() {  }

/// \brief Add finite element matrix to left-hand side
/// @param [in] coef Value to multiply by the added matrix
    virtual void LHS(real_t coef=1.) { }

/// \brief Add body source term to right-hand side
/// @param [in] f Vector containing the source given function at mesh nodes
    virtual void BodyRHS(const Vect<real_t>& f) { }

/// \brief Add boundary source term to right-hand side
/// @param [in] h Vector containing the source given function at mesh nodes
    virtual void BoundaryRHS(const Vect<real_t>& h) { }

/// \brief Define Source right-hand side of the equation
/// @param f [in] Vector containing source values at nodes
    virtual void setSource(const Vect<real_t>& f) { }

/** \brief Build and solve the linear system of equations using an iterative method.
 *  \details The matrix is preconditioned by the diagonal ILU method.
 *  The linear system is solved either by the Conjugate Gradient method if the matrix is symmetric
 *  positive definite (<tt>eps=-1</tt>) or the GMRES method if not. The solution is stored in the vector
 *  <tt>u</tt> given in the constructor.
 *  @return Number of performed iterations. Note that the maximal number
 *  is 1000 and the tolerance is 1.e-8
 */
    int run()
    {
       build();
       int ret = this->solveLinearSystem(*_b,_uu);
       _u->insertBC(*_theMesh,_uu,*_bc);
       return ret;
    }


/** \brief Build global matrix and right-hand side.
 *  \details The problem matrix and right-hand side are the ones used in the constructor.
 *  They are updated in this member function.
 */
    void build()
    {
       *_A = 0;
       _b = new Vect<real_t>(_theMesh->getNbEq());
       MESH_EL {
          set(theElement);
          LHS();
          if (_bf_given)
             BodyRHS(*_bf);
          if (_bc_given)
             this->updateBC(*_bc);
          this->ElementAssembly(_A);
          this->ElementAssembly(*_b);
       }
       if (_sf_given) {
          MESH_SD {
             set(theSide);
             BoundaryRHS(*_sf);
             this->SideAssembly(*_b);
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
       MESH_EL {
          set(theElement);
          buildEigen();
          e.Assembly(TheElement,eA0.get(),eA1.get());
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
