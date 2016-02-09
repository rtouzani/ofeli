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

      Definition of abstract class 'Equa_Laplace' for the Laplace equation

  ==============================================================================*/


#ifndef __EQUA_LAPLACE_H
#define __EQUA_LAPLACE_H

#include "equations/Equation.h"

namespace OFELI {

/*! \defgroup Laplace Laplace equation
 *  \brief Gathers equation classes for the Laplace equation
 */

/*! \file Equa_Laplace.h
 *  \brief Definition file for class Equa_Laplace.
 */

/*! \class Equa_Laplace
 *  \ingroup Laplace
 * \brief Abstract class for classes about the Laplace equation.
 *
 * \b Template \b Arguments:
 *
 * \arg \b T_   : data type (double, float, ...)
 * \arg \b NEN_ : Number of element nodes
 * \arg \b NEE_ : Number of element equations
 * \arg \b NSN_ : Number of side nodes
 * \arg \b NSE_ : Number of side equations
 *
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
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_uu;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_u;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::_bc;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eA0;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eA1;
   using Equation<T_,NEN_,NEE_,NSN_,NSE_>::eRHS;

/// \brief Default constructor.
/// \details Constructs an empty equation.
    Equa_Laplace()
    {
       _terms = 0;
    }

/// \brief Destructor
    virtual ~Equa_Laplace() {  }

/// \brief Solve the equation
/*    int run()
    {
       int ret = 0;
       _b = new Vect<real_t>(_theMesh->getNbEq());
       build();
       ret = this->SolveLinearSystem(_A,*_b,_uu);
       _u->insertBC(*_theMesh,_uu,*_bc);
       delete _b;
       return ret;
       }*/

    virtual void build() { }

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
    void set(const Element *el) { }
};

} /* namespace OFELI */

#endif
