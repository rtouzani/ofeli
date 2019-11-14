/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

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

                       Functions for Finite Element Assembly

  ==============================================================================*/

/*! \file Assembly.h
 *  \brief A set of template functions for assembly purposes.
 *
 * \tparam <E_> Entity type to assemble (Element or Side classes)
 * \tparam <T_> Data type (double, float, complex<double>, ...)
 * \tparam <N_> Constant size of local vector
 */

#ifndef __ASSEMBLY_H
#define __ASSEMBLY_H

#include "OFELI_Config.h"
#include "mesh/Node.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/LocalVect.h"
#include "linear_algebra/LocalMatrix.h"
#include "linear_algebra/SkMatrix.h"
#include "linear_algebra/SkSMatrix.h"
#include "linear_algebra/SpMatrix.h"


namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/** \addtogroup VectMat
 *  @{
 */

/** \fn void element_assembly(const E_ &e, const LocalVect<T_,N_> &be, Vect<T_> &b)
 *  \ingroup Equation
 *  \brief Assemble local vector into global vector
 *
 *  @param [in] e Reference to local entity (Element or Side)
 *  @param [in] be Local vector
 *  @param [in,out] b Global vector
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_,class E_>
void element_assembly(const E_&               e,
                      const LocalVect<T_,N_>& be,
                      Vect<T_>&               b);


/** \fn void element_assembly(const E_ &e, const LocalMatrix<T_,N_,N_> &ae, Vect<T_> &b)
 *  \ingroup Equation
 *  \brief Assemble diagonal local vector into global vector
 *  @param [in] e Reference to local entity (Element or Side)
 *  @param [in] ae Local matrix
 *  @param [in,out] b Global vector
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_,class E_>
void element_assembly(const E_&                    e,
                      const LocalMatrix<T_,N_,N_>& ae,
                      Vect<T_>&                    b);


/** \fn void element_assembly(const E_ &e, const LocalMatrix<T_,N_,N_> &ae, Matrix<T_> *A)
 *  \ingroup Equation
 *  \brief Assemble local matrix into global matrix.
 *  \details This function is to be called with an abstract pointer to matrix (class Matrix)
 *  @param [in] e Reference to local entity (Element or Side)
 *  @param [in] ae Local matrix
 *  @param [in,out] A Pointer to global matrix
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_,class E_>
void element_assembly(const E_&                    e,
                      const LocalMatrix<T_,N_,N_>& ae,
                      Matrix<T_>*                  A);


/** \fn void element_assembly(const E_ &e, const LocalMatrix<T_,N_,N_> &ae, SkMatrix<T_> &A)
 *  \ingroup Equation
 *  \brief Assemble local matrix into global skyline matrix
 *  @param [in] e Reference to local entity (Element or Side)
 *  @param [in] ae Local matrix
 *  @param [in,out] A Global matrix
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_,class E_>
void element_assembly(const E_&                    e,
                      const LocalMatrix<T_,N_,N_>& ae,
                      SkMatrix<T_>&                A);


/** \fn void element_assembly(const E_ &e, const LocalMatrix<T_,N_,N_> &ae, SkSMatrix<T_> &A)
 *  \ingroup Equation
 *  \brief Assemble local matrix into global symmetric skyline matrix
 *  @param [in] e Reference to local entity (Element or Side)
 *  @param [in] ae Local matrix
 *  @param [in,out] A Global matrix
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_,class E_>
void element_assembly(const E_&                    e,
                      const LocalMatrix<T_,N_,N_>& ae,
                      SkSMatrix<T_>&               A);


/** \fn void element_assembly(const E_ &e, const LocalMatrix<T_,N_,N_> &ae, SpMatrix<T_> &A)
 *  \ingroup Equation
 *  \brief Assemble local matrix into global sparse matrix
 *  @param [in] e Reference to local entity (Element or Side)
 *  @param [in] ae Local matrix
 *  @param [in,out] A Global matrix
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_,class E_>
void element_assembly(const E_&                    e,
                      const LocalMatrix<T_,N_,N_>& ae,
                      SpMatrix<T_>&                A);


#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** \fn void side_assembly(const E_ &e, const LocalVect<T_,N_> &be, Vect<T_> &b)
 *  \ingroup Equation
 *  \brief Assemble local vector into global vector by selecting one DOF.
 *  @param [in] e Reference to local entity (Element or Side)
 *  @param [in] be Local vector
 *  @param [in,out] b Global vector
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_,class E_>
void side_assembly(const E_&               e,
                   const LocalVect<T_,N_>& be,
                   Vect<T_>&               b);


/** \fn void side_assembly(const E_ &e, const LocalMatrix<T_,N_,N_> &ae, Vect<T_> &b)
 *  \ingroup Equation
 *  \brief Assemble local diagonal matrix into global vector by selecting one DOF.
 *  @param [in] e Reference to local entity (Element or Side)
 *  @param [in] ae Local matrix 
 *  @param [in,out] b Global vector
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_,class E_>
void side_assembly(const E_&                    e,
                   const LocalMatrix<T_,N_,N_>& ae,
                   Vect<T_>&                    b);


/** \fn void side_assembly(const E_ &e, const LocalMatrix<T_,N_,N_> &ae, SkSMatrix<T_> &A)
 *  \ingroup Equation
 *  \brief Assemble local matrix into global matrix by selecting one DOF.
 *  @param [in] e Reference to local entity (Element or Side)
 *  @param [in] ae Local vector
 *  @param [in,out] A Global vector
 */
template<class T_,size_t N_,class E_>
void side_assembly(const E_&                    e,
                   const LocalMatrix<T_,N_,N_>& ae,
                   SkSMatrix<T_>&               A);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/** \fn void side_assembly(const Element &e, const LocalMatrix<T_,N_,N_> &ae, SpMatrix<T_> &A)
 *  \brief %Side assembly of local matrix into global matrix (as instance of class SpMatrix).
 *  \ingroup Equation
 *  @param [in] e Reference to local Element
 *  @param [in] ae Local matrix
 *  @param [in,out] A Global matrix
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_>
void side_assembly(const Element&               e,
                   const LocalMatrix<T_,N_,N_>& ae,
                   SpMatrix<T_>&                A);


#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_,size_t N_>
void side_assembly(const Element&               e,
                   const LocalMatrix<T_,N_,N_>& ae,
                   Matrix<T_>*                  A);


template<class T_, size_t N_>
void Assembly(const Element&                   el,
              const LocalMatrix<real_t,N_,N_>& Ae,
              SpMatrix<T_>&                    A);


template<class T_, size_t N_>
void Assembly(const Element&              el,
              const LocalVect<real_t,N_>& be,
              Vect<T_>&                   b);
#endif  /* DOXYGEN_SHOULD_SKIP_THIS */


/** \fn void side_assembly(const Element &e, const LocalMatrix<T_,N_,N_> &ae, SkSMatrix<T_> &A)
 *  \ingroup Equation
 *  \brief %Side assembly of local matrix into global matrix (as instance of class SkSMatrix).
 *  @param [in] e Reference to local Element
 *  @param [in] ae Local matrix
 *  @param [in,out] A Global matrix
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_>
void side_assembly(const Element&               e,
                   const LocalMatrix<T_,N_,N_>& ae,
                   SkSMatrix<T_>&               A);


/** \fn void side_assembly(const Element &e, const LocalMatrix<T_,N_,N_> &ae, SkMatrix<T_> &A)
 *  \ingroup Equation
 *  \brief Side assembly of local matrix into global matrix (as instance of class SkMatrix).
 *  @param [in] e Reference to local Element
 *  @param [in] ae Local matrix
 *  @param [in,out] A Global matrix
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_>
void side_assembly(const Element&               e,
                   const LocalMatrix<T_,N_,N_>& ae,
                   SkMatrix<T_>&                A);


/** \fn void side_assembly(const Element &e, const LocalVect<T_,N_> &be, Vect<T_> &b)
 *  \ingroup Equation
 *  \brief %Side assembly of local vector into global vector.
 *  @param [in] e Reference to local Element
 *  @param [in] be Local vector
 *  @param [in,out] b Global vector
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_>
void side_assembly(const Element&          e,
                   const LocalVect<T_,N_>& be,
                   Vect<T_>&               b);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class E_,class T_>
void element_assembly(const E_& e,
                      const T_* ae,
                      const T_* be,
                      Vect<T_>& A,
                      Vect<T_>& b);


template<class E_>
void Element_Assembly(const E_&     e,
                      const real_t* be,
                      Vect<real_t>& b);


template<class E_,class T_>
void element_assembly(const E_& e,
                      const T_* be,
                      Vect<T_>& b);


template<class E_,class T_>
void element_assembly(const E_&   e,
                      const T_*   ae,
                      const T_*   be,
                      Matrix<T_>* A,
                      Vect<T_>&   b);


template<class E_,class T_>
void element_assembly(const E_&   e,
                      const T_*   ae,
                      Matrix<T_>* A);


template<class E_,class T_>
void element_assembly(const E_& e,
                      const T_* ae,
                      Vect<T_>* A);


template<class T_>
void update_bc(const Element&  el,
               const Vect<T_>& bc,
               const T_*       eA,
               T_*             eb);


template<class T_>
void update_bc(const Side&     sd,
               const Vect<T_>& bc,
               const T_*       sA,
               T_*             sb);


template<class T_>
void update_bc_diag(const Element&  el,
                    const Vect<T_>& bc,
                    const T_*       eA,
                    T_*             eb);


template<class T_>
void update_bc_diag(const Side&     sd,
                    const Vect<T_>& bc,
                    const T_*       sA,
                    T_*             sb);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*  @}  */

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
