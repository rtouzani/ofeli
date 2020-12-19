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

                      Definition of Preconditioner class

  ==============================================================================*/

#ifndef __PREC_H
#define __PREC_H

#include "OFELI_Config.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/SpMatrix.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Prec.h
 *  \brief Definition file for preconditioning classes.
 */


#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_,class M_>
int inv_diag(size_t      n,
             const M_&   A,
             vector<T_>& diag);

template<class T_>
int inv_diag(size_t            n,
             const Matrix<T_>* A,
             vector<T_>&       diag);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! \class Prec
 *  \ingroup Solver
 *  \brief To set a preconditioner.
 *  \details The preconditioner type is chosen in the constructor
 *
 * \tparam <T_> Data type (real_t, float, complex<real_t>, ...)
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

template<class T_> class Prec
{

 public:

/// \brief Default constructor.
    Prec();

/** \brief Constructor that chooses preconditioner.
 *  @param [in] type Preconditioner type:
 *  <ul>
 *    <li><tt>IDENT_PREC</tt>: Identity preconditioner (No preconditioning)
 *    <li><tt>DIAG_PREC</tt>: Diagonal preconditioner
 *    <li><tt>DILU_PREC</tt>: Diagonal Incomplete factorization preconditioner
 *    <li><tt>ILU_PREC</tt>: Incomplete factorization preconditioner
 *    <li><tt>SSOR_PREC</tt>: SSOR (Symmetric Successive Over Relaxation) preconditioner
 *  </ul>
 */
    Prec(int type);

/** \brief Constructor using matrix of the linear system to precondition
 *  @param [in] A %Matrix to precondition
 *  @param [in] type Preconditioner type:
 *  <ul>
 *    <li><tt>IDENT_PREC</tt>: Identity preconditioner (No preconditioning)
 *    <li><tt>DIAG_PREC</tt>: Diagonal preconditioner
 *    <li><tt>DILU_PREC</tt>: Diagonal Incomplete factorization preconditioner
 *    <li><tt>ILU_PREC</tt>: Incomplete factorization preconditioner
 *    <li><tt>SSOR_PREC</tt>: SSOR (Symmetric Successive Over Relaxation) preconditioner
 *  </ul>
 */
    Prec(const SpMatrix<T_>& A,
         int                 type=DIAG_PREC);

/** \brief Constructor using matrix of the linear system to precondition
 *  @param [in] A Pointer to abstract Matrix class to precondition
 *  @param [in] type Preconditioner type:
 *  <ul>
 *    <li><tt>IDENT_PREC</tt>: Identity preconditioner (No preconditioning)
 *    <li><tt>DIAG_PREC</tt>: Diagonal preconditioner
 *    <li><tt>DILU_PREC</tt>: Diagonal Incomplete factorization preconditioner
 *    <li><tt>ILU_PREC</tt>: Incomplete factorization preconditioner
 *    <li><tt>SSOR_PREC</tt>: SSOR (Symmetric Successive Over Relaxation) preconditioner
 *  </ul>
 */
    Prec(const Matrix<T_>* A,
         int               type=DIAG_PREC);

/// \brief Destructor
    ~Prec();

/** \brief Define preconditioner type
 *  @param [in] type Preconditioner type:
    <ul>
       <li> <tt>IDENT_PREC</tt>: Identity preconditioner (No preconditioning)
       <li> <tt>DIAG_PREC</tt>: Diagonal preconditioner
       <li> <tt>DILU_PREC</tt>: Diagonal Incomplete factorization preconditioner
       <li> <tt>ILU_PREC</tt>: Incomplete factorization preconditioner
       <li> <tt>SSOR_PREC</tt>: SSOR (Symmetric Successive Over Relaxation) preconditioner
    </ul>
 */
    void setType(int type);

/// \brief Define pointer to matrix for preconditioning (if this one is abstract)
/// @param [in] A %Matrix to precondition
    void setMatrix(const Matrix<T_>* A);

/// \brief Define the matrix for preconditioning
/// @param [in] A %Matrix to precondition (instance of class SpMatrix)
    void setMatrix(const SpMatrix<T_>& A);

/// \brief Solve a linear system with preconditioning matrix.
/// @param [in,out] x Right-hand side on input and solution on output.
    void solve(Vect<T_>& x) const;

/** \brief Solve a linear system with preconditioning matrix.
 *  @param [in] b Right-hand side
 *  @param [out] x Solution vector
 */
    void solve(const Vect<T_>& b,
               Vect<T_>&       x) const;

/// \brief Solve a linear system with transposed preconditioning matrix.
/// @param [in,out] x Right-hand side in input and solution in output.
    void TransSolve(Vect<T_>& x) const;

/** \brief Solve a linear system with transposed preconditioning matrix.
 *  @param [in] b Right-hand side vector
 *  @param [out] x Solution vector
 */
    void TransSolve(const Vect<T_>& b,
                    Vect<T_>&       x) const;

/// Return i-th pivot of preconditioning matrix
    T_ & getPivot(size_t i) const;

 private:
   vector<T_>         _pivot;
   const SpMatrix<T_> *_a;
   SpMatrix<T_>       _aa;
   size_t             _size, _length;
   vector<size_t>     _id, _row_ptr, _col_ind;
   int                _type;
   int inv_diag();
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
