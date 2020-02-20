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

  ==============================================================================*/

#ifndef __SCHUR_H
#define __SCHUR_H

#include "OFELI_Config.h"
#include "mesh/Mesh.h"
#include "linear_algebra/Vect.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */


/** \fn void Schur(SkMatrix<T_>& A, SpMatrix<T_>& U, SpMatrix<T_>& L, SpMatrix<T_>& D, Vect<T_>& b, Vect<T_>& c)
 *  \ingroup Solver
 *  \brief Solve a linear system of equations with a 2x2-block matrix.
 *  \details The linear system is of the form
 *
 *  @verbatim
        | A  U | |x|   |b|
        |      | | | = | |
        | L  D | |y|   |c|
    @endverbatim
 *
 * @param[in] A Instance of class SkMatrix class for the first diagonal block. The matrix must
 * be invertible and factorizable (Do not use SpMatrix class)
 * where <tt>A</tt>, <tt>U</tt>, <tt>L</tt>, <tt>D</tt> are instances of matrix classes,
 * @param [in] U Instance of class SpMatrix for the upper triangle block. The matrix can
 * be rectangular
 * @param [in] L Instance of class SpMatrix for the lower triangle block. The matrix can
 * be rectangular
 * @param [in] D Instance of class SpMatrix for the second diagonal block. The matrix must
 * be factorizable (Do not use SpMatrix class)
 * @param [in,out] b Vector (Instance of class Vect) that contains the first block of
 * right-hand side on input and the first block of the solution on output. <tt>b</tt> must
 * have the same size as the dimension of <tt>A</tt>.
 * @param [in,out] c Vect instance that contains the second block of
 * right-hand side on output and the first block of the solution on output. <tt>c</tt> must
 * have the same size as the dimension of <tt>D</tt>.
 *
 * Template Argument:
 *
 * \tparam <T_> data type (real_t, float, ...)
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
template <class T_>
void Schur(SkMatrix<T_>& A,
           SpMatrix<T_>& U,
           SpMatrix<T_>& L,
           SpMatrix<T_>& D,
           Vect<T_>&     b,
           Vect<T_>&     c)
{
   size_t nr=U.getNbRows(), nc=U.getNbColumns();
   int ret = A.Factor();
   Vect<T_> bb(b);
   A.Solve(b);
   c -= L*b;
   for (size_t j=1; j<=nc; j++) {
      Vect<T_> v = U.Column(j);
      A.Solve(v);
      Vect<T_> w = L*v;
      for (size_t i=1; i<=nc; i++)
         D.Add(i,j,-w(i));
   }
   ret = D.Factor();
   D.Solve(c);
   bb -= U*c;
   A.Solve(bb);
   b = bb;
}


#if defined(USE_PETSC)
/** \fn void Schur(PETScMatrix<T_>& A, PETScMatrix<T_>& U, PETScMatrix<T_>& L, PETScMatrix<T_>& D, PETScVect<T_>& b,
                   PETScVect<T_>& c)
 *  \ingroup Solver
 *  \brief Solve a linear system of equations with a 2x2-block matrix.
 *  \details The linear system is of the form
 *
 *  @verbatim
        | A  U | |x|   |b|
        |      | | | = | |
        | L  D | |y|   |c|
    @endverbatim
 *
 * @param[in] A Instance of class SkMatrix class for the first diagonal block. The matrix must
 * be invertible and factorizable (Do not use SpMatrix class)
 * where <tt>A</tt>, <tt>U</tt>, <tt>L</tt>, <tt>D</tt> are instances of matrix classes,
 * @param [in] U Instance of class PETScMatrix for the upper triangle block. The matrix can
 * be rectangular
 * @param [in] L Instance of class PETScMatrix for the lower triangle block. The matrix can
 * be rectangular
 * @param [in] D Instance of class PETScMatrix for the second diagonal block. The matrix must
 * be factorizable (Do not use SpMatrix class)
 * @param [in,out] b Vector (Instance of class PETScVect) that contains the first block of
 * right-hand side on input and the first block of the solution on output. <tt>b</tt> must
 * have the same size as the dimension of <tt>A</tt>.
 * @param [in,out] c PETScVect instance that contains the second block of
 * right-hand side on output and the first block of the solution on output. <tt>c</tt> must
 * have the same size as the dimension of <tt>D</tt>.
 *
 * Template Argument:
 *
 * \tparam <T_> data type (real_t, float, ...)
 */
template <class T_>
void Schur(PETScMatrix<T_>& A,
           PETScMatrix<T_>& U,
           PETScMatrix<T_>& L,
           PETScMatrix<T_>& D,
           PETScVect<T_>&   b,
           PETScVect<T_>&   c)
{
   MatGetSchurComplement(A,is0,is0,is1,is1,MAT_INITIAL_MATRIX,&S,PETSC_FALSE,MAT_IGNORE_MATRIX,NULL);
   size_t nr=U.getNbRows(), nc=U.getNbColumns();
   int ret = A.Factor();
   Vect<T_> bb(b);
   A.Solve(b);
   c -= L*b;
   for (size_t j=1; j<=nc; j++) {
      Vect<T_> v = U.Column(j);
      A.Solve(v);
      PETScVect<T_> w = L*v;
      for (size_t i=1; i<=nc; i++)
         D.Add(i,j,-w(i));
   }
   ret = D.Factor();
   D.Solve(c);
   bb -= U*c;
   A.Solve(bb);
   b = bb;
}
#endif

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
