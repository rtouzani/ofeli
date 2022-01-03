/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani

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

      Definition of Template Class LocalMatrix for dense small size matrices

  ==============================================================================*/


#ifndef __LOCAL_MATRIX_H
#define __LOCAL_MATRIX_H

#include <valarray>
#include "OFELI_Config.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file LocalMatrix.h
 *  \brief Definition file for class LocalMatrix.
 */

/*! \class LocalMatrix
 *  \ingroup VectMat
 *  \brief Handles small size matrices like element matrices, with a priori known size.
 *
 * \details The template class LocalMatrix treats small size matrices. Typically, this class is
 * recommended to store element and side arrays.\n
 * Internally, no dynamic storage is used.
 *
 * \tparam T_  Data type (double, float, complex<double>, ...)
 * \tparam NR_ number of rows of matrix
 * \tparam NC_ number of columns of matrix
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class Element;
template<class T_> class SkMatrix;
template<class T_> class SkSMatrix;
template<class T_> class SpMatrix;
template<class T_,size_t N_> class LocalVect;

template<class T_, size_t NR_, size_t NC_> class LocalMatrix
{

 public:

/// \brief Default constructor
/// \details Constructs a matrix with <tt>0</tt> rows and <tt>0</tt> columns
    LocalMatrix();

/// \brief Copy constructor
    LocalMatrix(const LocalMatrix<T_,NR_,NC_>& m);

/** \brief Constructor of a local matrix associated to element from a SpMatrix
 *  @param [in] el Pointer to Element
 *  @param [in] a Global matrix as instance of class SpMatrix.
 */
    LocalMatrix(Element*            el,
                const SpMatrix<T_>& a);

/** \brief Constructor of a local matrix associated to element from a SkMatrix
 *  @param [in] el Pointer to Element
 *  @param [in] a Global matrix as instance of class SkMatrix.
 */
    LocalMatrix(Element*            el,
                const SkMatrix<T_>& a);

/** \brief Constructor of a local matrix associated to element from a SkSMatrix
 *  @param [in] el Pointer to Element
 *  @param [in] a Global matrix as instance of class SkSMatrix.
 */
    LocalMatrix(Element*             el,
                const SkSMatrix<T_>& a);

/// Destructor
    ~LocalMatrix();

/// \brief Operator <tt>()</tt> (Non constant version)
/// \details Returns entry at row <tt>i</tt> and column <tt>j</tt>.
    T_& operator()(size_t i,
                   size_t j);

/// \brief Operator <tt>()</tt> (Constant version)
/// \details Returns entry at row <tt>i</tt> and column <tt>j</tt>.
    T_ operator()(size_t i,
                  size_t j) const;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void set(int opt);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Initialize matrix as element matrix from global SpMatrix
 *  @param [in] el Pointer to Element
 *  @param [in] a Global matrix as instance of class SpMatrix.
 *  This function is called by its corresponding constructor.
 */
    void Localize(Element*            el,
                  const SpMatrix<T_>& a);

/** \brief Initialize matrix as element matrix from global SkMatrix
 *  @param [in] el Pointer to Element
 *  @param [in] a Global matrix as instance of class SkMatrix.
 *  This function is called by its corresponding constructor.
 */
    void Localize(Element*            el,
                  const SkMatrix<T_>& a);

/** \brief Initialize matrix as element matrix from global SkSMatrix.
 *  @param [in] el Pointer to Element
 *  @param [in] a Global matrix as instance of class SkSMatrix.
 *  This function is called by its corresponding constructor.
 */
    void Localize(Element*             el,
                  const SkSMatrix<T_>& a);

/// \brief Operator <tt>=</tt>
/// \details Copy instance <tt>m</tt> into current instance.
    LocalMatrix<T_,NR_,NC_> & operator=(const LocalMatrix<T_,NR_,NC_>& m);

/// \brief Operator <tt>=</tt>
/// \details Assign matrix to identity times <tt>x</tt>
    LocalMatrix<T_,NR_,NC_> & operator=(const T_& x);

/// \brief Operator <tt>+=</tt>
/// \details Add <tt>m</tt> to current matrix.
    LocalMatrix<T_,NR_,NC_> & operator+=(const LocalMatrix<T_,NR_,NC_>& m);

/// \brief Operator <tt>-=</tt>
/// \details Subtract <tt>m</tt> from current matrix.
    LocalMatrix<T_,NR_,NC_> & operator-=(const LocalMatrix<T_,NR_,NC_>& m);

/// \brief Operator <tt>*</tt>
/// \details Return a Vect instance as product of current matrix by vector <tt>x</tt>.
    LocalVect<T_,NR_> operator*(LocalVect<T_,NC_>& x);

/// \brief Operator <tt>+=</tt>
/// \details Add constant <tt>x</tt> to current matrix entries.
    LocalMatrix<T_,NR_,NC_> & operator+=(const T_& x);

/// \brief Operator <tt>-=</tt>
/// \details Subtract <tt>x</tt> from current matrix entries.
    LocalMatrix<T_,NR_,NC_> & operator-=(const T_& x);

/// \brief Operator <tt>*=</tt>
/// \details Multiply matrix entries by constant <tt>x</tt>.
    LocalMatrix<T_,NR_,NC_> & operator*=(const T_& x);

/// \brief Operator <tt>/=</tt>
/// \details Divide by <tt>x</tt> current matrix entries.
    LocalMatrix<T_,NR_,NC_> & operator/=(const T_ &x);

/** \brief Multiply matrix by vector and add result to vector.
 *  @param [in] x Vector to multiply matrix by.
 *  @param [out] y Resulting vector (<tt>y += a * x</tt>)
 */
    void MultAdd(const LocalVect<T_,NC_>& x,
                 LocalVect<T_,NR_>&       y);

/** \brief Multiply matrix by scaled vector and add result to vector.
 *  @param [in] a Constant to premultiply by vector <tt>x</tt>.
 *  @param [in] x (Scaled) vector to multiply matrix by.
 *  @param [out] y Resulting vector (<tt>y += a * x</tt>)
 */
    void MultAddScal(const T_&                a,
                     const LocalVect<T_,NC_>& x,
                     LocalVect<T_,NR_>&       y);

/** \brief Multiply matrix by vector.
 *  @param [in] x Vector to multiply matrix by.
 *  @param [out] y Resulting vector.
 */
    void Mult(const LocalVect<T_,NC_>& x,
              LocalVect<T_,NR_>&       y);

/// \brief Symmetrize matrix
/// \details Fill upper triangle to form a symmetric matrix.
    void Symmetrize();

/** \brief Factorize matrix
 *  \details Performs a LU factorization.
 *  @return
 *  <ul>
 *     <li><tt>0</tt>: Factorization has ended normally,
 *     <li><tt>n</tt>: <tt>n</tt>-th pivot was zero.
 *  </ul>
 */
    int Factor();

/** \brief Forward and backsubstitute to solve a linear system.
 *  @param [in] b Right-hand side in input and solution vector in output.
 *  @return
 *  <ul>
 *     <li><tt>0</tt>: Solution was performed normally.
 *     <li><tt>n</tt>: <tt>n</tt>-th pivot is zero.
 *  </ul>
 *  @note %Matrix must have been factorized at first.
 */
    int solve(LocalVect<T_,NR_>& b);

/** \brief Factorize matrix and solve linear system.
 *  @param [in,out] b Right-hand side in input and solution vector in output.
 *  @return
 *     <tt>0</tt> if solution was performed normally.
 *     <tt>n</tt> if <tt>n</tt>-th pivot is zero.
 *  This function simply calls \b Factor() then \b Solve(b).
 */
    int FactorAndSolve(LocalVect<T_,NR_>& b);

/// \brief Calculate inverse of matrix
/// @param [out] A Inverse of matrix
    void Invert(LocalMatrix<T_,NR_,NC_>& A);

/** \brief Calculate inner product witrh respect to matrix.
 *  \details Returns the product <tt>x<sup>T</sup>Ay</tt> 
 *  @param [in] x Left vector
 *  @param [in] y Right vector
 *  @return Resulting product
 */
    T_ getInnerProduct(const LocalVect<T_,NC_>& x,
                       const LocalVect<T_,NR_>& y);

/// \brief Return pointer to matrix as a C-array.
    T_ *get();

 private:
    size_t            _length;
    std::valarray<T_> _a;
};


///////////////////////////////////////////////////////////////////////////////
//              A S S O C I A T E D    F U N C T I O N S                     //
///////////////////////////////////////////////////////////////////////////////

/** \fn LocalMatrix<T_,NR_,NC_> operator*(T_ a, const LocalMatrix<T_,NR_,NC_> &x)
    \brief Operator * (Multiply matrix <tt>x</tt> by scalar <tt>a</tt>)
    \ingroup VectMat
    \return <tt>a*x</tt>
 */
template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_> operator*(T_                             a,
                                  const LocalMatrix<T_,NR_,NC_>& x);

/** \fn LocalVect<T_,NR_,NC_> operator*(const LocalMatrix<T_,NR_,NC_> &x, const LocalVect<T_,NC_> &x)
    \brief Operator <tt>*</tt> (Multiply matrix <tt>A</tt> by vector <tt>x</tt>)
    \details This function performs a matrix-vector product and returns resulting vector as
    a reference to LocalVect instance
    \ingroup VectMat
    \return  <tt>A*x</tt>
 */
template<class T_,size_t NR_,size_t NC_>
LocalVect<T_,NR_> operator*(const LocalMatrix<T_,NR_,NC_>& A,
                            const LocalVect<T_,NC_>&       x);

/** \fn LocalMatrix<T_,NR_,NC_> operator/(T_ a, const LocalMatrix<T_,NR_,NC_> &x)
    \brief Operator <tt>/</tt> (Divide matrix <tt>x</tt> by scalar <tt>a</tt>)
    \ingroup VectMat
    \return <tt>x/a</tt>
 */
template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_> operator/(T_                             a,
                                  const LocalMatrix<T_,NR_,NC_>& x);

/** \fn LocalMatrix<T_,NR_,NC_> operator+(const LocalMatrix<T_,NR_,NC_> &x, const LocalMatrix<T_,NR_,NC_> &y)
    \brief Operator <tt>+</tt> (Add matrix x to y)
    \ingroup VectMat
    \return <tt>x+y</tt>
 */
template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_> operator+(const LocalMatrix<T_,NR_,NC_>& x,
                                  const LocalMatrix<T_,NR_,NC_>& y);

/** \fn LocalMatrix<T_,NR_,NC_> operator-(const LocalMatrix<T_,NR_,NC_> &x, const LocalMatrix<T_,NR_,NC_> &y)
    \brief Operator <tt>-</tt> (Subtract matrix y from x)
    \ingroup VectMat
    \return  <tt>x-y</tt>
 */
template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_> operator-(const LocalMatrix<T_,NR_,NC_>& x,
                                  const LocalMatrix<T_,NR_,NC_>& y);

/// \fn ostream& operator<<(ostream &s, const LocalMatrix<T_,NR_,NC_> &A)
/// \brief Output vector in output stream
/// \ingroup VectMat
template<class T_, size_t NR_, size_t NC_>
ostream& operator<<(ostream&                       s,
                    const LocalMatrix<T_,NR_,NC_>& A);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
