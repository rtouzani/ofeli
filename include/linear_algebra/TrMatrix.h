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

                   Template Class Matrix for tridiagonal matrices

  ==============================================================================*/


#ifndef _TR_MATRIX_H
#define _TR_MATRIX_H

#include "linear_algebra/Matrix.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file TrMatrix.h
 *  \brief Definition file for class TrMatrix.
 */

class Mesh;


/*! \class TrMatrix
 *  \ingroup VectMat
 * \brief To handle tridiagonal matrices.
 * \details 
 * This class enables storing and manipulating tridiagonal matrices.
 * The template parameter is the type of matrix entries.
 * Any matrix entry can be accessed by the () operator: For instance,
 * if \c A is an instance of this class, \c A(i,j) stands for the entry
 * at the i-th row and j-th column, \c i and \c j starting from 1.
 * If \c is difference from \c i-1, \c i or \c i+1, the returned value is \c 0.
 * Entries of \c A can be assigned a value by the same operator. Only nonzero 
 * entries can be assigned.
 *
 * \tparam T_ Data type (double, float, complex<double>, ...)
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */


template<class T_>
class TrMatrix : public Matrix<T_>
{

 public:

   using Matrix<T_>::operator();
   using Matrix<T_>::_size;
   using Matrix<T_>::_nb_cols;
   using Matrix<T_>::_nb_rows;
   using Matrix<T_>::_length;
   using Matrix<T_>::_a;
   using Matrix<T_>::_zero;
   using Matrix<T_>::_temp;

/// \brief Default constructor.
/// \details Initialize a zero dimension tridiagonal matrix
   TrMatrix();

/// \brief Constructor for a tridiagonal matrix with <tt>size</tt> rows.
    TrMatrix(size_t size);

/// \brief Copy Constructor
    TrMatrix(const TrMatrix& m);

/// \brief Destructor
    ~TrMatrix();

/// \brief Define matrix as identity matrix
    void Identity();

/// \brief Define matrix as a diagonal one
    void Diagonal();

/// \brief Define matrix as a diagona one and assign value \c a to all
/// diagonal entries
    void Diagonal(const T_& a);

/** \brief Define matrix as the one of 3-point finite difference discretization
 *  of the second derivative
 *  @param [in] h mesh size
 */    
    void Laplace1D(real_t h);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setGraph(const Vect<RC>& I,
                  int             opt=1);

    void setMesh(Mesh&  mesh,
                 size_t dof=0);

    void setMesh(Mesh&  mesh,
                 size_t dof,
                 size_t type);

    virtual void setMesh(size_t dof,
                         Mesh&  mesh,
                         int    code=0);

    virtual void setMesh(size_t dof,
                         size_t nb_eq,
                         Mesh&  mesh);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Set size (number of rows) of matrix.
/// @param [in] size Number of rows and columns.
    void setSize(size_t size);


/// \brief Multiply matrix by vector <tt>x</tt> and add result to <tt>y</tt>.
    void MultAdd(const Vect<T_>& x, Vect<T_>& y) const;

/// \brief Multiply matrix by vector <tt>a*x</tt> and add result to <tt>y</tt>.
    void MultAdd(T_              a,
                 const Vect<T_>& x,
                 Vect<T_>&       y) const;

/// \brief Multiply matrix by vector <tt>x</tt> and save result in <tt>y</tt>.
    void Mult(const Vect<T_>& x,
              Vect<T_>&       y) const;

/// \brief Multiply transpose of matrix by vector <tt>x</tt> and save result in <tt>y</tt>.
    void TMult(const Vect<T_>& x,
               Vect<T_>&       y) const;

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                  a,
              const TrMatrix<T_>& m);
   
/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                a,
              const Matrix<T_>* m);

/// \brief Assign constant <tt>val</tt> to an entry <tt>(i,j)</tt> of the matrix.
    void set(size_t    i,
             size_t    j,
             const T_& val);

/// \brief Add constant <tt>val</tt> value to an entry <tt>(i,j)</tt> of the matrix.
    void add(size_t    i,
             size_t    j,
             const T_& val);

/** \brief Operator () (Constant version).
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    T_ operator()(size_t i,
                  size_t j) const;

/** \brief Operator () (Non constant version).
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    T_ & operator()(size_t i,
                    size_t j);

/// \brief Operator =.
/// \details Copy matrix <tt>m</tt> to current matrix instance.
    TrMatrix<T_> & operator=(const TrMatrix<T_>& m);

/// \brief Operator =
/// Assign matrix to identity times <tt>x</tt>.
    TrMatrix<T_> & operator=(const T_& x);

/// \brief Operator *=.
/// \details Premultiply matrix entries by constant value <tt>x</tt>.
    TrMatrix<T_> & operator*=(const T_& x);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int Factor();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Solve a linear system with current matrix (forward and back substitution).
 *  @param [in,out] b Vect instance that contains right-hand side on input and solution on output.
 *  @param [in] fact Ununsed argument
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 * 
 *  \b Warning: %Matrix is modified after this function.
 */
    int solve(Vect<T_>& b,
              bool      fact=true);

/** \brief Solve a linear system with current matrix (forward and back substitution).
 *  @param [in] b Vect instance that contains right-hand side.
 *  @param [out] x Vect instance that contains solution.
 *  @param [in] fact Unused argument
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 * 
 *  \b Warning: %Matrix is modified after this function.
 */
    int solve(const Vect<T_>& b,
              Vect<T_>&       x,
              bool            fact=false);

/// \brief Return C-Array.
    T_ *get() const;

/// \brief Return entry <tt>(i,j)</tt> of matrix
    T_ get(size_t i, size_t j) const;

};

///////////////////////////////////////////////////////////////////////////////
//                         I M P L E M E N T A T I O N                       //
///////////////////////////////////////////////////////////////////////////////

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<>
inline void TrMatrix<real_t>::Laplace1D(real_t h);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


///////////////////////////////////////////////////////////////////////////////
//                           ASSOCIATED  FUNCTIONS                           //
///////////////////////////////////////////////////////////////////////////////

/** \fn Vect<T_> operator*(const TrMatrix<T_> &A, const Vect<T_> &b)
 *  \brief Operator * (Multiply vector by matrix and return resulting vector
 *  \ingroup VectMat
 *  @param [in] A TrMatrix instance to multiply by vector
 *  @param [in] b Vect instance 
 *  \return Vect instance containing <tt>A*b</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_>
Vect<T_> operator*(const TrMatrix<T_>& A,
                   const Vect<T_>&     b);

/** \fn TrMatrix<T_> operator*(T_ a, const TrMatrix<T_> &A)
 *  \brief Operator * (Premultiplication of matrix by constant)
 *  \ingroup VectMat
 *  @return a*A
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_>
TrMatrix<T_> operator*(T_                  a,
                       const TrMatrix<T_>& A);


/** \fn ostream& operator<<(ostream& s, const TrMatrix<T_>& a)
 *  \ingroup VectMat
 *  \brief Output matrix in output stream
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_>
ostream& operator<<(ostream&            s,
                    const TrMatrix<T_>& A);

} /* namespace OFELI */

#endif
