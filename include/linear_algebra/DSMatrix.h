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

               Definition of class DSMatrix for symmetric matrices

  ==============================================================================*/


#ifndef __DSMATRIX_H
#define __DSMATRIX_H

#include "linear_algebra/Matrix.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

class Mesh;


/*! \file DSMatrix.h
 *  \brief Definition file for abstract class DSMatrix.
 */

/*! \class DSMatrix
 *  \ingroup VectMat
 * \brief To handle symmetric dense matrices.
 *
 * This class enables storing and manipulating symmetric dense matrices.
 *
 * \tparam T_ Data type (double, float, complex<double>, ...)
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */


template<class T_>
class DSMatrix : public Matrix<T_>
{
   using Matrix<T_>::_msize;
   using Matrix<T_>::_a;
   using Matrix<T_>::_diag;
   using Matrix<T_>::_theMesh;
   using Matrix<T_>::_is_diagonal;

 public:

    using Matrix<T_>::operator();

/// \brief Default constructor.
    DSMatrix();

/** \brief Constructor that for a symmetric matrix with given number of r‡qows.
 *  @param [in] dim Number of rows
 */
    DSMatrix(size_t dim);

/// \brief Copy Constructor
/// @param [in] m DSMatrix instance to copy
    DSMatrix(const DSMatrix<T_>& m);

/** \brief Constructor using mesh to initialize matrix.
 *  @param [in] mesh Mesh instance for which matrix graph is determined.
 *  @param [in] dof Option parameter, with default value <tt>0</tt>.\n
 *  <tt>dof=1</tt> means that only one degree of freedom for each node (or element or side)
 *  is taken to determine matrix structure. The value <tt>dof=0</tt> means that matrix
 *  structure is determined using all DOFs.
 *  @param [in] is_diagonal Boolean argument to say is the matrix is actually a diagonal matrix or not.
 */
    DSMatrix(Mesh&  mesh,
             size_t dof=0,
             int    is_diagonal=false);

/// \brief Destructor
    ~DSMatrix();

/// \brief Store diagonal entries in a separate internal vector
    void setDiag();

/// \brief Set size (number of rows) of matrix.
/// @param [in] dim Number of rows and columns.
    void setSize(size_t dim);

/** \brief Assign constant to entry <tt>(i,j)</tt> of the matrix.
 *  @param [in] i row index
 *  @param [in] j column index
 *  @param [in] val value to assign to <tt>a(i,j)</tt>
 */
    void set(size_t    i,
             size_t    j,
             const T_& val);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setGraph(const vector<RC>& I,
                  int               opt=1);

    void setMesh(Mesh&  mesh,
                 size_t dof=0);

    void setMesh(size_t dof,
                 Mesh&  mesh,
                 int    code=0);

    void setMesh(size_t dof,
                 size_t nb_eq,
                 Mesh&  mesh);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Get <tt>j</tt>-th column vector
 *  @param [in] j Index of column to extract
 *  @param [out] v Reference to Vect instance where the column is stored
 *  @remark Vector v does not need to be sized before. It is resized in the function
 */
    void getColumn(size_t    j,
                   Vect<T_>& v) const;

/** \brief Get <tt>j</tt>-th column vector
 *  @param [in] j Index of column to extract
 *  @return Vect instance where the column is stored
 *  @remark Vector v does not need to be sized before. It is resized in the function
 */
    Vect<T_> getColumn(size_t j) const;

/** \brief Get <tt>i</tt>-th row vector
 *  @param [in] i Index of row to extract
 *  @param [out] v Reference to Vect instance where the row is stored
 *  @remark Vector v does not need to be sized before. It is resized in the function
 */
    void getRow(size_t    i,
                Vect<T_>& v) const;

/** \brief Get <tt>i</tt>-th row vector
 *  @param [in] i Index of row to extract
 *  @return Vect instance where the row is stored
 *  @remark Vector v does not need to be sized before. It is resized in the function
 */
    Vect<T_> getRow(size_t i) const;

/** \brief Copy a given vector to a prescribed row in the matrix.
 *  @param [in] i row index to be assigned
 *  @param [in] v Vect instance to copy
 */
    void setRow(size_t          i,
                const Vect<T_>& v);

/** \brief Copy a given vector to a prescribed column in the matrix.
 *  @param [in] j column index to be assigned
 *  @param [in] v Vect instance to copy
 */
    void setColumn(size_t          j,
                   const Vect<T_>& v);

/// \brief Set matrix as diagonal and assign its diagonal entries as a constant
/// @param [in] a Value to assign to all diagonal entries
    void setDiag(const T_& a);

/// \brief Set matrix as diagonal and assign its diagonal entries
/// @param [in] d Vector entries to assign to matrix diagonal entries
    void setDiag(const vector<T_>& d);

/** \brief Add constant to an entry ofthe matrix.
 *  @param [in] i row index
 *  @param [in] j column index
 *  @param [in] val value to add to <tt>a(i,j)</tt>
 */
    void add(size_t    i,
             size_t    j,
             const T_& val);

/** \brief Return a value of a matrix entry
 *  @param [in] i Row index (starts at 1)
 *  @param [in] j Column index (starts at 1)
 */
    T_ at(size_t i,
          size_t j);

/** \brief Operator <tt>()</tt> (Constant version).
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    T_ operator()(size_t i,
                  size_t j) const;

/** \brief Operator <tt>()</tt> (Non constant version).
 *  @param [in] i Row index
 *  @param [in] j Column index
 *  @warning To modify a value of an entry of the matrix it is safer not to
 *  modify both lower and upper triangles. Otherwise, wrong values will be 
 *  assigned. If not sure, use the member functions set or add.
 */
    T_ & operator()(size_t i,
                    size_t j);

/// \brief Operator <tt>=</tt>
/// Copy matrix <tt>m</tt> to current matrix instance.
    DSMatrix<T_> & operator=(const DSMatrix<T_>& m);

/// \brief Operator =
/// Assign matrix to identity times <tt>x</tt>.
    DSMatrix<T_> & operator=(const T_& x);

/// \brief Operator +=.
/// \details Add constant value <tt>x</tt> to all matrix entries.
    DSMatrix & operator+=(const T_& x);

/// \brief Operator -=.
/// \details Subtract constant value <tt>x</tt> from to all matrix entries.
    DSMatrix & operator-=(const T_& x);

/** \brief Factorize matrix (<tt>LDL<sup>T</sup></tt>)
 *  @return
 *  <ul>
 *    <li><tt>0</tt>, if factorization was normally performed,
 *    <li><tt>n</tt>, if the <tt>n</tt>-th pivot is null.
 *  </ul>
 */
    int setLDLt();

/// \brief Multiply matrix by vector <tt>a*x</tt> and add result to <tt>y</tt>.
    void MultAdd(const Vect<T_>& x,
                 Vect<T_>&       y) const;

/** \brief Multiply matrix by vector <tt>a*x</tt> and add to <tt>y</tt>.
 *  @param [in] a Constant to multiply by matrix
 *  @param [in] x Vector to multiply by matrix
 *  @param [in,out] y Vector to add to the result. <tt>y</tt> contains on output the result.
 */
    void MultAdd(T_              a,
                 const Vect<T_>& x,
                 Vect<T_>&       y) const;

/// \brief Multiply matrix by vector <tt>x</tt> and save result in <tt>y</tt>.
    void Mult(const Vect<T_>& x,
              Vect<T_>&       y) const;

/** \brief Multiply transpose of matrix by vector <tt>x</tt> and add result in <tt>y</tt>.
 *  @param [in] x Vector to add to <tt>y</tt>
 *  @param [in,out] y on input, vector to add to. On output, result.
 */
    void TMult(const Vect<T_>& x,
               Vect<T_>&       y) const;

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                  a,
              const DSMatrix<T_>& m);
 
/// \brief Add <tt>val</tt> to entry <tt>i</tt>.
    void add(size_t    i,
             const T_& val);

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                a,
              const Matrix<T_>* m);

/** \brief Solve linear system.
 *  \details The matrix is factorized using the LDLt (Crout) decomposition. If this one is already factorized,
 *  no further factorization is performed. If the matrix has been modified the user has to refactorize it 
 *  using the function setLDLt.
 *  @param [in,out] b Vect instance that contains right-hand side on input and solution on output.
 *  @param [in] fact Set true if matrix is to be factorized (Default value), false if not
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 */
    int solve(Vect<T_>& b,
              bool      fact=true);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int Factor();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Solve linear system.
 *  \details The matrix is factorized using the LDLt (Crout) decomposition. If this one is already factorized,
 *  no further factorization is performed. If the matrix has been modified the user has to refactorize it 
 *  using the function setLDLt.
 *  @param [in] b Vect instance that contains right-hand side.
 *  @param [out] x Vect instance that contains solution
 *  @param [in] fact Set true if matrix is to be factorized (Default value), false if not
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 */
    int solve(const Vect<T_>& b,
              Vect<T_>&       x,
              bool            fact=true);

/** \brief Return matrix as C-Array.
 *  Matrix is stored row by row.
 *  Only lower triangle is stored.
 */
    const T_ *getArray();

/// \brief Return entry <tt>(i,j)</tt> of matrix
    T_ get(size_t i,
           size_t j) const;

 protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   void set_(size_t i,
             size_t j,
             T_     x);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

};

///////////////////////////////////////////////////////////////////////////////
//                 A S S O C I A T E D   F U N C T I O N S                   //
///////////////////////////////////////////////////////////////////////////////

/** \fn Vect<T_> operator*(const DSMatrix<T_> &A, const Vect<T_> &b)
 *  \brief Operator * (Multiply vector by matrix and return resulting vector
 *  \ingroup VectMat
 *  @param [in] A DSMatrix instance to multiply by vector
 *  @param [in] b Vect instance 
 *  \return Vect instance containing <tt>A*b</tt>
 */
template<class T_>
Vect<T_> operator*(const DSMatrix<T_>& A,
                   const Vect<T_>&     b);

/// \fn ostream& operator<<(ostream& s, const DSMatrix<T_> &a)
/// \ingroup VectMat
/// \brief Output matrix in output stream
template<class T_>
ostream& operator<<(ostream&            s,
                    const DSMatrix<T_>& a);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
