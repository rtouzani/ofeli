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

              Definition of Template Class DMatrix for dense matrices

  ==============================================================================*/


#ifndef __DMATRIX_H
#define __DMATRIX_H

#include "linear_algebra/Matrix_impl.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file DMatrix.h
 *  \brief Definition file for class DMatrix.
 */

/*! \class DMatrix
 *  \ingroup VectMat
 * \brief To handle dense matrices.
 *
 * This class enables storing and manipulating general dense matrices. 
 * Matrices can be square or rectangle ones.
 *
 * \tparam T_ Data type (double, float, complex<double>, ...)
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */

class Mesh;

template<class T_>
class DMatrix : public Matrix<T_>
{
   using Matrix<T_>::_nb_rows;
   using Matrix<T_>::_nb_cols;
   using Matrix<T_>::_size;
   using Matrix<T_>::_length;
   using Matrix<T_>::_a;
   using Matrix<T_>::_diag;
   using Matrix<T_>::_ch;
   using Matrix<T_>::_zero;
   using Matrix<T_>::_theMesh;
   using Matrix<T_>::_is_diagonal;


 public:

    using Matrix<T_>::operator();

/*----------------------------  BASIC OPERATORS  ----------------------------*/

/// \brief Default constructor.
/// \details Initializes a zero-dimension matrix.
    DMatrix();

/// \brief Constructor for a matrix with <tt>nr</tt> rows and <tt>nr</tt> columns.
/// \details %Matrix entries are set to <tt>0</tt>.
    DMatrix(size_t nr);

/// \brief Constructor for a matrix with <tt>nr</tt> rows and <tt>nc</tt> columns.
/// \details Matrix entries are set to 0.
    DMatrix(size_t nr,
            size_t nc);

/** \brief Constructor that uses a Vect instance. The class uses the memory
 *  space occupied by this vector.
 *  @param [in] v Vector to copy
 */
    DMatrix(Vect<T_>& v);

/// \brief Copy Constructor.
/// @param [in] m %Matrix to copy
    DMatrix(const DMatrix<T_>& m);

/** \brief Constructor using mesh to initialize structure of matrix.
 *  @param [in] mesh Mesh instance for which matrix graph is determined.
 *  @param [in] dof Option parameter, with default value <tt>0</tt>.\n
 *  <tt>dof=1</tt> means that only one degree of freedom for each node (or element or side)
 *  is taken to determine matrix structure. The value <tt>dof=0</tt> means that matrix
 *  structure is determined using all DOFs.
 *  @param [in] is_diagonal Boolean argument to say is the matrix is actually a diagonal matrix or not.
 */
    DMatrix(Mesh&  mesh,
            size_t dof=0,
            int    is_diagonal=false);

/// \brief Destructor.
    ~DMatrix();

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setMesh(Mesh&  mesh,
                 size_t dof=0);

    void setMesh(size_t dof,
                 Mesh&  mesh,
                 int    code=0);

    void setMesh(size_t dof,
                 size_t nb_eq,
                 Mesh&  mesh);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Store diagonal entries in a separate internal vector
    void setDiag();

/// \brief Set matrix as diagonal and assign its diagonal entries as a constant
/// @param [in] a Value to assign to all diagonal entries
    void setDiag(const T_& a);

/// \brief Set matrix as diagonal and assign its diagonal entries
/// @param [in] d Vector entries to assign to matrix diagonal entries
    void setDiag(const vector<T_>& d);

/// \brief Set size (number of rows) of matrix.
/// @param [in] size Number of rows and columns.
    void setSize(size_t size);

/** \brief Set size (number of rows and columns) of matrix.
 *  @param [in] nr Number of rows.
 *  @param [in] nc Number of columns.
 */
    void setSize(size_t nr,
                 size_t nc);

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

/** \brief Assign a constant value to an entry of the matrix.
 *  @param [in] i row index of matrix
 *  @param [in] j column index of matrix
 *  @param [in] val Value to assign to <tt>a(i,j)</tt>.
 */
    void set(size_t    i,
             size_t    j,
             const T_& val);

/** \brief Set matrix to 0 and reset factorization parameter
 *  @warning This function must be used if after a factorization, the matrix has
 *  modified
 */
    void reset();

/** \brief Copy a given vector to a prescribed row in the matrix.
 *  @param [in] i row index to be assigned
 *  @param [in] v Vect instance to copy
 */
    void setRow(size_t          i,
                const Vect<T_>& v);

/** \brief Copy a given vector to a prescribed column in the matrix.
 *  @param [in] i column index to be assigned
 *  @param [in] v Vect instance to copy
 */
    void setColumn(size_t          j,
                   const Vect<T_>& v);

/** \brief Multiply matrix by vector <tt>a*x</tt> and add result to <tt>y</tt>.
 *  @param [in] a constant to multiply by
 *  @param [in] x Vector to multiply by <tt>a</tt>
 *  @param [in,out] y on input, vector to add to. On output, result.
 */
    void MultAdd(T_              a,
                 const Vect<T_>& x,
                 Vect<T_>&       y) const;

/** \brief Multiply matrix by vector <tt>x</tt> and add result to <tt>y</tt>.
 *  @param [in] x Vector to add to <tt>y</tt>
 *  @param [in,out] y on input, vector to add to. On output, result.
 */
    void MultAdd(const Vect<T_>& x,
                 Vect<T_>&       y) const;

/** \brief Multiply matrix by vector <tt>x</tt> and save result in <tt>y</tt>.
 *  @param [in] x Vector to add to <tt>y</tt>
 *  @param [out] y Result.
 */
    void Mult(const Vect<T_>& x,
              Vect<T_>&       y) const;

/** \brief Multiply transpose of matrix by vector <tt>x</tt> and add result in <tt>y</tt>.
 *  @param [in] x Vector to add to <tt>y</tt>
 *  @param [in,out] y on input, vector to add to. On output, result.
 */
    void TMult(const Vect<T_>& x,
               Vect<T_>&       y) const;

/** \brief Add constant <tt>val</tt> to entry <tt>(i,j)</tt> of the matrix.
 *  @param [in] i row index
 *  @param [in] j column index
 *  @param [in] val Constant to add
 */
    void add(size_t    i,
             size_t    j,
             const T_& val);

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                 a,
              const DMatrix<T_>& m);

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                a,
              const Matrix<T_>* m);

/** \brief Construct a QR factorization of the matrix
 *  \details This function constructs the QR decomposition using the Householder method.
 *  The upper triangular matrix R is returned in the upper triangle of the current matrix,
 *  except for the diagonal elements of R which are stored in an internal vector.
 *  The orthogonal matrix Q is represented as a product of n-1 Householder matrices
 *  Q1 . . . Qn-1, where Qj = 1 - uj.uj /cj . The i-th component of uj is zero for i = 1, ..., j-1
 *  while the nonzero components are returned in a[i][j] for i = j, ..., n. 
 *  @return \c 0 if the decomposition was successful, \c k is the \c k-th row is singular
 *  @remark The matrix can be square or rectangle
 */
    int setQR();

/** \brief Construct a QR factorization of the transpose of the matrix
 *  \details This function constructs the QR decomposition using the Householder method.
 *  The upper triangular matrix R is returned in the upper triangle of the current matrix,
 *  except for the diagonal elements of R which are stored in an internal vector.
 *  The orthogonal matrix Q is represented as a product of n-1 Householder matrices
 *  Q1 . . . Qn-1, where Qj = 1 - uj.uj /cj . The i-th component of uj is zero for i = 1, ..., j-1
 *  while the nonzero components are returned in a[i][j] for i = j, ..., n. 
 *  @return \c 0 if the decomposition was successful, \c k is the \c k-th row is singular
 *  @remark The matrix can be square or rectangle
 */
    int setTransQR();

/** \brief Solve a linear system by QR decomposition
 *  \details This function constructs the QR decomposition, if this was not already done
 *  by using the member function QR and solves the linear system 
 *  @param [in] b Right-hand side vector
 *  @param [out] x Solution vector. Must have been sized before using this function.
 *  @return The same value as returned by the function QR
 */
    int solveQR(const Vect<T_>& b,
                Vect<T_>&       x);

/** \brief Solve a transpose linear system by QR decomposition
 *  \details This function constructs the QR decomposition, if this was not already done
 *  by using the member function QR and solves the linear system 
 *  @param [in] b Right-hand side vector
 *  @param [out] x Solution vector. Must have been sized before using this function.
 *  @return The same value as returned by the function QR
 */
    int solveTransQR(const Vect<T_>& b,
                     Vect<T_>&       x);

/** \brief Operator () (Constant version).
 *  Return <tt>a(i,j)</tt>
 *  @param [in] i row index
 *  @param [in] j column index
 */
    T_ operator()(size_t i,
                  size_t j) const;

/** \brief Operator () (Non constant version).
 *  Return <tt>a(i,j)</tt>
 *  @param [in] i row index
 *  @param [in] j column index
 */
    T_ & operator()(size_t i,
                    size_t j);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setGraph(const Vect<RC>& I,
                  int             opt=1);

    int Factor();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Factorize the matrix (LU factorization)
 *  \details LU factorization of the matrix is realized. Note that since this
 *  is an in place factorization, the contents of the matrix are modified.
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if factorization was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 *  @remark A flag in this class indicates after factorization that this one
 *  has been realized, so that, if the member function solve is called after this
 *  no further factorization is done.
 */
    int setLU();

/** \brief Factorize the transpose of the matrix (LU factorization)
 *  \details LU factorization of the transpose of the matrix is realized. Note that since this
 *  is an in place factorization, the contents of the matrix are modified.
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if factorization was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 *  @remark A flag in this class indicates after factorization that this one
 *  has been realized, so that, if the member function solve is called after this
 *  no further factorization is done.
 */
    int setTransLU();

/** \brief Solve linear system.
 *  \details The linear system having the current instance as a matrix is solved by using the LU decomposition.
 *  Solution is thus realized after a factorization step and a forward/backward substitution step.
 *  The factorization step is realized only if this was not already done.\n
 *  Note that this function modifies the matrix contents is a factorization is performed. Naturally, if the
 *  the matrix has been modified after using this function, the user has to refactorize it using the function setLU.
 *  This is because the class has no non-expensive way to detect if the matrix has been modified.
 *  The function setLU realizes the factorization step only.
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

/** \brief Solve the transpose linear system.
 *  \details The linear system having the current instance as a transpose matrix is solved by using the LU decomposition.
 *  Solution is thus realized after a factorization step and a forward/backward substitution step.
 *  The factorization step is realized only if this was not already done.\n
 *  Note that this function modifies the matrix contents is a factorization is performed. Naturally, if the
 *  the matrix has been modified after using this function, the user has to refactorize it using the function setLU.
 *  This is because the class has no non-expensive way to detect if the matrix has been modified.
 *  The function setLU realizes the factorization step only.
 *  @param [in,out] b Vect instance that contains right-hand side on input and solution on output.
 *  @param [in] fact Set true if matrix is to be factorized (Default value), false if not
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 */
    int solveTrans(Vect<T_>& b,
                   bool      fact=true);

/** \brief Solve linear system.
 *  \details The linear system having the current instance as a matrix is solved by using the LU decomposition.
 *  Solution is thus realized after a factorization step and a forward/backward substitution step.
 *  The factorization step is realized only if this was not already done.\n
 *  Note that this function modifies the matrix contents is a factorization is performed. Naturally, if the
 *  the matrix has been modified after using this function, the user has to refactorize it using the function setLU.
 *  This is because the class has no non-expensive way to detect if the matrix has been modified.
 *  The function setLU realizes the factorization step only.
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

/** \brief Solve the transpose linear system.
 *  \details The linear system having the current instance as a transpose matrix is solved by 
 *  using the LU decomposition.
 *  Solution is thus realized after a factorization step and a forward/backward substitution step.
 *  The factorization step is realized only if this was not already done.\n
 *  Note that this function modifies the matrix contents is a factorization is performed. 
 *  Naturally, if the
 *  the matrix has been modified after using this function, the user has to refactorize 
 *  it using the function setLU.
 *  This is because the class has no non-expensive way to detect if the matrix has been modified.
 *  The function setLU realizes the factorization step only.
 *  @param [in] b Vect instance that contains right-hand side.
 *  @param [out] x Vect instance that contains solution
 *  @param [in] fact Set true if matrix is to be factorized (Default value), false if not
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 */
    int solveTrans(const Vect<T_>& b,
                   Vect<T_>&       x,
                   bool            fact=true);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int solveLU(const Vect<T_>& b,
                Vect<T_>&       x,
                bool            fact=true);

    int solveLU(Vect<T_>& b,
                bool      fact=true);

    int solveTransLU(const Vect<T_>& b,
                     Vect<T_>&       x,
                     bool            fact=true);

    int solveTransLU(Vect<T_>& b,
                     bool      fact=true);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Operator <tt>=</tt>
/// \details Copy matrix <tt>m</tt> to current matrix instance.
    DMatrix & operator=(DMatrix<T_>& m);

/// \brief Operator +=
/// \details Add matrix <tt>m</tt> to current matrix instance.
    DMatrix & operator+=(const DMatrix<T_>& m);

/// \brief Operator -=
/// \details Subtract matrix <tt>m</tt> from current matrix instance.
    DMatrix & operator-=(const DMatrix<T_>& m);

/// \brief Operator <tt>=</tt>
/// \details Assign matrix to identity times <tt>x</tt>
    DMatrix & operator=(const T_& x);

/// \brief Operator <tt>*=</tt>
/// \details Premultiply matrix entries by constant value <tt>x</tt>.
    DMatrix & operator*=(const T_& x);

/// \brief Operator <tt>+=</tt>
/// \details Add constant value <tt>x</tt> to matrix entries
    DMatrix & operator+=(const T_& x);

/// \brief Operator <tt>-=</tt>
/// \details Subtract constant value <tt>x</tt> from matrix entries.
    DMatrix & operator-=(const T_& x);

/// \brief Return matrix as C-Array.
/// \details %Matrix is stored row by row.
    T_ *getArray() const;

/// \brief Return entry <tt>(i,j)</tt> of matrix
    T_ get(size_t i,
           size_t j) const;

 private:
    vector<T_> _qr_c, _qr_d;
    int        _qr_set;
};

///////////////////////////////////////////////////////////////////////////////
//                           ASSOCIATED  FUNCTIONS                           //
///////////////////////////////////////////////////////////////////////////////


/** \fn Vect<T_> operator*(const DMatrix<T_> &A, const Vect<T_> &b)
 *  \brief Operator * (Multiply vector by matrix and return resulting vector
 *  \ingroup VectMat
 *  @param [in] A DMatrix instance to multiply by vector
 *  @param [in] b Vect instance 
 *  \return Vect instance containing <tt>A*b</tt>
 */
template<class T_>
Vect<T_> operator*(const DMatrix<T_>& A,
                   const Vect<T_>&    b);

/// \fn ostream& operator<<(ostream& s, const DMatrix<T_>& a)
/// \ingroup VectMat
/// \brief Output matrix in output stream
template<class T_>
ostream& operator<<(ostream&           s,
                    const DMatrix<T_>& a);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
