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

          Definition of class 'SkSMatrix' for Symmetric Skyline Matrix

  ==============================================================================*/


#ifndef __SKSMATRIX_H
#define __SKSMATRIX_H

#include "linear_algebra/Matrix.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file SkSMatrix.h
 *  \brief Definition file for class SkSMatrix.
 */

/*! \class SkSMatrix
 *  \ingroup VectMat
 *  \brief To handle symmetric matrices in skyline storage format.
 *
 * \details This template class allows storing and manipulating a symmetric matrix in skyline storage format.
 *
 * The matrix entries are stored column by column as in the following example:
 * 
 * @verbatim
       /                       \ 
       | a0   a1    0   0   a7 |
       |      a2   a3   0   a8 |
       | ...       a4  a5   a9 |
       |               a6  a10 |
       |                   a11 |
       \                       /
   @endverbatim
 *
 * \tparam T_ Data type (double, float, complex<double>, ...)
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */

template<class T_> class SkSMatrix;
class Mesh;

template<class T_>
class SkSMatrix : public Matrix<T_>
{
   using Matrix<T_>::_nb_rows;
   using Matrix<T_>::_nb_cols;
   using Matrix<T_>::_size;
   using Matrix<T_>::_length;
   using Matrix<T_>::_zero;
   using Matrix<T_>::_dof_type;
   using Matrix<T_>::_temp;
   using Matrix<T_>::_a;
   using Matrix<T_>::_diag;
   using Matrix<T_>::_ch;
   using Matrix<T_>::_is_diagonal;
   using Matrix<T_>::_theMesh;

public:

    using Matrix<T_>::operator();

/// \brief Default constructor.
/// \details Initializes a zero-dimension matrix
    SkSMatrix();

/** \brief Constructor that initializes a dense symmetric matrix.
 *  \details Normally, for a dense matrix this is not the right class.
 *  @param [in] size Number of matrix rows (and columns).
 *  @param [in] is_diagonal Boolean to select if the matrix is diagonal or not [Default: false]
 */
    SkSMatrix(size_t size,
              int    is_diagonal=false);

/** \brief Constructor using mesh to initialize skyline structure of matrix.
 *  @param [in] mesh Mesh instance for which matrix graph is determined.
 *  @param [in] dof Option parameter, with default value <tt>0</tt>.\n
 *  <tt>dof=1</tt> means that only one degree of freedom for each node (or element or side)
 *  is taken to determine matrix structure. The value <tt>dof=0</tt> means that matrix
 *  structure is determined using all DOFs.
 *  @param [in] is_diagonal Boolean argument to say is the matrix is actually a diagonal matrix or not.
 */
    SkSMatrix(Mesh&  mesh,
              size_t dof=0,
              int    is_diagonal=false);

/** \brief Constructor that initializes skyline structure of matrix using vector of column height.
 *  @param [in] ColHt Vect instance that contains rows lengths of matrix.
 */
    SkSMatrix(const Vect<size_t>& ColHt);

/** \brief Constructor for a square matrix using non zero row and column indices
 *  @param [in] I Vector containing row indices
 *  @param [in] J Vector containing column indices
 *  @param [in] opt Flag indicating if vectors I and J are cleaned and ordered
 *  (opt=1) or not (opt=0).\n In the latter case, these vectors can contain
 *  the same contents more than once and are not necessarily ordered.
 */
    SkSMatrix(const Vect<size_t>& I,
              const Vect<size_t>& J,
              int                 opt=1);

/** \brief Constructor for a square matrix using non zero row and column indices
 *  @param [in] I Vector containing row indices
 *  @param [in] J Vector containing column indices
 *  @param [in] a Vector containing matrix entries in the same order than
 *  the one given by <tt>I</tt> and <tt>J</tt>
 *  @param [in] opt Flag indicating if vectors <tt>I</tt> and <tt>J</tt> are cleaned and ordered
 *  (<tt>opt=1</tt>) or not (<tt>opt=0</tt>).\n In the latter case, these vectors can contain
 *  the same contents more than once and are not necessarily ordered              
 */
    SkSMatrix(const Vect<size_t>& I,
              const Vect<size_t>& J,
              const Vect<T_>&     a,
              int                 opt=1);

/// \brief Copy Constructor
    SkSMatrix(const SkSMatrix<T_>& m);

/// \brief Destructor
    ~SkSMatrix();

/** \brief Determine mesh graph and initialize matrix.
 *  \details This member function is called by constructor with the same arguments
 *  @param [in] mesh Mesh instance for which matrix graph is determined.
 *  @param [in] dof Option parameter, with default value <tt>0</tt>.\n
 *  <tt>dof=1</tt> means that only one degree of freedom for each node (or element or side)
 *  is taken to determine matrix structure. The value <tt>dof=0</tt> means that matrix
 *  structure is determined using all DOFs.
 */
    void setMesh(Mesh&  mesh,
                 size_t dof=0);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setGraph(const Vect<RC>& I,
                  int             opt=1);

/** \brief Determine mesh graph and initialize matrix.
 *  \details This member function is called by constructor with the same arguments
 *  @param [in] mesh Mesh instance for which matrix graph is determined.
 *  @param [in] dof_type Type of support of dof. To be given among the enumerated
 *  values: <tt>NODE_DOF</tt>, <tt>ELEMENT_DOF</tt>, <tt>SIDE_DOF</tt> 
 *  @param [in] dof1 Label of first degree of freedom for which numbering is performed.
 *  @param [in] dof2 Label of second degree of freedom for which numbering is performed.
 */
    void setGraph(Mesh&  mesh,
                  int    dof_type,
                  size_t dof1,
                  size_t dof2);

    void setMesh(size_t dof, 
                 Mesh&  mesh,
                 int    code=0);

    void setMesh(size_t dof,
                 size_t nb_eq,
                 Mesh&  mesh);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Determine matrix structure.
/// \details This member function calculates matrix structure using Mesh instance <tt>mesh</tt>.
    void setSkyline(Mesh& mesh);

/// \brief Store diagonal entries in a separate internal vector
    void setDiag();

/** \brief Assign a value to an entry ofthe matrix.
 *  @param [in] i Row index
 *  @param [in] j Column index
 *  @param [in] val Value to assign to <tt>a(i,j)</tt>
 */
    void set(size_t    i,
             size_t    j,
             const T_& val);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void SSet(size_t    i,
              size_t    j,
              const T_& val);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                   a,
              const SkSMatrix<T_>& m);

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m Pointer to %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                a,
              const Matrix<T_>* m);

/** \brief Multiply matrix by vector <tt>x</tt> and add to <tt>y</tt>.
 *  @param [in] x Vector to multiply by matrix
 *  @param [in,out] y Vector to add to the result. <tt>y</tt> contains on output the result.
 */
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

/** \brief Multiply matrix by vector <tt>x</tt> and save in <tt>y</tt>
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector that contains on output the result.
 */
    void Mult(const Vect<T_>& x,
              Vect<T_>&       y) const;

/** \brief Multiply transpose of matrix by vector x and save in y
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector that contains on output the result.
 */
    void TMult(const Vect<T_>& x,
               Vect<T_>&       y) const;

/** \brief Add a constant to an entry of the matrix
 *  @param [in] i Row index
 *  @param [in] j Column index
 *  @param [in] val Constant value to add to <tt>a(i,j)</tt>
 */
    void add(size_t    i,
             size_t    j,
             const T_& val);

/// \brief Return column height.
/// \details Column height at entry <tt>i</tt> is returned.
    size_t getColHeight(size_t i) const;

/// \brief Get <tt>j</tt>-th column vector.
    Vect<T_> getColumn(size_t j) const;

/// \brief Get <tt>i</tt>-th row vector.
    Vect<T_> getRow(size_t i) const;

/** \brief Operator () (Non constant version).
 *  @param [in] i Row index
 *  @param [in] j Column index
 *  @warning To modify a value of an entry of the matrix it is safer not to
 *  modify both lower and upper triangles. Otherwise, wrong values will be 
 *  assigned. If not sure, use the member functions set or add.
 */
    T_ & operator()(size_t i,
                    size_t j);

/** \brief Operator () (Constant version).
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    T_ operator()(size_t i,
                  size_t j) const;

/// \brief Operator =.
/// \details Copy matrix <tt>m</tt> to current matrix instance.
    SkSMatrix<T_> & operator=(const SkSMatrix<T_>& m);

/** \brief Operator =.
 *  \details define the matrix as a diagonal one with all diagonal entries equal
 *  to <tt>x</tt>.
 */
    SkSMatrix<T_> & operator=(const T_& x);

/// \brief Operator +=.
/// \details Add matrix <tt>m</tt> to current matrix instance.
    SkSMatrix<T_> & operator+=(const SkSMatrix<T_>& m);

/// \brief Operator *=.
/// \details Premultiply matrix entries by constant value <tt>x</tt>.
    SkSMatrix<T_> & operator*=(const T_& x);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int Factor();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Factorize matrix (LDLt (Crout) factorization).
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if factorization was normally performed
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null
 *  </ul>
 */
    int setLDLt();

/** \brief Solve a linear system using the LDLt (Crout) factorization
 *  \details This function solves a linear system. The LDLt factorization is 
 *  performed if this was not already done using the function setLU.
 *  @param [in] b Vect instance that contains right-hand side
 *  @param [out] x Vect instance that contains solution
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null
 *  </ul>
 *  Solution is performed only is factorization has previouly been invoked.
 */
    int solveLDLt(const Vect<T_>& b,
                  Vect<T_>&       x);

/** \brief Solve linear system.
 *  \details The linear system having the current instance as a matrix is solved by using the LDLt decomposition.
 *  Solution is thus realized after a factorization step and a forward/backward substitution step.
 *  The factorization step is realized only if this was not already done.\n
 *  Note that this function modifies the matrix contents is a factorization is performed. Naturally, if the
 *  the matrix has been modified after using this function, the user has to refactorize it using the function setLU.
 *  This is because the class has no non-expensive way to detect if the matrix has been modified.
 *  The function setLDLt realizes the factorization step only.
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

/** \brief Solve linear system.
 *  \details The linear system having the current instance as a matrix is solved by using the LDLt decomposition.
 *  Solution is thus realized after a factorization step and a forward/backward substitution step.
 *  The factorization step is realized only if this was not already done.\n
 *  Note that this function modifies the matrix contents is a factorization is performed. Naturally, if the
 *  the matrix has been modified after using this function, the user has to refactorize it using the function setLDLt.
 *  This is because the class has no non-expensive way to detect if the matrix has been modified.
 *  The function setLDLt realizes the factorization step only.
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

/// \brief Return C-Array
/// \details Skyline of matrix is stored row by row.
    T_ *get() const;

/// \brief Assign a value to the i-th entry of C-array containing matrix
    void set(size_t i,
             T_     x);

/// \brief Return entry <tt>(i,j)</tt> of matrix if this one is stored, <tt>0</tt> else
    T_ get(size_t i,
           size_t j) const;

 private:

   int _dof;
};


///////////////////////////////////////////////////////////////////////////////
//                           ASSOCIATED  FUNCTIONS                           //
///////////////////////////////////////////////////////////////////////////////


/** \fn Vect<T_> operator*(const SkSMatrix<T_> &A, const Vect<T_> &b)
 *  \brief Operator * (Multiply vector by matrix and return resulting vector
 *  \ingroup VectMat
 *  @param [in] A SkSMatrix instance to multiply by vector
 *  @param [in] b Vect instance 
 *  \return Vect instance containing <tt>A*b</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_>
Vect<T_> operator*(const SkSMatrix<T_>& A,
                   const Vect<T_>&      b);

/** \fn ostream & operator<<(ostream& s, const SkSMatrix<T_> &a)
 *  \ingroup VectMat
 *  \brief Output matrix in output stream
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_>
ostream& operator<<(ostream&             s,
                    const SkSMatrix<T_>& a);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
