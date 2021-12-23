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

                     Template Class Matrix for banded matrices

  ==============================================================================*/


#ifndef _B_MATRIX_H
#define _B_MATRIX_H

#include "linear_algebra/Matrix.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file BMatrix.h
 *  \brief Definition file for class BMatrix.
 */

class Mesh;

/*! \class BMatrix
 *  \ingroup VectMat
 * \brief To handle band matrices.
 *
 * \details This class enables storing and manipulating band matrices.
 * The matrix can have different numbers of lower and upper co-diagonals
 *
 * \tparam T_ Data type (double, float, complex<double>, ...)
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */


template<class T_>
class BMatrix : public Matrix<T_>
{

 public:

    using Matrix<T_>::operator();

/// \brief Default constructor.
/// \details Initialize a zero dimension band matrix
    BMatrix();

/** \brief Constructor that for a band matrix with given size and bandwidth.
 *  \details Assign 0 to all matrix entries.
 *  @param [in] size Number of rows and columns
 *  @param [in] ld Number of lower co-diagonals (must be > 0)
 *  @param [in] ud Number of upper co-diagonals (must be > 0)
 */
    BMatrix(size_t size,
            int    ld,
            int    ud);

/// \brief Copy Constructor
    BMatrix(const BMatrix& m);

/// \brief Destructor
    ~BMatrix();

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setGraph(const Vect<RC>& I,
                  int             opt=1);

    void setMesh(Mesh&  mesh,
                 size_t dof=0);

    void setMesh(Mesh&  mesh,
                 size_t dof,
                 size_t type=0);

    void setMesh(size_t dof, 
                 Mesh&  mesh,
                 int    code=0);

    void setMesh(size_t dof,
                 size_t nb_eq,
                 Mesh&  mesh);

    void setDiag();

    int Factor();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Set size (number of rows) and storage of matrix
 *  @param [in] size Number of rows and columns
 *  @param [in] ld Number of lower co-diagonals (must be > 0)
 *  @param [in] ud Number of upper co-diagonals (must be > 0)
 */
    void setSize(size_t size,
                 int    ld,
                 int    ud);

/// \brief Multiply matrix by vector <tt>x</tt> and add result to <tt>y</tt>
    void MultAdd(const Vect<T_>& x,
                 Vect<T_>&       y) const;

/// \brief Multiply matrix by vector <tt>a*x</tt> and add result to <tt>y</tt>
    void MultAdd(T_              a,
                 const Vect<T_>& x,
                 Vect<T_>&       y) const;

/// \brief Multiply matrix by vector <tt>x</tt> and save result in <tt>y</tt>
    void Mult(const Vect<T_>& x,
              Vect<T_>&       y) const;

/// \brief Multiply transpose of matrix by vector <tt>x</tt> and save result in <tt>y</tt>
    void TMult(const Vect<T_>& x,
               Vect<T_>&       y) const;

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] x %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                 a,
              const BMatrix<T_>& x);
   
/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] x %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                a,
              const Matrix<T_>* x);

/// \brief Add constant <tt>val</tt> to an entry <tt>(i,j)</tt> of the matrix.
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
    BMatrix<T_> & operator=(const BMatrix<T_>& m);

/// \brief Operator =
/// Assign matrix to identity times <tt>x</tt>.
    BMatrix<T_> & operator=(const T_& x);

/// \brief Operator *=.
/// \details Premultiply matrix entries by constant value <tt>x</tt>
    BMatrix<T_> & operator*=(const T_& x);

/// \brief Operator +=.
/// \details Add constant <tt>x</tt> to matrix entries.
    BMatrix<T_> & operator+=(const T_& x);

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

/** \brief Solve linear system.
 *  \details The linear system having the current instance as a matrix is solved by using the LU decomposition.
 *  Solution is thus realized after a factorization step and a forward/backward substitution step.
 *  The factorization step is realized only if this was not already done.\n
 *  Note that this function modifies the matrix contents is a factorization is performed. Naturally, if the
 *  the matrix has been modified after using this function, the user has to refactorize it using the function setLU.
 *  This is because the class has no non-expensive way to detect if the matrix has been modified.
 *  The function setLU realizes the factorization step only.
 *  @param [in,out] b Vect instance that contains right-hand side on input and solution on output.
 *  @param [in] fact Unused argument
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 */
    int solve(Vect<T_>& b,
              bool     fact=false);

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
 *  @param [in] fact Unused argument
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 */
    int solve(const Vect<T_>& b,
              Vect<T_>&       x,
              bool            fact=false);

/// \brief Return C-Array.
    T_* get() const;

/// \brief Return entry <tt>(i,j)</tt> of matrix
    T_ get(size_t i, size_t j) const;

 private:
   int                 _ld, _ud;
   vector<vector<T_> > _a;

   using Matrix<T_>::_size;
   using Matrix<T_>::_nb_cols;
   using Matrix<T_>::_nb_rows;
   using Matrix<T_>::_length;
   using Matrix<T_>::_diag;
   using Matrix<T_>::_zero;
   using Matrix<T_>::_temp;
   using Matrix<T_>::_ch;
};

///////////////////////////////////////////////////////////////////////////////
//                           ASSOCIATED  FUNCTIONS                           //
///////////////////////////////////////////////////////////////////////////////
 

/** \fn Vect<T_> operator*(const BMatrix<T_> &A, const Vect<T_> &b)
 *  \brief Operator * (Multiply vector by matrix and return resulting vector
 *  \ingroup VectMat
 *  @param [in] A BMatrix instance to multiply by vector
 *  @param [in] b Vect instance 
 *  \return Vect instance containing <tt>A*b</tt>
 */
template<class T_>
Vect<T_> operator*(const BMatrix<T_>& A,
                   const Vect<T_>&    b);


/** \fn BMatrix<T_> operator*(T_ a, const BMatrix<T_> &A)
 *  \brief Operator * (Premultiplication of matrix by constant)
 *  \ingroup VectMat
 *  @return a*A
 */
template<class T_>
BMatrix<T_> operator*(T_                 a,
                      const BMatrix<T_>& A);


/** \fn ostream& operator<<(ostream& s, const BMatrix<T_> &a)
 *  \ingroup VectMat
 *  \brief Output matrix in output stream
 */
template<class T_>
ostream& operator<<(ostream&           s,
                    const BMatrix<T_>& a);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
