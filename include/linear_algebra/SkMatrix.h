/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2025 Rachid Touzani

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

            Definition of class 'SkMatrix' for general skyline matrix

  ==============================================================================*/


#ifndef __SKMATRIX_H
#define __SKMATRIX_H

#include "linear_algebra/Matrix.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file SkMatrix.h
 *  \brief Definition file for class SkMatrix.
 */

class Mesh;
//template class Vect<size_t>;
 
/*! \class SkMatrix
 *  \ingroup VectMat
 *  \brief To handle square matrices in skyline storage format.
 *
 *  \details This template class allows storing and manipulating a matrix in
 *  skyline storage format.
 *
 * The matrix entries are stored in 2 vectors column by column as in the following 
 * example:
 * 
 * @verbatim
 /                    \        /                       \ 
 | l0            .    |        | u0   u1    0   0   u7 |
 | l1  l2        .    |        |      u2   u3   0   u8 |
 |  0  l3  l4    .    |        | ...       u4  u5   u9 |
 |  0   0  l5  l6     |        |               u6  u10 |
 | l7  l8  l9 l10 l11 |        |                   u11 |
 \                    /        \                       /
 @endverbatim
 * 
 * \tparam T_ Data type (double, float, complex<double>, ...)
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */

template<class T_>
class SkMatrix : public Matrix<T_>
{

 public:

   using Matrix<T_>::_msize;
   using Matrix<T_>::_temp;
   using Matrix<T_>::_a;
   using Matrix<T_>::_aU;
   using Matrix<T_>::_diag;
   using Matrix<T_>::_dof_type;
   using Matrix<T_>::_is_diagonal;
   using Matrix<T_>::_theMesh;
   using Matrix<T_>::operator();

/// \brief Default constructor.
/// \details Initializes a zero-dimension matrix
    SkMatrix();

/** \brief Constructor that initializes a dense symmetric matrix.
 *  \details Normally, for a dense matrix this is not the right class.
 *  @param [in] size Number of matrix rows (and columns).
 *  @param [in] is_diagonal Boolean to select if the matrix is diagonal or not [Default: false]
 */
    SkMatrix(size_t size,
             int    is_diagonal=false);

/** \brief Constructor using mesh to initialize skyline structure of matrix.
 *  @param [in] mesh Mesh instance for which matrix graph is determined.
 *  @param [in] dof Option parameter, with default value <tt>0</tt>.\n
 *  <tt>dof=1</tt> means that only one degree of freedom for each node (or element or side)
 *  is taken to determine matrix structure. The value <tt>dof=0</tt> means that matrix
 *  structure is determined using all DOFs.
 *  @param [in] is_diagonal Boolean argument to say is the matrix is actually a diagonal matrix or not.
 */
    SkMatrix(Mesh&  mesh,
             size_t dof=0,
             int    is_diagonal=false);

/** \brief Constructor that initializes skyline structure of matrix using vector of column heights.
 *  @param [in] ColHt Vect instance that contains rows lengths of matrix.
 */
    SkMatrix(const vector<size_t> &ColHt);

/// \brief Copy Constructor
    SkMatrix(const SkMatrix<T_>& m);

/// \brief Destructor
    ~SkMatrix();

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
    void setGraph(const vector<RC>& I,
                  int               opt=1);

    void setMesh(size_t dof,
                 Mesh&  mesh,
                 int    code=0);

    void setMesh(size_t dof,
                 size_t nb_eq,
                 Mesh&  mesh);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Determine matrix structure.
 *  \details This member function calculates matrix structure using a Mesh instance.
 *  @param [in] mesh Mesh instance
 */
    void setSkyline(Mesh& mesh);

/// \brief Store diagonal entries in a separate internal vector.
    void setDiag();

/** \brief Choose DOF to activate.
 *  \details This function is available only if variable <tt>dof</tt> is equal to 1 in the constructor
 *  @param [in] i Index of the DOF
 */
    void setDOF(size_t i);

/** \brief Assign a value to an entry ofthe matrix.
 *  @param [in] i Row index (starting at <tt>i=1</tt>)
 *  @param [in] j Column index (starting at <tt>i=1</tt>)
 *  @param [in] val Value to assign to entry <tt>a(i,j)</tt>
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
    void Axpy(T_                  a,
              const SkMatrix<T_>& m);

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m %Matrix by which <tt>a</tt> is multiplied. The result is added
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

/** \brief Multiply transpose of matrix by vector <tt>x</tt> and add to <tt>y</tt>.
 *  @param [in] x Vector to multiply by matrix
 *  @param [in,out] y Vector to add to the result. <tt>y</tt> contains on output the result.
 */
    void TMultAdd(const Vect<T_>& x,
                  Vect<T_>&       y) const;

/** \brief Multiply matrix by a vector and add to another one.
 *  @param [in] a Constant to multiply by matrix
 *  @param [in] x Vector to multiply by matrix
 *  @param [in,out] y Vector to add to the result. <tt>y</tt> contains on output the result.
 */
    void MultAdd(T_              a,
                 const Vect<T_>& x,
                 Vect<T_>&       y) const;

/** \brief Multiply matrix by vector <tt>x</tt> and save in <tt>y</tt>.
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector that contains on output the result.
 */
    void Mult(const Vect<T_>& x,
              Vect<T_>&       y) const;

/** \brief Multiply transpose of matrix by vector <tt>x</tt> and save in <tt>y</tt>.
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector that contains on output the result.
 */
    void TMult(const Vect<T_>& x,
               Vect<T_>&       y) const;

/** \brief Add a constant value to an entry ofthe matrix.
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

/** \brief Return a value of a matrix entry
 *  @param [in] i Row index (starts at 1)
 *  @param [in] j Column index (starts at 1)
 */
    T_ at(size_t i,
          size_t j);

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

/** \brief Impose an essential boundary condition.
 *  \details This member function modifies diagonal terms in matrix and terms
 *  in vector that correspond to degrees of freedom with nonzero code
 *  in order to impose a boundary condition.
 *  It can be modified by member function \b setPenal(..).
 *  @param [in] mesh Mesh instance from which information is extracted.
 *  @param [in] b Vect instance that contains right-hand side.
 *  @param [in] u Vect instance that conatins imposed valued at DOFs where they are to be imposed.
 *  @param [in] flag Parameter to determine whether only the right-hand side is to
 *  be modified (<tt>dof>0</tt>)\n
 *  or both matrix and right-hand side (<tt>dof=0</tt>, default value).
 */
    void DiagPrescribe(Mesh&           mesh,
                       Vect<T_>&       b,
                       const Vect<T_>& u,
                       int             flag=0);

/** \brief Impose an essential boundary condition using the Mesh instance provided by the constructor.
 *  \details This member function modifies diagonal terms in matrix and terms
 *  in vector that correspond to degrees of freedom with nonzero code
 *  in order to impose a boundary condition.
 *  It can be modified by member function \b setPenal(..).
 *  @param [in] b Vect instance that contains right-hand side.
 *  @param [in] u Vect instance that conatins imposed valued at DOFs where they are to be imposed.
 *  @param [in] flag Parameter to determine whether only the right-hand side is to
 *  be modified (<tt>dof>0</tt>)\n
 *  or both matrix and right-hand side (<tt>dof=0</tt>, default value).
 */
    void DiagPrescribe(Vect<T_>&       b,
                       const Vect<T_>& u,
                       int             flag=0);

/// \brief Operator =.
/// \details Copy matrix <tt>m</tt> to current matrix instance.
    SkMatrix<T_> & operator=(const SkMatrix<T_>& m);

/** \brief Operator =.
 *  \details define the matrix as a diagonal one with all diagonal entries equal
 *  to <tt>x</tt>.
 */
    SkMatrix<T_> & operator=(const T_& x);

/// \brief Operator +=.
/// \details Add matrix <tt>m</tt> to current matrix instance.
    SkMatrix<T_> & operator+=(const SkMatrix<T_>& m);

/// \brief Operator +=.
/// \details Add constant value <tt>x</tt> to matrix entries.
    SkMatrix<T_> & operator+=(const T_& x);

/// \brief Operator *=.
/// \details Premultiply matrix entries by constant value <tt>x</tt>.
    SkMatrix<T_> & operator*=(const T_& x);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
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

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int solveLU(const Vect<T_>& b,
                Vect<T_>&       x);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Return C-Array.
/// \details Skyline of matrix is stored row by row.
    T_ *get() const;

/// \brief Return entry <tt>(i,j)</tt> of matrix if this one is stored, 0 else
    T_ get(size_t i,
           size_t j) const;

/// \brief Add <tt>val</tt> to entry <tt>i</tt>.
    void add(size_t    i,
             const T_& val);

 private:
    int _dof;
};

///////////////////////////////////////////////////////////////////////////////
//                           ASSOCIATED  FUNCTIONS                           //
///////////////////////////////////////////////////////////////////////////////


/** \fn Vect<T_> operator*(const SkMatrix<T_> &A, const Vect<T_> &b)
 *  \brief Operator * (Multiply vector by matrix and return resulting vector
 *  \ingroup VectMat
 *  @param [in] A SkMatrix instance to multiply by vector
 *  @param [in] b Vect instance 
 *  \return Vect instance containing <tt>A*b</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_>
Vect<T_> operator*(const SkMatrix<T_>& A,
                   const Vect<T_>&     b);


/** \fn ostream & operator<<(ostream &s, const SkMatrix<T_> &a)
 *  \ingroup VectMat
 *  \brief Output matrix in output stream
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_>
ostream& operator<<(ostream&            s,
                    const SkMatrix<T_>& a);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
