/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2017 Rachid Touzani

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

#include "linear_algebra/Matrix.h"
#include "util/util.h"


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
   using Matrix<T_>::_fact;
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

/// \brief Destructor.
    ~DMatrix() { }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setMesh(Mesh&  mesh,
                 size_t dof=0);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setMesh(size_t dof,
                 Mesh&  mesh,
                 int    code=0);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
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
             const T_& val)
    { _a[_nb_cols*(i-1)+j-1] = val; }

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
    void setColumn(size_t          i,
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
    int Factor() { return setLU(); }
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
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 */
    int solve(Vect<T_>& b);

/** \brief Solve the transpose linear system.
 *  \details The linear system having the current instance as a transpose matrix is solved by using the LU decomposition.
 *  Solution is thus realized after a factorization step and a forward/backward substitution step.
 *  The factorization step is realized only if this was not already done.\n
 *  Note that this function modifies the matrix contents is a factorization is performed. Naturally, if the
 *  the matrix has been modified after using this function, the user has to refactorize it using the function setLU.
 *  This is because the class has no non-expensive way to detect if the matrix has been modified.
 *  The function setLU realizes the factorization step only.
 *  @param [in,out] b Vect instance that contains right-hand side on input and solution on output.
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 */
    int solveTrans(Vect<T_>& b);

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
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 */
    int solve(const Vect<T_>& b,
              Vect<T_>&       x);

/** \brief Solve the transpose linear system.
 *  \details The linear system having the current instance as a transpose matrix is solved by using the LU decomposition.
 *  Solution is thus realized after a factorization step and a forward/backward substitution step.
 *  The factorization step is realized only if this was not already done.\n
 *  Note that this function modifies the matrix contents is a factorization is performed. Naturally, if the
 *  the matrix has been modified after using this function, the user has to refactorize it using the function setLU.
 *  This is because the class has no non-expensive way to detect if the matrix has been modified.
 *  The function setLU realizes the factorization step only.
 *  @param [in] b Vect instance that contains right-hand side.
 *  @param [out] x Vect instance that contains solution
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 */
    int solveTrans(const Vect<T_>& b,
                   Vect<T_>&       x);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int Solve(Vect<T_>& b) { return solve(b); }
    int Solve(const Vect<T_>& b, Vect<T_>& x) { return solve(b,x); }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int solveLU(const Vect<T_>& b,
                Vect<T_>&       x);
    int solveLU(Vect<T_>& b);
    int solveTransLU(const Vect<T_>& b,
                     Vect<T_>&       x);
    int solveTransLU(Vect<T_>& b);
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
    T_ *getArray() const { return _a; }

/// \brief Return entry <tt>(i,j)</tt> of matrix
    T_ get(size_t i,
           size_t j) const;

 private:
    vector<T_> _qr_c, _qr_d;
    int        _qr_set;
};


///////////////////////////////////////////////////////////////////////////////
//                         I M P L E M E N T A T I O N                       //
///////////////////////////////////////////////////////////////////////////////


#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
DMatrix<T_>::DMatrix() : _qr_set(0)
{
   _length = _size = _nb_rows = _nb_cols = 0;
   _fact = false;
   _is_diagonal = false;
}


template<class T_>
DMatrix<T_>::DMatrix(size_t nr) : _qr_set(0)
{
   setSize(nr);
   _fact = false;
   _is_diagonal = false;
}


template<class T_>
DMatrix<T_>::DMatrix(size_t nr,
                     size_t nc) : _qr_set(0)
{
   setSize(nr,nc);
   _fact = false;
   _diag.resize(nr);
   _is_diagonal = false;
}


template<class T_>
DMatrix<T_>::DMatrix(Vect<T_>& v) : _qr_set(0)
{
   _fact = false;
   _is_diagonal = true;
   _length = v.size();
   _nb_rows = _nb_cols = _size = sqrt(real_t(_length));
   _a = v;
}


template<class T_>
DMatrix<T_>::DMatrix(const DMatrix<T_>& m)
{
   setSize(m._nb_rows,m._nb_cols);
   _a = m._a;
   _diag = m._diag;
   _fact = m._fact;
   _is_diagonal = m._is_diagonal;
   _qr_set = m._qr_set;
}


template<class T_>
void DMatrix<T_>::setMesh(Mesh&  mesh,
                          size_t dof)
{
   Matrix<T_>::init_set_mesh(mesh,dof);
   _length = _size*_size;
   _diag.resize(_size);
   _a.resize(_length);
   _zero = static_cast<T_>(0);
   _fact = false;
}


template<class T_>
void DMatrix<T_>::setMesh(size_t dof,
                          Mesh&  mesh,
                          int    code)
{
// This is just to avoid warning on unused variable
   dof  = 0;
   code = 0;
   _theMesh = &mesh;
   if (mesh.getDim()==0) { }
}


template<class T_>
void DMatrix<T_>::setMesh(size_t dof,
                          size_t nb_eq,
                          Mesh&  mesh)
{
// This is just to avoid warning on unused variable
   dof = 0;
   nb_eq = 0;
   if (mesh.getDim()==0) { }
   _theMesh = &mesh;
}


template<class T_>
void DMatrix<T_>::setDiag()
{
   for (size_t i=0; i<_size; i++)
      _diag.set(i+1,_a[_nb_cols*i+i]);
}


template<class T_>
void DMatrix<T_>::setDiag(const T_& a)
{
   _is_diagonal = true;
   for (size_t i=0; i<_size; i++)
      _diag.set(i+1,a);
   for (size_t i=0; i<_nb_rows; i++)
      _a[_nb_cols*i+i] = a;
}


template<class T_>
void DMatrix<T_>::setDiag(const vector<T_>& d)
{
   _is_diagonal = true;
   for (size_t i=0; i<_size; i++)
      _diag.set(i+1,d[i]);
   for (size_t i=0; i<_nb_rows; i++)
      _a[_nb_cols*i+i] = d[i];
}


template<class T_>
void DMatrix<T_>::setSize(size_t size)
{
   _nb_rows = _nb_cols = _size = size;
   _length = _nb_rows*_nb_cols;
   _a.resize(_length);
   _diag.resize(_size);
}


template<class T_>
void DMatrix<T_>::setSize(size_t nr,
                          size_t nc)
{
   _nb_rows = nr;
   _nb_cols = nc;
   _size = 0;
   if (_nb_rows==_nb_cols)
      _size = _nb_rows;
   _length = lsize_t(_nb_rows*_nb_cols);
   _a.resize(_length);
}


template<class T_>
void DMatrix<T_>::getColumn(size_t    j,
                            Vect<T_>& v) const
{
   v.resize(_nb_rows);
   for (size_t i=0; i<_nb_rows; i++)
      v[i] = _a[_nb_cols*i+j-1];
}


template<class T_>
Vect<T_> DMatrix<T_>::getColumn(size_t j) const
{
   Vect<T_> v(_nb_rows);
   for (size_t i=0; i<_nb_rows; i++)
      v[i] = _a[_nb_cols*i+j-1];
   return v;
}


template<class T_>
void DMatrix<T_>::setColumn(       size_t    j,
                             const Vect<T_>& v)
{
   for (size_t i=0; i<_nb_rows; i++)
      _a[_nb_cols*i+j-1] = v[i];
}


template<class T_>
void DMatrix<T_>::getRow(size_t    i,
                         Vect<T_>& v) const
{
   v.resize(_nb_cols);
   for (size_t j=0; j<_nb_cols; j++)
      v[j] = _a[_nb_cols*(i-1)+j];
}


template<class T_>
Vect<T_> DMatrix<T_>::getRow(size_t i) const
{
   Vect<T_> v(_nb_cols);
   for (size_t j=0; j<_nb_cols; j++)
      v[j] = _a[_nb_cols*(i-1)+j];
   return v;
}


template<class T_>
void DMatrix<T_>::setRow(      size_t    i,
                         const Vect<T_>& v)
{
   for (size_t j=0; j<_nb_cols; j++)
      _a[_nb_cols*(i-1)+j] = v[j];
}


template<class T_>
void DMatrix<T_>::MultAdd(T_              a,
                          const Vect<T_>& x,
                          Vect<T_>&       y) const
{
  for (size_t i=0; i<_nb_rows; i++)
     for (size_t j=0; j<_nb_cols; j++)
        y.add(i+1,a*_a[_nb_cols*i+j]*x[j]);
}


template<class T_>
void DMatrix<T_>::MultAdd(const Vect<T_>& x,
                          Vect<T_>&       y) const
{
   T_ s;
   for (size_t i=0; i<_nb_rows; i++) {
      s = 0;
      for (size_t j=0; j<_nb_cols; j++)
         s += _a[_nb_cols*i+j] * x[j];
      y.set(i+1,s);
   }
}


template<class T_>
void DMatrix<T_>::Mult(const Vect<T_>& x,
                       Vect<T_>&       y) const
{
   y = T_(0);
   MultAdd(x,y);
}


template<class T_>
void DMatrix<T_>::TMult(const Vect<T_>& x,
                              Vect<T_>& y) const
{
   for (size_t i=0; i<_nb_rows; i++)
      for (size_t j=0; j<_nb_cols; j++)
         y.add(i+1,_a[_nb_cols*(j-1)+i-1]*x[j]);
}


template<class T_>
void DMatrix<T_>::add(      size_t  i,
                            size_t  j,
                      const T_&     val)
{
   _a[_nb_cols*(i-1)+j-1] += val;
}


template<class T_>
T_ DMatrix<T_>::operator()(size_t i,
                           size_t j) const
{
#ifdef _BOUNDS
   assert(i>0);
   assert(i<=_nb_rows);
   assert(j>0);
   assert(j<=_nb_cols);
#endif
   return _a[_nb_cols*(i-1)+j-1];
}


template<class T_>
T_ & DMatrix<T_>::operator()(size_t i,
                             size_t j)
{
#ifdef _BOUNDS
   assert(i>0);
   assert(i<=_nb_rows);
   assert(j>0);
   assert(j<=_nb_cols);
#endif
   return _a[_nb_cols*(i-1)+j-1];
}


template<class T_>
int DMatrix<T_>::setLU()
{
   for (size_t i=1; i<_size; i++) {
      for (size_t j=1; j<=i; j++) {
         try {
            if (Abs(_a[_nb_rows*(j-1)+j-1]) < OFELI_EPSMCH)
               THROW_RT("setLU(): The " + itos(int(i)) + "-th pivot is null.");
         }
         CATCH_EXIT("DMatrix");
         _a[_nb_rows*i+j-1] /= _a[_nb_rows*(j-1)+j-1];
         for (size_t k=0; k<j; k++)
            _a[_nb_rows*i+j] -= _a[_nb_rows*i+k]*_a[_nb_rows*k+j];
      }
      for (size_t j=i+1; j<_size; j++)
         for (size_t k=0; k<i; k++)
            _a[_nb_rows*i+j] -= _a[_nb_rows*i+k]*_a[_nb_rows*k+j];
   }
   _fact = true;
   return 0;
}


template<class T_>
int DMatrix<T_>::setTransLU()
{
   for (size_t i=1; i<_size; i++) {
      for (size_t j=1; j<=i; j++) {
         try {
            if (Abs(_a[_nb_rows*(i-1)+i-1]) < OFELI_EPSMCH)
               THROW_RT("setTLU(): The " + itos(int(i)) + "-th pivot is null.");
         }
         CATCH_EXIT("DMatrix");
         _a[_nb_rows*j+i-1] /= _a[_nb_rows*(i-1)+i-1];
         for (size_t k=0; k<j; k++)
            _a[_nb_rows*j+i] -= _a[_nb_rows*k+i]*_a[_nb_rows*j+k];
      }
      for (size_t j=i+1; j<_size; j++)
         for (size_t k=0; k<i; k++)
            _a[_nb_rows*j+i] -= _a[_nb_rows*k+i]*_a[_nb_rows*j+k];
   }
   _fact = true;
   return 0;
}


template<class T_>
int DMatrix<T_>::solveTrans(const Vect<T_>& b,
                                  Vect<T_>& x)
{
   int ret = 0;
   if (_nb_rows != _nb_cols)
      ret = solveTransQR(b,x);
   else
      ret = solveTransLU(b,x);
   return ret;
}


template<class T_>
int DMatrix<T_>::solve(const Vect<T_>& b,
                             Vect<T_>& x)
{
   int ret = 0;
   if (_nb_rows != _nb_cols)
      ret = solveQR(b,x);
   else
      ret = solveLU(b,x);
   return ret;
}


template<class T_>
int DMatrix<T_>::solve(Vect<T_>& b)
{
   int ret = 0;
   Vect<T_> x(b.size());
   ret = solve(b,x);
   b = x;
   return ret;
}


template<class T_>
int DMatrix<T_>::solveTrans(Vect<T_>& b)
{
   int ret = 0;
   Vect<T_> x(b.size());
   ret = solveTrans(b,x);
   b = x;
   return ret;
}


template<class T_>
int DMatrix<T_>::solveLU(Vect<T_>& b)
{
   int ret = 0;
   if (!_fact)
      ret = setLU();
   for (size_t i=0; i<_size; i++) {
      T_ s=0;
      for (size_t j=0; j<i; j++)
         s += _a[_nb_cols*i+j] * b[j];
      b[i] -= s;
   }
   for (int ii=int(_size)-1; ii>-1; ii--) {
      try {
         if (Abs(_a[_nb_cols*ii+ii])<OFELI_EPSMCH)
            THROW_RT("solveLU(Vect<T_>): The " + itos(ii+1) + "-th pivot is null.");
      }
      CATCH_EXIT("DMatrix");
      b[ii] /= _a[_nb_cols*ii+ii];
      for (size_t j=0; j<size_t(ii); j++)
         b[j] -= b[ii] * _a[_nb_cols*j+ii];
   }
   return ret;
}


template<class T_>
int DMatrix<T_>::solveTransLU(Vect<T_>& b)
{
   int ret = 0;
   if (!_fact)
      ret = setTransLU();
   for (size_t i=0; i<_size; i++) {
      T_ s=0;
      for (size_t j=0; j<i; j++)
         s += _a[_nb_cols*j+i] * b[j];
      b[i] -= s;
   }
   for (int ii=int(_size)-1; ii>-1; ii--) {
      try {
         if (Abs(_a[_nb_cols*ii+ii])<OFELI_EPSMCH)
            THROW_RT("solveLU(Vect<real_t>): The " + itos(ii+1) + "-th pivot is null.");
      }
      CATCH_EXIT("DMatrix");
      b[ii] /= _a[_nb_cols*ii+ii];
      for (size_t j=0; j<size_t(ii); j++)
         b[j] -= b[ii] * _a[_nb_cols*ii+j];
   }
   return ret;
}


template<class T_>
int DMatrix<T_>::solveLU(const Vect<T_>& b,
                               Vect<T_>& x)
{
   x = b;
   return solveLU(x);
}


template<class T_>
int DMatrix<T_>::solveTransLU(const Vect<T_>& b,
                                    Vect<T_>& x)
{
   x = b;
   return solveTrans(x);
}


template<class T_>
DMatrix<T_>& DMatrix<T_>::operator=(DMatrix<T_>& m)
{
   _a = m._a;
   return *this;
}


template<class T_>
DMatrix<T_>& DMatrix<T_>::operator+=(const DMatrix<T_>& m)
{
   for (size_t i=0; i<_length; i++)
      _a[i] += m._a[i];
   return *this;
}


template<class T_>
DMatrix<T_>& DMatrix<T_>::operator-=(const DMatrix<T_>& m)
{
   for (size_t i=0; i<_length; i++)
      _a[i] -= -m._a[i];
   return *this;
}


template<class T_>
DMatrix<T_>& DMatrix<T_>::operator=(const T_& x)
{
   try {
      if (_nb_rows!=_nb_cols)
         THROW_RT("operator=(T_): Operator is valid for square matrices only.");
   }
   CATCH_EXIT("DMatrix");
   _fact = false;
   Clear(_a);
   for (size_t i=1; i<=_size; i++) {
      set(i,i,x);
      _diag[i-1] = x;
   }
   return *this;
}


template<class T_>
DMatrix<T_> & DMatrix<T_>::operator*=(const T_& x)
{
   for (size_t i=0; i<_length; i++)
      _a[i] *= x;
   return *this;
}


template<class T_>
DMatrix<T_> & DMatrix<T_>::operator+=(const T_& x)
{
   for (size_t i=0; i<_length; i++)
      _a[i] += x;
   return *this;
}


template<class T_>
DMatrix<T_> & DMatrix<T_>::operator-=(const T_& x)
{
   for (size_t i=0; i<_length; i++)
      _a[i] -= x;
   return *this;
}


template<class T_>
T_ DMatrix<T_>::get(size_t i,
                    size_t j) const 
{
   return _a[_nb_cols*(i-1)+j-1]; 
}


template<class T_>
void DMatrix<T_>::Axpy(      T_           a,
                       const DMatrix<T_>& m)
{
   Axpy(a,m._a,_a);
}


template<class T_>
void DMatrix<T_>::Axpy(      T_          a,
                       const Matrix<T_>* m)
{
   for (size_t i=0; i<_length; i++)
      _a[i] += a * m->_a[i];
}


template<class T_>
int DMatrix<T_>::setQR()
{
   return -1;
}


template<class T_>
int DMatrix<T_>::setTransQR()
{
   return -1;
}


template<>
inline int DMatrix<real_t>::setQR()
{
   size_t n=std::min(_nb_rows,_nb_cols);
   _qr_c.resize(n); Clear(_qr_c);
   _qr_d.resize(n); Clear(_qr_d);
   real_t scale = 0;
   for (size_t k=0; k<n; k++) {
      for (size_t i=k; i<_nb_rows; i++)
         scale = fmax(scale,fabs(_a[_nb_cols*i+k]));
      if (scale==0.0) {
         _qr_c[k] = _qr_d[k] = 0;
         return int(k)+1;
      }
      else {
         for (size_t i=k; i<_nb_rows; i++)
            _a[_nb_cols*i+k] /= scale;
         real_t s = 0;
         for (size_t i=k; i<_nb_rows; i++)
            s += _a[_nb_cols*i+k]*_a[_nb_cols*i+k];
         real_t sigma = sqrt(s)*Sgn(_a[_nb_cols*k+k]);
         _a[_nb_cols*k+k] += sigma;
         _qr_c[k] = sigma*_a[_nb_cols*k+k];
         _qr_d[k] = -scale*sigma;
         for (size_t j=k+1; j<_nb_cols; j++) {
            real_t s = 0;
            for (size_t i=k; i<_nb_rows; i++)
               s += _a[_nb_cols*i+k]*_a[_nb_cols*i+j];
            real_t tau = s/_qr_c[k];
            for (size_t i=k; i<_nb_rows; i++)
               _a[_nb_cols*i+j] -= tau*_a[_nb_cols*i+k];
         }
      }
   }
   if (_qr_d[n-1]==0)
      return int(n);
   _qr_set = 1;
   return 0;
}


template<>
inline int DMatrix<real_t>::setTransQR()
{
   size_t n=std::min(_nb_rows,_nb_cols);
   _qr_c.resize(n); Clear(_qr_c);
   _qr_d.resize(n); Clear(_qr_d);
   real_t scale = 0;
   for (size_t k=0; k<n; k++) {
      for (size_t i=k; i<_nb_rows; i++)
         scale = fmax(scale,fabs(_a[_nb_cols*k+i]));
      if (scale==0.0) {
         _qr_c[k] = _qr_d[k] = 0;
         return int(k)+1;
      }
      else {
         for (size_t i=k; i<_nb_rows; i++)
            _a[_nb_cols*k+i] /= scale;
         real_t s = 0;
         for (size_t i=k; i<_nb_rows; i++)
            s += _a[_nb_cols*k+i]*_a[_nb_cols*k+i];
         real_t sigma = sqrt(s)*Sgn(_a[_nb_cols*k+k]);
         _a[_nb_cols*k+k] += sigma;
         _qr_c[k] = sigma*_a[_nb_cols*k+k];
         _qr_d[k] = -scale*sigma;
         for (size_t j=k+1; j<_nb_cols; j++) {
            real_t s = 0;
            for (size_t i=k; i<_nb_rows; i++)
               s += _a[_nb_cols*k+i]*_a[_nb_cols*j+i];
            real_t tau = s/_qr_c[k];
            for (size_t i=k; i<_nb_rows; i++)
               _a[_nb_cols*j+i] -= tau*_a[_nb_cols*k+i];
         }
      }
   }
   if (_qr_d[n-1]==0)
      return int(n);
   _qr_set = 1;
   return 0;
}


template<class T_>
int DMatrix<T_>::solveQR(const Vect<T_>& b,
                         Vect<T_>&       x)
{
   return -1;
}


template<class T_>
int DMatrix<T_>::solveTransQR(const Vect<T_>& b,
                              Vect<T_>&      x)
{
   return -1;
}


template<>
inline int DMatrix<real_t>::solveQR(const Vect<real_t>& b,
                                    Vect<real_t>&       x)
{
   int ret = 0;
   Vect<real_t> c(b);
   if (_qr_set==0)
      ret = setQR();
   size_t n=std::min(_nb_rows,_nb_cols);
   for (size_t j=0; j<n; j++) {
      real_t s=0;
      for (size_t i=j; i<_nb_rows; i++)
         s += _a[_nb_cols*i+j] * c[i];
      real_t t=s/_qr_c[j];
      for (size_t i=j; i<_nb_rows; i++)
         c[i] -= t * _a[_nb_cols*i+j];
   }

   x[n-1] = c[n-1]/_qr_d[n-1];
   for (int i=int(n)-2; i>=0; i--) {
      real_t s=0;
      for (size_t j=i; j<n; j++)
         s += _a[_nb_cols*i+j]*x[j];
      x[i] = (c[i]-s)/_qr_d[i];
   }
   return ret;
}


template<>
inline int DMatrix<real_t>::solveTransQR(const Vect<real_t>& b,
                                         Vect<real_t>&       x)
{
   int ret = 0;
   Vect<real_t> c(b);
   if (_qr_set==0)
      ret = setQR();
   size_t n=std::min(_nb_rows,_nb_cols);
   for (size_t j=0; j<n; j++) {
      real_t s=0;
      for (size_t i=j; i<_nb_rows; i++)
         s += _a[_nb_cols*j+i] * c[i];
      real_t t=s/_qr_c[j];
      for (size_t i=j; i<_nb_rows; i++)
         c[i] -= t * _a[_nb_cols*j+i];
   }

   x[n-1] = c[n-1]/_qr_d[n-1];
   for (int i=int(n)-2; i>=0; i--) {
      real_t s=0;
      for (size_t j=i; j<n; j++)
         s += _a[_nb_cols*j+i]*x[j];
      x[i] = (c[i]-s)/_qr_d[i];
   }
   return ret;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


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
                   const Vect<T_>&    b)
{
   Vect<T_> v(A.getNbRows());
   A.Mult(b,v);
   return v;
}

/// \fn ostream& operator<<(ostream& s, const DMatrix<T_>& a)
/// \ingroup VectMat
/// \brief Output matrix in output stream
template<class T_>
ostream& operator<<(      ostream&     s,
                    const DMatrix<T_>& a)
{
   s.setf(ios::right|ios::scientific);
   s << endl;
   for (size_t i=1; i<=a.getNbRows(); i++) {
      s << "\nRow  " << setw(6) << i << endl;
      for (size_t j=1; j<=a.getNbColumns(); j++)
         s << "  " << setprecision(8) << std::setfill(' ')
           << setw(18) << a(i,j);
      s << endl;
   }
   return s;
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
