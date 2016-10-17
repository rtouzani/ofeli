/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzani

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

template<class T_> class DSMatrix;
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
 */


template<class T_>
class DSMatrix : public Matrix<T_>
{
   using Matrix<T_>::_nb_rows;
   using Matrix<T_>::_nb_cols;
   using Matrix<T_>::_size;
   using Matrix<T_>::_length;
   using Matrix<T_>::_a;
   using Matrix<T_>::_diag;
   using Matrix<T_>::_zero;
   using Matrix<T_>::_ch;
   using Matrix<T_>::_fact;
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

/// \brief Destructor
    ~DSMatrix() { }

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
    void set(      size_t  i,
                   size_t  j,
             const T_&     val);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setMesh(Mesh&  mesh,
                 size_t dof=0)
       { Matrix<T_>::init_set_mesh(mesh,dof); }

    void setMesh(size_t dof,
                 Mesh&  mesh,
                 int    code=0);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
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
    void setRow(      size_t    i,
                const Vect<T_>& v);

/** \brief Copy a given vector to a prescribed column in the matrix.
 *  @param [in] i column index to be assigned
 *  @param [in] v Vect instance to copy
 */
    void setColumn(      size_t    i,
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
    void add(      size_t  i,
                   size_t  j,
             const T_&     val);

/** \brief Operator <tt>()</tt> (Constant version).
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    T_ operator()(size_t i,
                  size_t j) const;

/** \brief Operator <tt>()</tt> (Non constant version).
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    T_ & operator()(size_t i,
                    size_t j);

/// \brief Operator <tt>=</tt>
/// Copy matrix <tt>m</tt> to current matrix instance.
    DSMatrix<T_> & operator=(const DSMatrix<T_>& m);

/// \brief Operator =
/// Assign matrix to identity times <tt>x</tt>.
    DSMatrix<T_> & operator=(const T_& x);

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
                       Vect<T_>& y) const;

/** \brief Multiply matrix by vector <tt>a*x</tt> and add to <tt>y</tt>.
 *  @param [in] a Constant to multiply by matrix
 *  @param [in] x Vector to multiply by matrix
 *  @param [in,out] y Vector to add to the result. <tt>y</tt> contains on output the result.
 */
    void MultAdd(      T_        a,
                 const Vect<T_>& x,
                       Vect<T_>& y) const;

/// \brief Multiply matrix by vector <tt>x</tt> and save result in <tt>y</tt>.
    void Mult(const Vect<T_>& x,
                    Vect<T_>& y) const;

/** \brief Multiply transpose of matrix by vector <tt>x</tt> and add result in <tt>y</tt>.
 *  @param [in] x Vector to add to <tt>y</tt>
 *  @param [in,out] y on input, vector to add to. On output, result.
 */
    void TMult(const Vect<T_>& x,
                     Vect<T_>& y) const;
   
/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(      T_            a,
              const DSMatrix<T_>& m);
   
/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(      T_          a,
              const Matrix<T_>* m);

/** \brief Solve linear system.
 *  \details The matrix is factorized using the LDLt (Crout) decomposition. If this one is already factorized,
 *  no further factorization is performed. If the matrix has been modified the user has to refactorize it 
 *  using the function setLDLt.
 *  @param [in,out] b Vect instance that contains right-hand side on input and solution on output.
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 */
    int solve(Vect<T_>& b);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int Solve(Vect<T_>& b) { return solve(b); }
    int Solve(const Vect<T_>& b, Vect<T_>& x) { return solve(b,x); }
    int Factor() { return setLDLt(); }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Solve linear system.
 *  \details The matrix is factorized using the LDLt (Crout) decomposition. If this one is already factorized,
 *  no further factorization is performed. If the matrix has been modified the user has to refactorize it 
 *  using the function setLDLt.
 *  @param [in] b Vect instance that contains right-hand side.
 *  @param [out] x Vect instance that contains solution
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 */
    int solve(const Vect<T_>& b,
                    Vect<T_>& x);

/** \brief Solve a linear system using the LDLt (Crout) factorization
 *  \details This function solves a linear system. The LDLt factorization is 
 *  performed if this was not already done using the function setLU.
 *  @param [in] b Vect instance that contains right-hand side
 *  @param [out] x Vect instance that contains solution
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null
 *  </ul>
 */
    int setLDLt(const Vect<T_>& b,
                      Vect<T_>& x);

/** \brief Return matrix as C-Array.
 *  Matrix is stored row by row.
 *  Only lower triangle is stored.
 */
    T_ *getArray() const { return _a; }

/// \brief Return entry <tt>(i,j)</tt> of matrix
    T_ get(size_t i,
           size_t j) const;

 protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   void set_(size_t i,
             size_t j,
             T_     x)
      { _a[(i-1)*(i-2)/2+j-1] = x; }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

};

///////////////////////////////////////////////////////////////////////////////
//                         I M P L E M E N T A T I O N                       //
///////////////////////////////////////////////////////////////////////////////


template<class T_>
DSMatrix<T_>::DSMatrix()
{
   _nb_rows = _length = _size = 0;
   _fact = false;
   _is_diagonal = false;
}


template<class T_>
DSMatrix<T_>::DSMatrix(size_t dim)
{
   _fact = false;
   _is_diagonal = false;
   setSize(dim);
}


template<class T_>
DSMatrix<T_>::DSMatrix(const DSMatrix<T_>& m)
{
   _size = _nb_rows = _nb_cols = m._size;
   _length = m._length;
   _fact = m._fact;
   _is_diagonal = false;
   _a.resize(_length);
   _a = m._a;
}


template<class T_>
void DSMatrix<T_>::setDiag()
{
   for (size_t i=0; i<_size; i++)
      _diag[i] = _a[_nb_cols*i+i];
}


template<class T_>
void DSMatrix<T_>::setSize(size_t dim)
{
   _size = _nb_rows = _nb_cols = dim;
   _length = _size*(_size+1)/2;
   _a.resize(_length);
   _diag.resize(_size);
}


template<class T_>
void DSMatrix<T_>::set(      size_t  i,
                             size_t  j,
                       const T_&     val)
{
   if (i>=j)
      _a[i*(i-1)/2+j-1] = val;
}


template<class T_>
void DSMatrix<T_>::setDiag(const T_& a)
{
   _is_diagonal = true;
   for (size_t i=0; i<_size; i++)
      _diag.set(i+1,a);
   for (size_t i=0; i<_nb_rows; i++)
      _a[i*(i-1)/2+i-1] = a;
}


template<class T_>
void DSMatrix<T_>::setDiag(const vector<T_>& d)
{
   _is_diagonal = true;
   for (size_t i=0; i<_size; i++)
      _diag.set(i+1,d[i]);
   for (size_t i=0; i<_nb_rows; i++)
      _a[i*(i-1)/2+i-1] = d[i];
}


template<class T_>
void DSMatrix<T_>::add(      size_t  i,
                             size_t  j,
                       const T_&     val)
{
   if (i>=j)
      _a[i*(i-1)/2+j-1] += val;
}


template<class T_>
T_ DSMatrix<T_>::operator()(size_t i,
                            size_t j) const
{
   if (i<j)
      return _a[j*(j-1)/2+i-1];
   else
      return _a[i*(i-1)/2+j-1];
}


template<class T_>
T_ & DSMatrix<T_>::operator()(size_t i,
                              size_t j)
{
   if (i>=j)
      return _a[i*(i-1)/2+j-1];
   else
      return _a[j*(j-1)/2+i-1];
}


template<class T_>
DSMatrix<T_> & DSMatrix<T_>::operator=(const DSMatrix<T_>& m)
{
   _a = m._a;
   return *this;
}


template<class T_>
DSMatrix<T_> & DSMatrix<T_>::operator=(const T_& x)
{
   _fact = false;
   Clear(_a);
   for (size_t i=1; i<=_size; i++) {
      _diag[i-1] = x;
      set(i,i,x);
   }
   return *this;
}


template<class T_>
int DSMatrix<T_>::setLDLt()
{
   int err = 0;
   T_ s, pivot;
   for (size_t i=0; i<_size; i++) {
      for (size_t j=0; j<i; j++)
         for (size_t k=0; k<j; k++)
            _a[(i+1)*i/2+j] -= _a[(i+1)*i/2+k]*_a[(j+1)*j/2+k];
      pivot = _a[(i+1)*i/2+i];
      for (size_t j=0; j<i; j++) {
         s = _a[(i+1)*i/2+j]*_a[(j+1)*j/2+j];
         pivot -= s*_a[(i+1)*i/2+j];
         _a[(i+1)*i/2+j] = s;
      }
      try {
         if (Abs(pivot) < OFELI_EPSMCH)
            THROW_RT("setLDLt(): The " + itos(int(i)+1) + "-th pivot is null.");
      }
      CATCH_EXIT("DSMatrix");
      if (pivot<0)
         err = -1;
      _a[(i+1)*i/2+i] = 1./pivot;
   }
   _fact = true;
   return err;
}


template<class T_>
void DSMatrix<T_>::getColumn(size_t    j,
                             Vect<T_>& v) const
{
   v.setSize(_size);
   for (size_t i=0; i<j-1; i++)
      v[i] = _a[_size*(j-1)/2+i];
   for (size_t i=j-1; i<_size; i++)
      v[i] = _a[_size*i/2+j-1];
}


template<class T_>
Vect<T_> DSMatrix<T_>::getColumn(size_t j) const
{
   Vect<T_> v(_size);
   for (size_t i=0; i<j-1; i++)
      v[i] = _a[_size*(j-1)/2+i];
   for (size_t i=j-1; i<_size; i++)
      v[i] = _a[_size*i/2+j-1];
   return v;
}


template<class T_>
void DSMatrix<T_>::setColumn(size_t          j,
                             const Vect<T_>& v)
{
   for (size_t i=0; i<j-1; i++)
      _a[_size*(j-1)/2+i] = v[i];
   for (size_t i=j-1; i<_size; i++)
      _a[_size*i/2+j-1] = v[i];
}


template<class T_>
void DSMatrix<T_>::getRow(size_t    i,
                          Vect<T_>& v) const
{
   v.resize(_size);
   for (size_t j=0; j<i; j++)
      v[j] = _a[_size*(i-1)/2+j];
   for (size_t j=i; j<_nb_cols; j++)
      v[j] = _a[_size*j/2+i-1];
}


template<class T_>
Vect<T_> DSMatrix<T_>::getRow(size_t i) const
{
   Vect<T_> v(_size);
   for (size_t j=0; j<i; j++)
      v[j] = _a[_size*(i-1)/2+j];
   for (size_t j=i; j<_size; j++)
      v[j] = _a[_size*j/2+i-1];
   return v;
}


template<class T_>
void DSMatrix<T_>::setRow(size_t          i,
                          const Vect<T_>& v)
{
   for (size_t j=0; j<i; j++)
      _a[_size*(i-1)/2+j] = v[j];
   for (size_t j=i; j<_size; j++)
      _a[_size*j/2+i-1] = v[j];
}


template<class T_>
void DSMatrix<T_>::MultAdd(const Vect<T_>& x,
                           Vect<T_>&       y) const
{
   for (size_t i=1; i<=_nb_rows; i++) {
      for (size_t j=1; j<=i; j++)
         y.add(i,_a[(i-1)*i/2+j-1]*x(j));
      for (size_t k=i+1; k<=_nb_rows; k++)
         y.add(i,_a[(k-1)*k/2+i-1]*x(k));
   }
}


template<class T_>
void DSMatrix<T_>::MultAdd(T_              a,
                           const Vect<T_>& x,
                           Vect<T_>&       y) const
{
   for (size_t i=1; i<=_nb_rows; i++) {
      for (size_t j=1; j<=i; j++)
         y.add(i,a*_a[(i-1)*i/2+j-1]*x(j));
      for (size_t k=i+1; k<=_nb_rows; k++)
         y.add(i,a*_a[(k-1)*k/2+i-1]*x(k));
   }
}


template<class T_>
void DSMatrix<T_>::Mult(const Vect<T_>& x,
                        Vect<T_>&       y) const
{
   y = 0;
   MultAdd(x,y);
}


template<class T_>
void DSMatrix<T_>::TMult(const Vect<T_>& x,
                         Vect<T_>&       y) const
{
   for (size_t i=1; i<=_nb_rows; i++)
      for (size_t j=1; j<=_nb_cols; j++)
         y.add(i,(*this)(i,j)*x(j));
}


#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
void DSMatrix<T_>::setMesh(size_t dof,
                           Mesh&  mesh,
                           int    code)
{
// This is just to avoid warning on unused variable
   dof  = 0;
   code = 0;
   _theMesh = &mesh;
   if (mesh.getDim()==0) { }
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
void DSMatrix<T_>::setMesh(size_t dof,
                           size_t nb_eq,
                           Mesh&  mesh)
{
// This is just to avoid warning on unused variable
   dof = 0;
   nb_eq = 0;
   if (mesh.getDim()==0) { }
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


template<class T_>
int DSMatrix<T_>::solve(Vect<T_>& b)
{
   int ret = 0;
   if (!_fact)
      ret = setLDLt();
   int i, j;
   for (i=0; i<int(_size); i++) {
      T_ s = 0;
      for (j=0; j<i; j++)
         s += _a[(i+1)*i/2+j] * b[j];
      b.add(i+1,-s);
   }
   for (i=0; i<int(_size); i++)
      b.set(i+1,b[i]*_a[(i+1)*i/2+i]);
   for (i=int(_size)-1; i>-1; i--)
      for (j=0; j<i; j++)
         b.add(j+1,-b[i]*_a[(i+1)*i/2+j]);
   return ret;
}


template<class T_>
int DSMatrix<T_>::solve(const Vect<T_>& b,
                        Vect<T_>&       x)
{
   x = b;
   return solve(x);
}


template<class T_>
T_ DSMatrix<T_>::get(size_t i,
                     size_t j) const
{
   if (i>=j)
      return _a[i*(i-1)/2+j-1];
   else
      return _a[j*(j-1)/2+i-1];
}
   
   
template<class T_>
void DSMatrix<T_>::Axpy(T_                  a,
                        const DSMatrix<T_>& m)
{
   _a += a * m._a;
}
   
   
template<class T_>
void DSMatrix<T_>::Axpy(T_                a,
                        const Matrix<T_>* m)
{
   for (size_t i=0; i<_length; i++)
      _a[i] += a * m->_a[i];
}


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
                   const Vect<T_>&     b)
{
   Vect<T_> v(A.getNbRows());
   A.Mult(b,v);
   return v;
}

/// \fn ostream& operator<<(ostream& s, const DSMatrix<T_> &a)
/// \ingroup VectMat
/// \brief Output matrix in output stream
template<class T_>
ostream& operator<<(ostream&            s,
                    const DSMatrix<T_>& a)
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
