/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2018 Rachid Touzani

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

template<class T_> class TrMatrix;
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
   using Matrix<T_>::_fact;
   using Matrix<T_>::_zero;
   using Matrix<T_>::_temp;

/// \brief Default constructor.
/// \details Initialize a zero dimension tridiagonal matrix
    TrMatrix() { _size = 0; }

/// \brief Constructor for a tridiagonal matrix with <tt>size</tt> rows.
    TrMatrix(size_t size)
    {
       _size = _nb_rows = _nb_cols = size;
       _length = 3*_size;
       _a.resize(_length,static_cast<T_>(0));
    }

/// \brief Copy Constructor
    TrMatrix(const TrMatrix& m)
    {
       _size = m._size;
       _length = m._length;
       _a.resize(_length);
       _a = m._a;
    }

/// \brief Destructor
    ~TrMatrix() { }

/// \brief Define matrix as identity matrix
    void Identity() { Diagonal(1); }

/// \brief Define matrix as a diagonal one
    void Diagonal()
    {
       for (size_t i=1; i<=_size; i++)
          _a[3*i-2] = static_cast<T_>(0);
    }

/// \brief Define matrix as a diagonal one
/// with diagonal entries equal to <tt>a</tt>
    void Diagonal(const T_& a)
    {
       _length = _nb_rows;
       for (size_t i=1; i<=_nb_rows; ++i) {
          _a[3*i-3] = static_cast<T_>(0);
          _a[3*i-2] = a;
          _a[3*i-1] = static_cast<T_>(0);
       }
    }

/** \brief Sets the matrix as the one for the Laplace equation in 1-D
 *  \details The matrix is initialized as the one resulting from P<sub>1</sub>
 *  finite element discretization of the classical elliptic operator
 *     -u'' = f
 *  with homogeneous Dirichlet boundary conditions
 *  \remark This function is available for real valued matrices only.
 *  @param [in] h %Mesh size (assumed constant)
 */
    void Laplace1D(real_t h);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setMesh(Mesh&  mesh,
                 size_t dof=0)
    {
       throw OFELIException("In TrMatrix::setMesh(Mesh,size_t): This member "
                            "function is not valid for class TrMatrix");
    }

    void setMesh(Mesh&  mesh,
                 size_t dof,
                 size_t type)
    {
       throw OFELIException("In TrMatrix::setMesh(Mesh,size_t,size_t): "
                            "This member function is not valid for class TrMatrix");
    }

    virtual void setMesh(size_t dof,
                         Mesh&  mesh,
                         int    code=0)
    {
       throw OFELIException("In TrMatrix::setMesh(Mesh,Mesh,int): "
                            "This member function is not valid for class TrMatrix");
    }

    virtual void setMesh(size_t dof,
                         size_t nb_eq,
                         Mesh&  mesh)
    {
       dof = 0;
       nb_eq = 0;
       if (mesh.getDim()==0) { }
    }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Set size (number of rows) of matrix.
/// @param [in] size Number of rows and columns.
    void setSize(size_t size)
    {
       _nb_rows = _nb_cols = _size = size;
       _length = 3*_nb_rows;
       _a.resize(_length);
    }


/// \brief Multiply matrix by vector <tt>x</tt> and add result to <tt>y</tt>.
    void MultAdd(const Vect<T_>& x, Vect<T_>& y) const
    {
       y(1) += _a[1]*x(1) + _a[2]*x(2);
       for (size_t i=2; i<_size; ++i)
          y(i) += _a[3*i-3]*x(i-1) + _a[3*i-2]*x(i) + _a[3*i-1]*x(i+1);
       y(_size) += _a[3*_size-3]*x(_size-1) + _a[3*_size-2]*x(_size);
    }

/// \brief Multiply matrix by vector <tt>a*x</tt> and add result to <tt>y</tt>.
    void MultAdd(T_              a,
                 const Vect<T_>& x,
                 Vect<T_>&       y) const
    {
       y(1) += a*(_a[1]*x(1) + _a[2]*x(2));
       for (size_t i=2; i<_size; ++i)
          y(i) += a * (_a[3*i-3]*x(i-1) + _a[3*i-2]*x(i) + _a[3*i-1]*x(i+1));
       y(_size) += a * (_a[3*_size-3]*x(_size-1) + _a[3*_size-2]*x(_size));
    }

/// \brief Multiply matrix by vector <tt>x</tt> and save result in <tt>y</tt>.
    void Mult(const Vect<T_>& x,
              Vect<T_>&       y) const
    {
       y = static_cast<T_>(0);
       MultAdd(x,y);
    }

/// \brief Multiply transpose of matrix by vector <tt>x</tt> and save result in <tt>y</tt>.
    void TMult(const Vect<T_>& x,
               Vect<T_>&       y) const
    {
       y(1) = _a[1]*x(1) + _a[3]*x(2);
       for (size_t i=2; i<_size; ++i)
          y(i) = _a[3*i-4]*x(i-1) + _a[3*i-2]*x(i) + _a[3*i]*x(i+1);
       y(_size) = _a[3*_size-4]*x(_size-1) + _a[3*_size-2]*x(_size);
    }

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                  a,
              const TrMatrix<T_>& m)
    {
       _a += a * m._a;
    }
   
/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                a,
              const Matrix<T_>* m)
    {
       for (size_t i=0; i<_length; i++)
          _a[i] += a * m->_a[i];
    }

/// \brief Assign constant <tt>val</tt> to an entry <tt>(i,j)</tt> of the matrix.
    void set(size_t    i,
             size_t    j,
             const T_& val)
    {
       if (i==j)
          _a[3*i-2] = val;
       else if (i==j+1)
          _a[3*i-3] = val;
       else if (i==j-1)
          _a[3*i-1] = val;
}

/// \brief Add constant <tt>val</tt> value to an entry <tt>(i,j)</tt> of the matrix.
    void add(size_t    i,
             size_t    j,
             const T_& val)
    {
       if (i==j)
          _a[3*i-2] += val;
       else if (i==j+1)
          _a[3*i-3] += val;
       else if (i==j-1)
          _a[3*i-1] += val;
}

/** \brief Operator () (Constant version).
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    T_ operator()(size_t i,
                  size_t j) const
    {
       if (i==j)
          return _a[3*i-2];
       else if (i==j+1)
          return _a[3*i-3];
       else if (i==j-1)
          return _a[3*i-1];
       else
          return _zero;
}

/** \brief Operator () (Non constant version).
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    T_ & operator()(size_t i,
                    size_t j)
    {
       if (i==j)
          return _a[3*i-2];
       else if (i==j+1)
          return _a[3*i-3];
       else if (i==j-1)
          return _a[3*i-1];
       else {
          throw OFELIException("In TrMatrix::Operator(): Index pair (" +
                               itos(int(i)) + "," + itos(int(j)) +
                               ") is not compatible with tridiagonal structure.");
       }
       return _zero;
}

/// \brief Operator =.
/// \details Copy matrix <tt>m</tt> to current matrix instance.
    TrMatrix<T_> & operator=(const TrMatrix<T_>& m)
    {
       _a = m._a;
       return *this;
    }

/// \brief Operator =
/// Assign matrix to identity times <tt>x</tt>.
    TrMatrix<T_> & operator=(const T_& x)
    {
       Clear(_a);
       for (size_t i=1; i<=_size; ++i)
          _a[3*i-2] = x;
       return *this;
    }

/// \brief Operator *=.
/// \details Premultiply matrix entries by constant value <tt>x</tt>.
    TrMatrix<T_> & operator*=(const T_& x)
    {
       for (size_t i=1; i<=_length; i++)
          _a[i-1] *= x;
       return *this;
    }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int Factor()
    {
       throw OFELIException("In TrMatrix::Factor(): Matrix cannot be fatorized separately "
                            "Call Solve(b) directly.");
       return 1;
    }

    int Solve(Vect<T_>& b) { return solve(b); }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Solve a linear system with current matrix (forward and back substitution).
 *  @param [in,out] b Vect instance that contains right-hand side on input and solution on output.
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 * 
 *  \b Warning: %Matrix is modified after this function.
 */
    int solve(Vect<T_>& b)
    {
       int ret = 0;
       T_ p=0;
       for (size_t i=1; i<_size; ++i) {
          ret = int(i);
          if (Abs(_a[3*i-2]) < OFELI_EPSMCH)
             throw OFELIException("In TrMatrix::solve(Vect<T_>): The " + itos(int(i))
                                  + "-th pivot is null."); 
          else
             p = _a[3*i]/_a [3*i-2];
          _a[3*i+1] -= p*_a[3*i-1];
          b.add(i+1,-p*b(i));
       }
       if (Abs(_a[3*_size-2]) < OFELI_EPSMCH)
          throw OFELIException("In TrMatrix::solve(Vect<T_>): The " + itos(_size)
                               + "-th pivot is null.");
       else
          b(_size) /= _a[3*_size-2];
       for (int j=int(_size)-1; j>0; j--)
          b(j) = (b(j)-_a[3*j-1]*b(j+1))/_a[3*j-2];
       return ret;
    }

/** \brief Solve a linear system with current matrix (forward and back substitution).
 *  @param [in] b Vect instance that contains right-hand side.
 *  @param [out] x Vect instance that contains solution.
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 * 
 *  \b Warning: %Matrix is modified after this function.
 */
    int solve(const Vect<T_>& b,
              Vect<T_>&       x)
    {
       x = b;
       return solve(x);
    }

/// \brief Return C-Array.
    T_ *get() const { return &_a[0]; }

/// \brief Return entry <tt>(i,j)</tt> of matrix
    T_ get(size_t i, size_t j) const
    {
       if (i==j)
          return _a[3*i-2];
       else if (i==j+1)
          return _a[3*i-3];
       else if (i==j-1)
          return _a[3*i-1];
       else
          return _zero;
    }

};

///////////////////////////////////////////////////////////////////////////////
//                         I M P L E M E N T A T I O N                       //
///////////////////////////////////////////////////////////////////////////////

template<>
inline void TrMatrix<real_t>::Laplace1D(real_t h)
{
   for (size_t i=1; i<=_nb_rows; ++i) {
      _a[3*i-3] = -1.0/h;
      _a[3*i-2] =  2.0/h;
      _a[3*i-1] = -1.0/h;
   }
}

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
                   const Vect<T_>&     b)
{
#if defined (USE_EIGEN)
   return Vect<T_>(Matrix<T_,Eigen::Dynamic,1>(x)*b);
#else
   size_t n = A.getSize();
   Vect<T_> v(n);
   v(1) = A(1,1)*b(1) + A(1,2)*b(2);
   for (size_t i=2; i<n; ++i)
      v(i) = A(i,i-1)*b(i-1) + A(i,i)*b(i) + A(i,i+1)*b(i+1);
   v(n) = A(n,n-1)*b(n-1) + A(n,n)*b(n);
   return v;
#endif
}


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
                       const TrMatrix<T_>& A)
{
   TrMatrix<T_> v(A);
   for (size_t i=0; i<A.getLength(); ++i)
      v[i] *= a;
   return v;
}


/** \fn ostream& operator<<(ostream& s, const TrMatrix<T_>& a)
 *  \ingroup VectMat
 *  \brief Output matrix in output stream
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_>
ostream& operator<<(ostream&            s,
                    const TrMatrix<T_>& A)
{
   s.setf(ios::right|ios::scientific);
   s << endl;
   for (size_t i=1; i<=A.getNbRows(); i++) {
      s << "\nRow  " << setw(6) << i << endl;
      for (size_t j=1; j<=A.getNbRows(); j++)
         s << "  " << setprecision(8) << std::setfill(' ') << setw(18) << A(i,j);
      s << endl;
   }
   return s;
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
