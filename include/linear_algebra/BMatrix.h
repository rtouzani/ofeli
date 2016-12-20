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

template<class T_> class BMatrix;
class Mesh;

/*! \class BMatrix
 *  \ingroup VectMat
 * \brief To handle band matrices.
 *
 * \details This class enables storing and manipulating band matrices.
 * The matrix can have different numbers of lower and upper co-diagonals
 *
 * \tparam T_ Data type (double, float, complex<double>, ...)
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
    ~BMatrix() { }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setMesh(Mesh&  mesh,
                 size_t dof=0);
    void setMesh(Mesh&  mesh,
                 size_t dof=0,
                 size_t type=0);
    void setMesh(size_t dof, 
                 Mesh&  mesh,
                 int    code=0);
    void setMesh(size_t dof,
                 size_t nb_eq,
                 Mesh&  mesh);
    void setDiag() { }
    int Factor() { return setLU(); }
    int Solve(Vect<T_>& b) { return solve(b); }
    int Solve(const Vect<T_>& b, Vect<T_>& x) { return solve(b,x); }
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
                       Vect<T_>& y) const;

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
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 */
    int solve(Vect<T_>& b);

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
                    Vect<T_>& x);

/// \brief Return C-Array.
    T_* get() const { return _a; }

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
   using Matrix<T_>::_fact;
   using Matrix<T_>::_zero;
   using Matrix<T_>::_temp;
   using Matrix<T_>::_ch;
};

///////////////////////////////////////////////////////////////////////////////
//                         I M P L E M E N T A T I O N                       //
///////////////////////////////////////////////////////////////////////////////

template<class T_>
BMatrix<T_>::BMatrix()
{
   _size = _ld = _ud = 0;
   _fact = false;
}


template<class T_>
BMatrix<T_>::BMatrix(size_t size,
                     int    ld,
                     int    ud)
{
   try {
       if (size<3 || ld<0 || ud<0 || ud+ld+1>int(size))
          THROW_RT("BMatrix(size_t,int,int): Illegal arguments.");
   }
   CATCH("BMatrix");
   _fact = false;
   setSize(size,ld,ud);
}


template<class T_>
BMatrix<T_>::BMatrix(const BMatrix& m)
{
   _size = m._size;
   _ld = m._ld;
   _ud = m._ud;
   _length = m._length;
   _a = m._a;
   _fact = m._fact;
}


#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
void BMatrix<T_>::setMesh(Mesh&  mesh,
                          size_t dof)
{
   try {
      THROW_RT("setMesh(Mesh,size_t): This member function is not valid for class BMatrix");
   }
   CATCH("BMatrix");
}


template<class T_>
void BMatrix<T_>::setMesh(Mesh&  mesh,
                          size_t dof,
                          size_t type)
{
   try {
      THROW_RT("setMesh(Mesh,size_t,size_t): This member function is not valid for class BMatrix");
   }
   CATCH("BMatrix");
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
void BMatrix<T_>::setMesh(size_t dof,
                          Mesh&  mesh,
                          int    code)
{
   try {
      THROW_RT("setMesh(size_t,Mesh,int): This member function is not valid for class BMatrix");
   }
   CATCH("BMatrix");
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
void BMatrix<T_>::setMesh(size_t dof,
                          size_t nb_eq,
                          Mesh&  mesh)
{
   try {
      THROW_RT("setMesh(size_t,size_t,Mesh): This member function is not valid for class BMatrix");
   }
   CATCH("BMatrix");
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

template<class T_>
void BMatrix<T_>::setSize(size_t size,
                          int    ld,
                          int    ud)
{
   _nb_rows = _nb_cols = _size = size;
   _length = (ld+ud+1) * _nb_rows;
   _ld = ld; _ud = ud;
   _a.resize(_size);
   for (int i=0; i<int(_size); i++)
      _a[i].resize(_ld+_ud+1,static_cast<T_>(0));
   _ch.resize(_size);
   _ch[0] = 0;
   for (int i=1; i<int(_size); i++)
      _ch[i] = i+1;
   _fact = false;
}


template<class T_>
void BMatrix<T_>::MultAdd(const Vect<T_>& x,
                          Vect<T_>&       y) const
{
  for (int i=0; i<int(_size); ++i)
     for (int j=std::max(i-_ld,0); j<std::min(i+_ud+_ld,int(_size)); ++j)
        y[i] += _a[i][j+_ld-i]*x[j];
}


template<class T_>
void BMatrix<T_>::MultAdd(T_              a,
                          const Vect<T_>& x,
                          Vect<T_>&       y) const
{
  for (int i=0; i<int(_size); ++i)
     for (int j=std::max(i-_ld,0); j<std::min(i+_ud+_ld,int(_size)); ++j)
        y[i] += a*_a[i][j+_ld-i]*x[j];
}


template<class T_>
void BMatrix<T_>::Mult(const Vect<T_>& x,
                       Vect<T_>&       y) const
{
   y = static_cast<T_>(0);
   MultAdd(x,y);
}


template<class T_>
void BMatrix<T_>::TMult(const Vect<T_>& x,
                        Vect<T_>&       y) const
{
}


template<class T_>
void BMatrix<T_>::Axpy(T_                 a,
                       const BMatrix<T_>& x)
{
}


template<class T_>
void BMatrix<T_>::Axpy(T_                a,
                       const Matrix<T_>* x)
{
}


template<class T_>
void BMatrix<T_>::set(size_t    i,
                      size_t    j,
                      const T_& val)
{
   if (_ld+int(j)-int(i)>=0 && _ld+int(j-i)<=_ud+_ld)
      _a[i-1][_ld+j-i] = val;
}


template<class T_>
void BMatrix<T_>::add(size_t    i,
                      size_t    j,
                      const T_& val)
{
   if (_ld+int(j)-int(i)>=0 && _ld+int(j)-int(i)<=_ud+_ld)
      _a[i-1][_ld+j-i] += val;
}


template<class T_>
T_ BMatrix<T_>::operator()(size_t i,
                           size_t j) const
{
   if (_ld+int(j)-int(i)>=0 && _ld+int(j)-int(i)<=_ud+_ld)
      return _a[i-1][_ld+j-i];
   else
      return _zero;
}


template<class T_>
T_ & BMatrix<T_>::operator()(size_t i,
                             size_t j)
{
   try {
      if (_ld+int(j)-int(i)>=0 && _ld+int(j)-int(i)<=_ud+_ld)
         return _a[i-1][_ld+j-i];
      else
         THROW_RT("operator(): Illegal pair of indices "+itos(int(i))+","+itos(int(j))+").");
   }
   CATCH("BMatrix");
   return _zero;
}


template<class T_>
BMatrix<T_> & BMatrix<T_>::operator=(const BMatrix<T_>& m)
{
   _a = m._a;
   return *this;
}


template<class T_>
BMatrix<T_> & BMatrix<T_>::operator=(const T_& x)
{
   for (size_t i=0; i<_size; ++i)
      for (int j=0; j<_ld+_ud+1; ++j)
         _a[i][j] = 0;
   for (size_t i=1; i<=_size; ++i)
      (*this)(i,i) = x;
   return *this;
}


template<class T_>
BMatrix<T_> & BMatrix<T_>::operator*=(const T_& x)
{
   for (size_t i=0; i<_size; ++i)
      for (int j=0; j<_ld+_ud+1; ++j)
         _a[i][j] *= x;
   return *this;
}


template<class T_>
BMatrix<T_> & BMatrix<T_>::operator+=(const T_& x)
{
   for (size_t i=0; i<_size; ++i)
      for (int j=0; j<_ld+_ud+1; ++j)
         _a[i][j] += x;
   return *this;
}


template<class T_>
int BMatrix<T_>::setLU()
{
   for (int i=0; i<int(_size)-1; i++) {
      int kend=std::min(_ld+1,int(_size)-i), kjend=std::min(_ud+1,int(_size)-i);
      try {
         if (Abs(_a[i][_ld]) < OFELI_EPSMCH)
            THROW_RT("Factor(): The " + itos(int(i)+1) + "-th pivot is null.");
      }
      CATCH_EXIT("BMatrix");
      for (int k=1; k!=kend; k++) {
         _a[k+i][_ld-k] /= _a[i][_ld];
         for (int j=1; j!=kjend; j++)
            _a[k+i][j+_ld-k] -= _a[k+i][_ld-k] * _a[i][j+_ld];
      }
   }
   _fact = true;
   return 0;
}


template<class T_>
int BMatrix<T_>::solve(Vect<T_>& b)
{
   int ret = 0;
   try {
      if (_size<3 || _ld<0 || _ud<0 || _ud+_ld+1>int(_size))
         THROW_RT("solve(Vect<double>): Illegal arguments.");
   }
   CATCH("BMatrix");
   if (_fact == false)
      ret = setLU();

   for (int i=0; i<int(_size)-1; i++) {
      int kend = std::min(_ld+1,int(_size)-i);
      for (int k=1; k!=kend; k++)
         b[k+i] -= _a[k+i][_ld-k] * b[i];
   }

   for (int i=int(_size)-1; i>=0; i--) {
      int kend = std::min(_ud+1,int(_size)-i);
      for (int k=1; k<kend; k++)
         b[i] -= _a[i][k+_ld] * b[i+k];
      b[i] /= _a[i][_ld];
  }
  return ret;
}


template<class T_>
int BMatrix<T_>::solve(const Vect<T_>& b,
                       Vect<T_>&       x)
{
   x = b;
   return solve(x);
}


template<class T_>
T_ BMatrix<T_>::get(size_t i,
                    size_t j) const
{
   if (_ld+int(j)-int(i)>=0 && _ld+int(j)-int(i)<=_ld+_ud)
      return _a[i-1][_ld+j-i];
   else
      return _zero;
}

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
                   const Vect<T_>&    b)
{
   Vect<T_> v(b.size());
   A.Mult(b,v);
}


/** \fn BMatrix<T_> operator*(T_ a, const BMatrix<T_> &A)
 *  \brief Operator * (Premultiplication of matrix by constant)
 *  \ingroup VectMat
 *  @return a*A
 */
template<class T_>
BMatrix<T_> operator*(T_                 a,
                      const BMatrix<T_>& A)
{
   BMatrix<T_> v(A);
   for (size_t i=0; i<A.getLength(); ++i)
      v[i] *= a;
   return v;
}


/** \fn ostream& operator<<(ostream& s, const BMatrix<T_> &a)
 *  \ingroup VectMat
 *  \brief Output matrix in output stream
 */
template<class T_>
ostream& operator<<(ostream&           s,
                    const BMatrix<T_>& a)
{
   s.setf(ios::right|ios::scientific);
   s << endl;
   for (size_t i=1; i<=a.getNbRows(); i++) {
      s << "\nRow  " << setw(6) << i << endl;
      for (size_t j=1; j<=a.getNbColumns(); j++)
         s << "  " << setprecision(8) << std::setfill(' ') << setw(18) << a(i,j);
      s << endl;
   }
   return s;
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
