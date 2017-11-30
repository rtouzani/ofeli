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

      Definition of Template Class LocalMatrix for dense small size matrices

  ==============================================================================*/


#ifndef __LOCAL_MATRIX_H
#define __LOCAL_MATRIX_H


#include "OFELI_Config.h"
#include "linear_algebra/LocalVect.h"
#include "linear_algebra/SkMatrix.h"
#include "linear_algebra/SkSMatrix.h"
#include "linear_algebra/SpMatrix.h"
#include "OFELIException.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file LocalMatrix.h
 *  \brief Definition file for class LocalMatrix.
 */

/*! \class LocalMatrix
 *  \ingroup VectMat
 *  \brief Handles small size matrices like element matrices, with a priori known size.
 *
 * \details The template class LocalMatrix treats small size matrices. Typically, this class is
 * recommended to store element and side arrays.\n
 * Internally, no dynamic storage is used.
 *
 * \tparam T_  Data type (double, float, complex<double>, ...)
 * \tparam NR_ number of rows of matrix
 * \tparam NC_ number of columns of matrix
 *
 */

template<class T_, size_t NC_> class LocalVect;
class Element;

template<class T_, size_t NR_, size_t NC_> class LocalMatrix
{

 public:

/// \brief Default constructor
/// \details Constructs a matrix with <tt>0</tt> rows and <tt>0</tt> columns
    LocalMatrix();

/// \brief Copy constructor
    LocalMatrix(const LocalMatrix<T_,NR_,NC_>& m);

/** \brief Constructor of a local matrix associated to element from a SpMatrix
 *  @param [in] el Pointer to Element
 *  @param [in] a Global matrix as instance of class SpMatrix.
 */
    LocalMatrix(Element*            el,
                const SpMatrix<T_>& a);

/** \brief Constructor of a local matrix associated to element from a SkMatrix
 *  @param [in] el Pointer to Element
 *  @param [in] a Global matrix as instance of class SkMatrix.
 */
    LocalMatrix(Element*            el,
                const SkMatrix<T_>& a);

/** \brief Constructor of a local matrix associated to element from a SkSMatrix
 *  @param [in] el Pointer to Element
 *  @param [in] a Global matrix as instance of class SkSMatrix.
 */
    LocalMatrix(Element*             el,
                const SkSMatrix<T_>& a);

/// Destructor
    ~LocalMatrix() { }

/// \brief Operator <tt>()</tt> (Non constant version)
/// \details Returns entry at row <tt>i</tt> and column <tt>j</tt>.
    T_& operator()(size_t i,
                   size_t j)
    {
       return _a[(i-1)*NC_+j-1];
    }

/// \brief Operator <tt>()</tt> (Constant version)
/// \details Returns entry at row <tt>i</tt> and column <tt>j</tt>.
    T_ operator()(size_t i,
                  size_t j) const
    {
       return _a[(i-1)*NC_+j-1];
    }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void set(int opt);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Initialize matrix as element matrix from global SpMatrix
 *  @param [in] el Pointer to Element
 *  @param [in] a Global matrix as instance of class SpMatrix.
 *  This function is called by its corresponding constructor.
 */
    void Localize(Element*            el,
                  const SpMatrix<T_>& a);

/** \brief Initialize matrix as element matrix from global SkMatrix
 *  @param [in] el Pointer to Element
 *  @param [in] a Global matrix as instance of class SkMatrix.
 *  This function is called by its corresponding constructor.
 */
    void Localize(Element*            el,
                  const SkMatrix<T_>& a);

/** \brief Initialize matrix as element matrix from global SkSMatrix.
 *  @param [in] el Pointer to Element
 *  @param [in] a Global matrix as instance of class SkSMatrix.
 *  This function is called by its corresponding constructor.
 */
    void Localize(Element*             el,
                  const SkSMatrix<T_>& a);

/// \brief Operator <tt>=</tt>
/// \details Copy instance <tt>m</tt> into current instance.
    LocalMatrix<T_,NR_,NC_> & operator=(const LocalMatrix<T_,NR_,NC_>& m);

/// \brief Operator <tt>=</tt>
/// \details Assign matrix to identity times <tt>x</tt>
    LocalMatrix<T_,NR_,NC_> & operator=(const T_& x);

/// \brief Operator <tt>+=</tt>
/// \details Add <tt>m</tt> to current matrix.
    LocalMatrix<T_,NR_,NC_> & operator+=(const LocalMatrix<T_,NR_,NC_>& m);

/// \brief Operator <tt>-=</tt>
/// \details Subtract <tt>m</tt> from current matrix.
    LocalMatrix<T_,NR_,NC_> & operator-=(const LocalMatrix<T_,NR_,NC_>& m);

/// \brief Operator <tt>*</tt>
/// \details Return a Vect instance as product of current matrix by vector <tt>x</tt>.
    LocalVect<T_,NR_> operator*(LocalVect<T_,NC_>& x);

/// \brief Operator <tt>+=</tt>
/// \details Add constant <tt>x</tt> to current matrix entries.
    LocalMatrix<T_,NR_,NC_> & operator+=(const T_& x);

/// \brief Operator <tt>-=</tt>
/// \details Subtract <tt>x</tt> from current matrix entries.
    LocalMatrix<T_,NR_,NC_> & operator-=(const T_& x);

/// \brief Operator <tt>*=</tt>
/// \details Multiply matrix entries by constant <tt>x</tt>.
    LocalMatrix<T_,NR_,NC_> & operator*=(const T_& x);

/// \brief Operator <tt>/=</tt>
/// \details Divide by <tt>x</tt> current matrix entries.
    LocalMatrix<T_,NR_,NC_> & operator/=(const T_ &x);

/** \brief Multiply matrix by vector and add result to vector.
 *  @param [in] x Vector to multiply matrix by.
 *  @param [out] y Resulting vector (<tt>y += a * x</tt>)
 */
    void MultAdd(const LocalVect<T_,NC_>& x,
                 LocalVect<T_,NR_>&       y);

/** \brief Multiply matrix by scaled vector and add result to vector.
 *  @param [in] a Constant to premultiply by vector <tt>x</tt>.
 *  @param [in] x (Scaled) vector to multiply matrix by.
 *  @param [out] y Resulting vector (<tt>y += a * x</tt>)
 */
    void MultAddScal(const T_&                a,
                     const LocalVect<T_,NC_>& x,
                     LocalVect<T_,NR_>&       y);

/** \brief Multiply matrix by vector.
 *  @param [in] x Vector to multiply matrix by.
 *  @param [out] y Resulting vector.
 */
    void Mult(const LocalVect<T_,NC_>& x,
              LocalVect<T_,NR_>&       y);

/// \brief Symmetrize matrix
/// \details Fill upper triangle to form a symmetric matrix.
    void Symmetrize();

/** \brief Factorize matrix
 *  \details Performs a LU factorization.
 *  @return
 *  <ul>
 *     <li><tt>0</tt>: Factorization has ended normally,
 *     <li><tt>n</tt>: <tt>n</tt>-th pivot was zero.
 *  </ul>
 */
    int Factor();

/** \brief Forward and backsubstitute to solve a linear system.
 *  @param [in] b Right-hand side in input and solution vector in output.
 *  @return
 *  <ul>
 *     <li><tt>0</tt>: Solution was performed normally.
 *     <li><tt>n</tt>: <tt>n</tt>-th pivot is zero.
 *  </ul>
 *  @note %Matrix must have been factorized at first.
 */
    int Solve(LocalVect<T_,NR_>& b);

/** \brief Factorize matrix and solve linear system.
 *  @param [in,out] b Right-hand side in input and solution vector in output.
 *  @return
 *     <tt>0</tt> if solution was performed normally.
 *     <tt>n</tt> if <tt>n</tt>-th pivot is zero.
 *  This function simply calls \b Factor() then \b Solve(b).
 */
    int FactorAndSolve(LocalVect<T_,NR_>& b);

/// \brief Calculate inverse of matrix
/// @param [out] A Inverse of matrix
    void Invert(LocalMatrix<T_,NR_,NC_>& A);

/** \brief Calculate inner product witrh respect to matrix.
 *  \details Returns the product <tt>x<sup>T</sup>Ay</tt> 
 *  @param [in] x Left vector
 *  @param [in] y Right vector
 *  @return Resulting product
 */
    T_ getInnerProduct(const LocalVect<T_,NC_>& x,
                       const LocalVect<T_,NR_>& y);

/// \brief Return pointer to matrix as a C-array.
    T_ *get() { return _a; }

 private:
    size_t _length;
    T_ _a[NR_*NC_];
};


///////////////////////////////////////////////////////////////////////////////
//                         I M P L E M E N T A T I O N                       //
///////////////////////////////////////////////////////////////////////////////

template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_>::LocalMatrix()
{
   _length = NR_ * NC_;
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_>::LocalMatrix(const LocalMatrix<T_,NR_,NC_>& m)
{
   _length = m._length;
   for (size_t k=0; k<_length; k++)
      _a[k] = m._a[k];
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_>::LocalMatrix(Element*            el,
                                     const SpMatrix<T_>& a)
{
   Localize(el,a);
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_>::LocalMatrix(Element*            el,
                                     const SkMatrix<T_>& a)
{
   Localize(el,a);
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_>::LocalMatrix(Element*             el,
                                     const SkSMatrix<T_>& a)
{
   Localize(el,a);
}


#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_,size_t NR_,size_t NC_>
void LocalMatrix<T_,NR_,NC_>::set(int opt)
{
   if (opt==IDENTITY) {
      if (NR_!=NC_)
         throw OFELIException("In LocalMatrix::set(opt): This argument is valid for square matrices only.");
      for (size_t i=0; i<NR_*NC_; i++)
         _a[i] = 0;
      for (size_t j=0; j<NR_; j++)
         _a[j*(NC_+1)] = 1;
   }
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


template<class T_,size_t NR_,size_t NC_>
void LocalMatrix<T_,NR_,NC_>::Localize(Element*            el,
                                       const SpMatrix<T_>& a)
{
   for (size_t n=1; n<=el->getNbNodes(); n++) {
      Node *nd = (*el)(n);
      for (size_t i=1; i<=nd->getNbDOF(); i++)
         for (size_t j=1; j<=nd->getNbDOF(); j++)
            (*this)(i,j) = a(nd->getDOF(i),nd->getDOF(j));
   }
}


template<class T_,size_t NR_,size_t NC_>
void LocalMatrix<T_,NR_,NC_>::Localize(Element*            el,
                                       const SkMatrix<T_>& a)
{
   for (size_t n=1; n<=el->getNbNodes(); n++) {
      Node *nd = (*el)(n);
      for (size_t i=1; i<=nd->getNbDOF(); i++)
         for (size_t j=1; j<=nd->getNbDOF(); j++)
            (*this)(i,j) = a(nd->getDOF(i),nd->getDOF(j));
   }
}


template<class T_,size_t NR_,size_t NC_>
void LocalMatrix<T_,NR_,NC_>::Localize(Element*             el,
                                       const SkSMatrix<T_>& a)
{
   for (size_t n=1; n<=el->getNbNodes(); n++) {
      Node *nd = (*el)(n);
      for (size_t i=1; i<=nd->getNbDOF(); i++)
         for (size_t j=1; j<=nd->getNbDOF(); j++)
            (*this)(i,j) = a(nd->getDOF(i),nd->getDOF(j));
   }
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_> & LocalMatrix<T_,NR_,NC_>::operator=(const LocalMatrix<T_,NR_,NC_>& m)
{
   _length = m._length;
   for (size_t k=0; k<_length; k++)
      _a[k] = m._a[k];
   return *this;
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_> & LocalMatrix<T_,NR_,NC_>::operator=(const T_& x)
{
   for (size_t k=0; k<_length; k++)
      _a[k] = 0;
   for (size_t k=1; k<=NR_; k++)
      (*this)(k,k) = x;
   return *this;
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_> & LocalMatrix<T_,NR_,NC_>::operator+=(const LocalMatrix<T_,NR_,NC_>& m)
{
   for (size_t k=0; k<_length; k++)
      _a[k] += m._a[k];
   return *this;
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_> & LocalMatrix<T_,NR_,NC_>::operator-=(const LocalMatrix<T_,NR_,NC_>& m)
{
   for (size_t k=0; k<_length; k++)
      _a[k] -= m._a[k];
   return *this;
}


template<class T_,size_t NR_,size_t NC_>
LocalVect<T_,NR_> LocalMatrix<T_,NR_,NC_>::operator*(LocalVect<T_,NC_>& x)
{
   LocalVect<T_,NR_> v;
   size_t k=0;
   for (size_t i=0; i<NR_; i++)
      for (size_t j=0; j<NC_; j++)
         v[i] += _a[k++]*x[j];
   return v;
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_> & LocalMatrix<T_,NR_,NC_>::operator+=(const T_& x)
{
   for (size_t k=0; k<_length; k++)
      _a[k] += x;
   return *this;
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_> & LocalMatrix<T_,NR_,NC_>::operator-=(const T_& x)
{
   for (size_t k=0; k<_length; k++)
      _a[k] -= x;
   return *this;
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_> & LocalMatrix<T_,NR_,NC_>::operator*=(const T_& x)
{
   for (size_t k=0; k<_length; k++)
      _a[k] *= x;
   return *this;
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_> & LocalMatrix<T_,NR_,NC_>::operator/=(const T_& x)
{
   for (size_t k=0; k<_length; k++)
      _a[k] /= x;
   return *this;
}


template<class T_,size_t NR_,size_t NC_>
void LocalMatrix<T_,NR_,NC_>::MultAdd(const LocalVect<T_,NC_>& x,
                                      LocalVect<T_,NR_>&       y)
{
   size_t k=0;
   for (size_t i=0; i<NR_; i++)
      for (size_t j=0; j<NC_; j++)
         y[i] += _a[k++] * x[j];
}


template<class T_,size_t NR_,size_t NC_>
void LocalMatrix<T_,NR_,NC_>::MultAddScal(const T_&                a,
                                          const LocalVect<T_,NC_>& x,
                                          LocalVect<T_,NR_>&       y)
{
   size_t k=0;
   for (size_t i=0; i<NR_; i++)
      for (size_t j=0; j<NC_; j++)
         y[i] += a * _a[k++] * x[j];
}


template<class T_,size_t NR_,size_t NC_>
void LocalMatrix<T_,NR_,NC_>::Mult(const LocalVect<T_,NC_>& x,
                                   LocalVect<T_,NR_>&       y)
{
   y = 0;
   MultAdd(x,y);
}


template<class T_,size_t NR_,size_t NC_>
void LocalMatrix<T_,NR_,NC_>::Symmetrize()
{
   for (size_t i=0; i<NR_; i++)
      for (size_t j=0; j<i; j++)
         _a[i*NC_+j] = _a[j*NC_+i];
}


template<class T_,size_t NR_,size_t NC_>
int LocalMatrix<T_,NR_,NC_>::Factor()
{
   register size_t j, k;
   if (NR_!=NC_)
      throw OFELIException("In LocalMatrix::Factor(): Can't factor a rectangle matrix.");
   for (size_t i=1; i<NR_; ++i) {
      for (j=1; j<=i; j++) {
         if (Abs(_a[NR_*(j-1)+j-1]) < OFELI_EPSMCH)
            throw OFELIException("In LocalMatrix::Factor(): The "+itos(i)+"-th pivot is too small.");
         _a[NR_*i+j-1] /= _a[NR_*(j-1)+j-1];
         for (k=0; k<j; k++)
            _a[NR_*i+j] -= _a[NR_*i+k]*_a[NR_*k+j];
      }
      for (j=i+1; j<NR_; ++j)
         for (k=0; k<i; ++k)
            _a[NR_*i+j] -= _a[NR_*i+k]*_a[NR_*k+j];
   }
   return 0;
}


template<class T_,size_t NR_,size_t NC_>
int LocalMatrix<T_,NR_,NC_>::Solve(LocalVect<T_,NR_>& b)
{
   register int i, j;
   if (NR_!=NC_)
      throw OFELIException("In LocalMatrix::Factor(): Can't solve with a rectangle matrix.");
   for (i=0; i<int(NR_); i++) {
      T_ s = 0;
      for (j=0; j<i; j++)
         s += _a[NR_*i+j] * b[j];
      b[i] -= s;
   }
   for (i=NR_-1; i>-1; i--) {
      if (Abs(_a[NR_*i+i]) < OFELI_EPSMCH)
         throw OFELIException("In LocalMatrix::Solve(b): The "+itos(i+1)+"-th diagonal entry is too small.");
      b[i] /= _a[NR_*i+i];
      for (j=0; j<i; j++)
         b[j] -= b[i] * _a[NR_*j+i];
   }
   return 0;
}


template<class T_,size_t NR_,size_t NC_>
int LocalMatrix<T_,NR_,NC_>::FactorAndSolve(LocalVect<T_,NR_>& b)
{
   Factor();
   Solve(b);
   return 0;
}


template<class T_,size_t NR_,size_t NC_>
void LocalMatrix<T_,NR_,NC_>::Invert(LocalMatrix<T_,NR_,NC_>& A)
{
   if (NR_!=NC_)
      throw OFELIException("In LocalMatrix::Invert(A): This function is valid for square matrices only.");
   LocalVect<T_,NR_> b;
   Factor();
   for (size_t i=1; i<=NR_; i++) {
      b = 0;
      b(i) = 1;
      Solve(b);
      for (size_t j=1; j<=NR_; j++)
         A(j,i) = b(j);
   }
}


template<class T_,size_t NR_,size_t NC_>
T_ LocalMatrix<T_,NR_,NC_>::getInnerProduct(const LocalVect<T_,NC_>& x,
                                            const LocalVect<T_,NR_>& y)
{
   double s = 0;
   for (size_t i=1; i<=NC_; i++)
      for (size_t j=1; i<=NR_; i++)
         s += _a[(i-1)*NC_+j-1]*x(i)*y(j);
   return s;
}


///////////////////////////////////////////////////////////////////////////////
//              A S S O C I A T E D    F U N C T I O N S                     //
///////////////////////////////////////////////////////////////////////////////

/** \fn LocalMatrix<T_,NR_,NC_> operator*(T_ a, const LocalMatrix<T_,NR_,NC_> &x)
    \brief Operator * (Multiply matrix <tt>x</tt> by scalar <tt>a</tt>)
    \ingroup VectMat
    \return <tt>a*x</tt>
 */
template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_> operator*(T_                             a,
                                  const LocalMatrix<T_,NR_,NC_>& x)
{
   LocalMatrix<T_,NR_,NC_> z;
   for (size_t i=1; i<=NR_; ++i)
      for (size_t j=1; j<=NC_; ++j)
         z(i,j) = a*x(i,j);
   return z;
}


/** \fn LocalMatrix<T_,NR_,NC_> operator*(const LocalVect<T_,NC_> &x, LocalVect<T_,NR_> &x)
    \brief Operator <tt>*</tt> (Multiply matrix <tt>A</tt> by vector <tt>x</tt>)
    \ingroup VectMat
    \return  <tt>A*x</tt>
 */
template<class T_,size_t NR_,size_t NC_>
LocalVect<T_,NR_> operator*(const LocalMatrix<T_,NR_,NC_>& A,
                            const LocalVect<T_,NC_>&       x)
{
   LocalVect<T_,NR_> b;
   for (size_t i=1; i<=NR_; ++i) {
      T_ s = T_(0);
      for (size_t j=1; j<=NC_; ++j)
        s += A(i,j)*x(j);
      b(i) = s;
   }
   return b;
}


/** \fn LocalMatrix<T_,NR_,NC_> operator/(T_ a, const LocalMatrix<T_,NR_,NC_> &x)
    \brief Operator <tt>/</tt> (Divide matrix <tt>x</tt> by scalar <tt>a</tt>)
    \ingroup VectMat
    \return <tt>x/a</tt>
 */
template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_> operator/(T_                             a,
                                  const LocalMatrix<T_,NR_,NC_>& x)
{
   LocalMatrix<T_,NR_,NC_> z;
   for (size_t i=1; i<=NR_; ++i)
      for (size_t j=1; j<=NC_; ++j)
         z(i,j) = x(i,j)/a;
   return z;
}


/** \fn LocalMatrix<T_,NR_,NC_> operator+(const LocalMatrix<T_,NR_,NC_> &x, const LocalMatrix<T_,NR_,NC_> &y)
    \brief Operator <tt>+</tt> (Add matrix x to y)
    \ingroup VectMat
    \return <tt>x+y</tt>
 */
template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_> operator+(const LocalMatrix<T_,NR_,NC_>& x,
                                  const LocalMatrix<T_,NR_,NC_>& y)
{
   LocalMatrix<T_,NR_,NC_> z;
   for (size_t i=1; i<=NR_; ++i)
      for (size_t j=1; j<=NC_; ++j)
         z(i,j) = x(i,j) + y(i,j);
   return z;
}


/** \fn LocalMatrix<T_,NR_,NC_> operator-(const LocalMatrix<T_,NR_,NC_> &x, const LocalMatrix<T_,NR_,NC_> &y)
    \brief Operator <tt>-</tt> (Subtract matrix y from x)
    \ingroup VectMat
    \return  <tt>x-y</tt>
 */
template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_> operator-(const LocalMatrix<T_,NR_,NC_>& x,
                                  const LocalMatrix<T_,NR_,NC_>& y)
{
   LocalMatrix<T_,NR_,NC_> z;
   for (size_t i=0; i<NR_; ++i)
      for (size_t j=1; j<=NC_; ++j)
         z(i,j) = x(i,j) - y(i,j);
   return z;
}


/// \fn ostream& operator<<(ostream &s, const LocalMatrix<T_,NR_,NC_> &a)
/// \brief Output vector in output stream
/// \ingroup VectMat
template<class T_, size_t NR_, size_t NC_>
ostream& operator<<(ostream&                       s,
                    const LocalMatrix<T_,NR_,NC_>& a)
{
   s.width(6);
   s.setf(ios::right|ios::scientific);
   s << endl;
   for (size_t i=1; i<=NR_; i++) {
      s << "\nRow  " << setw(6) << i << endl;
      for (size_t j=1; j<=NC_; j++)
         s << "  " << setprecision(8) << std::setfill(' ') << setw(18) << a(i,j);
      s << endl;
   }
   return s;
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
