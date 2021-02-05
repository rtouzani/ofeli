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

                     Implementation of class LocalMatrix

  ==============================================================================*/


#ifndef __LOCAL_MATRIX_IMPL_H
#define __LOCAL_MATRIX_IMPL_H

#include "linear_algebra/LocalMatrix.h"
#include "linear_algebra/LocalVect.h"
#include "linear_algebra/SkMatrix.h"
#include "linear_algebra/SkSMatrix.h"
#include "linear_algebra/SpMatrix.h"
#include "mesh/Element.h"
#include "util/util.h"
#include "OFELIException.h"

namespace OFELI {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

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
{ Localize(el,a); }


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_>::LocalMatrix(Element*            el,
                                     const SkMatrix<T_>& a)
{ Localize(el,a); }


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_>::LocalMatrix(Element*             el,
                                     const SkSMatrix<T_>& a)
{ Localize(el,a); }


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_>::~LocalMatrix()
{ }


template<class T_,size_t NR_,size_t NC_>
T_& LocalMatrix<T_,NR_,NC_>::operator()(size_t i,
                                        size_t j)
{
   return _a[(i-1)*NC_+j-1];
}


template<class T_,size_t NR_,size_t NC_>
T_ LocalMatrix<T_,NR_,NC_>::operator()(size_t i,
                                       size_t j) const
{
   return _a[(i-1)*NC_+j-1];
}


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


template<class T_,size_t NR_,size_t NC_>
void LocalMatrix<T_,NR_,NC_>::Localize(Element*            el,
                                       const SpMatrix<T_>& a)
{
   for (size_t n=1; n<=el->getNbNodes(); n++) {
      Node *nd = (*el)(n);
      for (size_t i=1; i<=nd->getNbDOF(); i++)
         for (size_t j=1; j<=nd->getNbDOF(); j++)
            _a[(i-1)*NC_+j-1] = a(nd->getDOF(i),nd->getDOF(j));
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
            _a[(i-1)*NC_+j-1] = a(nd->getDOF(i),nd->getDOF(j));
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
            _a[(i-1)*NC_+j-1] = a(nd->getDOF(i),nd->getDOF(j));
   }
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_>& LocalMatrix<T_,NR_,NC_>::operator=(const LocalMatrix<T_,NR_,NC_>& m)
{
   _length = m._length;
   for (size_t k=0; k<_length; k++)
      _a[k] = m._a[k];
   return *this;
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_>& LocalMatrix<T_,NR_,NC_>::operator=(const T_& x)
{
   for (size_t k=0; k<_length; k++)
      _a[k] = 0;
   for (size_t k=1; k<=NR_; k++)
      _a[(k-1)*NC_+k-1] = x;
   return *this;
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_>& LocalMatrix<T_,NR_,NC_>::operator+=(const LocalMatrix<T_,NR_,NC_>& m)
{
   for (size_t k=0; k<_length; k++)
      _a[k] += m._a[k];
   return *this;
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_>& LocalMatrix<T_,NR_,NC_>::operator-=(const LocalMatrix<T_,NR_,NC_>& m)
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
LocalMatrix<T_,NR_,NC_>& LocalMatrix<T_,NR_,NC_>::operator+=(const T_& x)
{
   for (size_t k=0; k<_length; k++)
      _a[k] += x;
   return *this;
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_>& LocalMatrix<T_,NR_,NC_>::operator-=(const T_& x)
{
   for (size_t k=0; k<_length; k++)
      _a[k] -= x;
   return *this;
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_>& LocalMatrix<T_,NR_,NC_>::operator*=(const T_& x)
{
   for (size_t k=0; k<_length; k++)
      _a[k] *= x;
   return *this;
}


template<class T_,size_t NR_,size_t NC_>
LocalMatrix<T_,NR_,NC_>& LocalMatrix<T_,NR_,NC_>::operator/=(const T_ &x)
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
   size_t j, k;
   if (NR_!=NC_)
      throw OFELIException("In LocalMatrix::Factor(): Can't factor a rectangle matrix.");
   for (size_t i=1; i<NR_; ++i) {
      for (j=1; j<=i; j++) {
         if (Abs(_a[NR_*(j-1)+j-1]) < OFELI_EPSMCH)
            throw OFELIException("In LocalMatrix::Factor(): The "+to_string(i)+"-th pivot is too small.");
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
int LocalMatrix<T_,NR_,NC_>::solve(LocalVect<T_,NR_>& b)
{
   int i, j;
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
         throw OFELIException("In LocalMatrix::Solve(b): The "+to_string(i+1)+"-th diagonal entry is too small.");
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
   solve(b);
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
      solve(b);
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


template<class T_,size_t NR_,size_t NC_>
T_* LocalMatrix<T_,NR_,NC_>::get()
{ return _a; }


///////////////////////////////////////////////////////////////////////////////
//              A S S O C I A T E D    F U N C T I O N S                     //
///////////////////////////////////////////////////////////////////////////////


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


template<class T_, size_t NR_, size_t NC_>
ostream& operator<<(ostream&                       s,
                    const LocalMatrix<T_,NR_,NC_>& A)
{
   s.width(6);
   s.setf(ios::right|ios::scientific);
   s << endl;
   for (size_t i=1; i<=NR_; i++) {
      s << "\nRow  " << setw(6) << i << endl;
      for (size_t j=1; j<=NC_; j++)
         s << "  " << setprecision(8) << std::setfill(' ') << setw(18) << A(i,j);
      s << endl;
   }
   return s;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */

#endif
