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

                      Implementation of class TrMatrix

  ==============================================================================*/


#ifndef _TR_MATRIX_IMPL_H
#define _TR_MATRIX_IMPL_H

#include "linear_algebra/TrMatrix.h"
#include "mesh/Mesh.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template<class T_>
TrMatrix<T_>::TrMatrix()
{
   _msize.mt = TRIDIAGONAL;
}


template<class T_>
TrMatrix<T_>::TrMatrix(size_t size)
{
   _msize.mt = TRIDIAGONAL;
   _msize.size = _msize.nb_rows = _msize.nb_cols = size;
   _msize.length = 3*_msize.size;
   _a.resize(_msize.length,static_cast<T_>(0));
}


template<class T_>
TrMatrix<T_>::TrMatrix(const TrMatrix& m)
{
   _msize = m._msize;
   _a.resize(_msize.length);
   _a = m._a;
}


template<class T_>
TrMatrix<T_>::~TrMatrix()
{ }


template<class T_>
void TrMatrix<T_>::Identity()
{
   Diagonal(1);
}


template<class T_>
void TrMatrix<T_>::Diagonal()
{
   for (size_t i=1; i<=_msize.size; i++)
      _a[3*i-2] = static_cast<T_>(0);
}


template<class T_>
void TrMatrix<T_>::Diagonal(const T_& a)
{
   _msize.length = _msize.nb_rows;
   for (size_t i=1; i<=_msize.nb_rows; ++i) {
      _a[3*i-3] = static_cast<T_>(0);
      _a[3*i-2] = a;
      _a[3*i-1] = static_cast<T_>(0);
   }
}


template<class T_>
void TrMatrix<T_>::setGraph(const vector<RC>& I,
                            int               opt)
{ }


template<class T_>
void TrMatrix<T_>::setMesh(Mesh&  mesh,
                           size_t dof)
{
   throw OFELIException("In TrMatrix::setMesh(Mesh,size_t): This member "
                        "function is not valid for class TrMatrix");
}


template<class T_>
void TrMatrix<T_>::setMesh(Mesh&  mesh,
                           size_t dof,
                           size_t type)
{
   throw OFELIException("In TrMatrix::setMesh(Mesh,size_t,size_t): "
                        "This member function is not valid for class TrMatrix");
}

 
template<class T_>
void TrMatrix<T_>::setMesh(size_t dof,
                           Mesh&  mesh,
                           int    code)
{
   throw OFELIException("In TrMatrix::setMesh(Mesh,Mesh,int): "
                        "This member function is not valid for class TrMatrix");
}


template<class T_>
void TrMatrix<T_>::setMesh(size_t dof,
                           size_t nb_eq,
                           Mesh&  mesh)
{
   dof = 0;
   nb_eq = 0;
   if (mesh.getDim()==0) { }
}


template<class T_>
void TrMatrix<T_>::setSize(size_t size)
{
   _msize.nb_rows = _msize.nb_cols = _msize.size = size;
   _msize.length = 3*_msize.nb_rows;
   _a.resize(_msize.length);
}


template<class T_>
void TrMatrix<T_>::MultAdd(const Vect<T_>& x, Vect<T_>& y) const
{
   y(1) += _a[1]*x(1) + _a[2]*x(2);
   for (size_t i=2; i<_msize.size; ++i)
      y(i) += _a[3*i-3]*x(i-1) + _a[3*i-2]*x(i) + _a[3*i-1]*x(i+1);
   y(_msize.size) += _a[3*_msize.size-3]*x(_msize.size-1) + _a[3*_msize.size-2]*x(_msize.size);
}


template<class T_>
void TrMatrix<T_>::MultAdd(T_              a,
                           const Vect<T_>& x,
                           Vect<T_>&       y) const
{
   y(1) += a*(_a[1]*x(1) + _a[2]*x(2));
   for (size_t i=2; i<_msize.size; ++i)
      y(i) += a * (_a[3*i-3]*x(i-1) + _a[3*i-2]*x(i) + _a[3*i-1]*x(i+1));
   y(_msize.size) += a * (_a[3*_msize.size-3]*x(_msize.size-1) + _a[3*_msize.size-2]*x(_msize.size));
}


template<class T_>
void TrMatrix<T_>::Mult(const Vect<T_>& x,
                        Vect<T_>&       y) const
{
   y = static_cast<T_>(0);
   MultAdd(x,y);
}


template<class T_>
void TrMatrix<T_>::TMult(const Vect<T_>& x,
                         Vect<T_>&       y) const
{
   y(1) = _a[1]*x(1) + _a[3]*x(2);
   for (size_t i=2; i<_msize.size; ++i)
      y(i) = _a[3*i-4]*x(i-1) + _a[3*i-2]*x(i) + _a[3*i]*x(i+1);
   y(_msize.size) = _a[3*_msize.size-4]*x(_msize.size-1) + _a[3*_msize.size-2]*x(_msize.size);
}


template<class T_>
void TrMatrix<T_>::Axpy(T_                  a,
                        const TrMatrix<T_>& m)
{
   _a += a * m._a;
}
   

template<class T_>
void TrMatrix<T_>::Axpy(T_                a,
                        const Matrix<T_>* m)
{
   for (size_t i=0; i<_msize.length; i++)
      _a[i] += a * m->_a[i];
}


template<class T_>
void TrMatrix<T_>::set(size_t    i,
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


template<class T_>
void TrMatrix<T_>::add(size_t    i,
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


template<class T_>
void TrMatrix<T_>::add(size_t    i,
                       const T_& val)
{
   _a[i-1] += val;
}


template<class T_>
T_ TrMatrix<T_>::at(size_t i,
                    size_t j)
{
   if (i==j)
      return _a[3*i-2];
   else if (i==j+1)
      return _a[3*i-3];
   else if (i==j-1)
      return _a[3*i-1];
   else
      return static_cast<T_>(0);
}


template<class T_>
T_ TrMatrix<T_>::operator()(size_t i,
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


template<class T_>
T_& TrMatrix<T_>::operator()(size_t i,
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
                            to_string(i) + "," + to_string(j) +
                           ") not compatible with tridiagonal storage.");
   }
   return _zero;
}

 
template<class T_>
TrMatrix<T_>& TrMatrix<T_>::operator=(const TrMatrix<T_>& m)
{
   _a = m._a;
   return *this;
}


template<class T_>
TrMatrix<T_>& TrMatrix<T_>::operator=(const T_& x)
{
   clear(_a);
   for (size_t i=1; i<=_msize.size; ++i)
      _a[3*i-2] = x;
   return *this;
}


template<class T_>
TrMatrix<T_>& TrMatrix<T_>::operator*=(const T_& x)
{
   for (size_t i=1; i<=_msize.length; i++)
      _a[i-1] *= x;
   return *this;
}


template<class T_>
int TrMatrix<T_>::Factor()
{
   throw OFELIException("In TrMatrix::Factor(): Matrix cannot be fatorized separately "
                        "Call Solve(b,fact) directly.");
   return 1;
}


template<class T_>
int TrMatrix<T_>::solve(Vect<T_>& b,
                        bool      fact)
{
   int ret = 0;
   T_ p = 0;
   for (size_t i=1; i<_msize.size; ++i) {
      ret = int(i);
      if (Abs(_a[3*i-2]) < OFELI_EPSMCH)
         throw OFELIException("In TrMatrix::solve(Vect<T_>): " + to_string(int(i))
                              + "-th pivot is null."); 
      else
         p = _a[3*i]/_a [3*i-2];
      _a[3*i+1] -= p*_a[3*i-1];
      b.add(i+1,-p*b(i));
   }
   if (Abs(_a[3*_msize.size-2]) < OFELI_EPSMCH)
      throw OFELIException("In TrMatrix::solve(Vect<T_>): " + to_string(_msize.size)
                           + "-th pivot is null.");
   else
      b(_msize.size) /= _a[3*_msize.size-2];
   for (int j=int(_msize.size)-1; j>0; j--)
      b(j) = (b(j)-_a[3*j-1]*b(j+1))/_a[3*j-2];
   return ret;
}

template<class T_>
int TrMatrix<T_>::solve(const Vect<T_>& b,
                        Vect<T_>&       x,
                        bool            fact)
{
   x = b;
   return solve(x,fact);
}


template<class T_>
T_* TrMatrix<T_>::get() const
{
   return &_a[0];
}


template<class T_>
T_ TrMatrix<T_>::get(size_t i, size_t j) const
{
   if (i==j)
      return _a[3*i-2];
   else if (i==j+1)
      return _a[3*i-3];
   else if (i==j-1)
      return _a[3*i-1];
   else
      return static_cast<T_>(0);
}


template<>
inline void TrMatrix<real_t>::Laplace1D(real_t h)
{
   for (size_t i=1; i<=_msize.size; ++i) {
      _a[3*i-3] = -1.0/h;
      _a[3*i-2] =  2.0/h;
      _a[3*i-1] = -1.0/h;
   }
}

///////////////////////////////////////////////////////////////////////////////
//                           ASSOCIATED  FUNCTIONS                           //
///////////////////////////////////////////////////////////////////////////////

template<class T_>
Vect<T_> operator*(const TrMatrix<T_>& A,
                   const Vect<T_>&     b)
{
   size_t n = A.getSize();
   Vect<T_> v(n);
   v(1) = A(1,1)*b(1) + A(1,2)*b(2);
   for (size_t i=2; i<n; ++i)
      v(i) = A(i,i-1)*b(i-1) + A(i,i)*b(i) + A(i,i+1)*b(i+1);
   v(n) = A(n,n-1)*b(n-1) + A(n,n)*b(n);
   return v;
}


template<class T_>
TrMatrix<T_> operator*(T_                  a,
                       const TrMatrix<T_>& A)
{
   TrMatrix<T_> v(A);
   for (size_t i=0; i<A.getLength(); ++i)
      v[i] *= a;
   return v;
}


template<class T_>
ostream& operator<<(ostream&            s,
                    const TrMatrix<T_>& A)
{
   s.setf(ios::right|ios::scientific);
   s << endl;
   for (size_t i=1; i<=A.getNbRows(); ++i) {
      s << "\nRow  " << setw(6) << i << endl;
      for (size_t j=1; j<=A.getNbRows(); ++j)
         s << "  " << setprecision(8) << std::setfill(' ') << setw(18) << A(i,j);
      s << endl;
   }
   return s;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */

#endif