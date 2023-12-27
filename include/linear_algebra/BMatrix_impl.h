/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2024 Rachid Touzani

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

                     Implementation of class BMatrix

  ==============================================================================*/


#ifndef _B_MATRIX_IMPL_H
#define _B_MATRIX_IMPL_H

#include "linear_algebra/BMatrix.h"
#include "mesh/Mesh.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template<class T_>
BMatrix<T_>::BMatrix()
{
   _size = _ld = _ud = 0;
}


template<class T_>
BMatrix<T_>::BMatrix(size_t size,
                     int    ld,
                     int    ud)
{
   if (size<3 || ld<0 || ud<0 || ud+ld+1>int(size))
      throw OFELIException("In BMatrix::BMatrix(size_t,int,int): Illegal arguments.");
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
}


template<class T_>
BMatrix<T_>::~BMatrix()
{ }


template<class T_>
void BMatrix<T_>::setGraph(const Vect<RC>& I,
                           int             opt)
{ }


template<class T_>
void BMatrix<T_>::setMesh(Mesh&  mesh,
                          size_t dof)
{
   throw OFELIException("In BMatrix::setMesh(Mesh,size_t): "
                        "This member function is not valid for class BMatrix");
}


template<class T_>
void BMatrix<T_>::setMesh(Mesh&  mesh,
                          size_t dof,
                          size_t type)
{
   throw OFELIException("In SkSMatrix::setMesh(Mesh,size_t,size_t): "
                        "This member function is not valid for class BMatrix");
}


template<class T_>
void BMatrix<T_>::setMesh(size_t dof, 
                          Mesh&  mesh,
                          int    code)
{
   throw OFELIException("In BMatrix::setMesh(size_t,Mesh,int): "
                        "This member function is not valid for class BMatrix");
}


template<class T_>
void BMatrix<T_>::setMesh(size_t dof,
                          size_t nb_eq,
                          Mesh&  mesh)
{
   throw OFELIException("In SkSMatrix::setMesh(size_t,size_t,Mesh): "
                        "This member function is not valid for class BMatrix");
}


template<class T_>
void BMatrix<T_>::setDiag()
{ }


template<class T_>
int BMatrix<T_>::Factor()
{
   return setLU();
}

 
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
void BMatrix<T_>:: MultAdd(T_              a,
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
{ }
   

template<class T_>
void BMatrix<T_>::Axpy(T_                 a,
                       const BMatrix<T_>& x)
{ }


template<class T_>
void BMatrix<T_>::Axpy(T_                a,
                       const Matrix<T_>* x)
{ }


template<class T_>
void BMatrix<T_>:: set(size_t    i,
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
T_& BMatrix<T_>::operator()(size_t i,
                                 size_t j)
{
   if (_ld+int(j)-int(i)>=0 && _ld+int(j)-int(i)<=_ud+_ld)
      return _a[i-1][_ld+j-i];
   else
      throw OFELIException("In BMatrix::operator(): Illegal pair of indices " +
                           to_string(i)+","+to_string(j)+").");
   return _zero;
}


template<class T_>
BMatrix<T_>& BMatrix<T_>::operator=(const BMatrix<T_>& m)
{
   _a = m._a;
   return *this;
}


template<class T_>
BMatrix<T_>& BMatrix<T_>::operator=(const T_& x)
{
   for (size_t i=0; i<_size; ++i)
      for (int j=0; j<_ld+_ud+1; ++j)
         _a[i][j] = 0;
   for (size_t i=1; i<=_size; ++i)
      (*this)(i,i) = x;
   return *this;
}


template<class T_>
BMatrix<T_>& BMatrix<T_>::operator*=(const T_& x)
{
   for (size_t i=0; i<_size; ++i)
      for (int j=0; j<_ld+_ud+1; ++j)
         _a[i][j] *= x;
   return *this;
}


template<class T_>
BMatrix<T_>& BMatrix<T_>::operator+=(const T_& x)
{
   for (size_t i=0; i<_size; ++i)
      for (int j=0; j<_ld+_ud+1; ++j)
         _a[i][j] += x;
   return *this;
}


template<class T_>
void BMatrix<T_>::add(size_t    i,
                      const T_& val)
{
   _a[i-1][0] += val;
}


template<class T_>
int BMatrix<T_>::setLU()
{
   for (int i=0; i<int(_size)-1; i++) {
      int kend=std::min(_ld+1,int(_size)-i), kjend=std::min(_ud+1,int(_size)-i);
      if (Abs(_a[i][_ld]) < OFELI_EPSMCH)
         throw OFELIException("In BMatrix::Factor(): " + to_string(i+1) + "-th pivot is null.");
      for (int k=1; k!=kend; k++) {
         _a[k+i][_ld-k] /= _a[i][_ld];
         for (int j=1; j!=kjend; j++)
            _a[k+i][j+_ld-k] -= _a[k+i][_ld-k] * _a[i][j+_ld];
      }
   }
   return 0;
}


template<class T_>
int BMatrix<T_>::solve(Vect<T_>& b,
                       bool      fact)
{
   if (_size<3 || _ld<0 || _ud<0 || _ud+_ld+1>int(_size))
      throw OFELIException("In BMatrix::solve(Vect<double>): Illegal arguments.");
   int ret = setLU();
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
                       Vect<T_>&       x,
                       bool            fact)
{
   x = b;
   return solve(x,fact);
}


template<class T_>
T_* BMatrix<T_>::get() const
{
   return _a;
}


template<class T_>
T_ BMatrix<T_>::get(size_t i, size_t j) const
{
   if (_ld+int(j)-int(i)>=0 && _ld+int(j)-int(i)<=_ld+_ud)
      return _a[i-1][_ld+j-i];
   else
      return _zero;
}


///////////////////////////////////////////////////////////////////////////////
//                           ASSOCIATED  FUNCTIONS                           //
///////////////////////////////////////////////////////////////////////////////
 
template<class T_>
Vect<T_> operator*(const BMatrix<T_>& A,
                   const Vect<T_>&    b)
{
   Vect<T_> v(b.size());
   A.Mult(b,v);
}


template<class T_>
BMatrix<T_> operator*(T_                 a,
                      const BMatrix<T_>& A)
{
   BMatrix<T_> v(A);
   for (size_t i=0; i<A.getLength(); ++i)
      v[i] *= a;
   return v;
}


template<class T_>
ostream& operator<<(ostream&           s,
                    const BMatrix<T_>& A)
{
   s.setf(ios::right|ios::scientific);
   s << endl;
   for (size_t i=1; i<=A.getNbRows(); i++) {
      s << "\nRow  " << setw(6) << i << endl;
      for (size_t j=1; j<=A.getNbColumns(); j++)
         s << "  " << setprecision(8) << std::setfill(' ') << setw(18) << A(i,j);
      s << endl;
   }
   return s;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */

#endif
