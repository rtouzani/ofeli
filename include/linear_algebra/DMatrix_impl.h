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

                     Implementation of class DMatrix

  ==============================================================================*/


#ifndef __DMATRIX_IMPL_H
#define __DMATRIX_IMPL_H

#include "linear_algebra/DMatrix.h"
#include "util/util.h"

namespace OFELI {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template<class T_>
DMatrix<T_>::DMatrix() : _qr_set(0)
{
   _length = _size = _nb_rows = _nb_cols = 0;
   _is_diagonal = false;
}


template<class T_>
DMatrix<T_>::DMatrix(size_t nr) : _qr_set(0)
{
   setSize(nr);
   _is_diagonal = false;
}


template<class T_>
DMatrix<T_>::DMatrix(size_t nr,
                     size_t nc) : _qr_set(0)
{
   setSize(nr,nc);
   _diag.resize(nr);
   _is_diagonal = false;
}


template<class T_>
DMatrix<T_>::DMatrix(Vect<T_>& v) : _qr_set(0)
{
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
   _is_diagonal = m._is_diagonal;
   _qr_set = m._qr_set;
}


template<class T_>
DMatrix<T_>::DMatrix(Mesh&  mesh,
                     size_t dof,
                     int    is_diagonal)
{
   _is_diagonal = is_diagonal;
   setMesh(mesh,dof);
}


template<class T_>
DMatrix<T_>::~DMatrix()
{ }


template<class T_>
void DMatrix<T_>::setMesh(Mesh&  mesh,
                          size_t dof)
{
   Matrix<T_>::init_set_mesh(mesh,dof);
   _length = _size*_size;
   _diag.resize(_size);
   _a.resize(_length);
   _zero = static_cast<T_>(0);
}


template<class T_>
void DMatrix<T_>::setMesh(size_t dof,
                          Mesh&  mesh,
                          int    code)
{
// This is just to avoid warning on unused variable
   dof = 0;
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
void DMatrix<T_>::set(size_t    i,
                      size_t    j,
                      const T_& val)
{
   _a[_nb_cols*(i-1)+j-1] = val;
}


template<class T_>
void DMatrix<T_>::reset()
{
   _a.clear();
}


template<class T_>
void DMatrix<T_>::setRow(size_t          i,
                         const Vect<T_>& v)
{
   for (size_t j=0; j<_nb_cols; j++)
      _a[_nb_cols*(i-1)+j] = v[j];
}


template<class T_>
void DMatrix<T_>::setColumn(size_t          j,
                            const Vect<T_>& v)
{
   for (size_t i=0; i<_nb_rows; i++)
      _a[_nb_cols*i+j-1] = v[i];
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
                        Vect<T_>&       y) const
{
   for (size_t i=0; i<_nb_rows; i++)
      for (size_t j=0; j<_nb_cols; j++)
         y.add(i+1,_a[_nb_cols*(j-1)+i-1]*x[j]);
}


template<class T_>
void DMatrix<T_>::add(size_t    i,
                      size_t    j,
                      const T_& val)
{
   _a[_nb_cols*(i-1)+j-1] += val;
}


template<class T_>
void DMatrix<T_>::Axpy(T_                 a,
                       const DMatrix<T_>& m)
{
   Axpy(a,m._a,_a);
}


template<class T_>
void DMatrix<T_>::Axpy(T_                a,
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
int DMatrix<T_>::solveQR(const Vect<T_>& b,
                         Vect<T_>&       x)
{
   return -1;
}


template<class T_>
int DMatrix<T_>::solveTransQR(const Vect<T_>& b,
                              Vect<T_>&       x)
{
   return -1;
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
T_& DMatrix<T_>::operator()(size_t i,
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
void DMatrix<T_>::setGraph(const Vect<RC>& I,
                           int             opt)
{ }


template<class T_>
int DMatrix<T_>::Factor()
{
   return setLU();
}


template<class T_>
int DMatrix<T_>::setLU()
{
   for (size_t i=1; i<_size; i++) {
      for (size_t j=1; j<=i; j++) {
         if (Abs(_a[_nb_rows*(j-1)+j-1]) < OFELI_EPSMCH)
            throw OFELIException("In DMatrix::setLU(): " + to_string(int(i)) + "-th pivot is null.");
         _a[_nb_rows*i+j-1] /= _a[_nb_rows*(j-1)+j-1];
         for (size_t k=0; k<j; k++)
            _a[_nb_rows*i+j] -= _a[_nb_rows*i+k]*_a[_nb_rows*k+j];
      }
      for (size_t j=i+1; j<_size; j++)
         for (size_t k=0; k<i; k++)
            _a[_nb_rows*i+j] -= _a[_nb_rows*i+k]*_a[_nb_rows*k+j];
   }
   return 0;
}


template<class T_>
int DMatrix<T_>::setTransLU()
{
   for (size_t i=1; i<_size; i++) {
      for (size_t j=1; j<=i; j++) {
         if (Abs(_a[_nb_rows*(i-1)+i-1]) < OFELI_EPSMCH)
            throw OFELIException("In DMatrix::setTLU(): " + to_string(i) + "-th pivot is null.");
         _a[_nb_rows*j+i-1] /= _a[_nb_rows*(i-1)+i-1];
         for (size_t k=0; k<j; k++)
            _a[_nb_rows*j+i] -= _a[_nb_rows*k+i]*_a[_nb_rows*j+k];
      }
      for (size_t j=i+1; j<_size; j++)
         for (size_t k=0; k<i; k++)
            _a[_nb_rows*j+i] -= _a[_nb_rows*k+i]*_a[_nb_rows*j+k];
   }
   return 0;
}


template<class T_>
int DMatrix<T_>::solve(Vect<T_>& b,
                       bool      fact)
{
   int ret = 0;
   Vect<T_> x(b.size());
   ret = solve(b,x);
   b = x;
   return ret;
}


template<class T_>
int DMatrix<T_>::solveTrans(Vect<T_>& b,
                            bool      fact)
{
   int ret = 0;
   Vect<T_> x(b.size());
   ret = solveTrans(b,x);
   b = x;
   return ret;
}


template<class T_>
int DMatrix<T_>::solve(const Vect<T_>& b,
                       Vect<T_>&       x,
                       bool            fact)
{
   int ret = 0;
   if (_nb_rows != _nb_cols) {
      setQR();
      ret = solveQR(b,x);
   }
   else
      ret = solveLU(b,x,fact);
   return ret;
}


template<class T_>
int DMatrix<T_>::solveTrans(const Vect<T_>& b,
                            Vect<T_>&       x,
                            bool            fact)
{
   int ret = 0;
   if (_nb_rows != _nb_cols)
      ret = solveTransQR(b,x);
   else
     ret = solveTransLU(b,x,fact);
   return ret;
}


template<class T_>
int DMatrix<T_>::solveLU(const Vect<T_>& b,
                         Vect<T_>&       x,
                         bool            fact)
{
   x = b;
   return solveLU(x,fact);
}


template<class T_>
int DMatrix<T_>::solveLU(Vect<T_>& b,
                         bool      fact)
{
   int ret = 0;
   if (fact)
      ret = setLU();
   for (size_t i=0; i<_size; i++) {
      T_ s = 0;
      for (size_t j=0; j<i; j++)
         s += _a[_nb_cols*i+j] * b[j];
      b[i] -= s;
   }
   for (int ii=int(_size)-1; ii>-1; ii--) {
      if (Abs(_a[_nb_cols*ii+ii])<OFELI_EPSMCH)
         throw OFELIException("In DMatrix::solveLU(Vect<T_>): " + to_string(ii+1) + "-th pivot is null.");
      b[ii] /= _a[_nb_cols*ii+ii];
      for (size_t j=0; j<size_t(ii); j++)
         b[j] -= b[ii] * _a[_nb_cols*j+ii];
   }
   return ret;
}


template<class T_>
int DMatrix<T_>::solveTransLU(const Vect<T_>& b,
                              Vect<T_>&       x,
                              bool            fact)
{
   x = b;
   return solveTrans(x,fact);
}


template<class T_>
int DMatrix<T_>::solveTransLU(Vect<T_>& b,
                              bool      fact)
{
   int ret = 0;
   if (fact)
      ret = setTransLU();
   for (size_t i=0; i<_size; i++) {
      T_ s=0;
      for (size_t j=0; j<i; j++)
         s += _a[_nb_cols*j+i] * b[j];
      b[i] -= s;
   }
   for (int ii=int(_size)-1; ii>-1; ii--) {
      if (Abs(_a[_nb_cols*ii+ii])<OFELI_EPSMCH)
         throw OFELIException("In DMatrix::solveTransLU(Vect<real_t>): " + to_string(ii+1)
                              + "-th pivot is null.");
      b[ii] /= _a[_nb_cols*ii+ii];
      for (size_t j=0; j<size_t(ii); j++)
         b[j] -= b[ii] * _a[_nb_cols*ii+j];
   }
   return ret;
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
   if (_nb_rows!=_nb_cols)
      throw OFELIException("In DMatrix::operator=(T_): Operator is valid "
                           "for square matrices only.");
   clear(_a);
   for (size_t i=1; i<=_size; i++) {
      set(i,i,x);
      _diag[i-1] = x;
   }
   return *this;
}


template<class T_>
DMatrix<T_>& DMatrix<T_>::operator*=(const T_& x)
{
   for (size_t i=0; i<_length; i++)
      _a[i] *= x;
   return *this;
}


template<class T_>
DMatrix<T_>& DMatrix<T_>::operator+=(const T_& x)
{
   for (size_t i=0; i<_length; i++)
      _a[i] += x;
   return *this;
}


template<class T_>
DMatrix<T_>& DMatrix<T_>::operator-=(const T_& x)
{
   for (size_t i=0; i<_length; i++)
      _a[i] -= x;
   return *this;
}


template<class T_>
T_* DMatrix<T_>::getArray() const
{
   return _a;
}


template<class T_>
T_ DMatrix<T_>::get(size_t i,
                    size_t j) const
{
   return _a[_nb_cols*(i-1)+j-1];
}


template<>
inline int DMatrix<real_t>::setQR()
{
   size_t n = std::min(_nb_rows,_nb_cols);
   _qr_c.resize(n); OFELI::clear(_qr_c);
   _qr_d.resize(n); OFELI::clear(_qr_d);
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
   _qr_c.resize(n); OFELI::clear(_qr_c);
   _qr_d.resize(n); OFELI::clear(_qr_d);
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


///////////////////////////////////////////////////////////////////////////////
//                           ASSOCIATED  FUNCTIONS                           //
///////////////////////////////////////////////////////////////////////////////


template<class T_>
Vect<T_> operator*(const DMatrix<T_>& A,
                   const Vect<T_>&    b)
{
   Vect<T_> v(A.getNbRows());
   A.Mult(b,v);
   return v;
}


template<class T_>
ostream& operator<<(ostream&           s,
                    const DMatrix<T_>& A)
{
   for (size_t i=1; i<=A.getNbRows(); i++) {
      s << "\nRow  " << setw(6) << i << endl;
      for (size_t j=1; j<=A.getNbColumns(); j++)
         s << "  " << setprecision(8) << std::setfill(' ')
           << setw(18) << A(i,j);
      s << endl;
   }
   return s;
}


#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */

#endif
