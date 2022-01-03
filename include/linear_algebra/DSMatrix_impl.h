/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani

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


#ifndef __DSMATRIX_IMPL_H
#define __DSMATRIX_IMPL_H

#include "linear_algebra/DSMatrix.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/Matrix_impl.h"

namespace OFELI {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template<class T_>
DSMatrix<T_>::DSMatrix()
{
   _nb_rows = _length = _size = 0;
   _is_diagonal = false;
}


template<class T_>
DSMatrix<T_>::DSMatrix(size_t dim)
{
   _is_diagonal = false;
   setSize(dim);
}


template<class T_>
DSMatrix<T_>::DSMatrix(const DSMatrix<T_>& m)
{
   _size = _nb_rows = _nb_cols = m._size;
   _length = m._length;
   _is_diagonal = false;
   _a.resize(_length);
   _a = m._a;
}


template<class T_>
DSMatrix<T_>::DSMatrix(Mesh&  mesh,
                       size_t dof,
                       int    is_diagonal)
{
   _is_diagonal = is_diagonal;
   setMesh(mesh,dof);
}


template<class T_>
DSMatrix<T_>::~DSMatrix()
{ }


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
void DSMatrix<T_>::set(size_t    i,
                       size_t    j,
                       const T_& val)
{
   if (i>=j)
      _a[i*(i-1)/2+j-1] = val;
}


template<class T_>
void DSMatrix<T_>::setGraph(const Vect<RC>& I,
                            int             opt)
{ }


template<class T_>
void DSMatrix<T_>::setMesh(Mesh&  mesh,
                           size_t dof)
{
   Matrix<T_>::init_set_mesh(mesh,dof);
   _length = _size*_size;
   _diag.resize(_size);
   _a.resize(_length);
   _zero = static_cast<T_>(0);
}


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
void DSMatrix<T_>::setColumn(size_t          j,
                             const Vect<T_>& v)
{
   for (size_t i=0; i<j-1; i++)
      _a[_size*(j-1)/2+i] = v[i];
   for (size_t i=j-1; i<_size; i++)
      _a[_size*i/2+j-1] = v[i];
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
void DSMatrix<T_>::add(size_t    i,
                       size_t    j,
                       const T_& val)
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
T_& DSMatrix<T_>::operator()(size_t i,
                             size_t j)
{
   if (i>=j)
      return _a[i*(i-1)/2+j-1];
   else
      return _a[j*(j-1)/2+i-1];
}


template<class T_>
DSMatrix<T_>& DSMatrix<T_>::operator=(const DSMatrix<T_>& m)
{
   _a = m._a;
   return *this;
}


template<class T_>
DSMatrix<T_>& DSMatrix<T_>::operator=(const T_& x)
{
   clear(_a);
   for (size_t i=1; i<=_size; i++) {
      _diag[i-1] = x;
      set(i,i,x);
   }
   return *this;
}


template<class T_>
DSMatrix<T_>& DSMatrix<T_>::operator+=(const T_& x)
{
   for (size_t i=1; i<=_length; ++i)
      _a.add(i,x);
   return *this;
}


template<class T_>
DSMatrix<T_>& DSMatrix<T_>::operator-=(const T_& x)
{
   for (size_t i=1; i<=_length; ++i)
      _a.add(i,-x);
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
      if (Abs(pivot) < OFELI_EPSMCH)
         throw OFELIException("In DSMatrix::setLDLt(): " + to_string(i+1) + "-th pivot is null.");
      _a[(i+1)*i/2+i] = 1./pivot;
   }
   return err;
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
   y = _zero;
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


template<class T_>
int DSMatrix<T_>::solve(Vect<T_>& b,
                        bool      fact)
{
   int ret = 0;
   if (fact)
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
int DSMatrix<T_>::Factor()
{
   return setLDLt();
}


template<class T_>
int DSMatrix<T_>::solve(const Vect<T_>& b,
                        Vect<T_>&       x,
                        bool            fact)
{
   x = b;
   return solve(x,fact);
}


template<class T_>
T_* DSMatrix<T_>::getArray() const
{
   return _a;
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


///////////////////////////////////////////////////////////////////////////////
//                 A S S O C I A T E D   F U N C T I O N S                   //
///////////////////////////////////////////////////////////////////////////////

template<class T_>
Vect<T_> operator*(const DSMatrix<T_>& A,
                   const Vect<T_>&     b)
{
   Vect<T_> v(A.getNbRows());
   A.Mult(b,v);
   return v;
}


template<class T_>
ostream& operator<<(ostream&            s,
                    const DSMatrix<T_>& A)
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
