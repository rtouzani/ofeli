/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

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

                         Implementation of class SkSMatrix

  ==============================================================================*/


#ifndef __SKSMATRIX_IMPL_H
#define __SKSMATRIX_IMPL_H

#include "linear_algebra/SkSMatrix.h"
#include "linear_algebra/GraphOfMatrix.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/Matrix_impl.h"

namespace OFELI {

#ifndef DOXYGEN_SHOULD_SKIP_THIS


template<class T_>
SkSMatrix<T_>::SkSMatrix() : _dof(0)
{
   _is_diagonal = false;
   _dof_type = NODE_DOF;
}


template<class T_>
SkSMatrix<T_>::SkSMatrix(size_t size,
                         int    is_diagonal) : _dof(0)
{
   _zero = 0;
   _is_diagonal = is_diagonal;
   _dof_type = NODE_DOF;
   _nb_rows = _nb_cols = _size = size;
   _ch.resize(size);
   _diag.resize(_size);
   _ch[0] = 0;
   if (_is_diagonal) {
      for (size_t i=1; i<_size; i++)
         _ch[i] = _ch[i-1] + 1;
   }
   else {
      for (size_t i=1; i<_size; i++)
         _ch[i] = _ch[i-1] + i + 1;
   }
   _length = _ch[_size-1] + 1;
   _a.resize(_length);
}


template<class T_>
SkSMatrix<T_>::SkSMatrix(Mesh&  mesh,
                         size_t dof,
                         int    is_diagonal)
{
   _is_diagonal = is_diagonal;
   setMesh(mesh,dof);
}


template<class T_>
SkSMatrix<T_>::SkSMatrix(const Vect<size_t>& ColHt) : _dof(0)
{
   _is_diagonal = false;
   _dof_type = NODE_DOF;
   _size = ColHt.size();
   _zero = 0;
   _ch.resize(_size);
   _ch[0] = 0;
   for (size_t i=1; i<_size; ++i)
      _ch[i] = _ch[i-1] + ColHt[i];
   _length = _ch[_size-1] + 1;
   _a.resize(_length);
   _diag.resize(_size);
}


template<class T_>
SkSMatrix<T_>::SkSMatrix(const Vect<size_t>& I,
                         const Vect<size_t>& J,
                         int                 opt) : _dof(0)
{
   _is_diagonal = false;
   _size = 0;
   _dof_type = NODE_DOF;
   size_t n = I.size();
   std::vector<RC> pp(n);
   for (size_t i=0; i<n; i++) {
      _size = std::max(_size,I[i]);
      pp[i] = RC(I[i]-1,J[i]-1);
   }
   _ch.resize(_size);
   _nb_rows = _nb_cols = _size;
   if (opt==0) {
      sort(pp.begin(),pp.end());
      vector<RC>::iterator new_end = std::unique(pp.begin(),pp.end());
      pp.erase(new_end,pp.end());
   }
   for (size_t i=0; i<n; i++) {
      if (I[i]>J[i])
         _ch[I[i]-1] = std::max(static_cast<unsigned>(abs(int(I[i])-int(J[i]))),_ch[I[i]-1]);
   }
   _ch[0] = 0;
   for (size_t i=1; i<_size; ++i)
      _ch[i] += _ch[i-1] + 1;
   _length = _ch[_size-1]+1;
   _a.resize(_length);
   _diag.resize(_size);
}


template<class T_>
SkSMatrix<T_>::SkSMatrix(const Vect<size_t>& I,
                         const Vect<size_t>& J,
                         const Vect<T_>&     a,
                         int                 opt) : _dof(0)
{
   _is_diagonal = false;
   _dof_type = NODE_DOF;
   size_t i,  n=I.size();
   std::vector<RC> pp(n);
   _size = 0;
   for (i=0; i<n; i++) {
      _size = std::max(_size,I[i]);
      pp[i] = RC(I[i]-1,J[i]-1);
   }
   _ch.resize(_size,0);
   _nb_rows = _nb_cols = _size;
   if (opt==0) {
      sort(pp.begin(),pp.end());
      vector<RC>::iterator new_end = std::unique(pp.begin(),pp.end());
      pp.erase(new_end,pp.end());
   }
   for (i=0; i<n; i++)
      if (I[i]>J[i])
         _ch[I[i]-1] = std::max(static_cast<unsigned>(abs(int(I[i])-int(J[i]))),_ch[I[i]-1]);
   _ch[0] = 0;
   for (i=1; i<_size; ++i)
      _ch[i] += _ch[i-1] + 1;
   _length = _ch[_size-1]+1;
   _a.resize(_length);
   size_t k=0;
   for (i=0; i<n; i++)
      set(I[i],J[i],a[k++]);
   _diag.resize(_size);
}


template<class T_>
SkSMatrix<T_>::SkSMatrix(const SkSMatrix<T_>& m) : _dof(m._dof)
{
   _length = m._length;
   _size = m._size;
   _ch.resize(_size);
   _ch = m._ch;
   _a.resize(_length);
   _a = m._a;
   _diag.resize(_size);
   _diag = m._diag;
   _zero = T_(0);
   _theMesh = m._theMesh;
   _is_diagonal = m._is_diagonal;
}


template<class T_>
SkSMatrix<T_>::~SkSMatrix()
{
}


template<class T_>
void SkSMatrix<T_>::setMesh(Mesh&  mesh,
                            size_t dof)
{
   _dof_type = mesh.getDOFSupport();
   Matrix<T_>::init_set_mesh(mesh,dof);
   if (_dof_type==NODE_DOF) {
      if (dof)
         _length = NodeSkyline(*_theMesh,_ch,dof);
      else
         _length = NodeSkyline(*_theMesh,_ch);
   }
   else if (_dof_type==SIDE_DOF) {
      if (dof)
         _length = SideSkyline(*_theMesh,_ch,dof);
      else
         _length = SideSkyline(*_theMesh,_ch);
   }
   else if (_dof_type==ELEMENT_DOF) {
      if (dof)
         _length = ElementSkyline(*_theMesh,_ch,dof);
      else
         _length = ElementSkyline(*_theMesh,_ch);
   }
   else
      ;
   _diag.resize(_size);
   _a.resize(_length);
   _ch[0] = 0;
   for (size_t i=1; i<_size; i++)
      _ch[i] += _ch[i-1];
   _dof = 0;
   _nb_rows = _nb_cols = _size;
}


template<class T_>
void SkSMatrix<T_>::setGraph(const Vect<RC>& I,
                             int             opt)
{ }


template<class T_>
void SkSMatrix<T_>::setGraph(Mesh&  mesh,
                             int    dof_type,
                             size_t dof1,
                             size_t dof2)
{
   _theMesh = &mesh;
   _theMesh->selectDOF(dof_type,dof1,dof2);
   if (dof_type==NODE_DOF)
      _length = NodeSkyline(mesh,_ch,dof1,dof2);
   else if (dof_type==SIDE_DOF) {
      mesh.getAllSides();
      _length = SideSkyline(mesh,_ch,dof1,dof2);
   }
   else if (dof_type==ELEMENT_DOF)
      _length = ElementSkyline(mesh,_ch);
   else;
   _size = _ch.size();
   _diag.resize(_size);
   _a.resize(_length);
   _ch[0] = 0;
   for (size_t i=1; i<_size; i++)
      _ch[i] += _ch[i-1];
   _dof = 0;
   _nb_rows = _nb_cols = _size;
}


template<class T_>
void SkSMatrix<T_>::setMesh(size_t dof, 
                            Mesh&  mesh,
                            int    code)
{
   dof = 0;
   if (mesh.getDim()==0) { }
   code = 0;
   _theMesh = &mesh;
}


template<class T_>
void SkSMatrix<T_>::setMesh(size_t dof,
                            size_t nb_eq,
                            Mesh&  mesh)
{
   dof = 0;
   nb_eq = 0;
   if (mesh.getDim()==0) { }
   _theMesh = &mesh;
}


template<class T_>
void SkSMatrix<T_>::setSkyline(Mesh& mesh)
{
   _zero = 0;
   int set_sides = mesh.SidesAreDOF();
   _size = mesh.getNbEq();
   if (_dof)
      _size = mesh.getNbNodes();
   if (set_sides) {
      if (_dof) {
         _size = mesh.getNbSides();
         _length = SideSkyline(mesh,_ch,_dof);
      }
      _length = SideSkyline(mesh,_ch);
   }
   else {
      if (_dof) {
         _size = mesh.getNbNodes();
         _length = NodeSkyline(mesh,_ch,_dof);
      }
      _length = NodeSkyline(mesh,_ch);
   }
   _diag.resize(_size);
   _ch[0] = 0;
   for (size_t i=1; i<_size; i++)
      _ch[i] += _ch[i-1];
   _a.resize(_length);
   _theMesh = &mesh;
}


template<class T_>
void SkSMatrix<T_>::setDiag()
{
   for (size_t i=0; i<_size; i++)
      _diag[i] = _a[_ch[i]];
}


template<class T_>
void SkSMatrix<T_>::set(size_t    i,
                        size_t    j,
                        const T_& val)
{
   if (i>=j)
      _a[_ch[i-1]+j-i] = val;
   else
      throw OFELIException("In SkSMatrix::set(i,j,x): Index pair (" + itos(int(i)) +
                           "," + itos(int(j)) + ") not compatible with skyline symmeric storage.");
}


template<class T_>
void SkSMatrix<T_>::SSet(size_t    i,
                         size_t    j,
                         const T_& val)
{
   _a[_ch[i-1]+j-i] = val;
}


template<class T_>
void SkSMatrix<T_>::Axpy(T_                   a,
                         const SkSMatrix<T_>& m)
{
   _a += a * m._a;
}


template<class T_>
void SkSMatrix<T_>::Axpy(T_                a,
                         const Matrix<T_>* m)
{
   for (size_t i=0; i<_length; i++)
      _a[i] += a * m->_a[i];
}


template<class T_>
void SkSMatrix<T_>::MultAdd(const Vect<T_>& x,
                            Vect<T_>&       y) const
{
   for (size_t i=0; i<_size; i++) {
      for (size_t j=0; j<_size; j++)
         y[i] += operator()(i+1,j+1)*x[j];
   }
}


template<class T_>
void SkSMatrix<T_>::MultAdd(T_              a,
                            const Vect<T_>& x,
                            Vect<T_>&       y) const
{
   for (size_t i=0; i<_size; i++)
      for (size_t j=0; j<_size; j++)
         y[i] += a * operator()(i+1,j+1)*x[j];
}


template<class T_>
void SkSMatrix<T_>::Mult(const Vect<T_>& x,
                         Vect<T_>&       y) const
{
   y = static_cast<T_>(0);
   MultAdd(x,y);
}


template<class T_>
void SkSMatrix<T_>::TMult(const Vect<T_>& x,
               Vect<T_>&       y) const
{
   Mult(x,y);
}


template<class T_>
void SkSMatrix<T_>::add(size_t    i,
                        size_t    j,
                        const T_& val)
{
   if (i>=j)
      _a[_ch[i-1]+j-i] += val;
}


template<class T_>
size_t SkSMatrix<T_>::getColHeight(size_t i) const
{
   if (i==1)
      return 1;
   else
      return _ch[i-1]-_ch[i-2];
}


template<class T_>
Vect<T_> SkSMatrix<T_>::getColumn(size_t j) const
{
   Vect<T_> v(_nb_rows);
   for (size_t i=1; i<=_nb_rows; i++)
      v(i) = (*this)(i,j);
   return v;
}


template<class T_>
Vect<T_> SkSMatrix<T_>::getRow(size_t i) const
{
    return getColumn(i);
}


template<class T_>
T_& SkSMatrix<T_>::operator()(size_t i,
                              size_t j)
{
   int k=0, l=0;
   if (i>1)
      k = int(j-i+_ch[i-1]-_ch[i-2]-1);
   if (j>1)
      l = int(i-j+_ch[j-1]-_ch[j-2]-1);
   if (k>=0 && i>=j)
      return _a[_ch[i-1]+j-i];
   else if (l>=0 && i<j)
      return _a[_ch[j-1]+i-j];
   else
      return _zero;
}


template<class T_>
T_ SkSMatrix<T_>::operator()(size_t i,
                             size_t j) const
{
   int k=0, l=0;
   if (i>1)
      k = int(j-i+_ch[i-1]-_ch[i-2]-1);
   if (j>1)
      l = int(i-j+_ch[j-1]-_ch[j-2]-1);
   if (k>=0 && i>=j)
      return _a[_ch[i-1]+j-i];
   else if (l>=0 && i<j)
      return _a[_ch[j-1]+i-j];
   else
      return _zero;
}


template<class T_>
SkSMatrix<T_>& SkSMatrix<T_>::operator=(const SkSMatrix<T_>& m)
{
   _a = m._a;
   return *this;
}


template<class T_>
SkSMatrix<T_>& SkSMatrix<T_>::operator=(const T_& x)
{
   for (size_t i=0; i<_length; i++)
      _a[i] = 0;
   for (size_t i=0; i<_nb_rows; i++) {
      _diag[i] = x;
      set(i+1,i+1,x);
   }
   return *this;
}


template<class T_>
SkSMatrix<T_>& SkSMatrix<T_>::operator+=(const SkSMatrix<T_>& m)
{
   for (size_t i=0; i<_length; i++)
      _a[i] += m._a[i];
   return *this;
}


template<class T_>
SkSMatrix<T_>& SkSMatrix<T_>::operator*=(const T_& x)
{
   for (size_t i=0; i<_length; i++)
      _a[i] *= x;
   return *this;
}


template<class T_>
int SkSMatrix<T_>::Factor()
{
   return setLDLt();
}


template<class T_>
int SkSMatrix<T_>::setLDLt()
{
   if (_is_diagonal)
      return 0;
   size_t di=0, dij;
   T_ s, pivot=_a[_ch[0]];
   if (Abs(pivot) < OFELI_EPSMCH)
      throw OFELIException("In SkSMatrix::setLDLt(): The first pivot is null.");
   else
      _a[_ch[0]] = T_(1.)/pivot;
   for (size_t i=1; i<_size; i++) {
      size_t dj = 0;
      for (size_t j=di=i+1+_ch[i-1]-_ch[i]; j<i; j++) {
         if (j>0)
            dj = j+1+_ch[j-1]-_ch[j];
         dij = std::max(di,dj);
         for (size_t k=0; k<j-dij; k++)
            _a[_ch[i]+j-i] -= _a[_ch[i]+dij+k-i]*_a[_ch[j]+dij+k-j];
      }
      pivot = _a[_ch[i]];
      for (size_t k=di; k<i; k++) {
         s = _a[_ch[i]+k-i]*_a[_ch[k]];
         pivot -= s*_a[_ch[i]+k-i];
         _a[_ch[i]+k-i] = s;
      }
      if (Abs(pivot) < OFELI_EPSMCH)
         throw OFELIException("In SkSMatrix::setLDLt(): " + itos(int(i)+1) + "-th pivot is null.");
      else
         _a[_ch[i]] = T_(1.)/pivot;
   }
   return 0;
}


template<class T_>
int SkSMatrix<T_>::solveLDLt(const Vect<T_>& b,
                             Vect<T_>&       x)
{
   int ret = setLDLt();
   x = b;
   if (_is_diagonal) {
      for (size_t i=0; i<_size; i++)
         x[i] /= _a[i];
      return 0;
   }
   for (size_t i=1; i<_size; i++) {
      size_t di = i+1+_ch[i-1]-_ch[i];
      T_ s = 0;
      for (size_t j=0; j<i-di; j++)
         s += _a[_ch[i]+di+j-i] * x[di+j];
      x[i] -= s;
   }
   for (size_t i=0; i<_size; i++)
      x[i] *= _a[_ch[i]];
   for (int k=int(_size-1); k>0; k--) {
      size_t di = k+1+_ch[k-1]-_ch[k];
      for (size_t j=0; j<k-di; j++)
         x[j+di] -= x[k] * _a[_ch[k]+di+j-k];
   }
   return ret;
}


template<class T_>
int SkSMatrix<T_>::solve(Vect<T_>& b,
                         bool      fact)
{
  int ret = 0;
   if (_is_diagonal) {
      for (size_t i=0; i<_size; i++) {
         if (Abs(_a[i]) < OFELI_EPSMCH)
            throw OFELIException("In SkSMatrix::solve(b): " + itos(i+1) + "-th diagonal is null.");
         b[i] /= _a[i];
      }
      return ret;
   }
   if (fact)
      ret = setLDLt();
   for (size_t i=1; i<_size; i++) {
      size_t di = i+1+_ch[i-1]-_ch[i];
      T_ s = 0;
      for (size_t j=0; j<i-di; j++)
         s += _a[_ch[i]+di+j-i] * b[di+j];
      b[i] -= s;
   }
   for (size_t i=0; i<_size; i++)
      b[i] *= _a[_ch[i]];
   for (int k=int(_size-1); k>0; k--) {
      size_t di = k+1+_ch[k-1]-_ch[k];
      for (size_t j=0; j<k-di; j++)
         b[j+di] -= b[k] * _a[_ch[k]+di+j-k];
   }
   return ret;
}


template<class T_>
int SkSMatrix<T_>::solve(const Vect<T_>& b,
                         Vect<T_>&       x,
                         bool            fact)
{
   x = b;
   return solve(x,fact);
}


template<class T_>
T_ *SkSMatrix<T_>::get() const
{
   return _a;
}


template<class T_>
void SkSMatrix<T_>::set(size_t i,
                        T_     x)
{
   _a[i] = x;
}


template<class T_>
T_ SkSMatrix<T_>::get(size_t i,
                      size_t j) const
{
   if (i>=j)
      return _a[_ch[i-1]+j-i];
   else
      return _a[_ch[j-1]+i-j];
}


///////////////////////////////////////////////////////////////////////////////
//                           ASSOCIATED  FUNCTIONS                           //
///////////////////////////////////////////////////////////////////////////////


template<class T_>
Vect<T_> operator*(const SkSMatrix<T_>& A,
                   const Vect<T_>&      b)
{
   Vect<T_> v(b.size());
   A.Mult(b,v);
}


template<class T_>
ostream& operator<<(ostream&             s,
                    const SkSMatrix<T_>& a)
{
   size_t nr=a.getNbRows();
   if (Verbosity<=1) {
      s << "Matrix dimension: " << nr << endl;
      return s;
   }
   if (Verbosity==2)
      nr = !(10<nr)?nr:10;
   else if (Verbosity==3)
      nr = !(50<nr)?nr:50;
   s.setf(ios::right|ios::scientific);
   s << endl;
   for (size_t i=1; i<=nr; i++) {
      s << "\nRow:  " << setw(6) << i << endl;
      for (size_t j=1; j<=nr; j++)
          s << "  " << setprecision(8) << std::setfill(' ') << setw(18) << a(i,j);
      s << endl;
   }
   return s;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */

#endif
