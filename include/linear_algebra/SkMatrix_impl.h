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

                   Implementation of class SkMatrix

  ==============================================================================*/


#ifndef __SKMATRIX_IMPL_H
#define __SKMATRIX_IMPL_H

#include "linear_algebra/SkMatrix.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/Matrix_impl.h"

namespace OFELI {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

class Mesh;


template<class T_>
SkMatrix<T_>::SkMatrix()
{
   _dof = 0;
   _is_diagonal = false;
}


template<class T_>
SkMatrix<T_>::SkMatrix(size_t size,
                       int    is_diagonal)
{
   _dof = 0;
   _zero = 0;
   _is_diagonal = is_diagonal;
   _dof_type = NODE_DOF;
   _nb_rows = _nb_cols = _size = size;
   _ch.resize(size);
   _diag.resize(_size);
   _ch[0] = 0;
   for (size_t i=1; i<_size; i++)
      _ch[i] = _ch[i-1] + i + 1;
   if (_is_diagonal) {
      for (size_t i=1; i<_size; i++)
         _ch[i] = _ch[i-1] + 1;
   }
   _length = _ch[_size-1] + 1;
   _a.resize(_length);
   _aU.resize(_length);
}


template<class T_>
SkMatrix<T_>::SkMatrix(Mesh&  mesh,
                       size_t dof,
                       int    is_diagonal)
{
   _is_diagonal = is_diagonal;
   setMesh(mesh,dof);
}


template<class T_>
SkMatrix<T_>::SkMatrix(const Vect<size_t> &ColHt) : _dof(0)
{
   _is_diagonal = false;
   _zero = 0;
   _size = ColHt.size();
   _ch.resize(_size,0);
   for (size_t i=1; i<_size; i++)
      _ch[i] = _ch[i-1] + ColHt[i];
   _length = _ch[_size-1] + 1;
   _a.resize(_length);
   _aU.resize(_length);
   _diag.resize(_size);
}


template<class T_>
SkMatrix<T_>::SkMatrix(const SkMatrix<T_>& m)
{
   _is_diagonal = m._is_diagonal;
   _size = m._size;
   _length = m._length;
   _ch.resize(_size);
   _ch = m._ch;
   _diag.resize(_size);
   _diag = m._diag;
   _a.resize(_length);
   _aU.resize(_length);
   _a = m._a;
   _aU = m._aU;
   _dof = m._dof;
   _zero = static_cast<T_>(0);
}


template<class T_>
SkMatrix<T_>::~SkMatrix()
{ }


template<class T_>
void SkMatrix<T_>::setMesh(Mesh&  mesh,
                           size_t dof)
{
   _dof_type = mesh.getDOFSupport();
   Matrix<T_>::init_set_mesh(mesh,dof);
   if (_dof_type==NODE_DOF) {
      if (dof)
         _length = NodeSkyline(mesh,_ch,dof);
      else
         _length = NodeSkyline(mesh,_ch);
   }
   else if (_dof_type==SIDE_DOF) {
      if (dof)
         _length = SideSkyline(mesh,_ch,dof);
      else
         _length = SideSkyline(mesh,_ch);
   }
   else if (_dof_type==ELEMENT_DOF) {
      if (dof)
         _length = ElementSkyline(mesh,_ch,dof);
      else
         _length = ElementSkyline(mesh,_ch);
   }
   else
      ;
   _diag.resize(_size);
   _a.resize(_length);
   _aU.resize(_length);
   _ch[0] = 0;
   for (size_t i=1; i<_size; i++)
       _ch[i] += _ch[i-1];
    _zero = T_(0);
}


template<class T_>
void SkMatrix<T_>::setGraph(const Vect<RC>& I,
                            int             opt)
{ }


template<class T_>
void SkMatrix<T_>::setMesh(size_t dof,
                           Mesh&  mesh,
                           int    code)
{
// This is just to avoid warning on unused variable
   dof = 0;
   if (mesh.getDim()==0) { }
   code = 0;
   _theMesh = &mesh;
}


template<class T_>
void SkMatrix<T_>::setMesh(size_t dof,
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
void SkMatrix<T_>::setSkyline(Mesh& mesh)
{
   _zero = 0;
   int set_sides = mesh.SidesAreDOF();
   _size = mesh.getNbEq();
   _theMesh = &mesh;
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
         _length = SideSkyline(mesh,_ch,_dof);
      }
      _length = SideSkyline(mesh,_ch);
   }
   _diag.resize(_size);
   _ch[0] = 0;
   for (size_t i=1; i<_size; i++)
      _ch[i] += _ch[i-1];
   _a.resize(_length);
   _aU.resize(_length);
}


template<class T_>
void SkMatrix<T_>::setDiag()
{
   for (size_t i=0; i<_size; i++)
      _diag[i] = _aU[_ch[i]];
}


template<class T_>
void SkMatrix<T_>::setDOF(size_t i)
{
   _dof = i;
}


template<class T_>
void SkMatrix<T_>::set(size_t    i,
                       size_t    j,
                       const T_& val)
{
   int k=0, l=0;
   if (i>1)
      k = int(j-i+_ch[i-1]-_ch[i-2]-1);
   if (j>1)
      l = int(i-j+_ch[j-1]-_ch[j-2]-1);
   if (k>=0 && i>j)
      _a[_ch[i-1]+j-i] = val;
   else if (l>=0 && i<=j)
      _aU[_ch[j-1]+i-j] = val;
   else
      throw OFELIException("In SkMatrix::Set(i,j,x): Index pair: ("+to_string(i)+"," +
                            to_string(j)+") not " + "compatible with skyline symmeric storage.");
}


template<class T_>
void SkMatrix<T_>::SSet(size_t    i,
                        size_t    j,
                        const T_& val)
{
   int k=0, l=0;
   if (i>1)
      k = j-i+_ch[i-1]-_ch[i-2]-1;
   if (j>1)
      l = i-j+_ch[j-1]-_ch[j-2]-1;
   if (k>=0 && i>j)
      _a[_ch[i-1]+j-i] = val;
   else if (l>=0 && i<=j)
      _aU[_ch[j-1]+i-j] = val;
}


template<class T_>
void SkMatrix<T_>::Axpy(T_                  a,
                        const SkMatrix<T_>& m)
{
   _a  += a * m._a;
   _aU += a * m._aU;
}


template<class T_>
void SkMatrix<T_>::Axpy(T_                a,
                        const Matrix<T_>* m)
{
   for (size_t i=0; i<_length; i++) {
      _a[i]  += a * m->_a[i];
      _aU[i] += a * m->_aU[i];
   }
}


template<class T_>
void SkMatrix<T_>::MultAdd(const Vect<T_>& x,
                           Vect<T_>&       y) const
{
   for (size_t i=0; i<_size; i++)
      for (size_t j=0; j<_size; j++)
         y[i] += operator()(i+1,j+1)*x[j];
}


template<class T_>
void SkMatrix<T_>::TMultAdd(const Vect<T_>& x,
                            Vect<T_>&       y) const
{
   cerr << "TMultAdd is not implemented for class SkMatrix" << endl;
}


template<class T_>
void SkMatrix<T_>::MultAdd(T_              a,
                           const Vect<T_>& x,
                           Vect<T_>&       y) const
{
   for (size_t i=0; i<_size; i++)
      for (size_t j=0; j<_size; j++)
         y[i] += a * operator()(i+1,j+1)*x[j];
}


template<class T_>
void SkMatrix<T_>::Mult(const Vect<T_>& x,
                        Vect<T_>&       y) const
{
   y = static_cast<T_>(0);
   MultAdd(x,y);
}


template<class T_>
void SkMatrix<T_>::TMult(const Vect<T_>& x,
                         Vect<T_>&       y) const
{
   y = static_cast<T_>(0);
   TMultAdd(x,y);
}


template<class T_>
void SkMatrix<T_>::add(size_t    i,
                       size_t    j,
                       const T_& val)
{
   if (i>j)
      _a[_ch[i-1]+j-i] += val;
   else if (i<=j)
      _aU[_ch[j-1]+i-j] += val;
}


template<class T_>
size_t SkMatrix<T_>::getColHeight(size_t i) const
{
   if (i==1)
      return 1;
   else
      return _ch[i-1]-_ch[i-2];
}


template<class T_>
T_ SkMatrix<T_>::operator()(size_t i,
                            size_t j) const
{
   int k=0, l=0;
   if (i>1)
      k = int(j-i+_ch[i-1]-_ch[i-2]-1);
   if (j>1)
      l = int(i-j+_ch[j-1]-_ch[j-2]-1);
   if (k>=0 && i>j)
      return _a[_ch[i-1]+j-i];
   else if (l>=0 && i<=j)
      return _aU[_ch[j-1]+i-j];
   else
      return _zero;
}


template<class T_>
T_& SkMatrix<T_>::operator()(size_t i,
                             size_t j)
{
   int k=0, l=0;
   if (i>1)
      k = int(j-i+_ch[i-1]-_ch[i-2]-1);
   if (j>1)
      l = int(i-j+_ch[j-1]-_ch[j-2]-1);
   if (k>=0 && i>j)
      return _a[_ch[i-1]+j-i];
   else if (l>=0 && i<=j)
      return _aU[_ch[j-1]+i-j];
   else
      throw OFELIException("In SkMatrix::Operator(): Index pair (" + to_string(i) + "," +
                            to_string(j) + ") not compatible with skyline structure");
   return _temp;
}


template<class T_>
void SkMatrix<T_>::DiagPrescribe(Mesh&           mesh,
                                 Vect<T_>&       b,
                                 const Vect<T_>& u,
                                 int             flag)
{
   real_t p=0;
   for (size_t l=0; l<_size; l++)
      p = std::max(p,_aU[_ch[l]]);
   node_loop(&mesh) {
      for (size_t i=1; i<=The_node.getNbDOF(); ++i) {
         if (The_node.getCode(i)>0) {
            size_t ii = The_node.getDOF(i)-1;
            for (size_t j=ii+1+_ch[ii-1]-_ch[ii]; j<=ii; j++) {
               b[ii] = p*u[ii];
               _a[_ch[ii]+j-ii] = _aU[_ch[ii]+j-ii] = 0;
            }
            _diag[ii] = _aU[_ch[ii]] = p;
         }
      }
   }
}


template<class T_>
void SkMatrix<T_>::DiagPrescribe(Vect<T_>&       b,
                                 const Vect<T_>& u,
                                 int             flag)
{
   real_t p=0;
   for (size_t l=0; l<_size; l++)
      p = std::max(p,_aU[_ch[l]]);
   MESH_ND {
      for (size_t i=1; i<=The_node.getNbDOF(); ++i) {
         if (The_node.getCode(i)>0) {
            size_t ii = The_node.getDOF(i)-1;
            for (size_t j=ii+1+_ch[ii-1]-_ch[ii]; j<=ii; j++) {
               b[ii] = p*u[ii];
               _a[_ch[ii]+j-ii] = _aU[_ch[ii]+j-ii] = 0;
            }
            _diag[ii] = _aU[_ch[ii]] = p;
         }
      }
   }
}


template<class T_>
SkMatrix<T_>& SkMatrix<T_>::operator=(const SkMatrix<T_>& m)
{
   _a = m._a;
   _aU = m._aU;
   return *this;
}


template<class T_>
SkMatrix<T_>& SkMatrix<T_>::operator=(const T_& x)
{
   for (size_t i=0; i<_length; i++)
      _a[i]  = _aU[i] = 0;
   for (size_t i=0; i<_nb_rows; i++) {
      _diag[i] = x;
      set(i+1,i+1,x);
   }
   return *this;
}


template<class T_>
SkMatrix<T_>& SkMatrix<T_>::operator+=(const SkMatrix<T_>& m)
{
   _a  += m._a;
   _aU += m._aU;
   return *this;
}


template<class T_>
SkMatrix<T_>& SkMatrix<T_>::operator+=(const T_& x)
{
   _a += x;
   _aU += x;
   return *this;
}


template<class T_>
SkMatrix<T_>& SkMatrix<T_>::operator*=(const T_& x)
{
   _a  *= x;
   _aU *= x;
   return *this;
}


template<class T_>
int SkMatrix<T_>::Factor()
{
   return setLU();
}


template<class T_>
int SkMatrix<T_>::setLU()
{
   if (_is_diagonal)
       return 0;
   size_t k, di, dij, i=0;
   if (Abs(_aU[_ch[0]]) < OFELI_EPSMCH)
      throw OFELIException("In SkMatrix::Factor(): First pivot is null.");
   for (i=1; i<_size; i++) {
      size_t dj = 0;
      for (size_t j=di=i+1+_ch[i-1]-_ch[i]; j<i; j++) {
         if (j>0)
            dj = j+1+_ch[j-1]-_ch[j];
         dij = std::max(di,dj);
         for (k=0; k<j-dij; k++)
            _a[_ch[i]+j-i] -= _a[_ch[i]+dij+k-i]*_aU[_ch[j]+dij+k-j];
         _a[_ch[i]+j-i] /= _aU[_ch[j]];
         for (k=0; k<j-dij; k++)
            _aU[_ch[i]+j-i] -= _a[_ch[j]+dij+k-j]*_aU[_ch[i]+dij+k-i];
      }
      for (k=0; k<i-di; k++)
         _aU[_ch[i]] -= _a[_ch[i]+k+di-i]*_aU[_ch[i]+k+di-i];
      if (Abs(_aU[_ch[i]]) < OFELI_EPSMCH)
         throw OFELIException("In SkMatrix::Factor(): " + to_string(i+1) + "-th pivot is null.");
   }
   return 0;
}


template<class T_>
int SkMatrix<T_>::solve(Vect<T_>& b,
                        bool      fact)
{
   int ret = 0;
   if (_is_diagonal) {
      for (size_t i=0; i<_size; i++) {
         if (Abs(_aU[i]) < OFELI_EPSMCH)
            throw OFELIException("In SkMatrix::solve(b): " + to_string(i+1) + "-th diagonal is null.");
         b[i] /= _aU[i];
      }
      return ret;
   }
   size_t di;
   if (fact)
      ret = setLU();
   size_t i, j;
   for (i=1; i<_size; i++) {
      di = i+1+_ch[i-1]-_ch[i];
      T_ s = 0;
      for (j=0; j<i-di; j++)
         s += _a[_ch[i]+di+j-i] * b[di+j];
      b[i] -= s;
   }
   for (int k=int(_size-1); k>0; k--) {
      if (Abs(_aU[_ch[k]]) < OFELI_EPSMCH)
         throw OFELIException("In SkMatrix::solve(b): " + to_string(k+1) + "-th pivot is null.");
      b[k] /= _aU[_ch[k]];
      di = k+1+_ch[k-1]-_ch[k];
      for (j=0; j<k-di; j++)
         b[j+di] -= b[k] * _aU[_ch[k]+di+j-k];
   }
   b[0] /= _aU[_ch[0]];
   return ret;
}


template<class T_>
int SkMatrix<T_>::solve(const Vect<T_>& b,
                        Vect<T_>&       x,
                        bool            fact)
{
       x = b;
       return solve(x,fact);
}


template<class T_>
int SkMatrix<T_>::solveLU(const Vect<T_>& b,
                          Vect<T_>&       x)
{
   int ret = setLU();
   x = b;
   if (_is_diagonal) {
      for (size_t i=0; i<_size; i++)
         x[i] /= _aU[i];
      return 0;
   }
   size_t di, i, j;
   for (i=1; i<_size; i++) {
      di = i+1+_ch[i-1]-_ch[i];
      T_ s = 0;
      for (j=0; j<i-di; j++)
         s += _a[_ch[i]+di+j-i] * x[di+j];
      x[i] -= s;
   }
   for (int k=int(_size-1); k>0; k--) {
      if (Abs(_aU[_ch[k]]) < OFELI_EPSMCH)
         throw OFELIException("In SkMatrix::solveLU(b,x): " + to_string(k+1) + "-th pivot is null.");
      x[k] /= _aU[_ch[k]];
      di = k+1+_ch[k-1]-_ch[k];
      for (j=0; j<k-di; j++)
         x[j+di] -= b[k] * _aU[_ch[k]+di+j-k];
   }
   x[0] /= _aU[_ch[0]];
   return ret;
}


template<class T_>
T_* SkMatrix<T_>::get() const
{
   return _a;
}


template<class T_>
T_ SkMatrix<T_>::get(size_t i,
                     size_t j) const
{
   if (i>j)
      return _a[_ch[i-1]+j-i];
   else
      return _aU[_ch[j-1]+i-j];
}


///////////////////////////////////////////////////////////////////////////////
//                           ASSOCIATED  FUNCTIONS                           //
///////////////////////////////////////////////////////////////////////////////


template<class T_>
Vect<T_> operator*(const SkMatrix<T_>& A,
                   const Vect<T_>&     b)
{
   Vect<T_> v(b.size());
   A.Mult(b,v);
}


template<class T_>
ostream& operator<<(ostream&            s,
                    const SkMatrix<T_>& A)
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
