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

                         Implementation of class SpMatrix

  ==============================================================================*/


#ifndef __SPMATRIX_IMPL_H
#define __SPMATRIX_IMPL_H

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include "linear_algebra/SpMatrix.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/GraphOfMatrix.h"
#include "solvers/LinearSolver.h"

#ifdef USE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Sparse>

typedef Eigen::Matrix<T_,Eigen::Dynamic,1> VectorX;
typedef SparseMatrix<T_>                   SpMat;
typedef Triplet<real_t>                    Tr;
#endif

namespace OFELI {

template<class T_> class LinearSolver;

template<class T_>
SpMatrix<T_>::SpMatrix()
             : _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
{
   _max_it = 1000;
   _toler = 1.e-8;
   _is_diagonal = 0;
}


template<class T_>
SpMatrix<T_>::SpMatrix(size_t nr,
                       size_t nc)
             : _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
{
   _is_diagonal = 0;
   _nb_rows = nr;
   _nb_cols = nc;
   _size = 0;
   if (_nb_rows==_nb_cols)
      _size = _nb_rows;
   _length = _nb_rows*_nb_cols;
#ifdef USE_EIGEN
   _A.resize(_nb_rows,_nb_cols);
#else
   _row_ptr.resize(_nb_rows+1);
   _col_ind.resize(_length);
   _row_ptr[0] = 0;
   for (size_t i=0; i<_nb_rows; i++)
      _row_ptr[i+1] = _row_ptr[i] + _nb_cols;
   size_t l=0;
   for (size_t i=0; i<_nb_rows; i++)
      for (size_t j=0; j<_nb_cols; j++)
         _col_ind[l++] = j+1;
   _a.resize(_length,T_(0.));
#endif
   _max_it = 1000;
   _toler = 1.e-8;
}


template<class T_>
SpMatrix<T_>::SpMatrix(size_t size,
                       int    is_diagonal)
             : _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
{
   _is_diagonal = is_diagonal;
   _size = _nb_rows = _nb_cols = size;
#ifdef USE_EIGEN
   _A.resize(_nb_rows,_nb_cols);
#else
   _length = _size*_size;
   _row_ptr.resize(_nb_rows+1);
   _col_ind.resize(_length);
   _row_ptr[0] = 0;
   for (size_t i=0; i<_nb_rows; i++)
      _row_ptr[i+1] = _row_ptr[i] + _nb_cols;
   size_t l = 0;
   for (size_t i=0; i<_nb_rows; i++)
      for (size_t j=0; j<_nb_cols; j++)
         _col_ind[l++] = j+1;
   _a.resize(_length,T_(0.));
#endif
   _max_it = 1000;
   _toler = 1.e-8;
}


template<class T_>
SpMatrix<T_>::SpMatrix(Mesh&  mesh,
                       size_t dof,
                       int    is_diagonal)
             : _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
{
   _is_diagonal = is_diagonal;
   setMesh(mesh,dof);
   _max_it = 1000;
   _toler = 1.e-8;
}


template<class T_>
SpMatrix<T_>::SpMatrix(size_t dof,
                       Mesh&  mesh,
                       int    code)
             : _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
{
   _is_diagonal = false;
   setMesh(dof,mesh,mesh.getDOFSupport());
   _max_it = 1000;
   _toler = 1.e-8;
}

#ifndef USE_EIGEN
template<class T_>
SpMatrix<T_>::SpMatrix(size_t dof,
                       size_t nb_eq,
                       Mesh&  mesh)
             : _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
{
   _is_diagonal = false;
   setMesh(dof,nb_eq,mesh);
   _max_it = 1000;
   _toler = 1.e-8;
}
#endif

#ifndef USE_EIGEN
template<class T_>
SpMatrix<T_>::SpMatrix(const Vect<RC>& I,
                       int             opt)
             : _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
{
   _is_diagonal = false;
   setGraph(I,opt);
   _max_it = 1000;
   _toler = 1.e-8;
}
#endif

#ifndef USE_EIGEN
template<class T_>
SpMatrix<T_>::SpMatrix(const Vect<RC>& I,
                       const Vect<T_>& a,
                       int             opt)
             : _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
{
   _is_diagonal = false;
   size_t n=I.size();
   _nb_rows = _nb_cols = 0;
   for (size_t i=0; i<n; i++) {
      _nb_rows = std::max(_nb_rows,I[i].first);
      _nb_cols = std::max(_nb_cols,I[i].second);
      _IJ.push_back(RC(I[i].first-1,I[i].second-1));
   }
   if (_nb_rows==_nb_cols)
      _size = _nb_rows;
   if (opt==0) {
      sort(_IJ.begin(),_IJ.end());
      vector<RC>::iterator new_end=unique(_IJ.begin(),_IJ.end());
      _IJ.erase(new_end,_IJ.end());
   }
   _row_ptr.resize(_size+1);
   _col_ind.resize(n);
   StoreGraph(_size,_IJ,_row_ptr,_col_ind);
   _length = _IJ.size();
   _a.resize(_length);
   for (size_t j=0; j<n; j++)
      _a[_row_ptr[I[j].first-1]+_col_index(I[j].first,I[j].second)] = a[j];
   _max_it = 1000;
   _toler = 1.e-8;
}
#endif

#ifndef USE_EIGEN
template<class T_>
SpMatrix<T_>::SpMatrix(size_t                nr,
                       size_t                nc,
                       const vector<size_t>& row_ptr,
                       const vector<size_t>& col_ind)
             : _type(0), _dof(0), _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
{
   _is_diagonal = false;
   _nb_rows = nr;
   _nb_cols = nc;
   _size = 0;
   _length = col_ind.size();
   _row_ptr.resize(_nb_rows+1);
   _row_ptr = row_ptr;
   _col_ind.resize(_length);
   _col_ind = col_ind;
   _a.resize(_length,T_(0.));
   _diag.resize(_size);
   _max_it = 1000;
   _toler = 1.e-8;
}
#endif

#ifndef USE_EIGEN
template<class T_>
SpMatrix<T_>::SpMatrix(size_t                nr,
                       size_t                nc,
                       const vector<size_t>& row_ptr,
                       const vector<size_t>& col_ind,
                       const vector<T_>&     a)
             : _type(0), _dof(0), _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
{
   _is_diagonal = false;
   _size = 0;
   _nb_rows = nr; _nb_cols = nc;
   _length = col_ind.size();
   _row_ptr.resize(_nb_rows+1);
   _row_ptr = row_ptr;
   _col_ind.resize(_length);
   _col_ind = col_ind;
   _a.resize(_length);
   _a = a;
   _diag.resize(_size);
   _max_it = 1000;
   _toler = 1.e-8;
}
#endif

#ifndef USE_EIGEN
template<class T_>
SpMatrix<T_>::SpMatrix(const vector<size_t>& row_ptr,
                       const vector<size_t>& col_ind)
             : _type(0), _dof(0), _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
{
   _is_diagonal = false;
   _size = _nb_rows = _nb_cols = row_ptr.size()-1;
   _length = col_ind.size();
   _row_ptr.resize(_size+1);
   _row_ptr = row_ptr;
   _col_ind.resize(_length);
   _col_ind = col_ind;
   _a.resize(_length,T_(0.));
   _diag.resize(_size);
   _max_it = 1000;
   _toler = 1.e-8;
}
#endif

#ifndef USE_EIGEN
template<class T_>
SpMatrix<T_>::SpMatrix(const vector<size_t>& row_ptr,
                       const vector<size_t>& col_ind,
                       const vector<T_>&     a)
             : _type(0), _dof(0), _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
{
   _is_diagonal = false;
   _size = _nb_rows = _nb_cols = row_ptr.size()-1;
   _length = col_ind.size();
   _row_ptr.resize(_size+1);
   _row_ptr = row_ptr;
   _col_ind.resize(_length);
   _col_ind = col_ind;
   _a.resize(_length);
   _a = a;
   _diag.resize(_size);
   _max_it = 1000;
   _toler = 1.e-8;
}
#endif


template<class T_>
SpMatrix<T_>::SpMatrix(const SpMatrix& m)
{
   _is_diagonal = m._is_diagonal;
   _dof = m._dof;
   _size = m._size;
   _is_dense = m._is_dense;
   _extended = m._extended;
#ifdef USE_EIGEN
   if (_nb_rows)
      _A.resize(_nb_rows,_nb_cols);
#else
   _length = m._length;
   _row_ptr.resize(_size+1);
   _col_ind.resize(_length);
   _col_ind = m._col_ind;
   _row_ptr = m._row_ptr;
   _a.resize(_length);
   _a = m._a;
#endif
   _diag.resize(_size);
   _diag = m._diag;
   _solver = -1;
   _prec = -1;
   _max_it = 1000;
   _toler = 1.e-8;
   _solver = -1;
   _prec = -1;
   _max_it = 1000;
   _toler = 1.e-8;
}


template<class T_>
SpMatrix<T_>::~SpMatrix()
{ }


template<class T_>
void SpMatrix<T_>::Dense()
{
   _is_dense = 1;
#ifdef USE_EIGEN
   _length = _nb_rows*_nb_cols;
   _A.reserve(_length);
   clear();
#endif
}


template<class T_>
void SpMatrix<T_>::Identity()
{
#ifdef USE_EIGEN
   _A.reserve(_length);
   for (size_t i=0; i<_nb_rows; i++)
      _A.coeffRef(i,i) = static_cast<T_>(1.);
#else
   for (size_t i=0; i<_nb_rows; ++i) {
      for (size_t j=0; j<_row_ptr[i+1]-_row_ptr[i]; ++j)
         _a[_row_ptr[i]+j-1] = static_cast<T_>(0.);
      _a[_row_ptr[i]+i-1] = static_cast<T_>(1.);
   }
#endif
}


template<class T_>
void SpMatrix<T_>::Diagonal()
{
   _is_dense = 0;
   _is_diagonal = 1;
#ifdef USE_EIGEN
   _length = _nb_rows;
   _A.reserve(_length);
   for (size_t i=0; i<_nb_rows; i++)
      _A.insert(i,i) = 0;
#else
   for (size_t i=0; i<_length; i++)
      _a[i] = static_cast<T_>(0);
#endif
}


template<class T_>
void SpMatrix<T_>::Diagonal(const T_& a)
{
#ifdef USE_EIGEN
   _A.reserve(_length);
   for (size_t i=0; i<_nb_rows; i++)
      _A.coeffRef(i,i) = a;
#else
   for (size_t i=0; i<_nb_rows; ++i) {
      for (size_t j=0; j<_row_ptr[i+1]-_row_ptr[i]; ++j)
         _a[_row_ptr[i]+j-1] = 0;
      _a[_row_ptr[i]+i-1] = a;
   }
#endif
}




template<>
inline void SpMatrix<real_t>::Laplace2D(size_t nx,
                                        size_t ny)
{
   _is_dense = 0;
   _is_diagonal = 0;
   _size = _nb_rows = _nb_cols = nx*ny;
#ifdef USE_EIGEN
   _length = 5*_size;
   _A.reserve(_length);
   for (size_t ii=0; ii<_size; ii++) {
      _A.insert(ii,ii) =  4.;
      size_t i=ii/ny, j=ii-i*ny;
      if (i>0)
         _A.insert(ii,ii-ny) = -1.;
      if (i<nx-1)
         _A.insert(ii,ii+ny) = -1.;
      if (j>0)
         _A.insert(ii,ii- 1) = -1.;
      if (j<ny-1)
         _A.insert(ii,ii+ 1) = -1.;
   }
#else
   _row_ptr.push_back(0);
   _length = 0;
   for (size_t ii=0; ii<_size; ii++) {
      size_t i=ii/ny, j=ii-i*ny;
      if (i>0) {
         _a.push_back(-1.);
         _col_ind.push_back(ii-ny+1);
         _length++;
      }
      if (j>0) {
         _a.push_back(-1.);
         _col_ind.push_back(ii);
         _length++;
      }
      _a.push_back(4.);
      _col_ind.push_back(ii+1);
      _length++;
      if (j<ny-1) {
         _a.push_back(-1.);
         _col_ind.push_back(ii+2);
         _length++;
      }
      if (i<nx-1) {
         _a.push_back(-1.);
         _col_ind.push_back(ii+ny+1);
         _length++;
      }
      _row_ptr.push_back(_length);
   }
#endif
}


template<class T_>
void SpMatrix<T_>::setMesh(Mesh&  mesh,
                           size_t dof)
{
   _dof_type = mesh.getDOFSupport();
   Matrix<T_>::init_set_mesh(mesh,dof);
   if (_dof_type==NODE_DOF) {
      if (_extended)
         _length = XGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
      else if (dof)
         _length = NodeGraphScal(mesh,_row_ptr,_col_ind,_IJ,_nbc);
      else
         _length = NodeGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
   }
   else if (_dof_type==SIDE_DOF)
      _length = SideGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
   else if (_dof_type==ELEMENT_DOF)
      _length = ElementGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
   else
      ;
#ifdef USE_EIGEN
   _A.resize(_size,_size);
   _A.reserve(_nbc);
   clear();
#else
   _a.resize(_length,static_cast<T_>(0));
#endif
   _diag.resize(_size);
}

template<class T_>
void SpMatrix<T_>::setMesh(size_t dof,
                           Mesh&  mesh,
                           int    code)
{
   _size = _nb_rows = _nb_cols = mesh.getNbEq();
   if (dof)
      _size = _nb_rows = _nb_cols = mesh.getNbNodes();
   if (code!=0)
      _length = XGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
   else
      _length = NodeGraphScal(mesh,_row_ptr,_col_ind,_IJ,_nbc);
#ifdef USE_EIGEN
   _A.resize(_size,_size);
   _A.reserve(_nbc);
   clear();
#else
   _a.resize(_length,T_(0.));
#endif
   _diag.resize(_size);
}


template<class T_>
void SpMatrix<T_>::setMesh(size_t dof, 
                           size_t nb_eq,
                           Mesh&  mesh)
{
   _type = 0;
   _dof = 0;
   _size = _nb_rows = _nb_cols = nb_eq;
   _length = NodeGraphScal(mesh,dof,nb_eq,_row_ptr,_col_ind,_IJ,_nbc);
#ifdef USE_EIGEN
   _A.resize(_size,_size);
   _A.reserve(_nbc);
   clear();
#else
   _a.resize(_length,T_(0.));
#endif
   _diag.resize(_size);
}


template<class T_>
void SpMatrix<T_>::setMesh(Mesh&  mesh,
                           size_t dof,
                           size_t type)
{
   _dof_type = mesh.getDOFSupport();
   Matrix<T_>::init_set_mesh(mesh,dof);
   if (_dof_type==NODE_DOF) {
      if (type && dof)
         _length = XGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
      else if (dof)
         _length = NodeGraphScal(mesh,_row_ptr,_col_ind,_IJ,_nbc);
      else if (type==0 && dof==0)
         _length = NodeGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
   }
   else if (_dof_type==SIDE_DOF)
      _length = SideGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
   else if (_dof_type==ELEMENT_DOF)
      _length = ElementGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
   else
      ;
#ifdef USE_EIGEN
   _A.resize(_size,_size);
   _A.reserve(_nbc);
   clear();
#else
   _a.resize(_length,T_(0.));
#endif
   _diag.resize(_size);
}


template<class T_>
void SpMatrix<T_>::setExtendedGraph()
{
   _extended = true;
}


template<class T_>
void SpMatrix<T_>::setOneDOF()
{
   _one_dof = true;
}


template<class T_>
void SpMatrix<T_>::setSides()
{
   _sides = true;
}


template<class T_>
void SpMatrix<T_>::setDiag()
{
#ifdef USE_EIGEN
   for (size_t i=0; i<_nb_rows; i++)
      _diag[i] = _A.coeff(i,i);
#else
   for (size_t i=0; i<_size; i++) {
      for (size_t j=0; j<_row_ptr[i+1]-_row_ptr[i]; j++) {
         if (i==j) {
            _diag[i] = _a[_row_ptr[i]+i-1];
            break;
         }
      }
   }
#endif
}


template<class T_>
void SpMatrix<T_>::DiagPrescribe(Mesh&           mesh,
                                 Vect<T_>&       b,
                                 const Vect<T_>& u)
{
   real_t p = _zero;
   for (size_t j=1; j<=_nb_rows; j++)
      p = std::max(p,Abs(get(j,j)));
#if !defined(USE_EIGEN)
   size_t k=0;
   node_loop(&mesh) {
      for (size_t i=1; i<=The_node.getNbDOF(); ++i) {
         size_t ii=The_node.getDOF(i)-1;
         for (size_t j=0; j<_row_ptr[ii+1]-_row_ptr[ii]; j++,k++)
            if (The_node.getCode(i)>0) {
               b[ii] = p*u[ii];
               _a[k] = 0;
               if (ii+1==_col_ind[k])
                  _a[k] = p;
            }
      }
   }
#endif
}


template<class T_>
void SpMatrix<T_>::DiagPrescribe(Vect<T_>&       b,
                                 const Vect<T_>& u)
{
#if !defined(USE_EIGEN)
   real_t p = 0;
   for (size_t j=1; j<=_nb_rows; j++)
      p = std::max(p,Abs(get(j,j)));
   size_t k=0;
   MESH_ND {
      for (size_t i=1; i<=The_node.getNbDOF(); ++i) {
        size_t ii=The_node.getDOF(i)-1;
         for (size_t j=0; j<_row_ptr[ii+1]-_row_ptr[ii]; j++,k++)
            if (The_node.getCode(i)>0) {
               b[ii] = p*u[ii];
               _a[k] = 0;
               if (ii+1==_col_ind[k])
                  _a[k] = p;
            }
      }
   }
#endif
}


template<class T_>
void SpMatrix<T_>::setSize(size_t size)
{
   _nb_rows = _nb_cols = _size = size;
#ifdef USE_EIGEN
   _A.resize(size,size);
#else
   _length = lsize_t(_nb_rows*_nb_cols);
   _row_ptr.resize(_nb_rows+1);
   _col_ind.resize(_length);
   _row_ptr[0] = 0;
   for (size_t i=1; i<=_nb_rows; i++)
      _row_ptr[i] = _row_ptr[i-1] + _nb_cols;
   size_t l = 0;
   for (size_t i=0; i<_nb_rows; i++)
      for (size_t j=0; j<_nb_cols; j++)
         _col_ind[l++] = j+1;
   _a.resize(_length,T_(0.));
#endif
}


template<class T_>
void SpMatrix<T_>::setSize(size_t nr,
                           size_t nc)
{
   _nb_rows = nr;
   _nb_cols = nc;
   _size = 0;
   if (_nb_rows==_nb_cols)
      _size = _nb_rows;
   _length = lsize_t(_nb_rows*_nb_cols);
#ifdef USE_EIGEN
   _A.resize(nr,nc);
#else
   _row_ptr.resize(_nb_rows+1);
   _col_ind.resize(_length);
   _row_ptr[0] = 1;
   for (size_t i=1; i<=_nb_rows; i++)
      _row_ptr[i] = _row_ptr[i-1] + _nb_cols;
   size_t l = 0;
   for (size_t i=0; i<_nb_rows; i++)
      for (size_t j=0; j<_nb_cols; j++)
         _col_ind[l++] = j + 1;
   _a.resize(_length,T_(0.));
#endif
}


template<class T_>
void SpMatrix<T_>::setGraph(const Vect<RC>& I,
                            int             opt)
{
   _length = I.size();
   _nb_rows = _nb_cols = 0;
   _IJ.resize(_length);
   for (size_t i=0; i<_length; i++) {
      _nb_rows = std::max(_nb_rows,I[i].first);
      _nb_cols = std::max(_nb_cols,I[i].second);
      _IJ[i] = RC(I[i].first-1,I[i].second-1);
   }
   _size = _nb_rows;
   if (opt==0) {
      sort(_IJ.begin(),_IJ.end());
      vector<RC>::iterator new_end = unique(_IJ.begin(),_IJ.end());
      _IJ.erase(new_end,_IJ.end());
   }
   _length = _IJ.size();
   _row_ptr.resize(_size+1);
   _col_ind.resize(_length);
   StoreGraph(_IJ,_row_ptr,_col_ind);
#ifdef USE_EIGEN
   for (size_t i=0; i<_row_ptr.size()-1; i++)
      _nbc.push_back(_row_ptr[i+1]-_row_ptr[i]);
   _A.reserve(_nbc);
   clear();
#else
   _a.resize(_length,0);
#endif
}


template<class T_>
Vect<T_> SpMatrix<T_>::getRow(size_t i) const
{
   Vect<T_> v(_nb_cols);
#ifdef USE_EIGEN
   size_t j=0;
   for (size_t k=0; k<_length; k++) {
      if (_IJ[k].first==i-1)
         v[j++] = _A.coeff(_IJ[k].first,_IJ[k].second);
   }
#else
   for (size_t j=1; j<=_nb_cols; j++)
      v(j) = get(i,j);
#endif
   return v;
}


template<class T_>
Vect<T_> SpMatrix<T_>::getColumn(size_t j) const
{
   Vect<T_> v(_nb_rows);
#ifdef USE_EIGEN
   size_t i=0;
   for (size_t k=0; k<_length; k++) {
      if (_IJ[k].second==j-1)
         v[i++] = _A.coeff(_IJ[k].first,_IJ[k].second);
   }
#else
   for (size_t i=1; i<=_nb_rows; i++)
      v(i) = get(i,j);
#endif
   return v;
}


template<class T_>
T_& SpMatrix<T_>::operator()(size_t i,
                             size_t j)
{
#ifdef USE_EIGEN
   return _A.coeffRef(i-1,j-1);
#else
   int k=_col_index(i,j);
   if (k<0)
      throw OFELIException("In SpMatrix::set(i,j,x): Index pair: (" + to_string(i) +
                           "," + to_string(j) + ") not compatible with sparse storage.");
   else
      return _a[_row_ptr[i-1]+k];
   return _temp;
#endif
}


template<class T_>
T_ SpMatrix<T_>::operator()(size_t i,
                            size_t j) const
{
#ifdef USE_EIGEN
   return _A.coeff(i-1,j-1);
#else
   int k=_col_index(i,j);
   if (k<0)
      return _zero;
   else
      return _a[_row_ptr[i-1]+k];
#endif
}


template<class T_>
T_ SpMatrix<T_>::operator()(size_t i) const { return _a[i-1]; }


template<class T_>
T_ SpMatrix<T_>::operator[](size_t i) const { return _a[i]; }


template<class T_>
Vect<T_> SpMatrix<T_>::operator*(const Vect<T_>& x) const
{
   Vect<T_> y(_nb_rows);
   y.clear();
#ifdef USE_EIGEN
   for (size_t i=0; i<_length; i++)
      y[_IJ[i].first] += _A.coeff(_IJ[i].first,_IJ[i].second)*x[_IJ[i].second];
#else
   Mult(x,y);
#endif
   return y;
}


template<class T_>
SpMatrix<T_>& SpMatrix<T_>::operator*=(const T_& a)
{
#ifdef USE_EIGEN
   for (size_t i=0; i<_length; i++)
      _A.coeffRef(_IJ[i].first,_IJ[i].second) *= a;
#else
   for (size_t k=0; k<_length; ++k)
      _a[k] *= a;
#endif
   return *this;
}


template<class T_>
void SpMatrix<T_>::getMesh(Mesh& mesh)
{
#if !defined(USE_EIGEN)
   if (_sides)
       SideGraph(mesh,_row_ptr,_col_ind);
   else {
       if (_dof)
          _size = _nb_rows = _nb_cols = mesh.getNbNodes();
       else
          _size = _nb_rows = _nb_cols = mesh.getNbEq();
       if (_type)
          XGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
//             XGraphScal(mesh,_row_ptr,_col_ind);
       else
          NodeGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
   }
   _a.resize(_length,T_(0.));
#endif
   _diag.resize(_size);
}


template<class T_>
void SpMatrix<T_>::Mult(const Vect<T_>& x,
                        Vect<T_>&       y) const
{
#ifdef USE_EIGEN
   for (size_t i=0; i<_length; i++)
      y[_IJ[i].first] = _A.coeff(_IJ[i].first,_IJ[i].second)*x[_IJ[i].second];
#else
   y = static_cast<T_>(0);
   MultAdd(x,y);
#endif
}


template<class T_>
void SpMatrix<T_>::MultAdd(const Vect<T_>& x,
                           Vect<T_>&       y) const
{
#ifdef USE_EIGEN
   for (size_t i=0; i<_length; i++)
      y[_IJ[i].first] += _A.coeff(_IJ[i].first,_IJ[i].second)*x[_IJ[i].second];
#else
   size_t l=0;
   for (size_t i=0; i<_nb_rows; ++i)
      for (size_t j=0; j<_row_ptr[i+1]-_row_ptr[i]; ++j)
         y[i] += _a[_row_ptr[i]+j] * x[_col_ind[l++]-1];
#endif
}


template<class T_>
void SpMatrix<T_>::MultAdd(T_              a,
                           const Vect<T_>& x,
                           Vect<T_>&       y) const
{
#ifdef USE_EIGEN
   for (size_t i=0; i<_length; i++)
      y[_IJ[i].first] += a*_A.coeff(_IJ[i].first,_IJ[i].second)*x[_IJ[i].second];
#else
   size_t l=0;
   for (size_t i=0; i<_nb_rows; i++)
      for (size_t j=0; j<_row_ptr[i+1]-_row_ptr[i]; j++)
         y[i] += a * _a[_row_ptr[i]+j] * x[_col_ind[l++]-1];
#endif
}


template<class T_>
void SpMatrix<T_>::TMult(const Vect<T_>& x,
                         Vect<T_>&       y) const
{
#ifdef USE_EIGEN
   for (size_t i=0; i<_length; i++)
      y[_IJ[i].first] += _A.coeff(_IJ[i].second,_IJ[i].first)*x[_IJ[i].second];
#else
   for (size_t i=0; i<_nb_rows; i++)
      for (size_t j=0; j<_row_ptr[i+1]-_row_ptr[i]; j++)
         y[_col_ind[j+_row_ptr[i]-1]] += _a[_row_ptr[i]+j] * x[i];
#endif
}


template<class T_>
void SpMatrix<T_>:: Axpy(T_                  a,
                         const SpMatrix<T_>& m)
{
#ifdef USE_EIGEN
   _A += a * m._A;
#else
   _a += a * m._a;
#endif
}


template<class T_>
void SpMatrix<T_>::Axpy(T_                a,
                        const Matrix<T_>* m)
{ }


template<class T_>
void SpMatrix<T_>::set(size_t    i,
                       size_t    j,
                       const T_& val)
{
#ifdef USE_EIGEN
   _A.coeffRef(i-1,j-1) = val;
#else
   int k=_col_index(i,j);
   if (k<0)
      throw OFELIException("In SpMatrix::set(i,j,x): Index pair (" + to_string(i) +
                           "," + to_string(j) + ") not compatible with sparse storage.");
   else
      _a[_row_ptr[i-1]+k] = val;
#endif
}


template<class T_>
void SpMatrix<T_>::add(size_t    i,
                       size_t    j,
                       const T_& val)
{
#ifdef USE_EIGEN
   _A.coeffRef(i-1,j-1) += val;
#else
   _a[_row_ptr[i-1]+_col_index(i,j)] += val;
#endif
}


template<class T_>
void SpMatrix<T_>::operator=(const T_& x)
{
#ifdef USE_EIGEN
   for (size_t i=0; i<_length; i++) {
      _A.coeffRef(_IJ[i].first,_IJ[i].second) = static_cast<T_>(0);
      if (_IJ[i].first==_IJ[i].second)
         _A.coeffRef(_IJ[i].first,_IJ[i].first) = x;	
}
#else
   for (size_t i=0; i<_length; ++i)
      _a[i] = x;
#endif
}


template<class T_>
size_t SpMatrix<T_>::getColInd(size_t i) const
{
   return _col_ind[i-1];
}


template<class T_>
size_t SpMatrix<T_>::getRowPtr(size_t i) const
{
   return _row_ptr[i-1];
}


template<class T_>
int SpMatrix<T_>::Factor()
{
   return -1;
}


#ifdef USE_EIGEN
template<class T_>
int SpMatrix<T_>::DILUFactorize(vector<size_t>& id,
                                    vector<T_>&     pivot) const;
#else
template<class T_>
int SpMatrix<T_>::DILUFactorize(vector<size_t>& id,
                                vector<T_>&     pivot) const
{
   id.resize(_size);
   pivot.resize(_size);
   size_t k=0;
   for (size_t i=0; i<_size; i++) {
      for (size_t j=_row_ptr[i]; j<_row_ptr[i+1]; j++, k++) {
         if (_col_ind[j]==i+1) {
            id[i] = k + 1;
            if (_a[k]==static_cast<T_>(0))
               throw OFELIException("In SpMatrix::DILUFactorize(...): Zero pivot detected in row "
                                    +to_string(i+1));
            pivot[i] = _a[k];
         }
      }
   }
   int found=0;
   T_ c=static_cast<T_>(0);
   for (size_t i=0; i<_size; ++i) {
      pivot[i] = static_cast<T_>(1)/pivot[i];
      for (size_t j=id[i]; j<_row_ptr[i+1]; ++j) {
         found = 0;
         size_t l=_col_ind[j]-1;
         for (k=_row_ptr[l]; k<id[l]-1; ++k) {
            if (_col_ind[k]==i+1)
               found = 1, c = _a[k];
         }
         if (found)
            pivot[l] -= c*pivot[i]*_a[j];
      }
   }
   return 0;
}
#endif


#ifdef USE_EIGEN
template<class T_>
int SpMatrix<T_>::ILUFactorize(const SpMatrix<T_>& A);
#else
template<class T_>
int SpMatrix<T_>::ILUFactorize(const SpMatrix<T_>& A)
{
   _nb_rows = A._nb_rows, _nb_cols = A._nb_cols;
   _size = A._size, _length = A._length;
   _row_ptr = A._row_ptr; _col_ind = A._col_ind;
   _lnnz = (_length-_nb_rows)/2, _unnz = (_length+_nb_rows)/2;
   _aL.resize(_lnnz); _aU.resize(_unnz);
   _l_col_ind.resize(_lnnz);
   _u_col_ind.resize(_unnz);
   _l_row_ptr.resize(_nb_rows+1);
   _u_row_ptr.resize(_nb_rows+1);
   _l_row_ptr[0] = _u_row_ptr[0] = 0;

   for (size_t i=0; i<_nb_rows; i++) {
      _l_row_ptr[i+1] = _l_row_ptr[i];
      _u_row_ptr[i+1] = _u_row_ptr[i];
      for (size_t j=A._row_ptr[i]; j<A._row_ptr[i+1]; j++) {
         if (A._col_ind[j]<i+1) {
            size_t k=_l_row_ptr[i+1]++;
            _aL[k] = A._a[j];
            _l_col_ind[k] = A._col_ind[j];
         }
         else {
            size_t k=_u_row_ptr[i+1]++;
            _aU[k] = A._a[j];
            _u_col_ind[k] = A._col_ind[j];
         }
      }
   }
   for (size_t i=1; i<_nb_rows; i++) {
      for (size_t j=_l_row_ptr[i]; j<_l_row_ptr[i+1]; j++) {
         size_t pn=_u_row_ptr[_col_ind[j]-1], qn=j+1, rn=_u_row_ptr[i];
         T_ p = (_aL[j]/=_aU[pn]);
         for (pn++; pn<_u_row_ptr[_l_col_ind[j]] && _u_col_ind[pn]<i+1; pn++) {
            while (qn<_l_row_ptr[i+1] && _l_col_ind[qn]<_u_col_ind[pn])
               qn++;
            if (qn<_l_row_ptr[i+1] && _u_col_ind[pn]==_l_col_ind[qn])
               _aL[qn] -= p*_aU[pn];
         }
         for (; pn<_u_row_ptr[_l_col_ind[j]]; pn++) {
            while (rn<_u_row_ptr[i+1] && _u_col_ind[rn]<_u_col_ind[pn])
               rn++;
            if (rn<_u_row_ptr[i+1] && _u_col_ind[pn]==_u_col_ind[rn])
               _aU[rn] -= p*_aU[pn];
         }
      }
   }

   size_t k=0, kl=0, ku=0;
   _a.resize(_length);
   for (size_t i=0; i<_nb_rows; i++) {
      for (size_t j=_row_ptr[i]; j<_row_ptr[i+1]; j++) {
         if (_col_ind[j]<i+1)
            _a[k++] = _aL[kl++];
         else
            _a[k++] = _aU[ku++];
      }
   }
   return 0;
}
#endif


#ifndef USE_EIGEN
template<class T_>
void SpMatrix<T_>::DILUSolve(const vector<size_t>& id,
                             const vector<T_>&     pivot,
                             const Vect<T_>&       b,
                             Vect<T_>&             x) const
{
   vector<T_> z(_size);
   for (size_t i=0; i<_size; i++) {
      T_ s = 0;
      for (size_t j=_row_ptr[i]; j<id[i]; ++j)
         s += _a[j] * z[_col_ind[j]-1];
      z[i] = pivot[i] * (b[i]-s);
   }
   for (size_t i=0; i<_size; ++i) {
      T_ s = 0;
      for (size_t j=id[_size-i-1]; j<_row_ptr[_size-i]; ++j)
         s += _a[j] * x(_col_ind[j]);
      x[_size-i-1] = z[_size-i-1] - pivot[_size-i-1] * s;
   }
}
#endif


#ifndef USE_EIGEN
template<class T_>
void SpMatrix<T_>::ILUSolve(const Vect<T_>& b,
                            Vect<T_>&       x) const
{
   for (size_t i=0; i<_size; i++) {
      T_ s=0;
      for (size_t j=_l_row_ptr[i]; j<_l_row_ptr[i+1]; j++)
         s += _aL[_l_col_ind[j]-1]*b[_l_col_ind[j]-1];
      x[i] = b[i] - s;
   }
   for (int i=int(_size)-1; i>=0; i--) {
      T_ s=0;
      for (size_t j=_l_row_ptr[i]; j<_l_row_ptr[i+1]; j++)
         s += _aU[_u_col_ind[j]-1]*x[_u_col_ind[j]-1];
      x[i] = (x[i]-s)/_aU[_u_col_ind[_l_row_ptr[i]]-1];
   }
}
#endif


#ifndef USE_EIGEN
template<class T_>
void SpMatrix<T_>::SSORSolve(const Vect<T_>& b,
                             Vect<T_>&       x) const
{
   size_t k=0;
   vector<size_t> id(_size);
   for (size_t i=0; i<_size; i++) {
      for (size_t j=_row_ptr[i]; j<_row_ptr[i+1]; j++, k++) {
         if (_col_ind[j]==i+1) {
            id[i] = k + 1;
            if (_a[k]==static_cast<T_>(0))
               throw OFELIException("In SpMatrix::SSORSolve(b,x):"
                                    " Zero pivot detected in row " + to_string(i+1));
         }
      }
   }
   Vect<T_> z(_size);
   for (size_t i=0; i<_size; i++) {
      T_ s = 0;
      for (size_t j=_row_ptr[i]; j<id[i]-1; j++)
         s += _a[j] * z(_col_ind[j]);
      z[i] = (b[i]-s)/_a[id[i]-1];
   }
   for (int i=int(_size)-1; i>=0; i--) {
      T_ s = 0;
      for (size_t j=id[i]; j<_row_ptr[i+1]; j++)
         s += _a[j] * x(_col_ind[j]);
      x[i] = z[i] - s/_a[id[i]-1];
   }
}
#endif


template<class T_>
int SpMatrix<T_>::solve(Vect<T_>& b,
                        bool      fact)
{
   Vect<T_> x(b.size());
   int ret = solve(b,x,fact);
   b = x;
   return ret;
}


template<class T_>
int SpMatrix<T_>::solve(const Vect<T_>& b,
                        Vect<T_>&       x,
                        bool            fact)
{
   if (_solver==DIRECT_SOLVER)
      throw OFELIException("In SpMatrix::solve(...): No solver provided.");
   LinearSolver<T_> ls(*this,b,x);
   return ls.solve(_solver,_prec);
}


template<class T_>
void SpMatrix<T_>::setSolver(Iteration      solver,
                             Preconditioner prec,
                             int            max_it,
                             real_t         toler)
{
   _solver = solver;
   _prec = prec;
   _max_it = max_it;
   _toler = toler;
}


template<class T_>
void SpMatrix<T_>::clear()
{
#ifdef USE_EIGEN
   for (size_t i=0; i<_length; i++)
      _A.coeffRef(_IJ[i].first,_IJ[i].second) = static_cast<T_>(0);
#else
   for (size_t i=0; i<_length; ++i)
      _a[i] = static_cast<T_>(0);
#endif
}


template<class T_>
T_ *SpMatrix<T_>::get() const
{
   return &_a[0];
}


template<class T_>
T_ SpMatrix<T_>::get(size_t i,
                     size_t j) const
{
#ifdef USE_EIGEN
   return _A.coeff(i-1,j-1);
#else
   int k=_col_index(i,j);
   if (k<0)
      return _zero;
   else
      return _a[_row_ptr[i-1]+k-1];
#endif
}

#ifdef USE_EIGEN
template<class T_>
SpMat& SpMatrix<T_>::getEigenMatrix()
{
   return _A;
}
#endif

#if !defined(USE_EIGEN)
template<class T_>
int SpMatrix<T_>::_col_index(size_t i, size_t j) const
{
   for (int k=0; k<int(_row_ptr[i]-_row_ptr[i-1]); ++k)
      if (_col_ind[_row_ptr[i-1]+k]==j)
         return k;
   return -1;
}
#endif


template<>
inline void SpMatrix<real_t>::Laplace1D(size_t n,
                                        real_t h)
{
   _is_dense = 0;
   _is_diagonal = 0;
   _size = _nb_rows = _nb_cols = n;
#ifdef USE_EIGEN
   _length = 3*_size - 2;
   _A.reserve(_length);
   _A.insert(0,0) =  2./h;
   _A.insert(0,1) = -1./h;
   for (size_t i=1; i<_nb_rows-1; i++) {
      _A.insert(i,i-1) = -1./h;
      _A.insert(i,i  ) =  2./h;
      _A.insert(i,i+1) = -1./h;
   }
   _A.insert(_nb_rows-1,_nb_rows-2) = -1./h;
   _A.insert(_nb_rows-1,_nb_rows-1) =  2./h;
#else
   _row_ptr.push_back(0);
   _row_ptr.push_back(2);
   _col_ind.push_back(1);
   _col_ind.push_back(2);
   _a.push_back( 2./h);
   _a.push_back(-1./h);
   for (size_t i=1; i<_size-1; i++) {
      _row_ptr.push_back(_row_ptr[i]+3);
      _col_ind.push_back(i);
      _col_ind.push_back(i+1);
      _col_ind.push_back(i+2);
      _a.push_back(-1./h);
      _a.push_back( 2./h);
      _a.push_back(-1./h);
   }
   _col_ind.push_back(_size-1);
   _col_ind.push_back(_size);
   _row_ptr.push_back(_row_ptr[_size-1]+2);
   _a.push_back(-1./h);
   _a.push_back( 2./h);
   _length = _row_ptr[_size];
#endif
}


template<class T_>
Vect<T_> operator*(const SpMatrix<T_>& A,
                   const Vect<T_>&     b)
{
   Vect<T_> v(b.size());
   A.Mult(b,v);
}


template<class T_>
ostream& operator<<(ostream&            s,
                    const SpMatrix<T_>& A)
{
#ifdef USE_EIGEN
   s << endl << A._A << endl;
#else
   s.setf(ios::right|ios::scientific);
   s << endl;
   size_t k = 0;
   for (size_t i=0; i<A.getNbRows(); ++i) {
      for (size_t j=A._row_ptr[i]; j<A._row_ptr[i+1]; ++j)
         s << "(" << setw(6) << i+1 << "," << setw(6) << A._col_ind[j] << "): "
           << setprecision(8) << std::setfill(' ') << setw(18) << A._a[k++] << endl;
   }
#endif
   return s;
}

} /* namespace OFELI */

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif

