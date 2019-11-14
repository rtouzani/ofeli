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

                     Implementation of class Matrix

  ==============================================================================*/


#ifndef __MATRIX_IMPL_H
#define __MATRIX_IMPL_H

#include <iostream>
#include <algorithm>

#include "mesh/Mesh.h"
#include "linear_algebra/GraphOfMatrix.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/LocalMatrix_impl.h"
#include "util/util.h"
#include "OFELIException.h"

namespace OFELI {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template<class T_>
Matrix<T_>::Matrix()
           : _nb_rows(0), _nb_cols(0), _size(0), _length(0), _zero(T_(0)),
             _penal(1.e20), _is_diagonal(false)
{ }


template<class T_>
Matrix<T_>::Matrix(const Matrix<T_> &m)
           : _nb_rows(m._nb_rows), _nb_cols(m._nb_cols), _size(m._size), _length(m._length),
             _zero(T_(0)), _penal(m._penal), _is_diagonal(m._is_diagonal)
{
   _ch.resize(_size);
   _diag.setSize(_size);
   _ch = m._ch;
   _diag = m._diag;
   _theMesh = m._theMesh;
}


template<class T_>
Matrix<T_>::~Matrix() { }


template<class T_>
void Matrix<T_>::reset() { }


template<class T_>
size_t Matrix<T_>::getNbRows() const { return _nb_rows; }


template<class T_>
size_t Matrix<T_>::getNbColumns() const { return _nb_cols; }


template<class T_>
void Matrix<T_>::setPenal(real_t p) { _penal = p; }


template<class T_>
void Matrix<T_>::setDiagonal()
{
   _size = _theMesh->getNbEq();
   _ch.resize(_size);
   _ch[0] = 0;
   for (size_t i=1; i<_size; i++)
      _ch[i] = i+1;
   _a.resize(_size);
   clear(_a);
   _dof = 0;
   _length = _nb_rows = _nb_cols = _size;
   _is_diagonal = true;
}


template<class T_>
T_ Matrix<T_>::getDiag(size_t k) const { return _diag[k-1]; }


template<class T_>
size_t Matrix<T_>::size() const { return _size; }


template<class T_>
void Matrix<T_>::setDiagonal(Mesh& mesh)
{
   init_set_mesh(mesh);
   _size = _theMesh->getNbEq();
   _ch.resize(_size);
   _ch[0] = 0;
   for (size_t i=1; i<_size; i++)
      _ch[i] = i+1;
   _a.resize(_size);
   clear(_a);
   _dof = 0;
   _length = _nb_rows = _nb_cols = _size;
   _is_diagonal = true;
}

 
template<class T_>
void Matrix<T_>::init_set_mesh(Mesh&  mesh,
                               size_t dof)
{
   _theMesh = &mesh;
   _zero = T_(0);
   _dof_type = 0;
   if (_theMesh->NodesAreDOF())
      _dof_type = NODE_DOF;
   else if (_theMesh->SidesAreDOF())
      _dof_type = SIDE_DOF;
   else if (_theMesh->ElementsAreDOF())
      _dof_type = ELEMENT_DOF;
   _dof = dof;
   _size = _nb_rows = _nb_cols = _theMesh->getNbEq();
   if (_dof_type==NODE_DOF)
      if (_dof)
         _size = _nb_rows = _nb_cols = _theMesh->getNbNodes();
      else
         _size = _nb_rows = _nb_cols = _theMesh->getNbEq();
   else if (_dof_type==SIDE_DOF)
      if (_dof)
         _size = _nb_rows = _nb_cols = _theMesh->getNbSides();
      else
         _size = _nb_rows = _nb_cols = _theMesh->getNbEq();
   else if (_dof_type==ELEMENT_DOF)
      _size = _nb_rows = _nb_cols = _theMesh->getNbElements();
   else;
}


template<class T_>
void Matrix<T_>::clear()
{
#ifndef USE_EIGEN
   for (size_t i=0; i<_length; i++)
      _a[i] = static_cast<T_>(0);
#endif
}


template<class T_>
void Matrix<T_>::Assembly(const Element& el,
                          T_*            a)
{
   size_t kk=0;
   if (_is_diagonal) {
      for (size_t i=1; i<=el.getNbNodes(); ++i) {
         Node *nd = el(i);
         size_t nb_dof = nd->getNbDOF();
         for (size_t k=1; k<=nb_dof; ++k) {
            size_t n=nb_dof*(nd->n()-1) + k;
            add(n,n,a[kk]);
            kk += nb_dof*el.getNbNodes() + 1;
         }
      }
      return;
   }
   for (size_t i=1; i<=el.getNbNodes(); ++i) {
      Node *nd1 = el(i);
      for (size_t k=1; k<=nd1->getNbDOF(); ++k) {
         for (size_t j=1; j<=el.getNbNodes(); ++j) {
            Node *nd2 = el(j);
            for (size_t l=1; l<=nd2->getNbDOF(); ++l, kk++) {
               if (nd1->getDOF(k) && nd2->getDOF(l))
                  add(nd1->getDOF(k),nd2->getDOF(l),a[kk]);
            }
         }
      }
   }
}

/*
template<class T_>
void Matrix<T_>::Assembly(const Element&     el,
                          const DMatrix<T_>& A)
{
   size_t i=1;
   for (size_t in=1; in<=el.getNbNodes(); ++in) {
      Node *nd1=el(in);
      for (size_t k=1; k<=nd1->getNbDOF(); ++k, i++) {
         size_t j=1;
         for (size_t jn=1; jn<=el.getNbNodes(); ++jn) {
            Node *nd2=el(jn);
            for (size_t l=1; l<=nd2->getNbDOF(); ++l, j++) {
               if (nd1->getDOF(k)!=0 && nd2->getDOF(l)!=0)
                  add(nd1->getDOF(k),nd2->getDOF(l),A(i,j));
            }
         }
      }
   }
   }*/


template<class T_>
void Matrix<T_>::DGAssembly(const Element& el,
                            T_*            a)
{
   size_t kk=0;
   if (_is_diagonal) {
      for (size_t i=1; i<=el.getNbDOF(); ++i) {
         Node *nd = el(i);
         size_t nb_dof = nd->getNbDOF();
         for (size_t k=1; k<=nb_dof; ++k) {
            size_t n=nb_dof*(nd->n()-1) + k;
            add(n,n,a[kk]);
            kk += nb_dof*el.getNbNodes() + 1;
         }
      }
      return;
   }
   for (size_t i=1; i<=el.getNbNodes(); ++i) {
      Node *nd1 = el(i);
      for (size_t k=1; k<=nd1->getNbDOF(); ++k) {
         for (size_t j=1; j<=el.getNbNodes(); ++j) {
            Node *nd2=el(j);
            for (size_t l=1; l<=nd2->getNbDOF(); ++l, kk++) {
               if (nd1->getDOF(k)!=0 && nd2->getDOF(l)!=0)
                  add(nd1->getDOF(k),nd2->getDOF(l),a[kk]);
            }
         }
      }
   }
}


template<class T_>
void Matrix<T_>::DGAssembly(const Element&                                               el,
                            const LocalMatrix<T_,MAX_NB_ELEMENT_DOF,MAX_NB_ELEMENT_DOF>& a)
{
   for (size_t i=1; i<=el.getNbDOF(); ++i) {
      for (size_t j=1; j<=el.getNbDOF(); ++j) {
         if (el.getDOF(i)!=0 && el.getDOF(j)!=0)
            add(el.getDOF(i),el.getDOF(j),a(i,j));
      }
   }
}


template<class T_>
void Matrix<T_>::DGAssembly(const Side&                                            sd,
                            const LocalMatrix<T_,MAX_NB_SIDE_DOF,MAX_NB_SIDE_DOF>& a)
{
   for (size_t i=1; i<=sd.getNbDOF(); ++i) {
      for (size_t j=1; j<=sd.getNbDOF(); ++j) {
         if (sd.getDOF(i)!=0 && sd.getDOF(j)!=0)
            add(sd.getDOF(i),sd.getDOF(j),a(i,j));
      }
   }
}


template<class T_>
void Matrix<T_>::Assembly(const Side& sd, 
                          T_*         a)
{
   size_t kk = 0;
   for (size_t in=1; in<=sd.getNbNodes(); ++in) {
      Node *nd1 = sd(in);
      for (size_t k=1; k<=nd1->getNbDOF(); ++k) {
         for (size_t jn=1; jn<=sd.getNbNodes(); ++jn) {
            Node *nd2=sd(jn);
            for (size_t l=1; l<=nd2->getNbDOF(); ++l, kk++) {
               if (nd1->getDOF(k)!=0 && nd2->getDOF(l)!=0)
                  add(nd1->getDOF(k),nd2->getDOF(l),a[kk]);
            }
         }
      }
   }
}

/*
template<class T_>
void Matrix<T_>::Assembly(const Side&        sd,
                          const DMatrix<T_>& a)
{
   size_t i=1;
   for (size_t in=1; in<=sd.getNbNodes(); ++in) {
      Node *nd1 = sd(in);
      for (size_t k=1; k<=nd1->getNbDOF(); ++k, i++) {
         size_t j=1;
         for (size_t jn=1; jn<=sd.getNbNodes(); ++jn) {
            Node *nd2=sd(jn);
            for (size_t l=1; l<=nd2->getNbDOF(); ++l, j++) {
               if (nd1->getDOF(k)!=0 && nd2->getDOF(l)!=0)
                  add(nd1->getDOF(k),nd2->getDOF(l),a(i,j));
            }
         }
      }
   }
   }*/


template<class T_>
void Matrix<T_>::Prescribe(Mesh&           mesh,
                           Vect<T_>&       b,
                           const Vect<T_>& u,
                           int             flag)
{
   MeshNodes(mesh) {
      for (size_t i=1; i<=theNode->getNbDOF(); ++i) {
         if (TheNode.getCode(i)>0) {
            size_t k=TheNode.getDOF(i);
            if (flag==0) {
               _diag[k-1] = get(k,k)*_penal;
               set(k,k,_diag[k-1]);
            }
            b.set(k,u(k)*_diag[k-1]);
         }
      }
   }
}


template<class T_>
void Matrix<T_>::Prescribe(Vect<T_>&       b,
                           const Vect<T_>& u,
                           int             flag)
{
   Prescribe(*_theMesh,b,u,flag);
}


template<class T_>
void Matrix<T_>::Prescribe(int             dof,
                           int             code,
                           Mesh&           mesh,
                           Vect<T_>&       b,
                           const Vect<T_>& u,
                           int             flag)
{
   MeshNodes(mesh) {
      if (theNode->getCode(dof)==code) {
         size_t k=theNode->getDOF(dof);
         if (flag==0) {
            _diag[k-1] = get(k,k)*_penal;
            set(k,k,_diag[k-1]);
         }
         b.set(k,u(k)*_diag[k-1]);
      }
   }
}


template<class T_>
void Matrix<T_>::Prescribe(int             dof,
                           int             code,
                           Vect<T_>&       b,
                           const Vect<T_>& u,
                           int             flag)
{
   Prescribe(dof,code,*_theMesh,b,u,flag);
}


template<class T_>
void Matrix<T_>::Prescribe(Mesh&     mesh,
                           Vect<T_>& b,
                           int       flag)
{
   MeshNodes(mesh) {
      for (size_t j=1; j<=theNode->getNbDOF(); ++j)
         if (theNode->getCode(j)>0) {
            size_t k=theNode->getDOF(j);
            if (!flag) {
               _diag[k-1] = get(k,k)*_penal;
               set(k,k,_diag[k-1]);
            }
            b.set(k,0);
         }
   }
}


template<class T_>
void Matrix<T_>::Prescribe(Vect<T_>& b,
                           int       flag)
{
   Prescribe(*_theMesh,b,flag);
}


template<class T_>
void Matrix<T_>::Prescribe(size_t          dof,
                           Mesh&           mesh,
                           Vect<T_>&       b,
                           const Vect<T_>& u,
                           int             flag)
{
   mesh_nodes(mesh) {
      if (The_node.getCode(dof)>0) {
         size_t k=node_label;
         if (!flag) {
            _diag[k-1] = get(k,k)*_penal;
            set(k,k,_diag[k-1]);
         }
         b.set(k,u(k)*_diag[k-1]);
      }
   }
}


template<class T_>
void Matrix<T_>::Prescribe1(Mesh&           mesh,
                            Vect<T_>&       b,
                            const Vect<T_>& u,
                            int             flag)
{
   mesh_nodes(mesh) {
      for (size_t i=1; i<=The_node.getNbDOF(); i++) {
         if (The_node.getCode(i)>0) {
            size_t k=The_node.getDOF(i);
            if (!flag)
               add(k,k,_penal);
            b.set(k,u(k)*_penal);
         }
      }
   }
}


template<class T_>
void Matrix<T_>::Prescribe1(Vect<T_>&       b,
                            const Vect<T_>& u,
                            int             flag)
{
   Prescribe1(*_theMesh,b,u,flag);
}


template<class T_>
void Matrix<T_>::Prescribe(Mesh& mesh)
{
   mesh_nodes(mesh) {
      for (size_t i=1; i<=The_node.getNbDOF(); i++) {
         if (The_node.getCode(i)>0) {
            size_t k=The_node.getDOF(i);
            set(k,k,get(k,k)*_penal);
         }
      }
   }
}


template<class T_>
void Matrix<T_>::Prescribe()
{
   Prescribe(*_theMesh);
}


template<class T_>
void Matrix<T_>::PrescribeSide(Mesh& mesh)
{
   mesh_sides(mesh) {
      for (size_t i=1; i<=The_side.getNbDOF(); i++) {
         if (The_side.getCode(i)>0) {
            size_t k = The_side.getDOF(i);
            set(k,k,get(k,k)*_penal);
         }
      }
   }
}


template<class T_>
void Matrix<T_>::PrescribeSide()
{
   PrescribeSide(*_theMesh);
}


template<class T_>
void Matrix<T_>::Constraint(const Mesh& mesh)
{
   Prescribe(mesh);
}


template<class T_>
void Matrix<T_>::Constraint()
{
   Prescribe();
}


template<class T_>
int Matrix<T_>::FactorAndSolve(Vect<T_>& b)
{
   Factor();
   return solve(b);
}


template<class T_>
int Matrix<T_>::FactorAndSolve(const Vect<T_>& b,
                               Vect<T_>&       x)
{
   int ret = Factor();
   solve(b,x,false);
   return ret;
}


template<class T_>
size_t Matrix<T_>::getLength() const
{
   return _length;
}


template<class T_>
int Matrix<T_>::isDiagonal() const
{
   return _is_diagonal;
}


template<class T_>
size_t Matrix<T_>::getColInd(size_t i) const
{
   i = 0;
   return 0;
}


template<class T_>
size_t Matrix<T_>::getRowPtr(size_t i) const
{
   i = 0;
   return 0;
}


template<class T_>
T_ Matrix<T_>::operator()(size_t i) const
{
   return _a[i-1];
}


template<class T_>
T_& Matrix<T_>::operator()(size_t i)
{
   return _a[i-1];
}


template<class T_>
T_& Matrix<T_>::operator[](size_t k)
{
   return _a[k];
}


template<class T_>
T_ Matrix<T_>::operator[](size_t k) const
{
   return _a[k];
}

template<class T_>
Matrix<T_>& Matrix<T_>::operator=(Matrix<T_>& m)
{
   _a = m._a;
   return *this;
}


template<class T_>
Matrix<T_>& Matrix<T_>::operator+=(const Matrix<T_>& m)
{
   for (size_t i=1; i<=_length; ++i)
      _a.add(i,m._a[i-1]);
   return *this;
}


template<class T_>
Matrix<T_>& Matrix<T_>::operator-=(const Matrix<T_>& m)
{
   for (size_t i=0; i<_length; ++i)
      _a[i] -= m._a[i];
   return *this;
}


template<class T_>
Matrix<T_>& Matrix<T_>::operator=(const T_ &x)
{
   for (size_t i=0; i<_length; ++i)
      _a[i] = x;
   return *this;
}


template<class T_>
Matrix<T_>& Matrix<T_>::operator*=(const T_& x)
{
   for (size_t i=0; i<_length; ++i)
      _a[i] *= x;
   return *this;
}


template<class T_>
Matrix<T_>& Matrix<T_>::operator+=(const T_& x)
{
   for (size_t i=0; i<_length; ++i)
      _a[i] += x;
   return *this;
}


template<class T_>
Matrix<T_>& Matrix<T_>::operator-=(const T_& x)
{
   for (size_t i=0; i<_length; ++i)
      _a[i] = -x;
   return *this;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */

#endif
