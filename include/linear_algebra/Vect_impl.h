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

                           Implementation of class Vect 

  ==============================================================================*/

#ifndef __VECT_IMPL_H
#define __VECT_IMPL_H

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include "util/macros.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/Point.h"
#include "linear_algebra/LocalVect.h"
#include "linear_algebra/LocalMatrix.h"
#include "mesh/Mesh.h"
#include "mesh/MeshUtil.h"
#include "mesh/Grid.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Tetra4.h"
#include "OFELIException.h"

using std::to_string;

/** \defgroup VectMat Vector and Matrix
 *  \brief Vector and matrix classes
 */

/*! \file Vect.h
 *  \brief Definition file for class Vect.
 */

#if defined (USE_PETSC)
template<class T_> class PETScVect;
#endif

namespace OFELI {

template<class T_>
Vect<T_>::Vect() :
#if !defined (USE_EIGEN)
          vector<T_>(),
#endif
          _dof_type(NODE_DOF), _size(0), _nx(0), _ny(1), _nz(1), _nb_dof(1),
          _dg_degree(-1), _grid(true), _with_mesh(false), _theMesh(nullptr),
          _name("#"), _time(0)
{
   for (size_t i=0; i<10; ++i)
      _with_regex[i] = false;
}


template<class T_>
Vect<T_>:: Vect(size_t n) :
#if !defined (USE_EIGEN)
           vector<T_>(n),
#endif
           _dof_type(NODE_DOF), _size(n), _nx(n), _ny(1), _nz(1), _nb_dof(1), _dg_degree(-1),
           _grid(true), _with_mesh(false), _theMesh(nullptr), _name("#"), _time(0)
{
#if defined (USE_EIGEN)
   _v.conservativeResize(_size);
   clear();
#else
   for (size_t i=0; i<_size; i++)
      (*this)[i] = 0;
#endif
   for (size_t i=0; i<10; ++i)
      _with_regex[i] = false;
}


template<class T_>
Vect<T_>::Vect(size_t nx,
               size_t ny) :
#if !defined (USE_EIGEN)
          vector<T_>(nx*ny),
#endif
          _dof_type(NODE_DOF), _size(nx*ny), _nx(nx), _ny(ny), _nz(1), _nb_dof(1), _dg_degree(-1),
          _grid(false), _with_mesh(false), _theMesh(nullptr), _name("#"), _time(0)
{
#if defined (USE_EIGEN)
   _v.conservativeResize(_size);
   clear();
#else
   clear();
#endif
   for (size_t i=0; i<10; ++i)
      _with_regex[i] = false;
}


template<class T_>
Vect<T_>::Vect(size_t nx,
               size_t ny,
               size_t nz) :
#if !defined (USE_EIGEN)
          vector<T_>(nx*ny*nz),
#endif
          _dof_type(NODE_DOF), _size(nx*ny*nz), _nx(nx), _ny(ny), _nz(nz),
          _nb_dof(1), _dg_degree(-1),
          _grid(false), _with_mesh(false), _theMesh(nullptr), _name("#"), _time(0)
{
#if defined (USE_EIGEN)
   _v.conservativeResize(_size);
   clear();
#else
   for (size_t i=0; i<_size; i++)
      (*this)[i] = 0;
#endif
   for (size_t i=0; i<10; ++i)
      _with_regex[i] = false;
}


template<class T_>
Vect<T_>::Vect(size_t n,
               T_*    x) :
#if !defined (USE_EIGEN)
           vector<T_>(n),
#endif
           _dof_type(NODE_DOF), _size(n), _nx(n), _ny(1), _nz(1),
           _nb_dof(1), _dg_degree(-1),
           _grid(false), _with_mesh(false), _theMesh(nullptr), _name("#"), _time(0)
{
#if defined (USE_EIGEN)
   _v.conservativeResize(_size);
#endif
   for (size_t i=1; i<=n; ++i)
      set(i,x[i-1]);
   for (size_t i=0; i<10; ++i)
      _with_regex[i] = false;
}


template<class T_>
Vect<T_>::Vect(Grid& g) :
#if !defined (USE_EIGEN)
           vector<T_>((g.getNx()+1)*(g.getNy()+1)*(g.getNz()+1)),
#endif
           _dof_type(NODE_DOF), _nb_dof(1), _dg_degree(-1),
           _grid(true), _with_mesh(false), _theMesh(nullptr), _name("#"), _time(0)
{
   setGrid(g);
   for (size_t i=0; i<10; ++i)
      _with_regex[i] = false;
}


template<class T_>
Vect<T_>::Vect(Mesh&      m,
               DOFSupport dof_type,
               int        nb_dof)
         : _dg_degree(-1), _grid(false), _with_mesh(true), _name("#"), _time(0)
{
   setMesh(m,dof_type,nb_dof);
   for (size_t i=0; i<10; ++i)
      _with_regex[i] = false;
}


template<class T_>
Vect<T_>::Vect(Mesh&      m,
               DOFSupport dof_type,
               string     name,
               int        nb_dof,
               real_t     t) :
#if !defined (USE_EIGEN)
          vector<T_>(),
#endif
          _dg_degree(-1), _grid(false), _with_mesh(true), _name(name), _time(t)
{
   setMesh(m,dof_type,nb_dof);
   for (size_t i=0; i<10; ++i)
      _with_regex[i] = false;
}


template<class T_>
Vect<T_>::Vect(const Element*  el,
               const Vect<T_>& v) :
#if !defined (USE_EIGEN)
          vector<T_>(),
#endif
          _nx(el->getNbNodes()), _ny(v._ny), _nz(1),
          _nb(el->getNbNodes()), _nb_dof(v._nb_dof), _dg_degree(-1),
          _grid(false), _with_mesh(false), _name(v._name), _time(v._time)
{
   setSize(_nx,_ny);
   for (size_t n=1; n<=el->getNbNodes(); ++n) {
      Node *nd=(*el)(n);
      for (size_t j=1; j<=nd->getNbDOF(); ++j)
         set(n,j,v(nd->n(),j));
   }
   for (size_t i=0; i<10; ++i)
      _with_regex[i] = false;
}


template<class T_>
Vect<T_>::Vect(const Side*     sd,
               const Vect<T_>& v) :
#if !defined (USE_EIGEN)
           vector<T_>(),
#endif
           _nx(sd->getNbNodes()), _ny(v._nb_dof), _nz(1),
           _nb(sd->getNbNodes()), _nb_dof(v._nb_dof), _dg_degree(-1),
           _grid(false), _with_mesh(false), _name(v._name), _time(v._time)
{
   setSize(_nx,_ny);
   size_t i=0;
   for (size_t n=1; n<=sd->getNbNodes(); ++n) {
      Node *nd=(*sd)(n);
      size_t k=nd->getFirstDOF()-1;
      for (size_t j=1; j<=nd->getNbDOF(); ++j)
         set(++i,v[k++]);
   }
   for (size_t i=0; i<10; ++i)
      _with_regex[i] = false;
}


template<class T_>
Vect<T_>::Vect(const Vect<T_>& v,
               const Vect<T_>& bc) :
#if !defined (USE_EIGEN)
           vector<T_>(bc.size()),
#endif
           _dof_type(v._dof_type), _size(v._nb*v._nb), _nx(v._nb), _ny(v._nb_dof), _nz(1),
           _nb(v._nb), _nb_dof(v._nb_dof),
           _dg_degree(v._dg_degree), _grid(v._grid), _with_mesh(v._with_mesh),
           _theMesh(v._theMesh), _name(v._name), _time(v._time)
{
   size_t i=1, n=0;
#if defined (USE_EIGEN)
   resize(bc.size());
#endif
   MESH_ND {
      for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
         set(i,bc[n++]);
         if (The_node.getCode(k) == 0)
            set(i,v[The_node.getDOF(k)-1]);
         i++;
      }
   }
   for (size_t i=0; i<10; ++i)
      _with_regex[i] = false;
}


template<class T_>
Vect<T_>::Vect(const Vect<T_>& v,
               size_t          nb_dof,
               size_t          first_dof)
         : _dof_type(v._dof_type), _size(v._size), _nx(v._nx), _ny(v._ny), _nz(v._nz),
           _nb(v._nb), _nb_dof(v._nb_dof),
           _dg_degree(v._dg_degree), _grid(v._grid), _with_mesh(v._with_mesh),
           _theMesh(v._theMesh), _name(v._name), _time(v._time)
{
   setSize(_nb,_nb_dof,1);
   for (size_t i=1; i<=_nb; i++)
      for (size_t j=1; j<=_nb_dof; j++)
         set(i,j,v(i,j+first_dof-1));
   for (size_t i=0; i<10; ++i) {
      _with_regex[i] = v._with_regex[i];
   }
}


template<class T_>
Vect<T_>::Vect(const Vect<T_>& v)
         : _dof_type(v._dof_type), _nb(v._nb), _nb_dof(v._nb_dof),
           _dg_degree(v._dg_degree), _grid(v._grid), _with_mesh(v._with_mesh),
           _theMesh(v._theMesh), _name(v._name), _time(v._time)
{
   setSize(v._nx,v._ny,v._nz);
   for (size_t i=1; i<=_size; i++)
      set(i,v[i-1]);
   for (size_t i=0; i<10; ++i)
      _with_regex[i] = v._with_regex[i];
}


template<class T_>
Vect<T_>::Vect(const Vect<T_>& v,
               size_t          n)
         : _nb_dof(v._nb_dof), _grid(v._grid), _with_mesh(v._with_mesh),
           _theMesh(v._theMesh), _name(v._name), _time(v._time)
{
   if (n==1) {
      setSize(v._nx,1,1);
      for (size_t i=1; i<=_nx; i++)
         set(i,v(i,1,1));
   }
   else if (n==2) {
      setSize(1,v._ny,1);
      for (size_t j=1; j<=_ny; j++)
         set(j,v(1,j,1));
   }
   else if (n==3) {
      setSize(1,1,v._nz);
      for (size_t k=1; k<=_nz; k++)
         set(k,v(1,1,k));
   }
   for (size_t i=0; i<10; ++i)
      _with_regex[i] = false;
}


template<class T_>
Vect<T_>::Vect(size_t          d,
               const Vect<T_>& v,
               const string&   name)
         : _dof_type(v._dof_type), _nb(v._nb), _with_mesh(v._with_mesh),
           _theMesh(v._theMesh), _name(name), _time(v._time)
{
   if (d<=0)
      throw OFELIException("In Vect::Vect(size_t,Vect<T_>,string): Illegal value of nb_dof = "+to_string(d));
   size_t nd=v.getNbDOF();
   vector<size_t> dof_list(nd);
   dof_select(d,dof_list);
   if (_nb_dof>nd)
      throw OFELIException("In Vect::Vect(size_t,Vect<T_>,string): Illegal value of dof = "+to_string(nd));
   _theMesh = &(v.getMesh());
   setSize(_nb,_nb_dof,1);
   for (size_t i=1; i<=_nb; i++) {
      for (size_t k=0; k<nd; k++) {
         if (dof_list[k]!=0)
            set(i,dof_list[k],v(i,k+1));
      }
   }
   for (size_t i=0; i<10; ++i)
      _with_regex[i] = false;
}


template <class T_>
void Vect<T_>::dof_select(size_t          d,
                          vector<size_t>& dof_list)
{
   size_t k, kk, m, j=d, nd=dof_list.size();
   for (k=0; k<nd; k++) {
      kk = size_t(pow(10.,real_t(nd-k-1)));
      m = j/kk;
      dof_list[k] = m;
      j -= m*kk;
   }
   _nb_dof = 0;
   for (k=0; k<nd; k++)
      if (dof_list[k]!=0)
         _nb_dof++;
}


#if defined (USE_EIGEN)
template<class T_>
Vect<T_>::Vect(const VectorX& v)
         : _size(v.size()), _nx(_size), _ny(1), _nz(1), _dof_type(NONE), _nb_dof(1),
           _grid(true), _with_mesh(false), _theMesh(nullptr), _name("#"), _time(0)
{
   for (size_t i=0; i<10; ++i)
      _with_regex[i] = false;
}
#endif


template<class T_>
Vect<T_>::~Vect() { }


template<class T_>
void Vect<T_>::set(const T_* v,
                   size_t    n)
{
   setSize(n);
   for (size_t i=1; i<=n; ++i)
      set(i,v[i-1]);
}


template<class T_>
void Vect<T_>::select(const Vect<T_>& v,
                      size_t          nb_dof,
                      size_t          first_dof)
{
   _size = nb_dof*v._nb;
   setSize(_size);
   size_t i=1;
   for (size_t n=1; n<=v._nb; ++n)
      for (size_t j=first_dof; j<=nb_dof+first_dof-1; j++)
         set(i++,v(n,j));
}


template<class T_>
void Vect<T_>::set(const string& exp,
                   size_t        dof)
{
   if (_theMesh==nullptr)
      throw OFELIException("In Vect::set(string,dof): No mesh defined");
   _theFct.set(exp,_var);
   if (_dof_type==NODE_DOF) {
      MESH_ND
         set(The_node.getNbDOF()*(node_label-1)+dof,_theFct(The_node.getCoord(),_time));
   }
   else if (_dof_type==SIDE_DOF) {
      MESH_SD
         set(The_side.getNbDOF()*(side_label-1)+dof,_theFct(The_side.getCenter(),_time));
   }
   else if (_dof_type==BOUNDARY_SIDE_DOF) {
      MESH_BD_SD
         set(The_side.getNbDOF()*(side_label-1)+dof,_theFct(The_side.getCenter(),_time));
   }
   else if (_dof_type==ELEMENT_DOF) {
      MESH_EL
         set(element_label,_theFct(The_element.getCenter(),_time));
   }
   else
      throw OFELIException("In Vect::set(string,size_t): Unknown vector type.");
}


template<class T_>
void Vect<T_>::set(const string&       exp,
                   const Vect<real_t>& x)
{
   _theFct.set(exp,_var);
   for (size_t i=0; i<x.size(); ++i)
      set(i+1,_theFct(x[i],0.,0.,_time));
}


template <class T_>
void Vect<T_>::set(Mesh&  ms,
                   const  string& exp,
                   size_t dof)
{
   setMesh(ms);
   if (_theMesh==nullptr)
      throw OFELIException("In Vect::set(ms,string,dof): No mesh defined");
   _theFct.set(exp,_var);
   if (_dof_type==NODE_DOF) {
      MESH_ND
         set(The_node.getNbDOF()*(node_label-1)+dof,_theFct(The_node.getCoord(),_time));
   }
   else
      throw OFELIException("In Vect::set(string,size_t): This member function is "
                           "for nodewise vectors only.");
}


template <class T_>
void Vect<T_>::set(const Vect<real_t>& x,
                   const string&       exp)
{
   setSize(x._nx,x._ny,x._nz);
   _theFct.set(exp,_var_xit);
   vector<real_t> xv(3);
   for (size_t i=0; i<_size; i++) {
      xv[0] = x[i], xv[1] = i+1, xv[2] = _time;
      set(i+1,_theFct(xv));
   }
}


template<class T_>
void Vect<T_>::setMesh(Mesh&      m,
                       DOFSupport dof_type,
                       size_t     nb_dof)
{
   _theMesh = &m;
   _with_mesh = true;
   _nb_dof = nb_dof;
   if (nb_dof==0)
      _nb_dof = _theMesh->getNbDOF()/_theMesh->getNbNodes();
   _dof_type = dof_type;
   if (dof_type==NODE_DOF)
      _nb = _theMesh->getNbNodes();
   else if (dof_type==SIDE_DOF) {
      _theMesh->getAllSides();
      _nb = _theMesh->getNbSides();
   }
   else if (dof_type==BOUNDARY_SIDE_DOF) {
      _theMesh->getAllSides();
      _nb = _theMesh->getNbSides();
   }
   else if (dof_type==ELEMENT_DOF)
      _nb = _theMesh->getNbElements();
   setSize(_nb,_nb_dof,1);
   clear();
}


template<class T_>
void Vect<T_>::setGrid(Grid& g)
{
   setSize(g.getNx()+1,g.getNy()+1,g.getNz()+1);
   clear();
}


template<class T_>
size_t Vect<T_>::size() const { return _size; }


template<class T_>
void Vect<T_>::setSize(size_t nx,
                       size_t ny,
                       size_t nz,
                       size_t nt)
{
   _nx = nx, _ny = ny, _nz = nz, _nt = nt;
   _size = _nx*_ny*_nz*_nt;
#if defined (USE_EIGEN)
   _v.conservativeResize(_size,1);
   clear();
#else
   vector<T_>::resize(_size);
   clear();
#endif
}


template<class T_>
void Vect<T_>::resize(size_t n) { setSize(n); }


template<class T_>
void Vect<T_>::resize(size_t n,
                      T_     v)
{
   _nx = n, _ny = _nz = 1;
   _size = _nx*_ny*_nz;
#if defined (USE_EIGEN)
   _v.conservativeResize(_size,1);
   for (size_t i=1; i<=_size; i++)
      set(i,v);
#else
   vector<T_>::resize(_size,v);
#endif
}


template<class T_>
void Vect<T_>::setDOFType(DOFSupport dof_type) { _dof_type = dof_type; }


template<class T_>
void Vect<T_>::setDG(int degree)
{
   if (_with_mesh==false)
      throw OFELIException("In Vect::setDG(int): To be used only if a Mesh instance "
                           "is associated to the vector");
   _dg_degree = degree;
   if (degree<0)
      return;
   _nb_dof = 0;
   _dof_type = ELEMENT_DOF;
   switch (_theMesh->getDim()) {

      case 1:
         _nb_dof = _dg_degree+1;
         break;

      case 2:
         if (_dg_degree<10)
            _nb_dof = (_dg_degree+1)*(_dg_degree+2)/2;
         else
            _nb_dof = (_dg_degree+1)*(_dg_degree+1);
         break;

      case 3:
         if (_dg_degree<10)
            _nb_dof = (_dg_degree+1)*(_dg_degree+2)/2;
         else if (_dg_degree<20)
            _nb_dof = (_dg_degree+1)*(_dg_degree+1)*(_dg_degree+1);
         else
            _nb_dof = (_dg_degree+1)*(_dg_degree+1)*(_dg_degree+2)/2;
         break;
   }
   setSize(_theMesh->getNbElements(),_nb_dof,1);
}


template<class T_>
bool Vect<T_>::isGrid() const { return _grid; }


template<class T_>
size_t Vect<T_>::getNbDOF() const { return _nb_dof; }


template<class T_>
size_t Vect<T_>::getNb() const { return size()/_nb_dof; }


template<class T_>
Mesh& Vect<T_>::getMesh() const { return *_theMesh; }


template<class T_>
bool Vect<T_>::WithMesh() const { return _with_mesh; }


template<class T_>
DOFSupport Vect<T_>::getDOFType() const { return _dof_type; }


template<class T_>
void Vect<T_>::setTime(real_t t) { _time = t; }


template<class T_>
real_t Vect<T_>::getTime() const { return _time; }


template<class T_>
void Vect<T_>::setName(string name) { _name = name; }


template<class T_>
string Vect<T_>::getName() const { return _name; }


template<>
inline real_t Vect<real_t>::getNorm1() const
{
   real_t s=0.;
   for (size_t i=0; i<size(); ++i)
      s += std::abs((*this)[i]);
   return s;
}


template<>
inline real_t Vect<real_t>::getNorm2() const
{
   real_t s=0;
   for (size_t i=0; i<size(); ++i)
      s += (*this)[i]*(*this)[i];
   s = std::sqrt(s);
   return s;
}


template<>
inline real_t Vect<complex_t>::getNorm2() const
{
   real_t s=0;
   for (size_t i=0; i<size(); ++i) {
      complex_t z = (*this)[i];
      s += z.real()*z.real() + z.imag()*z.imag();
   }
   return std::sqrt(s);
}


template<>
inline real_t Vect<real_t>::getNormMax() const
{
   real_t s=0;
   for (size_t i=0; i<size(); ++i) {
      real_t z = std::abs((*this)[i]);
      s = (z > s ? z : s);
   }
   return s;
}


template<>
inline real_t Vect<complex_t>::getNormMax() const
{
   real_t s=0;
   for (size_t i=0; i<size(); ++i) {
      real_t z = std::abs((*this)[i]);
      s = (z > s ? z : s);
   }
   return s;
}


template<>
inline real_t Vect<real_t>::Norm(NormType t) const
{
   if (t==NORM1)
      return getNorm1();
   else if (t==WNORM1)
      return getWNorm1();
   else if (t==NORM2)
      return getNorm2();
   else if (t==WNORM2)
      return getWNorm2();
   else if (t==NORM_MAX)
      return getNormMax();
   else
      return 0.;
}


template<class T_>
real_t Vect<T_>::getWNorm1() const { return getNorm1()/size(); }


template<class T_>
real_t Vect<T_>::getWNorm2() const { return getNorm2()/sqrt(real_t(size())); }


template<class T_>
T_ Vect<T_>::getMin() const
{
   T_ s = (*this)[0];
   for (size_t i=1; i<size(); ++i)
      s = (*this)[i] < s ? (*this)[i] : s;
   return s;
}


template<class T_>
T_ Vect<T_>::getMax() const
{
  T_ s = (*this)[0];
   for (size_t i=1; i<size(); ++i)
      s = (*this)[i] > s ? (*this)[i] : s;
   return s;
}


template<class T_>
size_t Vect<T_>::getNx() const { return _nx; }


template<class T_>
size_t Vect<T_>::getNy() const { return _ny; }


template<class T_>
size_t Vect<T_>::getNz() const { return _nz; }


template<class T_>
size_t Vect<T_>::getNt() const { return _nt; }


template <class T_>
void Vect<T_>::setIJKL(const string& exp)
{
   _theFct.set(exp,_var_ijkt);
   vector<real_t> xv(4);
   for (size_t i=1; i<=_nx; ++i) {
      for (size_t j=1; j<=_ny; ++j) {
         for (size_t k=1; k<=_nz; ++k) {
            for (size_t l=1; l<=_nt; ++l) {
               xv[0] = i, xv[1] = j, xv[2] = k, xv[3] = _time;
               set(i,j,k,l,_theFct(xv));
            }
	 }
      }
   }
}


template <class T_>
void Vect<T_>::setIJK(const string& exp)
{
   _theFct.set(exp,_var_ijkt);
   vector<real_t> xv(4);
   for (size_t i=1; i<=_nx; ++i) {
      for (size_t j=1; j<=_ny; ++j) {
         for (size_t k=1; k<=_nz; ++k) {
            xv[0] = i, xv[1] = j, xv[2] = k, xv[3] = _time;
            set(i,j,k,_theFct(xv));
         }
      }
   }
}


template<class T_>
void Vect<T_>::setNodeBC(Mesh&  m,
                         int    code,
                         T_     val,
                         size_t dof)
{
   if (m.getDOFSupport()==NODE_DOF) {
      node_loop(&m) {
         for (size_t i=1; i<=The_node.getNbDOF(); i++) {
            if (The_node.getCode(dof)==code && code!=0)
               set(node_label,dof,val);
         }
      }
   }
   if (m.getDOFSupport()==SIDE_DOF) {
      boundary_side_loop(&m) {
         for (size_t i=1; i<=The_side.getNbDOF(); i++) {
            if (The_side.getCode(dof)==code && code!=0)
               set(side_label,dof,val);
         }
      }
   }
}


template<class T_>
void Vect<T_>::setNodeBC(Mesh& m,
                         int   code,
                         T_    val)
{
   int c[10];
   if (m.getDOFSupport()==NODE_DOF) {
      node_loop(&m) {
         DOFCode(code,The_node.getNbDOF(),c);
         for (size_t i=1; i<=The_node.getNbDOF(); i++) {
            if (The_node.getCode(i)==c[i-1] && c[i-1]!=0)
               set(node_label,i,val);
         }
      }
   }
   if (m.getDOFSupport()==SIDE_DOF) {
      boundary_side_loop(&m) {
         DOFCode(code,The_side.getNbDOF(),c);
         for (size_t i=1; i<=The_side.getNbDOF(); ++i) {
            if (The_side.getCode(i)==c[i-1] && c[i-1]!=0)
               set(side_label,i,val);
         }
      }
   }
}


template<class T_>
void Vect<T_>::setSideBC(Mesh&  m,
                         int    code,
                         T_     val,
                         size_t dof)
{
   setSideBC(*_theMesh,code,val,dof);
}


template<class T_>
void Vect<T_>::setSideBC(Mesh& m,
                         int   code,
                         T_    val)
{
   setSideBC(*_theMesh,code,val);
}


template<class T_>
void Vect<T_>::setNodeBC(Mesh&         m,
                         int           code,
                         const string& exp,
                         size_t        dof)
{
   _theFct.set(exp,_var);
   node_loop(&m) {
      for (size_t i=1; i<=The_node.getNbDOF(); i++) {
         if (The_node.getCode(dof)==code && code!=0)
            set(node_label,dof,_theFct(The_node.getCoord(),_time));
      }
   }
}


template<class T_>
void Vect<T_>::setNodeBC(Mesh&         m,
                         int           code,
                         const string& exp)
{
  _theFct.set(exp,_var);
   int c[6];
   node_loop(&m) {
      DOFCode(code,The_node.getNbDOF(),c);
      for (size_t i=1; i<=The_node.getNbDOF(); i++) {
         if (The_node.getCode(i)==c[i-1] && c[i-1]!=0)
            set(node_label,i,_theFct(The_node.getCoord(),_time));
      }
   }
}


template<class T_>
void Vect<T_>::setSideBC(Mesh&         m,
                         int           code,
                         const string& exp,
                         size_t        dof)
{
  _theFct.set(exp,_var);
   side_loop(&m) {
      if (The_side.getCode(dof)==code && code!=0) {
         for (size_t i=1; i<=The_side.getNbDOF(); i++)
            set(side_label,dof,_theFct(The_side.getCenter(),_time));
      }
   }
}


template<class T_>
void Vect<T_>::setSideBC(Mesh&         m,
                         int           code,
                         const string& exp)
{
   _theFct.set(exp,_var);
   int c[6];
   side_loop(&m) {
      DOFCode(code,The_side.getNbDOF(),c);
      for (size_t i=1; i<=The_side.getNbDOF(); ++i) {
         if (The_side.getCode(i)==c[i-1] && c[i-1]!=0) {
            set(side_label,i,_theFct(The_side.getCenter(),_time));
         }
      }
   }
}


template<class T_>
void Vect<T_>::setNodeBC(int    code,
                         T_     val,
                         size_t dof)
{
   setNodeBC(*_theMesh,code,val,dof);
}


template<class T_>
void Vect<T_>::setNodeBC(int code,
                         T_  val)
{
   setNodeBC(*_theMesh,code,val);
}


template<class T_>
void Vect<T_>::setNodeBC(int           code,
                         const string& exp,
                         size_t        dof)
{
   setNodeBC(*_theMesh,code,exp,dof);
}


template<class T_>
void Vect<T_>::setNodeBC(int           code,
                         const string& exp)
{
   setNodeBC(*_theMesh,code,exp);
}


template<class T_>
void Vect<T_>::setSideBC(int           code,
                         const string& exp,
                         size_t        dof)
{
   setSideBC(*_theMesh,code,exp,dof);
}


template<class T_>
void Vect<T_>::setSideBC(int           code,
                         const string& exp)
{
   setSideBC(*_theMesh,code,exp);
}


template<class T_>
void Vect<T_>::setSideBC(int     code,
                         T_      val,
                         size_t  dof)
{
   setSideBC(*_theMesh,code,val,dof);
}


template<class T_>
void Vect<T_>::removeBC(const Mesh&     ms,
                        const Vect<T_>& v,
                        int             dof)
{
   if (dof==0) {
      size_t n = 1;
      node_loop(&ms) {
         for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
            if (The_node.getCode(k) == 0)
               set(The_node.getDOF(k),v(n));
            n++;
         }
      }
   }
   else {
      node_loop(&ms) {
         if (The_node.getCode(dof) == 0)
            set(The_node.getDOF(dof),v(node_label));
      }
   }
}


template<class T_>
void Vect<T_>::removeBC(const Vect<T_>& v,
                        int             dof)
{
   if (dof==0) {
      size_t n = 1;
      MESH_ND {
         for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
            if (The_node.getCode(k) == 0)
               set(The_node.getDOF(k),v(n));
            n++;
         }
      }
   }
   else {
      MESH_ND {
         if (The_node.getCode(dof) == 0)
            set(The_node.getDOF(dof),v(node_label));
      }
   }
}


template<class T_>
void Vect<T_>::transferBC(const Vect<T_>& bc,
                          int             dof)
{
   size_t i=1, k=1;
   if (dof==0) {
      if (_theMesh->NodesAreDOF()) {
         MESH_ND {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               if (The_node.getCode(k)>0)
                  set(i,bc(i));
               i++;
            }
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         MESH_SD {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               if (The_side.getCode(k)>0)
                  set(i,bc(i));
               i++;
            }
         }
      }
      else
         ;
   }
   else {
      if (_theMesh->NodesAreDOF()) {
         MESH_ND {
            if (The_node.getCode(dof)>0)
               set(i,bc(k));
            i++;
            k += The_node.getNbDOF();
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         MESH_SD {
            if (The_side.getCode(dof)>0)
               set(i,bc(k));
            i++, k += The_side.getNbDOF();
         }
      }
      else
         ;
   }
}


template<class T_>
void Vect<T_>::insertBC(Mesh&           m,
                        const Vect<T_>& v,
                        const Vect<T_>& bc,
                        int             dof)
{
   size_t i=1, j=1;
   if (dof==0) {
      if (m.NodesAreDOF()) {
         node_loop(&m) {
           for (size_t k=1; k<=The_node.getNbDOF(); ++k, i++) {
               if (The_node.getCode(k)==0)
                  set(i,v(j++));
               else
                  set(i,bc(i));
            }
         }
      }
      else if (m.SidesAreDOF()) {
         side_loop(&m) {
           for (size_t k=1; k<=The_side.getNbDOF(); ++k, i++) {
               if (The_side.getCode(k)>=0)
                  set(i,v(j++));
               else
                  set(i,bc(i));
            }
         }
      }
      else
         ;
   }
   else {
      if (m.NodesAreDOF()) {
         size_t k=dof;
         node_loop(&m) {
            if (The_node.getCode(dof)==0)
               set(i,v(The_node.getDOF(dof)));
            else
               set(i,bc(k));
            i++, k += The_node.getNbDOF();
         }
      }
      else if (m.SidesAreDOF()) {
         size_t k=dof;
         side_loop(&m) {
            if (The_side.getCode(dof)>=0)
               set(i,v(The_side.getDOF(dof)));
            else
               set(i,bc(k));
            i++, k += The_side.getNbDOF();
         }
      }
      else
         ;
   }
}


template<class T_>
void Vect<T_>::setRegex(int dof)
{
   if (_theMesh==nullptr)
      throw OFELIException("In Vect::setRegex(int): No mesh defined");
   _with_regex[dof-1] = true;
}


template<class T_>
bool Vect<T_>::withRegex(int dof) const
{
   return _with_regex[dof-1];
}


template<class T_>
void Vect<T_>::insertBC(Mesh&           m,
                        const Vect<T_>& v,
                        int             dof)
{
   size_t i=1;
   if (dof==0) {
      if (m.NodesAreDOF()) {
         node_loop(&m) {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               set(i,0);
               if (The_node.getCode(k)==0)
                  set(i,v(The_node.getDOF(k)));
               i++;
            }
         }
      }
      else if (m.SidesAreDOF()) {
         side_loop(&m) {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               set(i,0);
               if (The_side.getCode(k)==0)
                  set(i,v(The_side.getDOF(k)));
               i++;
            }
         }
      }
      else
         ;
   }
   else {
      if (m.NodesAreDOF()) {
         node_loop(&m) {
            set(i,0);
            if (The_node.getCode(dof)==0)
               set(i,v(The_node.getDOF(dof)));
            i++;
         }
      }
      else if (m.SidesAreDOF()) {
         side_loop(&m) {
            set(i,0);
            if (The_side.getCode(dof)==0)
               set(i,v(The_side.getDOF(dof)));
            i++;
         }
      }
      else
         ;
   }
}


template<class T_>
void Vect<T_>::insertBC(const Vect<T_>& v,
                        const Vect<T_>& bc,
                        int             dof)
{
   size_t i = 1;
   if (dof==0) {
      if (_theMesh->NodesAreDOF()) {
         MESH_ND {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               if (The_node.getCode(k)==0)
                  set(i,v(The_node.getDOF(k)));
               else
                  set(i,bc(i));
               i++;
            }
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         MESH_SD {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               if (The_side.getCode(k)>=0)
                  set(i,v(The_side.getDOF(k)));
               else
                  set(i,bc(i));
               i++;
            }
         }
      }
      else
         ;
   }
   else {
      if (_theMesh->NodesAreDOF()) {
         size_t k=dof;
         MESH_ND {
            if (The_node.getCode(dof)==0)
               set(i,v(The_node.getDOF(dof)));
            else
               set(i,bc(k));
            i++, k += The_node.getNbDOF();
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         size_t k=dof;
         MESH_SD {
            if (The_side.getCode(dof)>=0)
               set(i,v(The_side.getDOF(dof)));
            else
               set(i,bc(The_side.getDOF(dof)));
            i++, k += The_side.getNbDOF();
         }
      }
      else
         ;
   }
}


template<class T_>
void Vect<T_>::insertBC(const Vect<T_>& v,
                        int             dof)
{
   size_t i=1;
   if (dof==0) {
      if (_theMesh->NodesAreDOF()) {
         MESH_ND {
           for (size_t k=1; k<=The_node.getNbDOF(); ++k, i++) {
               set(i,0);
               if (The_node.getCode(k)==0)
                  set(i,v(The_node.getDOF(k)));
            }
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         MESH_SD {
           for (size_t k=1; k<=The_side.getNbDOF(); ++k, i++) {
               set(i,0);
               if (The_side.getCode(k)==0)
                  set(i,v(The_side.getDOF(k)));
            }
         }
      }
      else
         ;
   }
   else {
      if (_theMesh->NodesAreDOF()) {
         MESH_ND {
            set(i,0);
            if (The_node.getCode(dof)==0)
               set(i,v(The_node.getDOF(dof)));
            i++;
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         MESH_SD {
            set(i,0);
            if (The_side.getCode(dof)==0)
               set(i,v(The_side.getDOF(dof)));
            i++;
         }
      }
      else
         ;
   }
}


#if defined(USE_PETSC)
template<class T_>
void Vect<T_>::insertBC(Mesh&                m,
                        const PETScVect<T_>& v,
                        const Vect<T_>&      bc,
                        int                  dof)
{
   size_t i=1;
   if (dof==0) {
      if (m.NodesAreDOF()) {
         node_loop(&m) {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               if (The_node.getCode(k)==0)
                  set(i,v(The_node.getDOF(k)));
               else
                  set(i,bc(i));
               i++;
            }
         }
      }
      else if (m.SidesAreDOF()) {
         side_loop(&m) {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               if (The_side.getCode(k)>=0)
                  set(i,v(The_side.getDOF(k)));
               else
                  set(i,bc(i));
               i++;
            }
         }
      }
      else
         ;
   }
   else {
      if (m.NodesAreDOF()) {
         size_t k=dof;
         node_loop(&m) {
            if (The_node.getCode(dof)==0)
               set(i,v(The_node.getDOF(dof)));
            else
               set(i,bc(k));
            i++, k += The_node.getNbDOF();
         }
      }
      else if (m.SidesAreDOF()) {
         size_t k=dof;
         side_loop(&m) {
            if (The_side.getCode(dof)>=0)
               set(i,v(The_side.getDOF(dof)));
            else
               set(i,bc(k));
            i++, k += The_side.getNbDOF();
         }
      }
      else
         ;
   }
}


template<class T_>
void Vect<T_>::insertBC(Mesh&                m,
                        const PETScVect<T_>& v,
                        int                  dof)
{
   size_t i=1;
   if (dof==0) {
      if (m.NodesAreDOF()) {
         node_loop(&m) {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               set(i,0);
               if (The_node.getCode(k)==0)
                  set(i,v(The_node.getDOF(k)));
               i++;
            }
         }
      }
      else if (m.SidesAreDOF()) {
         side_loop(&m) {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               set(i,0);
               if (The_side.getCode(k)==0)
                  set(i,v(The_side.getDOF(k)));
               i++;
            }
         }
      }
      else
         ;
   }
   else {
      if (m.NodesAreDOF()) {
         node_loop(&m) {
            set(i,0);
            if (The_node.getCode(dof)==0)
               set(i,v(The_node.getDOF(dof)));
            i++;
         }
      }
      else if (m.SidesAreDOF()) {
         side_loop(&m) {
            set(i,0);
            if (The_side.getCode(dof)==0)
               set(i,v(The_side.getDOF(dof)));
            i++;
         }
      }
      else
         ;
   }
}


template<class T_>
void Vect<T_>::insertBC(const PETScVect<T_>& v,
                        const Vect<T_>&      bc,
                        int                  dof)
{
   size_t i = 1;
   if (dof==0) {
      if (_theMesh->NodesAreDOF()) {
         MESH_ND {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               if (The_node.getCode(k)==0)
                  set(i,v(The_node.getDOF(k)));
               else
                  set(i,bc(i));
               i++;
            }
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         MESH_SD {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               if (The_side.getCode(k)>=0)
                  set(i,v(The_side.getDOF(k)));
               else
                  set(i,bc(i));
               i++;
            }
         }
      }
      else
         ;
   }
   else {
      if (_theMesh->NodesAreDOF()) {
         size_t k=dof;
         MESH_ND {
            if (The_node.getCode(dof)==0)
               set(i,v(The_node.getDOF(dof)));
            else
               set(i,bc(k));
            i++, k += The_node.getNbDOF();
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         size_t k=dof;
         MESH_SD {
            if (The_side.getCode(dof)>=0)
               set(i,v(The_side.getDOF(dof)));
            else
               set(i,bc(The_side.getDOF(dof)));
            i++, k += The_side.getNbDOF();
         }
      }
      else
         ;
   }
}
#endif


template<class T_>
void Vect<T_>::Assembly(const Element&  el,
                        const Vect<T_>& b)
{
   size_t i=1;
   for (size_t n=1; n<=el.getNbNodes(); ++n, i++) {
      Node *nd=el(n);
      for (size_t k=1; k<=nd->getNbDOF(); ++k) {
         if (nd->getDOF(k))
            add(nd->getDOF(k),b(i));
      }
   }
}


template<class T_>
void Vect<T_>::Assembly(const Element& el,
                        const T_*      b)
{
   size_t i=0;
   for (size_t n=1; n<=el.getNbNodes(); ++n) {
      Node *nd = el(n);
      for (size_t k=1; k<=nd->getNbDOF(); ++k, ++i) {
         if (nd->getDOF(k))
            add(nd->getDOF(k),b[i]);
      }
   }
}


template<class T_>
void Vect<T_>::Assembly(const Side&     sd,
                        const Vect<T_>& b)
{
   size_t i=0;
   for (size_t n=1; n<=sd.getNbNodes(); ++n) {
      Node *nd = sd(n);
      for (size_t k=1; k<=nd->getNbDOF(); ++k, ++i) {
         if (nd->getDOF(k))
            add(nd->getDOF(k),b[i]);
      }
   }
}


template<class T_>
void Vect<T_>::Assembly(const Side& sd,
                        const T_*   b)
{
   size_t i=0;
   for (size_t n=1; n<=sd.getNbNodes(); ++n) {
      Node *nd = sd(n);
      for (size_t k=1; k<=nd->getNbDOF(); ++k, ++i) {
         if (nd->getDOF(k))
            add(nd->getDOF(k),b[i]);
      }
   }
}


template<class T_>
void Vect<T_>::DGAssembly(const Element&                          el,
                          const LocalVect<T_,MAX_NB_ELEMENT_DOF>& b)
{
   for (size_t i=1; i<=el.getNbDOF(); ++i) {
      if (el.getDOF(i)!=0)
         add(el.getDOF(i),b(i));
   }
}


template<class T_>
void Vect<T_>::DGAssembly(const Side&                          sd,
                          const LocalVect<T_,MAX_NB_SIDE_DOF>& b)
{
   for (size_t i=1; i<=sd.getNbDOF(); ++i) {
      if (sd.getDOF(i)!=0)
         add(sd.getDOF(i),b(i));
   }
}


template <class T_>
void Vect<T_>::getGradient(Vect<T_>& v)
{
   if (_theMesh==nullptr)
     throw OFELIException("In Vect::getGradient(Vect<>): No mesh defined for this vector.");
   T_ a;
   real_t b;
   Point<T_> aa;
   v.setMesh(*_theMesh,ELEMENT_DOF,_theMesh->getDim());
   v.setTime(_time);
   MESH_EL {
      if (The_element.getShape()==LINE) {
         a = (*this)(The_element(2)->n()) - (*this)(The_element(1)->n());
         b = The_element(2)->getCoord(1) - The_element(1)->getCoord(1);
         v.set(element_label,a/b);
      }
      else if (The_element.getShape()==TRIANGLE) {
         Triang3 t(the_element);
         vector<Point<real_t> > dsh = t.DSh();
         aa = (*this)(The_element(1)->n())*dsh[0] + 
              (*this)(The_element(2)->n())*dsh[1] + 
              (*this)(The_element(3)->n())*dsh[2];
         v.set(element_label,1,aa.x);
         v.set(element_label,2,aa.y);
      }
      else if (The_element.getShape()==TETRAHEDRON) {
         Tetra4 t(the_element);
         vector<Point<real_t> > dsh = t.DSh();
         aa = (*this)(The_element(1)->n())*dsh[0] + 
              (*this)(The_element(2)->n())*dsh[1] + 
              (*this)(The_element(3)->n())*dsh[2] + 
              (*this)(The_element(4)->n())*dsh[3];
         v.set(element_label,1,aa.x);
         v.set(element_label,2,aa.y);
         v.set(element_label,3,aa.z);
      }
      else
         throw OFELIException("In Vect::getGradient(): This function doesn't work for this element.");
   }
}


template <class T_>
void Vect<T_>::getGradient(Vect<Point<T_> >& v)
{
   if (_theMesh==nullptr)
      throw OFELIException("In Vect::getGradient(Vect<>): No mesh defined for this vector.");
   T_ a;
   real_t b;
   Point<T_> aa;
   v.setMesh(*_theMesh,ELEMENT_DOF,_theMesh->getDim());
   v.setTime(_time);
   MESH_EL {
      if (The_element.getShape()==LINE) {
         a = (*this)(The_element(2)->n()) - (*this)(The_element(1)->n());
         b = The_element(2)->getCoord(1) - The_element(1)->getCoord(1);
         v.set(element_label,a/b);
      }
      else if (the_element->getShape()==TRIANGLE) {
         Triang3 t(the_element);
         vector<Point<real_t> > dsh = t.DSh();
         aa = (*this)(The_element(1)->n())*dsh[0] + 
              (*this)(The_element(2)->n())*dsh[1] + 
              (*this)(The_element(3)->n())*dsh[2];
         v.set(element_label,aa);
      }
      else if (_theMesh->getShape()==TETRAHEDRON) {
         Tetra4 t(the_element);
         vector<Point<real_t> > dsh = t.DSh();
         aa = (*this)(The_element(1)->n())*dsh[0] +
              (*this)(The_element(2)->n())*dsh[1] +
              (*this)(The_element(3)->n())*dsh[2] +
              (*this)(The_element(4)->n())*dsh[3];
         v.set(element_label,aa);
      }
      else
         throw OFELIException("In Vect::getGradient(): This function doesn't work for this element.");
   }
}


template <class T_>
void Vect<T_>::getCurl(Vect<T_>& v)
{
   if (_theMesh==nullptr)
      throw OFELIException("In Vect::getCurl(Vect<>): No mesh defined for this vector.");
   T_ a;
   real_t b;
   Point<T_> du, dv, dw;
   v.setMesh(*_theMesh,ELEMENT_DOF,_theMesh->getDim());
   v.setTime(_time);
   MESH_EL {
      if (The_element.getShape()==LINE) {
         a = (*this)(The_element(2)) - (*this)(The_element(1));
         b = The_element(2)->getCoord(1) - The_element(1)->getCoord(1);
         v.set(element_label,a/b);
      }
      else if (The_element.getShape()==TRIANGLE) {
         Triang3 t(the_element);
         vector<Point<real_t> > dsh = t.DSh();
         du = (*this)(The_element(1)->n())*dsh[0] + 
              (*this)(The_element(2)->n())*dsh[1] + 
              (*this)(The_element(3)->n())*dsh[2];
         v.set(element_label,1, du.y);
         v.set(element_label,2,-du.x);
      }
      else if (The_element.getShape()==TETRAHEDRON) {
         Tetra4 t(the_element);
         vector<Point<real_t> > dsh = t.DSh();
         du = (*this)(The_element(1)->n(),1)*dsh[0] + 
              (*this)(The_element(2)->n(),1)*dsh[1] + 
              (*this)(The_element(3)->n(),1)*dsh[2] + 
              (*this)(The_element(4)->n(),1)*dsh[3];
         dv = (*this)(The_element(1)->n(),2)*dsh[0] + 
              (*this)(The_element(2)->n(),2)*dsh[1] + 
              (*this)(The_element(3)->n(),2)*dsh[2] + 
              (*this)(The_element(4)->n(),2)*dsh[3];
         dw = (*this)(The_element(1)->n(),3)*dsh[0] + 
              (*this)(The_element(2)->n(),3)*dsh[1] + 
              (*this)(The_element(3)->n(),3)*dsh[2] + 
              (*this)(The_element(4)->n(),3)*dsh[3];
         v.set(element_label,1,dw.y - dv.z);
         v.set(element_label,2,du.z - dw.x);
         v.set(element_label,3,dv.x - du.y);
      }
      else
         throw OFELIException("In Vect::getCurl(): This function doesn't work for this element.");
   }
}


template <class T_>
void Vect<T_>::getCurl(Vect<Point<T_> >& v)
{
   if (_theMesh==nullptr)
      throw OFELIException("In Vect::getCurl(Vect<>): No mesh defined for this vector.");
   T_ a;
   real_t b;
   Point<T_> du, dv, dw;
   v.setMesh(*_theMesh,ELEMENT_DOF,_theMesh->getDim());
   v.setTime(_time);
   MESH_EL {
      if (The_element.getShape()==LINE) {
         a = (*this)(The_element(2)->n()) - (*this)(The_element(1)->n());
         b = The_element(2)->getCoord(1) - The_element(1)->getCoord(1);
         v.set(element_label,a/b);
      }
      else if (The_element.getShape()==TRIANGLE) {
         Triang3 t(the_element);
         vector<Point<real_t> > dsh = t.DSh();
         du = (*this)(The_element(1)->n())*dsh[0] + 
              (*this)(The_element(2)->n())*dsh[1] + 
              (*this)(The_element(3)->n())*dsh[2];
         v.set(element_label,1, du.y);
         v.set(element_label,2,-du.x);
      }
      else if (The_element.getShape()==TETRAHEDRON) {
         Tetra4 t(the_element);
         vector<Point<real_t> > dsh = t.DSh();
         du = (*this)(The_element(1)->n(),1)*dsh[0] + (*this)(The_element(2)->n(),1)*dsh[1] + 
              (*this)(The_element(3)->n(),1)*dsh[2] + (*this)(The_element(4)->n(),1)*dsh[3];
         dv = (*this)(The_element(1)->n(),2)*dsh[0] + (*this)(The_element(2)->n(),2)*dsh[1] + 
              (*this)(The_element(3)->n(),2)*dsh[2] + (*this)(The_element(4)->n(),2)*dsh[2];
         dw = (*this)(The_element(1)->n(),3)*dsh[0] + (*this)(The_element(2)->n(),3)*dsh[1] + 
              (*this)(The_element(3)->n(),3)*dsh[2] + (*this)(The_element(4)->n(),3)*dsh[3];
         v.set(element_label,Point<T_>(dw.y-dv.z,du.z-dw.x,dv.x-du.y));
      }
      else
         throw OFELIException("In Vect::getCurl(): This function doesn't work for this element.");
   }
}


template <class T_>
void Vect<T_>::getSCurl(Vect<T_>& v)
{
   if (_theMesh==nullptr)
      throw OFELIException("In Vect::getSCurl(Vect<>): No mesh defined for this vector.");
   if (_theMesh->getDim()==1 || _theMesh->getDim()==3)
      throw OFELIException("In Vect::getSCurl(): This function is valid for 2-D only.");
   Point<T_> du, dv;
   v.setMesh(*_theMesh,ELEMENT_DOF,1);
   v.setTime(_time);
   MESH_EL {
      if (The_element.getShape()==TRIANGLE) {
         Triang3 t(the_element);
         vector<Point<real_t> > dsh = t.DSh();
         du = (*this)(The_element(1)->n(),1)*dsh[0] + (*this)(The_element(2)->n(),1)*dsh[1] + 
              (*this)(The_element(3)->n(),1)*dsh[2];
         dv = (*this)(The_element(1)->n(),2)*dsh[0] + (*this)(The_element(2)->n(),2)*dsh[1] + 
              (*this)(The_element(3)->n(),2)*dsh[2];
         v.set(element_label,dv.x - du.y);
      }
      else
         throw OFELIException("In Vect::getSCurl(): This function doesn't work for this element");
   }
}


template <class T_>
void Vect<T_>::getDivergence(Vect<T_>& v)
{
   if (_theMesh==nullptr)
      throw OFELIException("In Vect::getDivergence(Vect<>): No mesh defined for this vector.");
   T_ a;
   real_t b;
   Point<T_> du, dv, dw;
   v.setMesh(*_theMesh,ELEMENT_DOF,1);
   v.setTime(_time);
   MESH_EL {
      if (The_element.getShape()==LINE) {
         a = (*this)(The_element(2)->n()) - (*this)(The_element(1)->n());
         b = The_element(2)->getX() - The_element(1)->getX();
         v.set(element_label,a/b);
      }
      else if (The_element.getShape()==TRIANGLE) {
         Triang3 t(the_element);
         vector<Point<real_t> > dsh = t.DSh();
         du = (*this)(The_element(1)->n(),1)*dsh[0] + (*this)(The_element(2)->n(),1)*dsh[1] + 
              (*this)(The_element(3)->n(),1)*dsh[2];
         dv = (*this)(The_element(1)->n(),2)*dsh[0] + (*this)(The_element(2)->n(),2)*dsh[1] + 
              (*this)(The_element(3)->n(),2)*dsh[2];
         v.set(element_label,du.x + dv.y);
      }
      else if (The_element.getShape()==TETRAHEDRON) {
         Tetra4 t(the_element);
         vector<Point<real_t> > dsh = t.DSh();
         du = (*this)(The_element(1)->n(),1)*dsh[0] + (*this)(The_element(2)->n(),1)*dsh[1] + 
              (*this)(The_element(3)->n(),1)*dsh[1] + (*this)(The_element(4)->n(),1)*dsh[3];
         dv = (*this)(The_element(1)->n(),2)*dsh[0] + (*this)(The_element(2)->n(),2)*dsh[1] + 
              (*this)(The_element(3)->n(),2)*dsh[2] + (*this)(The_element(4)->n(),2)*dsh[3];
         dw = (*this)(The_element(1)->n(),3)*dsh[0] + (*this)(The_element(2)->n(),3)*dsh[1] + 
              (*this)(The_element(3)->n(),3)*dsh[2] + (*this)(The_element(4)->n(),3)*dsh[3];
         v.set(element_label,du.x + dv.y + dw.z);
      }
      else
         throw OFELIException("In Vect::getDivergence(): This function doesn't work for this element.");
   }
}


template<class T_>
real_t Vect<T_>::getAverage(const Element& el,
                            int            type) const
{
   switch (type) {

      case LINE2:
         return 0.5*((*this)(el(1)->n())+(*this)(el(2)->n()));

      case TRIANG3: 
         return OFELI_THIRD*((*this)(el(1)->n()) + (*this)(el(2)->n()) +
                             (*this)(el(3)->n()));

      case QUAD4:
         return 0.25*((*this)(el(1)->n()) + (*this)(el(2)->n()) +
                      (*this)(el(3)->n()) + (*this)(el(4)->n()));

      case TETRA4:
         return 0.25*((*this)(el(1)->n()) + (*this)(el(2)->n()) +
                      (*this)(el(3)->n()) + (*this)(el(4)->n()));

      case PENTA6:
         return OFELI_SIXTH*((*this)(el(1)->n()) + (*this)(el(2)->n()) +
                             (*this)(el(3)->n()) + (*this)(el(4)->n()) +
                             (*this)(el(5)->n()) + (*this)(el(6)->n()));

      case HEXA8:
         return 0.125*((*this)(el(1)->n()) + (*this)(el(2)->n()) +
                       (*this)(el(3)->n()) + (*this)(el(4)->n()) +
                       (*this)(el(5)->n()) + (*this)(el(6)->n()) +
                       (*this)(el(7)->n()) + (*this)(el(8)->n()));
   }
   return 0.;
}
   

template<class T_>
Vect<T_> &Vect<T_>::MultAdd(const Vect<T_>& x,
                            const T_&       a)
{
   for (size_t i=1; i<=size(); ++i)
      add(i,a*x(i));
   return *this;
}


template<class T_>
void Vect<T_>::Axpy(T_              a,
                    const Vect<T_>& x)
{
   for (size_t i=1; i<=size(); i++)
      add(i,a*x(i));
}


template<class T_>
void Vect<T_>::set(size_t i,
                   T_     val)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=size());
#endif
#if defined (USE_EIGEN)
   _v(i-1) = val;
#else
   (*this)[i-1] = val;
#endif
}


template<class T_>
void Vect<T_>::set(size_t i,
                   size_t j,
                   T_     val)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=_nx);
   assert(j>0 && j<=_ny);
#endif
#if defined (USE_EIGEN)
   _v(loc(i,j)) = val;
#else
   (*this)[loc(i,j)] = val;
#endif
}


template<class T_>
void Vect<T_>::set(size_t i,
                   size_t j,
                   size_t k,
                   T_     val)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=_nx);
   assert(j>0 && j<=_ny);
   assert(k>0 && k<=_nz);
#endif
#if defined (USE_EIGEN)
   _v(loc(i,j,k)) = val;
#else
   (*this)[loc(i,j,k)] = val;
#endif
}


template<class T_>
void Vect<T_>::add(size_t i,
                   T_     val)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=size());
#endif
#if defined (USE_EIGEN)
   _v(i-1) += val;
#else
   (*this)(i) += val;
#endif
}


template<class T_>
void Vect<T_>::add(size_t i,
                   size_t j,
                   T_     val)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=_nx);
   assert(j>0 && j<=_ny);
#endif
#if defined (USE_EIGEN)
   _v(loc(i,j)) += val;
#else
   (*this)(i,j) += val;
#endif
}


template<class T_>
void Vect<T_>::add(size_t i,
                   size_t j,
                   size_t k,
                   T_     val)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=_nx);
   assert(j>0 && j<=_ny);
   assert(k>0 && k<=_nz);
#endif
#if defined (USE_EIGEN)
   _v(loc(i,j)) = val;
#else
   (*this)[loc(i,j,k)] += val;
#endif
}


template<class T_>
void Vect<T_>::clear()
{
   for (size_t i=0; i<_size; i++)
#if defined (USE_EIGEN)
      _v[i] = static_cast<T_>(0);
#else
      (*this)[i] = static_cast<T_>(0);
#endif
}


class Fct;
template<>
inline void Vect<Fct>::clear()
{
}


struct fct;
template<>
inline void Vect<fct>::clear()
{
}


template<>
inline void Vect<string>::clear()
{
   for (size_t i=0; i<size(); ++i)
      (*this)[i] = " ";
}


#if defined (USE_EIGEN)
template<class T_>
T_ &Vect<T_>::operator[](size_t i)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=size());
#endif
   return _v[i];
}
#endif


#if defined (USE_EIGEN)
template<class T_>
T_ Vect<T_>::operator[](size_t i) const
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=size());
#endif
   return _v[i];
}
#endif


template<class T_>
T_ &Vect<T_>::operator()(size_t i)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=size());
#endif
#if defined (USE_EIGEN)
   return _v(i-1);
#else
   return (*this)[i-1];
#endif
}


template<class T_>
T_ Vect<T_>::operator()(size_t i) const
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=size());
#endif
#if defined (USE_EIGEN)
   return _v(i-1);
#else
   return (*this)[i-1];
#endif
}


template<class T_>
T_ &Vect<T_>::operator()(size_t i,
                         size_t j)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=_nx);
   assert(j>0 && j<=_ny);
#endif
#if defined (USE_EIGEN)
   return _v(loc(i,j));
#else
   return (*this)[loc(i,j)];
#endif
}


template<class T_>
T_ Vect<T_>::operator()(size_t i,
                        size_t j) const
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=_nx);
   assert(j>0 && j<=_ny);
#endif
#if defined (USE_EIGEN)
   return _v(loc(i,j));
#else
   return (*this)[loc(i,j)];
#endif
}


template<class T_>
T_ &Vect<T_>::operator()(size_t i,
                         size_t j,
                         size_t k)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=_nx);
   assert(j>0 && j<=_ny);
   assert(k>0 && k<=_nz);
#endif
#if defined (USE_EIGEN)
   return _v(loc(i,j,k));
#else
   return (*this)[loc(i,j,k)];
#endif
}


template<class T_>
T_ Vect<T_>::operator()(size_t i,
                        size_t j,
                        size_t k) const
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=_nx);
   assert(j>0 && j<=_ny);
   assert(k>0 && k<=_nz);
#endif
#if defined (USE_EIGEN)
   return _v(loc(i,j,k));
#else
   return (*this)[loc(i,j,k)];
#endif
}


template<class T_>
T_ &Vect<T_>::operator()(size_t i,
                         size_t j,
                         size_t k,
                         size_t l)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=_nx);
   assert(j>0 && j<=_ny);
   assert(k>0 && k<=_nz);
   assert(l>0 && l<=_nt);
#endif
#if defined (USE_EIGEN)
   return _v(loc(i,j,k,l));
#else
   return (*this)[loc(i,j,k,l)];
#endif
}


template<class T_>
T_ Vect<T_>::operator()(size_t i,
                        size_t j,
                        size_t k,
                        size_t l) const
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=_nx);
   assert(j>0 && j<=_ny);
   assert(k>0 && k<=_nz);
   assert(l>0 && l<=_nt);
#endif
#if defined (USE_EIGEN)
   return _v(loc(i,j,k,l));
#else
   return (*this)[loc(i,j,k,l)];
#endif
}


#if defined (USE_EIGEN)
template<class T_>
Vect<T_> &Vect<T_>::operator=(const VectorX& v)
{
   for (size_t i=0; i<_size; i++)
      (*this)[i] = v[i];
   return *this;
}
#endif


template<class T_>
Vect<T_> &Vect<T_>::operator=(const Vect<T_>& v)
{
   _theMesh = v._theMesh;
   _time = v._time;
   _name = v._name;
   _nb_dof = v._nb_dof;
   _nb = v._nb;
   _grid = v._grid;
   _nx = v._nx; _ny = v._ny; _nz = v._nz;
   _size = v._size;
#if defined (USE_EIGEN)
   setSize(_nx,_ny,_nz);
   for (size_t i=0; i<_size; i++)
      (*this)[i] = v[i];
#else
   for (size_t i=0; i<_size; i++)
      (*this)[i] = v[i];
#endif
   return *this;
}


template<class T_>
Vect<T_> &Vect<T_>::operator=(const T_& a)
{
   for (size_t i=1; i<=_size; ++i)
      set(i,a);
   return *this;
}


template<>
inline void Vect<real_t>::operator=(string s)
{
   if (_theMesh==NULL)
      throw OFELIException("In Vect::operator=(string): No mesh is defined");
   set(s);
}


template<class T_>
void Vect<T_>::setUniform(T_     vmin,
                          T_     vmax,
                          size_t n)
{
   setSize(n);
   for (size_t i=0; i<_nx; i++)
      (*this)[i] = T_(i)*(vmax-vmin)/T_(_nx-1);
}


template<class T_>
Vect<T_> &Vect<T_>::operator+=(const Vect<T_>& v)
{
#ifdef _OFELI_RANGE_CHECK
   assert(v.size() == _size);
#endif
   for (size_t i=1; i<=_size; ++i)
      add(i,v(i));
   return *this;
}


template<class T_>
Vect<T_> &Vect<T_>::operator+=(const T_& a)
{
   for (size_t i=1; i<=_size; i++)
      add(i,a);
   return *this;
}


template<class T_>
Vect<T_> &Vect<T_>::operator-=(const Vect<T_>& v)
{
#ifdef _OFELI_RANGE_CHECK
   assert(v.size() == _size);
#endif
   for (size_t i=1; i<=_size; ++i)
      add(i,-v(i));
   return *this;
}


template<class T_>
Vect<T_> &Vect<T_>::operator-=(const T_& a)
{
   for (size_t i=1; i<=_size; i++)
      add(i,-a);
   return *this;
}


template<class T_>
Vect<T_> &Vect<T_>::operator*=(const T_& a)
{
   for (size_t i=1; i<=_size; i++)
      set(i,(*this)(i)*a);
   return *this;
}


template<class T_>
Vect<T_> &Vect<T_>::operator/=(const T_& a)
{
   for (size_t i=1; i<=_size; i++)
      set(i,(*this)(i)/a);
   return *this;
}


template<class T_>
void Vect<T_>::push_back(const T_& v)
{
#if defined (USE_EIGEN)
       (*this)[_nx] = v;
#else
       vector<T_>::push_back(v);
#endif
       _nx++; _size++;
}


template<class T_>
const Mesh &Vect<T_>::getMeshPtr() const
{
   return *_theMesh;
}


template<class T_>
T_ Vect<T_>::operator,(const Vect<T_>& v) const
{
   T_ p = 0;
   for (size_t i=0; i<_size; i++)
      p += (*this)[i] * v[i];
   return p;
}


template<class T_>
Vect<complex_t> Vect<T_>::getFFT()
{
   void fft(vector<complex_t>& x);
   int logn = int(log(real_t(_size))/(log(2.0))+0.01);
   if (_size<2)
      throw OFELIException("In Vect<T_>::getFFT(v): Vector size is less than two.\n"
                           "Can't run FFT.");
   else if ((logn-1) > MAX_FFT_SIZE)
      throw OFELIException("In Vect<T_>::getFFT(v): FFT has too many points.");
   else if (int(pow(2.0,logn)) != _size)
      throw OFELIException("In Vect<T_>::getFFT(v): Vector size not a power of 2.\n"
                           "Can't run FFT.");

   Vect<complex_t> v(_size);
   for (size_t i=0; i<_size; i++)
      v[i] = (*this)[i];
 
// divide
   vector<complex_t> even, odd;
   OddEven(v,odd,even);

// conquer
   fft(even);
   fft(odd);

// combine
   for (size_t k=0; k<_size/2; ++k) {
      complex_t t = std::polar(1.0,-2*OFELI_PI*k/_size)*odd[k];
      v[k        ] = even[k] + t;
      v[k+_size/2] = even[k] - t;
   }
   return v;
}


template<class T_>
Vect<complex_t> Vect<T_>::getInvFFT()
{
   void fft(vector<complex_t>& x);
   int logn = int(log(real_t(_size))/(log(2.0))+0.01);
   if (_size<2)
      throw OFELIException("In Vect<T_>::getInvFFT(v): Vector size is less than two.\n"
                           "Can't run Inverse FFT.");

   else if ((logn-1) > MAX_FFT_SIZE)
      throw OFELIException("In Vect<T_>::getInvFFT(v): FFT has too many points.");

   else if (int(pow(2.0,logn)) != _size)
      throw OFELIException("In Vect<T_>::getInvFFT(v): Vector size is not a power of 2.\n"
                           "Can't run Inverse FFT.");

   Vect<complex_t> v(_size);
   for (size_t i=0; i<_size; i++)
      v[i] = complex_t(std::conj((*this)[i]));
 
// divide
   vector<complex_t> even, odd;
   OddEven(v,odd,even);

// conquer
   fft(even);
   fft(odd);

// combine
   for (size_t k=0; k<_size/2; ++k) {
      complex_t t = std::polar(1.0,-2*OFELI_PI*k/_size)*odd[k];
      v[k        ] = even[k] + t;
      v[k+_size/2] = even[k] - t;
   }

// conjugate and scale
   for (size_t i=0; i<_size; i++)
      v[i] = std::conj(v[i])*(1.0/_size);
   return v;
}

#if defined (USE_EIGEN)
template<class T_>
operator Vector<T_>::VectorX() const
{
   return _v;
}
#endif


///////////////////////////////////////////////////////////////////////////////
//                           ASSOCIATED  FUNCTIONS                           //
///////////////////////////////////////////////////////////////////////////////

template<class T_>
Vect<T_> operator+(const Vect<T_>& x,
                   const Vect<T_>& y)
{
#ifdef _OFELI_RANGE_CHECK
   assert(x.size() == y.size());
#endif
#if defined (USE_EIGEN)
   return Vect<T_>(Eigen::Matrix<T_,Eigen::Dynamic,1>(x)+Eigen::Matrix<T_,Eigen::Dynamic,1>(y));
#else
   Vect<T_> v(x);
   for (size_t i=0; i<x.size(); ++i)
      v.add(i+1,y[i]);
   return v;
#endif
}


template<class T_>
Vect<T_> operator-(const Vect<T_>& x,
                   const Vect<T_>& y)
{
#ifdef _OFELI_RANGE_CHECK
   assert(x.size()==y.size());
#endif
#if defined (USE_EIGEN)
   return Vect<T_>(Eigen::Matrix<T_,Eigen::Dynamic,1>(x)-Eigen::Matrix<T_,Eigen::Dynamic,1>(y));
#else
   Vect<T_> v(x);
   for (size_t i=0; i<x.size(); ++i)
      v.add(i+1,-y[i]);
   return v;
#endif
}


template<class T_>
Vect<T_> operator*(const T_&       a,
                   const Vect<T_>& x)
{
#if defined (USE_EIGEN)
   return a*x;
#else
   Vect<T_> v(x.size());
   for (size_t i=1; i<=x.size(); ++i)
      v.set(i,a*x(i));
   return v;
#endif
}


template<class T_>
Vect<T_> operator*(const Vect<T_>& x,
                   const T_&       a)
{
#if defined (USE_EIGEN)
   return Vect<T_>(VectorX(x)*a);
#else
   Vect<T_> v(x.size());
   for (size_t i=1; i<=x.size(); ++i)
      v.set(i,a*x(i));
   return v;
#endif
}


template<class T_>
Vect<T_> operator/(const Vect<T_>& x,
                   const T_&       a)
{
#if defined (USE_EIGEN)
   return Vect<T_>(Matrix<T_,Eigen::Dynamic,1>(x)/a);
#else
   Vect<T_> v(x);
   for (size_t i=1; i<=x.size(); ++i)
      v.set(i,x(i)/a);
   return v;
#endif
}


template<class T_>
T_ Dot(const Vect<T_>& x,
       const Vect<T_>& y)
{
#ifdef _OFELI_RANGE_CHECK
   assert(x.size() == y.size());
#endif
#if defined (USE_EIGEN)
   return Matrix<T_,Eigen::Dynamic,1>(x).dot(Matrix<T_,Eigen::Dynamic,1>(y));
#else
   T_ s=0;
   for (size_t i=0; i<x.size(); ++i)
      s += x[i]*y[i];
   return s;
#endif
}


inline real_t Discrepancy(Vect<real_t>&       x,
                          const Vect<real_t>& y,
                          int                 n,
                          int                 type)
{
   size_t s=x.size();
   real_t old=0., d=0.;
   if (n==0) {
      old = x.getNormMax();
      for (size_t i=0; i<s; i++) 
         if (d<std::abs(x[i]-y[i]))
            d = std::abs(x[i]-y[i]);
   }
   else if (n==1) {
      old = x.getWNorm1();
      for (size_t i=0; i<s; i++)
         d += std::abs(x[i]-y[i]);
      d /= s;
   }
   else if (n==2) {
      old = x.getWNorm2();
      for (size_t i=0; i<s; i++) {
         real_t z = std::abs(x[i]-y[i]);
         d += z*z;
      }
      d = std::sqrt(d/s);
   }
   else
      ;
   x = y;
   if (type==1 && old>0.)
      d /= old;
   return d;
}


inline real_t Discrepancy(Vect<complex_t>&       x,
                          const Vect<complex_t>& y,
                          int                    n,
                          int                    type)
{
   size_t s=x.size();
   real_t old=0., d=0.;
   if (n==0) {
      old = x.getNormMax();
      for (size_t i=0; i<s; i++) 
         if (d<std::abs(x[i]-y[i]))
            d = std::abs(x[i]-y[i]);
   }
   else if (n==1) {
      old = x.getWNorm1();
      for (size_t i=0; i<s; i++)
         d += std::abs(x[i]-y[i]);
      d /= s;
   }
   else if (n==2) {
      old = x.getWNorm2();
      for (size_t i=0; i<s; i++) {
         real_t z=std::abs(x[i]-y[i]);
         d += z*z;
      }
      d = std::sqrt(d/s);
   }
   else
      ;
   x = y;
   if (type==1 && old>0.)
      d /= old;
   return d;
}


inline void Modulus(const Vect<complex_t>& x,
                    Vect<real_t>&          y)
{
   y.setSize(x.size());
   for (size_t i=1; i<=x.size(); ++i) {
      complex_t z = x(i);
      y.set(i,std::sqrt(z.real()*z.real()+z.imag()*z.imag()));
   }
}


inline void Real(const Vect<complex_t>& x,
                 Vect<real_t>&          y)
{
   y.setSize(x.size());
   for (size_t i=1; i<=x.size(); ++i)
      y.set(i,x(i).real());
}


inline void Imag(const Vect<complex_t>& x,
                 Vect<real_t>&          y)
{
   y.setSize(x.size());
   for (size_t i=1; i<=x.size(); ++i)
      y.set(i,x(i).imag());
}


template<class T_>
istream& operator>>(istream&  s,
                    Vect<T_>& v)
{
   T_ z;
   for (size_t i=1; i<=v.size(); ++i) {
      s >> i >> z;
      v.set(i,z);
   }
   return s;
}


template<class T_>
ostream &operator<<(ostream&        s,
                    const Vect<T_>& v)
{
   int verb = Verbosity;
   Verbosity = 10;
   size_t nx=v.getNx(), ny=v.getNy(), nz=v.getNz();
   if (Verbosity==0) {
      s << "Vector size: " << v.size() << endl;
      return s;
   }
   if (Verbosity==1) {
      if (v.getName() != "#")
         s << v.getName() << " at time = " << v.getTime() << endl << endl;
      s << "Vector size: " << nx << "*" << ny << "*" << nz << endl;
      return s;
   }
   if (Verbosity==2)
      nx = !(10<nx)?nx:10, ny = !(10<ny)?ny:10, nz = !(10<nz)?nz:10;
   else if (Verbosity==3)
      nx = !(50<nx)?nx:50, ny = !(50<ny)?ny:50, nz = !(50<nz)?nz:50;
   else if (Verbosity==4)
      nx = !(100<nx)?nx:100, ny = !(100<ny)?ny:100, nz = !(100<nz)?nz:100;
   s.setf(ios::scientific);
   if ((nx<v.getNx() || ny<v.getNy() || nz<v.getNz()) && nx>0 && ny>0 && nz>0)
      cout << "Partial output of vector contents." << endl;
   if (nz>1) {
      for (size_t i=1; i<=nx; i++) {
         if (v.getNy()==1)
            s << setw(6) << i << "  " << setprecision(8) << setw(18) << v(i) << endl;
         else
            s << "\n[[ i = " << i << " ]]" << endl;
         for (size_t j=1; j<=ny; j++) {
            s << "\n[ j = " << j << " ]" << endl;
            for (size_t k=1; k<=nz; k++)
               s << setw(6) << setprecision(8) << setw(18) << v(i,j,k);
            s << endl;
         }
      }
   }
   else {
      s.setf(ios::scientific);
      for (size_t i=1; i<=nx; i++) {
         s << setw(6) << i << "   ";
         for (size_t j=1; j<=ny; j++)
            s << setprecision(8) << setw(18) << v(i,j);
         s << endl;
      }
      s << endl;
   }
   Verbosity = verb;
   return s;
}


template<class T_>
inline void OddEven(vector<T_>& x,
                    vector<T_>& odd,
                    vector<T_>& even)
{
   for (typename vector<T_>::iterator it=x.begin(); it!=x.end();) {
      even.push_back(*(it++));
      odd.push_back(*(it++));
   }
}


template<class T_>
inline void OddEven(Vect<T_>&   x,
                    vector<T_>& odd,
                    vector<T_>& even)
{
   for (typename Vect<T_>::iterator it=x.begin(); it!=x.end();) {
      even.push_back(*(it++));
      odd.push_back(*(it++));
   }
}


inline void fft(vector<complex_t>& x)
{
   size_t n = x.size();
   if (n <= 1)
      return;
 
// divide
   vector<complex_t> even, odd;
   OddEven(x,odd,even);
 
// conquer
   fft(even);
   fft(odd);
 
// combine
   for (size_t k=0; k<n/2; ++k) {
      complex_t t = std::polar(1.0,-2*OFELI_PI*k/n)*odd[k];
      x[k    ] = even[k] + t;
      x[k+n/2] = even[k] - t;
   }
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */

#endif
