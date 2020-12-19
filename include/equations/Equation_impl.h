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

                    Implementation of abstract class 'Equation'

  ==============================================================================*/


#ifndef __EQUATION_IMPL_H
#define __EQUATION_IMPL_H

#include "equations/Equation.h"
#include "mesh/Mesh.h"
#include "mesh/Side.h"
#include "mesh/Node.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/LocalVect_impl.h"
#include "linear_algebra/Matrix_impl.h"
#include "linear_algebra/SkMatrix_impl.h"
#include "linear_algebra/SkSMatrix_impl.h"
#include "linear_algebra/SpMatrix_impl.h"
#include "linear_algebra/BMatrix_impl.h"
#include "linear_algebra/TrMatrix_impl.h"
#include "OFELIException.h"

namespace OFELI {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
Equation<T_,NEN_,NEE_,NSN_,NSE_>::Equation()
                                 : AbsEqua<T_>()
{
   _nb_dof = NEE_/NEN_;
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
Equation<T_,NEN_,NEE_,NSN_,NSE_>::Equation(Mesh &mesh)
                                 : AbsEqua<T_>()
{
   if (NEE_/NEN_ != NSE_/NSN_)
      throw OFELIException("In Equation<>::Equation(Mesh): Numbers of element and "
                           "side nodes are inconsistent");
   _nb_dof = NEE_/NEN_;
   initEquation(mesh,0.,1.,0.1);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
Equation<T_,NEN_,NEE_,NSN_,NSE_>::Equation(Mesh&     mesh,
                                           Vect<T_>& u)
                                 : AbsEqua<T_>()
{
   if (NEE_/NEN_ != NSE_/NSN_)
      throw OFELIException("In Equation<>::Equation(Mesh&,Vect<>&): Numbers of "
                           "element and side nodes are inconsistent");
   _nb_dof = NEE_/NEN_;
   initEquation(mesh,0.,1.,0.1);
   this->setInput(SOLUTION,u);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
Equation<T_,NEN_,NEE_,NSN_,NSE_>::Equation(Mesh&     mesh,
                                           Vect<T_>& u,
                                           real_t&   init_time,
                                           real_t&   final_time,
                                           real_t&   time_step)
                                 : AbsEqua<T_>()
{
   if (NEE_/NEN_ != NSE_/NSN_)
      throw OFELIException("In Equation<>::Equation(Mesh&,...): Numbers of "
                           "element and side nodes are inconsistent");
   _nb_dof = NEE_/NEN_;
   initEquation(mesh,init_time,final_time,time_step);
   this->setInput(INITIAL_FIELD,u);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
Equation<T_,NEN_,NEE_,NSN_,NSE_>::~Equation()
{ }


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::initEquation(Mesh&  mesh,
                                                    real_t init_time,
                                                    real_t final_time,
                                                    real_t time_step)
{
   AbsEqua<T_>::setMesh(mesh);
   _TimeInt.delta = time_step;
   _TimeInt.time = _TimeInt.init = init_time;
   _TimeInt.final = final_time;
   AbsEqua<T_>::_b = new Vect<T_>(_nb_eq);
   AbsEqua<T_>::_uu.setSize(_nb_eq);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::getMatrix(const SpMatrix<T_>& A) const
{
   eMat.Localize(_theElement,A);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::getMatrix(const SkMatrix<T_>& A) const
{
   eMat.Localize(_theElement,A);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::getMatrix(const SkSMatrix<T_>& A) const
{
   eMat.Localize(_theElement,A);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::getSolution(const Vect<T_>& u) const
{
   _eu.Localize(_theElement,u);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::getRHS(const Vect<T_>& b) const
{
   eRHS.Localize(_theElement,b);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::updateBC(const Element&  el,
                                                const Vect<T_>& bc)
{
   size_t in=1;
   for (size_t i=1; i<=NEN_; ++i) {
      for (size_t k=1; k<=_nb_dof; ++k, ++in) {
         size_t jn=1;
         for (size_t j=1; j<=NEN_; ++j) {
            for (size_t l=1; l<=_nb_dof; ++l, ++jn) {
               if (el(j)->getCode(l) > 0)
                  eRHS(in) -= eMat(in,jn) * bc(el(j)->n(),l);
            }
         }
      }
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::updateBC(const Element& el,
                                                Vect<T_>*      bc)
{
   if (bc==nullptr)
      return;
   size_t in=1;
   for (size_t i=1; i<=NEN_; ++i) {
      for (size_t k=1; k<=_nb_dof; ++k, ++in) {
         size_t jn=1;
         for (size_t j=1; j<=NEN_; ++j) {
            for (size_t l=1; l<=_nb_dof; ++l, ++jn) {
               if (el(j)->getCode(l) > 0)
                  eRHS(in) -= eMat(in,jn) * (*bc)(el(j)->n(),l);
            }
         }
      }
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::DiagBC(DOFSupport dof_type,
                                              int        dof)
{
   size_t in, jn;
   if (dof_type==NODE_DOF) {
      if (dof) {
         for (size_t i=1; i<=NEN_; i++) {
            Node *nd1 = (*_theElement)(i);
            if (nd1->getCode(dof)) {
               for (size_t j=1; j<=NEN_; j++) {
                  Node *nd2 = (*_theElement)(j);
                  in = nd1->n(), jn = nd2->n();
                  if (in != jn)
                     eMat(in,jn) = T_(0.);
               }
            }
         }
      }
      else {
         size_t in=1;
         for (size_t i=1; i<=NEN_; i++) {
            for (size_t k=1; k<=_nb_dof; k++) {
               size_t jn = 1;
               for (size_t j=1; j<=NEN_; j++) {
                  for (size_t l=1; l<=_nb_dof; l++) {
                     if (in != jn && (*_theElement)(i)->getCode(k))
                        eMat(in,jn) = T_(0.);
                     jn++;
                  }
               }
               in++;
            }
         }
      }
   }

   else if (dof_type==SIDE_DOF) {
      if (dof) {
         for (size_t i=1; i<=NEN_; i++) {
            Side *sd1 = _theElement->getPtrSide(i);
            if (sd1->getCode(dof))
               for (size_t j=1; j<=NEN_; j++)
                  if (i != j)
                     eMat(i,j) = T_(0.);
         }
      }
      else {
         for (size_t i=1; i<=NEN_; i++) {
            for (size_t k=1; k<=_nb_dof; k++) {
               if (_theElement->getPtrSide(i)->getCode(k)) {
                  size_t in = (*_theElement)(i)->getDOF(k);
                  for (size_t j=1; j<=NEN_; j++) {
                     for (size_t l=1; l<=_nb_dof; l++) {
                        size_t jn = _theElement->getPtrSide(j)->getDOF(l);
                        if (in != jn)
                           eMat(in,jn) = T_(0.);
                     }
                  }
               }
            }
         }
      }
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::setBodyForce(const Vect<T_>& f)
{
   AbsEqua<T_>::_bf = &f;
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::setResidue()
{
   eMat.Mult(_eu,eRes);
   eRes -= eRHS;
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::LocalNodeVector(Vect<T_>& b)
{
   size_t k = 0;
   for (size_t n=1; n<=NEN_; ++n) {
      size_t nd = (*_theElement)(n)->n();
      for (size_t j=1; j<=_nb_dof; j++)
         _eu[k++] = b(nd,j);
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementNodeVector(const Vect<T_>&     b,
                                                         LocalVect<T_,NEE_>& be)
{
   size_t k = 0;
   for (size_t n=1; n<=NEN_; ++n) {
      size_t nd = (*_theElement)(n)->n();
      for (size_t i=0; i<_nb_dof; ++i)
         be[k++] = b[_nb_dof*(nd-1)+i];
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideNodeVector(const Vect<T_>&     b,
                                                      LocalVect<T_,NSE_>& bs)
{
   size_t k = 0;
   for (size_t n=1; n<=NSN_; ++n) {
      size_t nd = (*_theSide)(n)->n();
      for (size_t i=0; i<_nb_dof; ++i)
         bs[k++] = b[_nb_dof*(nd-1)+i];
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementNodeVectorSingleDOF(const Vect<T_>&     b,
                                                                  LocalVect<T_,NEN_>& be)
{
   for (size_t n=1; n<=NEN_; ++n)
      be(n) = b((*_theElement)(n)->n());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementNodeVector(const Vect<T_>&     b,
                                                         LocalVect<T_,NEN_>& be,
                                                         int                 dof)
{
   size_t k = 0;
   for (size_t n=1; n<=NEN_; ++n) {
      size_t nd = (*_theElement)(n)->n();
      for (size_t i=1; i<=_nb_dof; ++i)
         if (i==dof)
            be[k++] = b(_nb_dof*(nd-1)+i);
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementSideVector(const Vect<T_>&     b,
                                                         LocalVect<T_,NSE_>& be)
{
   size_t k = 0;
   Side *sd;
   for (size_t n=1; n<=NSE_; ++n) {
      sd = _theElement->getPtrSide(n);
      for (size_t j=1; j<=sd->getNbDOF(); j++)
         if (sd->getDOF(j))
            be[k++] = b(sd->getDOF(j));
   }
}

/*
template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideNodeVector(const Vect<T_>& b,
                                                      T_*             sb)
{
   size_t k = 0;
   for (size_t n=1; n<=NSN_; ++n) {
      size_t nd = (*_theSide)(n)->n();
      for (size_t i=0; i<_nb_dof; ++i)
         sb[k++] = b[_nb_dof*(nd-1)+i];
   }
   }*/


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementVector(const Vect<T_>& b,
                                                     DOFSupport      dof_type,
                                                     int             flag)
{
   size_t k=0;
   switch (dof_type) {

      case NODE_DOF:
         if (flag==0) {
            for (size_t n=1; n<=NEN_; ++n) {
               Node *nd = (*_theElement)(n);
               for (size_t j=1; j<=nd->getNbDOF(); j++)
                  if (nd->getDOF(j))
                     _eu[k++] = b(nd->getDOF(j));
            }
         }
         else {
            for (size_t n=1; n<=NEN_; ++n) {
               Node *nd = (*_theElement)(n);
               if (nd->getDOF(flag))
                  _eu[k++] = b(nd->getDOF(flag));
            }
         }
         break;

      case SIDE_DOF:
         if (flag==0) {
            for (size_t n=1; n<=NEN_; n++) {
               Side *sd=_theElement->getPtrSide(n);
               for (size_t j=1; j<=sd->getNbDOF(); j++)
                  if (sd->getDOF(j))
                     _eu[k++] = b(sd->getDOF(j));
            }
         }
         else {
            for (size_t n=1; n<=NEN_; n++) {
               Side *sd=_theElement->getPtrSide(n);
               if (sd->getDOF(flag))
                  _su[k++] = b(sd->getDOF(flag));
            }
         }

      default:
         break;
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideVector(const Vect<T_>& b,
                                                  T_*             sb)
{
   for (size_t i=0; i<_theSide->getNbDOF(); ++i)
      sb[i] = b(_theSide->getNbDOF()*(_theSide->n()-1)+i+1);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementNodeCoordinates()
{
   for (size_t n=1; n<=NEN_; n++)
      _x[n-1] = (*_theElement)(n)->getCoord();
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideNodeCoordinates()
{
   for (size_t n=1; n<=NSN_; n++)
      _x[n-1] = (*_theSide)(n)->getCoord();
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementAssembly(Matrix<T_>* A)
{
   A->Assembly(TheElement,eMat.get());
}


#if defined (USE_PETSC)
template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementAssembly(PETScMatrix<T_>& A)
{
   ElementAssembly(TheElement,A);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideAssembly(PETScMatrix<T_>& A)
{
   SideAssembly(TheSide,A);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementAssembly(PETScVect<T_>& b)
{
   ElementAssembly(TheElement,b);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideAssembly(PETScVect<T_>& b)
{
   SideAssembly(TheSide,b);
}
#endif


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementAssembly(Element&    el,
                                                       Matrix<T_>* A)
{
   A->Assembly(el,eMat.get());
}

template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementAssembly(Element&      el,
                                                       SkMatrix<T_>& A)
{
   A.Assembly(el,eMat.get());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementAssembly(Element&       el,
                                                       SkSMatrix<T_>& A)
{
   A.Assembly(el,eMat.get());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementAssembly(Element&      el,
                                                       SpMatrix<T_>& A)
{
   A.Assembly(el,eMat.get());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementAssembly(Element&     el,
                                                       BMatrix<T_>& A)
{
   A.Assembly(el,eMat.get());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementAssembly(Element&      el,
                                                       TrMatrix<T_>& A)
{
   A.Assembly(el,eMat.get());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideAssembly(Side&       sd,
                                                    Matrix<T_>* A)
{
   A->Assembly(sd,sMat.get());
}

template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideAssembly(Side&        sd,
                                                    BMatrix<T_>& A)
{
   A.Assembly(sd,sMat.get());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideAssembly(Side&         sd,
                                                    SkMatrix<T_>& A)
{
   A.Assembly(sd,sMat.get());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideAssembly(Side&          sd,
                                                    SkSMatrix<T_>& A)
{
   A.Assembly(sd,sMat.get());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideAssembly(Side&         sd,
                                                    SpMatrix<T_>& A)
{
   A.Assembly(sd,sMat.get());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementAssembly(Element&  el,
                                                       Vect<T_>& v)
{
   v.Assembly(el,eRHS.get());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideAssembly(Side&     sd,
                                                    Vect<T_>& v)
{
   v.Assembly(sd,sRHS.get());
}


#if defined(USE_PETSC)
template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::updateBC(const Element&       el,
                                                const PETScVect<T_>& bc)
{
   size_t in = 0;
   for (size_t i=1; i<=NEN_; ++i) {
      for (size_t k=1; k<=_nb_dof; ++k) {
         in++;
         size_t nn = el(i)->getFirstDOF() + k - 1;
         if (el(i)->getDOF(k) == 0) {
            size_t jn = 0;
            for (size_t j=1; j<=NEN_; ++j) {
               for (size_t l=1; l<=_nb_dof; ++l) {
                  jn++;
                  if (el(j)->getDOF(l) > 0)
                     eRHS(jn) -= eMat(jn,in) * bc(nn);
               }
            }
         }
      }
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::updateBC(const PETScVect<T_>& bc)
{
   updateBC(TheElement,bc);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementAssembly(Element&         el,
                                                       PETScMatrix<T_>& A)
{
   A.Assembly(el,eMat.get());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideAssembly(Side&            sd,
                                                    PETScMatrix<T_>& A)
{
   A.Assembly(sd,sMat.get());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementAssembly(Element&       el,
                                                       PETScVect<T_>& v)
{
   v.Assembly(el,eRHS.get());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideAssembly(Side&          sd,
                                                    PETScVect<T_>& v)
{
   v.Assembly(sd,sRHS.get());
}
#endif


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementAssembly(BMatrix<T_>& A)
{
   ElementAssembly(TheElement,A);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementAssembly(SkSMatrix<T_>& A)
{
   ElementAssembly(TheElement,A);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementAssembly(SkMatrix<T_>& A)
{
   ElementAssembly(TheElement,A);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementAssembly(SpMatrix<T_>& A)
{
   ElementAssembly(TheElement,A);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementAssembly(TrMatrix<T_>& A)
{
   ElementAssembly(TheElement,A);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::DGElementAssembly(Matrix<T_>* A)
{
   A->DGAssembly(TheElement,eMat.get());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::DGElementAssembly(SkSMatrix<T_>& A)
{
   A.DGAssembly(TheElement,eMat.get());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::DGElementAssembly(SkMatrix<T_>& A)
{
   A.DGAssembly(TheElement,eMat.get());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::DGElementAssembly(SpMatrix<T_>& A)
{
   A.DGAssembly(TheElement,eMat.get());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::DGElementAssembly(TrMatrix<T_>& A)
{
   A.DGAssembly(TheElement,eMat.get());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideAssembly(Matrix<T_>* A)
{
   SideAssembly(TheSide,A);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideAssembly(SkSMatrix<T_>& A)
{
   SideAssembly(TheSide,A);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideAssembly(SkMatrix<T_>& A)
{
   SideAssembly(TheSide,A);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideAssembly(SpMatrix<T_>& A)
{
   SideAssembly(TheSide,A);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementAssembly(Vect<T_>& v)
{
   ElementAssembly(TheElement,v);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideAssembly(Vect<T_>& v)
{
   SideAssembly(TheSide,v);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::AxbAssembly(const Element&  el,
                                                   const Vect<T_>& x,
                                                   Vect<T_>&       b)
{
   size_t ii=0, jj=0, ik, jl;
   for (size_t i=1; i<=NEN_; ++i) {
      Node *ndi = el(i);
      for (size_t k=1; k<=_nb_dof; ++k) {
         ii++, ik = ndi->getDOF(k);
         for (size_t j=1; j<=NEN_; ++j) {
            Node *ndj = el(j);
            for (size_t l=1; l<=_nb_dof; ++l) {
               jj++, jl = ndj->getDOF(l);
               if (ik)
                  b(ik) += eMat(ii,jj)*x(jl);
            }
         }
      }
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::AxbAssembly(const Side&     sd,
                                                   const Vect<T_>& x,
                                                   Vect<T_>&       b)
{
   size_t ii=0, jj=0, ik, jl;
   for (size_t i=1; i<=NSN_; ++i) {
      Node *ndi = sd(i);
      for (size_t k=1; k<=_nb_dof; ++k, ii++) {
         ik = ndi->getDOF(k);
         for (size_t j=1; j<=NSN_; ++j) {
            Node *ndj = sd(j);
            for (size_t l=1; l<=_nb_dof; ++l, jj++) {
               jl = ndj->getDOF(l);
               if (ik)
                  b(ik) += sMat(ii,jj)*x(jl);
            }
         }
      }
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
size_t Equation<T_,NEN_,NEE_,NSN_,NSE_>::getNbNodes() const
{
   if (_theElement)
      return NEN_;
   else
      return NSN_;
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
size_t Equation<T_,NEN_,NEE_,NSN_,NSE_>::getNbEq() const
{
   if (_theElement)
      return NEE_;
   else
      return NSE_;
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
real_t Equation<T_,NEN_,NEE_,NSN_,NSE_>::setMaterialProperty(const string& exp,
                                                             const string& prop)
{
   exprtk::expression<double> expression;
   exprtk::symbol_table<double> symbol_table;
   add_constants(symbol_table);
   symbol_table.add_variable("x",_x.x);
   symbol_table.add_variable("y",_x.y);
   symbol_table.add_variable("z",_x.z);
   expression.register_symbol_table(symbol_table);
   theParser.compile(exp,expression);
   return expression.value();
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */

#endif
