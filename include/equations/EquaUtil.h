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
‡   along with OFELI. If not, see <http://www.gnu.org/licenses/>.

  ==============================================================================
       Template utilities functions for finite element equations classes
  ==============================================================================*/


#ifndef __FEEQUA_UTIL_H
#define __FEEQUA_UTIL_H

#include <float.h>
#include <stdlib.h>
#include <math.h>

#include "OFELI_Config.h"
#include "Mesh.h"
#include "Element.h"
#include "Side.h"
#include "Vect.h"
#include "BCVect.h"
#include "SpMatrix.h"
#include "SkMatrix.h"
#include "SkSMatrix.h"
#include "LocalMatrix.h"
#include "LocalVect.h"
#include "Material.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
void getMatrix(const SpMatrix<T_>& A)
{
   _M.Localize(_theElement,A);
}


void getMatrix(const SkMatrix<T_>& A)
{
   _M.Localize(_theElement,A);
}


void getMatrix(const SkSMatrix<T_>& A)
{
   _M.Localize(_theElement,A);
}


void getSolution(const Vect<T_>& u)
{
   _P.Localize(_theElement,u);
}


void getRHS(const Vect<T_>& b)
{
   _B.Localize(_theElement,b);
}

// Calculate residual in element
void Residual()
{
   _M.Mult(_P,_R); _R -= _B;
}


// Set time value
inline void setTime(real_t t)
{
   _time = t;
}


// Return nb of element nodes
size_t getNbNodes() const
{
   if (_theElement)
      return NEN_;
   else return NSN_;
}


// Return nb of element equations
size_t getNbEq() const
{
   if (_theElement)
      return NEE_;
   else
      return NSE_;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Return element matrix.
/// \details %Matrix is returned as a C-array
    T_ *A() { if (theElement_) return _M.Val(); else return _sM.Val(); }

/// \brief Return element right-hand side.
/// \details Right-hand side is returned as a C-array 
    T_ *b() { return _B.Val(); }

/// \brief Return element previous vector.
/// \details This is the vector given in time dependent constructor. It is returned as a C-array.
    T_ *Prev() { return _P.Val(); }


#ifndef DOXYGEN_SHOULD_SKIP_THIS
    T_ *Res() { return _R.Val(); }
    LocalMatrix<T_,NEE_,NEE_> &EA() { return _M; }
    LocalMatrix<T_,NSE_,NSE_> &SA() { return _sM; }
    LocalVect<T_,NEE_> &Eb() { return _B; }
    LocalVect<T_,NEE_> &Ep() { return _P; }

template < class T_, class E_ >
void ElementVect(const E_&       e,
                 const Vect<T_>& b,
                 Vect<T_>&       p,
                 int             dof,
                 int             dof_type)
{
   size_t k = 0;
   if (dof_type==NODE_DOF) {
     Node *nd;
     if (dof) {
       for (size_t n=1; n<=e.getNbNodes(); n++)
          p[k++] = b(e.getNodeLabel(n));
     }
     else {
       for (size_t n=1; n<=e.getNbNodes(); n++) {
          for (size_t j=1; j<=nd->getNbDOF(); j++)
             p[k++] = b(e(n)->getDOF(j));
       }
     }
   }
   else if (dof_type==SIDE_DOF) {
     Side *sd;
     if (dof) {
       for (size_t n=1; n<=e.getNbNodes(); n++)
          p[k++] = b(e.getSideLabel(n));
     }
     else {
       for (size_t n=1; n<=NEN_; n++) {
          sd = e.getPtrSide(n);
          for (size_t j=1; j<=sd->getNbDOF(); j++)
             p[k++] = b(sd->getDOF(j));
       }
     }
   }
}


template < class T_, class E_ >
void SideVect(const E_&       e,
              const Vect<T_>& b,
              Vect<T_>&       p)
{
   Node *nd;
   size_t k = 0;
   for (size_t n=1; n<=e.getNbNodes(); n++) {
      nd = e.getPtrNode(n);
      for (size_t j=1; j<=nd->getNbDOF(); j++)
         p[k++] = b(nd->getDOF(j));
   }
}


template <class T_, class E_ >
void updateBC(const E_&       e,
              const Vect<T_>& bc)
{
   size_t in = 0;
   for (size_t i=1; i<=e.getNbNodes(); i++) {
      for (size_t k=1; k<=_nb_dof; k++) {
         in++;
         size_t nn = (*_theElement)(i)->getFirstDOF() + k - 1;
         if ((*_theElement)(i)->getDOF(k) == 0) {
            size_t jn = 0;
            for (size_t j=1; j<=NEN_; j++) {
               for (size_t l=1; l<=nb_dof_; l++) {
                  jn++;
                  if ((*_theElement)(j)->getDOF(l) > 0)
                     eRHS(jn) -= eMat(jn,in) * bc(nn);
               }
            }
         }
      }
   }
}


template <class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void DiagBC(int dof_type,
            int dof)
{
   if (dof_type==NODE_DOF) {
      if (dof) {
         for (size_t i=1; i<=NEN_; i++) {
            Node *nd1 = (*_theElement)(i);
            if (nd1->getCode(dof)) {
               size_t in = nd1->Label();
               for (size_t j=1; j<=NEN_; j++) {
                  Node *nd2 = (*_theElement)(j);
                  if (nd2->getCode(dof)) {
                     size_t jn = nd2->n();
                     if (in != jn)
                        eMat(in,jn) = 0;
                  }
               }
            }
         }
      }
      else {
         for (size_t i=1; i<=NEN_; i++) {
            for (size_t k=1; k<=_nb_dof; k++) {
               if ((*_theElement)(i)->getCode(k)) {
                  size_t in = (*_theElement)(i)->getDOF(k);
                  for (size_t j=1; j<=NEN_; j++) {
                     for (size_t l=1; l<=_nb_dof; l++) {
                        size_t jn = (*_theElement)(j)->getDOF(l);
                        if (in != jn)
                           eMat(in,jn) = 0;
                     }
                  }
               }
            }
         }
      }
   }

   else if (dof_type==SIDE_DOF) {
      if (dof) {
         for (size_t i=1; i<=NEN_; i++) {
            Side *sd1 = gettheElement->getPtrSide(i);
            if (sd1->Code(dof))
               for (size_t j=1; j<=NEN_; j++)
                  if (i != j)
                     eMat(i,j) = 0;
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
                           eMat(in,jn) = 0;
                     }
                  }
               }
            }
         }
      }
   }
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
