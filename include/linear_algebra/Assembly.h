/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2019 Rachid Touzani

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

                       Functions for Finite Element Assembly

  ==============================================================================*/

/*! \file Assembly.h
 *  \brief A set of template functions for assembly purposes.
 *
 * \tparam <E_> Entity type to assemble (Element or Side classes)
 * \tparam <T_> Data type (double, float, complex<double>, ...)
 * \tparam <N_> Constant size of local vector
 */

#ifndef __ASSEMBLY_H
#define __ASSEMBLY_H

#include "OFELI_Config.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/LocalVect.h"
#include "linear_algebra/LocalMatrix.h"
#include "linear_algebra/SkMatrix.h"
#include "linear_algebra/SkSMatrix.h"
#include "linear_algebra/SpMatrix.h"


namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/** \addtogroup VectMat
 *  @{
 */

/** \fn void element_assembly(const E_ &e, const LocalVect<T_,N_> &be, Vect<T_> &b)
 *  \ingroup Equation
 *  \brief Assemble local vector into global vector
 *
 *  @param [in] e Reference to local entity (Element or Side)
 *  @param [in] be Local vector
 *  @param [in,out] b Global vector
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_,class E_>
void element_assembly(const E_&               e,
                      const LocalVect<T_,N_>& be,
                      Vect<T_>&               b)
{
   size_t i=1;
   for (size_t n=1; n<=e.getNbNodes(); n++) {
      Node *nd=e.getPtrNode(n);
      for (size_t k=1; k<=nd->getNbDOF(); k++)
         b.add(nd->getDOF(k),be(i++));
   }
}


/** \fn void element_assembly(const E_ &e, const LocalMatrix<T_,N_,N_> &ae, Vect<T_> &b)
 *  \ingroup Equation
 *  \brief Assemble diagonal local vector into global vector
 *  @param [in] e Reference to local entity (Element or Side)
 *  @param [in] ae Local matrix
 *  @param [in,out] b Global vector
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_,class E_>
void element_assembly(const E_&                    e,
                      const LocalMatrix<T_,N_,N_>& ae,
                      Vect<T_>&                    b)
{
   size_t i=0;
   for (size_t n=1; n<=e.getNbNodes(); n++) {
      Node *nd = e.getPtrNode(n);
      for (size_t k=1; k<=nd->getNbDOF(); k++, i++)
         b.add(nd->getDOF(k),ae(i,i));
   }
}


/** \fn void element_assembly(const E_ &e, const LocalMatrix<T_,N_,N_> &ae, Matrix<T_> *A)
 *  \ingroup Equation
 *  \brief Assemble local matrix into global matrix.
 *  \details This function is to be called with an abstract pointer to matrix (class Matrix)
 *  @param [in] e Reference to local entity (Element or Side)
 *  @param [in] ae Local matrix
 *  @param [in,out] A Pointer to global matrix
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_,class E_>
void element_assembly(const E_&                    e,
                      const LocalMatrix<T_,N_,N_>& ae,
                      Matrix<T_>*                  A)
{
   size_t i=1;
   for (size_t in=1; in<=e.getNbNodes(); in++) {
      Node *nd1=e(in);
      for (size_t k=1; k<=nd1->getNbDOF(); k++) {
         size_t ii=nd1->getDOF(k);
         size_t j=1;
         for (size_t jn=1; jn<=e.getNbNodes(); jn++) {
            Node *nd2=e(jn);
            for (size_t l=1; l<=nd2->getNbDOF(); l++, j++)
               A->add(ii,nd2->getDOF(l),ae(i,j));
         }
         i++;
      }
   }
}


/** \fn void element_assembly(const E_ &e, const LocalMatrix<T_,N_,N_> &ae, SkMatrix<T_> &A)
 *  \ingroup Equation
 *  \brief Assemble local matrix into global skyline matrix
 *  @param [in] e Reference to local entity (Element or Side)
 *  @param [in] ae Local matrix
 *  @param [in,out] A Global matrix
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_,class E_>
void element_assembly(const E_&                    e,
                      const LocalMatrix<T_,N_,N_>& ae,
                      SkMatrix<T_>&                A)
{
   size_t i = 1;
   for (size_t in=1; in<=e.getNbNodes(); in++) {
      Node *nd1 = e.getPtrNode(in);
      for (size_t k=1; k<=nd1->getNbDOF(); k++) {
         size_t ii = nd1->getDOF(k);
         size_t j = 1;
         for (size_t jn=1; jn<=e.getNbNodes(); jn++) {
            Node *nd2 = e(jn);
            for (size_t l=1; l<=nd2->getNbDOF(); l++, j++)
               A.add(ii,nd2->getDOF(l),ae(i,j));
         }
         i++;
      }
   }
}


/** \fn void element_assembly(const E_ &e, const LocalMatrix<T_,N_,N_> &ae, SkSMatrix<T_> &A)
 *  \ingroup Equation
 *  \brief Assemble local matrix into global symmetric skyline matrix
 *  @param [in] e Reference to local entity (Element or Side)
 *  @param [in] ae Local matrix
 *  @param [in,out] A Global matrix
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_,class E_>
void element_assembly(const E_&                    e,
                      const LocalMatrix<T_,N_,N_>& ae,
                      SkSMatrix<T_>&               A)
{
   size_t i=1, ii, jj;
   for (size_t in=1; in<=e.getNbNodes(); in++) {
      Node *nd1 = e(in);
      for (size_t k=1; k<=nd1->getNbDOF(); k++) {
         if ((ii=nd1->getDOF(k)) > 0) {
            size_t j = 1;
            for (size_t jn=1; jn<=e.getNbNodes(); jn++) {
               Node *nd2 = e(jn);
               for (size_t l=1; l<=nd2->getNbDOF(); l++, j++) {
                  if ((jj=nd2->getDOF(l)) > 0)
                     A.add(ii,jj,ae(i,j));
               }
            }
         }
         i++;
      }
   }
}


/** \fn void element_assembly(const E_ &e, const LocalMatrix<T_,N_,N_> &ae, SpMatrix<T_> &A)
 *  \ingroup Equation
 *  \brief Assemble local matrix into global sparse matrix
 *  @param [in] e Reference to local entity (Element or Side)
 *  @param [in] ae Local matrix
 *  @param [in,out] A Global matrix
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_,class E_>
void element_assembly(const E_&                    e,
                      const LocalMatrix<T_,N_,N_>& ae,
                      SpMatrix<T_>&                A)
{
   size_t i = 1;
   for (size_t in=1; in<=e.getNbNodes(); in++) {
      Node *nd1 = e(in);
      for (size_t k=1; k<=nd1->getNbDOF(); k++) {
         size_t ii=nd1->getDOF(k), j=1;
         for (size_t jn=1; jn<=e.getNbNodes(); jn++) {
            Node *nd2 = e(jn);
            for (size_t l=1; l<=nd2->getNbDOF(); l++, j++)
               A.add(ii,nd2->getDOF(l),ae(i,j));
         }
         i++;
      }
   }
}


#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** \fn void side_assembly(const E_ &e, const LocalVect<T_,N_> &be, Vect<T_> &b)
 *  \ingroup Equation
 *  \brief Assemble local vector into global vector by selecting one DOF.
 *  @param [in] e Reference to local entity (Element or Side)
 *  @param [in] be Local vector
 *  @param [in,out] b Global vector
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_,class E_>
void side_assembly(const E_&               e,
                   const LocalVect<T_,N_>& be,
                   Vect<T_>&               b)
{
   size_t i = 1;
   for (size_t n=1; n<=N_; n++)
      b.add(e(n)->n(),be(i++));
}


/** \fn void side_assembly(const E_ &e, const LocalMatrix<T_,N_,N_> &ae, Vect<T_> &b)
 *  \ingroup Equation
 *  \brief Assemble local diagonal matrix into global vector by selecting one DOF.
 *  @param [in] e Reference to local entity (Element or Side)
 *  @param [in] ae Local matrix 
 *  @param [in,out] b Global vector
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_,class E_>
void side_assembly(const E_&                    e,
                   const LocalMatrix<T_,N_,N_>& ae,
                   Vect<T_>&                    b)
{
   for (size_t n=1; n<=N_; n++)
      b.add(e(n)->n(),ae(n,n));
}


/** \fn void side_assembly(const E_ &e, const LocalMatrix<T_,N_,N_> &ae, SkSMatrix<T_> &A)
 *  \ingroup Equation
 *  \brief Assemble local matrix into global matrix by selecting one DOF.
 *  @param [in] e Reference to local entity (Element or Side)
 *  @param [in] ae Local vector
 *  @param [in,out] A Global vector
 */
template<class T_,size_t N_,class E_>
void side_assembly(const E_&                    e,
                   const LocalMatrix<T_,N_,N_>& ae,
                   SkSMatrix<T_>&               A)
{
   size_t i = 1;
   for (size_t in=1; in<=N_; in++) {
      size_t ii=e(in)->n(), j=1;
      for (size_t jn=1; jn<=N_; jn++, j++) {
         size_t jj = e(jn)->n();
         if (ii >= jj)
            A.add(ii,jj,ae(i,j));
      }
      i++;
   }
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/** \fn void side_assembly(const Element &e, const LocalMatrix<T_,N_,N_> &ae, SpMatrix<T_> &A)
 *  \brief %Side assembly of local matrix into global matrix (as instance of class SpMatrix).
 *  \ingroup Equation
 *  @param [in] e Reference to local Element
 *  @param [in] ae Local matrix
 *  @param [in,out] A Global matrix
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_>
void side_assembly(const Element&               e,
                   const LocalMatrix<T_,N_,N_>& ae,
                   SpMatrix<T_>&                A)
{
   size_t i = 1;
   for (size_t n=1; n<=e.getNbSides(); n++) {
      Side *s1 = e.getPtrSide(n);
      for (size_t k=1; k<=s1->getNbDOF(); k++) {
         if (s1->getDOF(k)!=0) {
            size_t j = 1;
            for (size_t m=1; m<=e.getNbSides(); m++) {
               Side *s2 = e.getPtrSide(m);
               for (size_t l=1; l<=s2->getNbDOF(); l++, j++) {
                  if (s2->getDOF(l)!=0)
                     A.add(s1->getDOF(k),s2->getDOF(l),ae(i,j));
               }
            }
         }
         i++;
      }
   }
}


#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_,size_t N_>
void side_assembly(const Element&               e,
                   const LocalMatrix<T_,N_,N_>& ae,
                   Matrix<T_>*                  A)
{
   size_t i = 1;
   for (size_t n=1; n<=e.getNbSides(); n++) {
      Side *s1 = e.getPtrSide(n);
      for (size_t k=1; k<=s1->getNbDOF(); k++) {
         if (s1->getDOF(k)!=0) {
            size_t j = 1;
            for (size_t m=1; m<=e.getNbSides(); m++) {
               Side *s2 = e.getPtrSide(m);
               for (size_t l=1; l<=s2->getNbDOF(); l++, j++) {
                  if (s2->getDOF(l)!=0)
                     A->add(s1->getDOF(k),s2->getDOF(l),ae(i,j));
               }
            }
         }
         i++;
      }
   }
}


template<class T_, size_t N_>
void Assembly(const Element&                   el,
              const LocalMatrix<real_t,N_,N_>& Ae,
              SpMatrix<T_>&                    A)
{
   for (size_t i=1; i<=el.getNbNodes(); ++i) {
      Node *nd1 = el(i);
      for (size_t k=1; k<=nd1->getNbDOF(); ++k) {
         size_t ik = nd1->getNbDOF()*(i-1) + k;
         for (size_t j=1; j<=el.getNbNodes(); ++j) {
            Node *nd2 = el(j);
            for (size_t l=1; l<=nd2->getNbDOF(); ++l) {
               size_t jl = nd2->getNbDOF()*(j-1) + l;
               if (nd1->getDOF(k)!=0 && nd2->getDOF(l)!=0)
                  A.add(nd1->getDOF(k),nd2->getDOF(l),Ae(ik,jl));
            }
         }
      }
   }
}


template<class T_, size_t N_>
void Assembly(const Element&              el,
              const LocalVect<real_t,N_>& be,
              Vect<T_>&                   b)
{
   for (size_t i=1; i<=el.getNbNodes(); ++i) {
      Node *nd = el(i);
      for (size_t k=1; k<=nd->getNbDOF(); ++k) {
         size_t ik = nd->getNbDOF()*(i-1) + k;
         if (nd->getDOF(k)!=0)
            b.add(nd->getDOF(k),be(ik));
      }
   }
}

#endif  /* DOXYGEN_SHOULD_SKIP_THIS */


/** \fn void side_assembly(const Element &e, const LocalMatrix<T_,N_,N_> &ae, SkSMatrix<T_> &A)
 *  \ingroup Equation
 *  \brief %Side assembly of local matrix into global matrix (as instance of class SkSMatrix).
 *  @param [in] e Reference to local Element
 *  @param [in] ae Local matrix
 *  @param [in,out] A Global matrix
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_>
void side_assembly(const Element&               e,
                   const LocalMatrix<T_,N_,N_>& ae,
                   SkSMatrix<T_>&               A)
{
   size_t i=1, j, ii, jj;
   for (size_t n=1; n<=e.getNbSides(); n++) {
      Side *s1 = e.getPtrSide(n);
      for (size_t k=1; k<=s1->getNbDOF(); k++) {
         if ((ii=s1->getDOF(k)) > 0) {
            j = 1;
            for (size_t m=1; m<=e.getNbSides(); m++) {
               Side *s2 = e.getPtrSide(m);
               for (size_t l=1; l<=s2->getNbDOF(); l++, j++) {
                  if ((jj=s2->getDOF(l)) > 0)
                     A.add(ii,jj,ae(i,j));
               }
            }
         }
         i++;
      }
   }
}


/** \fn void side_assembly(const Element &e, const LocalMatrix<T_,N_,N_> &ae, SkMatrix<T_> &A)
 *  \ingroup Equation
 *  \brief Side assembly of local matrix into global matrix (as instance of class SkMatrix).
 *  @param [in] e Reference to local Element
 *  @param [in] ae Local matrix
 *  @param [in,out] A Global matrix
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_>
void side_assembly(const Element&               e,
                   const LocalMatrix<T_,N_,N_>& ae,
                   SkMatrix<T_>&                A)
{
   size_t i = 1;
   for (size_t n=1; n<=e.getNbSides(); n++) {
      Side *s1 = e.getPtrSide(n);
      for (size_t k=1; k<=s1->getNbDOF(); k++) {
         size_t ii = s1->getDOF(k);
         if (ii != 0) {
            size_t j = 1;
            for (size_t m=1; m<=e.getNbSides(); m++) {
               Side *s2 = e.getPtrSide(m);
               for (size_t l=1; l<=s2->getNbDOF(); l++, j++) {
                  size_t jj = s2->getDOF(l);
                  if (jj)
                     A.add(ii,jj,ae(i,j));
               }
            }
         }
         i++;
      }
   }
}


/** \fn void side_assembly(const Element &e, const LocalVect<T_,N_> &be, Vect<T_> &b)
 *  \ingroup Equation
 *  \brief %Side assembly of local vector into global vector.
 *  @param [in] e Reference to local Element
 *  @param [in] be Local vector
 *  @param [in,out] b Global vector
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_,size_t N_>
void side_assembly(const Element&          e,
                   const LocalVect<T_,N_>& be,
                   Vect<T_>&               b)
{
   size_t i = 1;
   for (size_t n=1; n<=e.getNbSides(); n++) {
      Side *sd = e.getPtrSide(n);
      for (size_t k=1; k<=sd->getNbDOF(); k++, i++) {
         if (sd->getDOF(k) != 0)
            b.add(sd->getDOF(k),be(i));
      }
   }
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class E_,class T_>
void element_assembly(const E_& e,
                      const T_* ae,
                      const T_* be,
                      Vect<T_>& A,
                      Vect<T_>& b)
{
   size_t i=0, j=0, nee=e.getNbNodes()*e.getNbDOF();
   for (size_t in=1; in<=e.getNbNodes(); in++) {
      Node *nd=e(in);
      for (size_t k=1; k<=e.getNbDOF(); k++, i++, j+=nee+1) {
         size_t ii=nd->getDOF(k);
         if (ii) {
            b.add(ii,be[i]);
            A.add(ii,ae[j]);
         }
      }
   }
}


template<class E_>
void Element_Assembly(const E_&     e,
                      const real_t* be,
                      Vect<real_t>& b)
{
   size_t i=0;
   for (size_t in=1; in<=e.getNbNodes(); in++) {
      Node *nd=e(in);
      for (size_t k=1; k<=e.getNbDOF(); k++, i++)
         b.add(nd->getDOF(k)*(nd->n())+k,be[i]);
   }
}


template<class E_,class T_>
void element_assembly(const E_& e,
                      const T_* be,
                      Vect<T_>& b)
{
   size_t i=0;
   for (size_t in=1; in<=e.getNbNodes(); in++) {
      Node *nd=e(in);
      for (size_t k=1; k<=e.getNbDOF(); k++, i++) {
         if (nd->getDOF(k))
            b.add(nd->getDOF(k),be[i]);
      }
   }
}


template<class E_,class T_>
void element_assembly(const E_&   e,
                      const T_*   ae,
                      const T_*   be,
                      Matrix<T_>* A,
                      Vect<T_>&   b)
{
   size_t i=0, j=0, jj=0;
   size_t nen=e.getNbNodes(), nb_dof=e.getNbDOF();
   for (size_t in=1; in<=nen; in++) {
      for (size_t k=1; k<=nb_dof; k++, j++) {
         size_t ii = e(in)->getDOF(k);
         if (ii>0)
            b.add(ii,be[j]);
         for (size_t jn=1; jn<=nen; jn++) {
            for (size_t l=1; l<=nb_dof; l++, i++) {
               if ((jj=e(jn)->getDOF(l)) && ii>0)
                  A->add(ii,jj,ae[i]);
            }
         }
      }
   }
}


template<class E_,class T_>
void element_assembly(const E_&   e,
                      const T_*   ae,
                      Matrix<T_>* A)
{
   size_t i=0;
   for (size_t in=1; in<=e.getNbNodes(); in++) {
      for (size_t k=1; k<=e(in)->getNbDOF(); k++) {
         size_t ii=e(in)->getDOF(k);
         for (size_t jn=1; jn<=e.getNbNodes(); jn++) {
            for (size_t l=1; l<=e(jn)->getNbDOF(); l++, i++) {
               size_t jj=e(jn)->getDOF(l);
               if (ii>0 && jj>0)
                  A->add(ii,jj,ae[i]);
            }
         }
      }
   }
}


template<class E_,class T_>
void element_assembly(const E_& e,
                      const T_* ae,
                      Vect<T_>* A)
{
   size_t i=0;
   for (size_t in=1; in<=e.getNbNodes(); in++) {
      for (size_t k=1; k<=e(in)->getNbDOF(); k++) {
         size_t ii=e(in)->getDOF(k);
         for (size_t jn=1; jn<=e.getNbNodes(); jn++) {
            for (size_t l=1; l<=e(jn)->getNbDOF(); l++, i++)
               if ((ii>0) && (ii==e(jn)->getDOF(l)))
                  A->add(ii,ae[i]);
         }
      }
   }
}


template<class T_>
void update_bc(const Element&  el,
               const Vect<T_>& bc,
               const T_*       eA,
               T_*             eb)
{
   size_t nen=el.getNbNodes(), in=0, ij=0;
   for (size_t i=1; i<=nen; i++) {
      for (size_t k=1; k<=el(i)->getNbDOF(); k++, in++) {
         size_t jn=0;
         for (size_t j=1; j<=nen; j++) {
            for (size_t l=1; l<=el(j)->getNbDOF(); l++, jn++, ij++) {
               if (el(j)->getCode(l)>0)
                  eb[in] -= eA[ij] * bc[jn];
            }
         }
      }
   }
}


template<class T_>
void update_bc(const Side&     sd,
               const Vect<T_>& bc,
               const T_*       sA,
               T_*             sb)
{
  if (sb==NULL)
     return;
   size_t nsn=sd.getNbNodes(), nb_dof=sd.getNbDOF(), in=0, ij=0;
   for (size_t i=1; i<=nsn; i++) {
      for (size_t k=1; k<=nb_dof; k++, in++) {
         size_t jn=0;
         for (size_t j=1; j<=nsn; j++) {
            for (size_t l=1; l<=nb_dof; l++, jn++, ij++) {
               if (sd(j)->getCode(l)>0)
                  sb[in] -= sA[ij] * bc[jn];
            }
         }
      }
   }
}


template<class T_>
void update_bc_diag(const Element&  el,
                    const Vect<T_>& bc,
                    const T_*       eA,
                    T_*             eb)
{
   size_t nen=el.getNbNodes(), nb_dof=el.getNbDOF(), in=0, ij=0;
   for (size_t i=1; i<=nen; i++) {
      for (size_t k=1; k<=nb_dof; k++, in+=nen*nb_dof+1, ij++) {
         if (el(i)->getCode(k)>0)
            eb[ij] -= eA[in] * bc[ij];
      }
   }
}


template<class T_>
void update_bc_diag(const Side&     sd,
                    const Vect<T_>& bc,
                    const T_*       sA,
                    T_*             sb)
{
  if (sA==NULL)
     return;
   size_t nsn=sd.getNbNodes(), nb_dof=sd.getNbDOF(), in=0, ij=0;
   for (size_t i=1; i<=nsn; i++) {
      for (size_t k=1; k<=nb_dof; k++, in+=nsn*nb_dof+1, ij++) {
         if (sd(i)->getCode(k)>0)
            sb[ij] -= sA[in] * bc[ij];
      }
   }
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*  @}  */

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
