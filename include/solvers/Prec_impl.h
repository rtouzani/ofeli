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

                    Implementation of Preconditioner classes

  ==============================================================================*/

#ifndef __PREC_IMPL_H
#define __PREC_IMPL_H

#include <iostream>
using std::ostream;
using std::cout;
using std::endl;

#include <iomanip>
using std::setw;

//#include "OFELI_Config.h"
#include "solvers/Prec.h"
#include "util/util.h"
#include "linear_algebra/Matrix_impl.h"
#include "linear_algebra/SpMatrix_impl.h"
#include "OFELIException.h"

namespace OFELI {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template<class T_> class SpMatrix;

template<class T_,class M_>
int inv_diag(size_t      n,
             const M_&   A,
             vector<T_>& diag)
{
   for (size_t i=1; i<n; ++i) {
      if (A(i,i)==static_cast<T_>(0))
         return int(i);
      if ((diag[i-1]=static_cast<T_>(1)/A(i,i)) == static_cast<T_>(0))
         return -int(i);
   }
   return 0;
}


template<class T_>
int inv_diag(size_t            n,
             const Matrix<T_>* A,
             vector<T_>&       diag)
{
   for (size_t i=1; i<=n; i++) {
      T_ d=(*A)(i,i);
      if (d==static_cast<T_>(0))
         return int(i);
      diag[i-1] = static_cast<T_>(1)/d;
      if (diag[i-1] == static_cast<T_>(0))
         return -int(i);
   }
   return 0;
}


template<class T_>
Prec<T_>::Prec()
         : _type(-1)
{ }


template<class T_>
Prec<T_>::Prec(int type)
         : _type(type)
{ }


template<class T_>
Prec<T_>::Prec(const SpMatrix<T_>& A,
               int                 type)
         : _type(type)
{
   setMatrix(A);
}


template<class T_>
Prec<T_>::Prec(const Matrix<T_>* A,
               int               type)
         : _type(type)
{
   setMatrix(A);
}


template<class T_>
Prec<T_>::~Prec()
{ }


template<class T_>
void Prec<T_>::setType(int type)
{
   _type = type;
}


template<class T_>
void Prec<T_>::setMatrix(const Matrix<T_>* A)
{
   if (_type==-1)
      throw OFELIException("In Prec::setMatrix(Matrix<T_> *): Choose preconditioner before "
                           "setting matrix !");
   _a = (SpMatrix<T_> *)A;
   _size = _a->size();
   _length = _a->getLength();

   switch (_type) {

      case IDENT_PREC:
         break;

      case DIAG_PREC:
          _pivot.resize(_size);
          inv_diag();
          break;

      case DILU_PREC:
          _a->DILUFactorize(_id,_pivot);
          break;

      case ILU_PREC:
          _aa.ILUFactorize(*_a);
          break;

      case SSOR_PREC:
          break;
   }
}


template<class T_>
void Prec<T_>::setMatrix(const SpMatrix<T_>& A)
{
   if (_type==-1)
       throw OFELIException("In Prec::setMatrix(SpMatrix<T_>): Choose preconditioner before "
                            "setting matrix !");
   _a = &A;
   _size = _a->size();
   _length = _a->getLength();

   switch (_type) {

      case IDENT_PREC:
         break;

      case DIAG_PREC:
         _pivot.resize(_size);
         inv_diag();
         break;

      case DILU_PREC:
         _a->DILUFactorize(_id,_pivot);
         break;

      case ILU_PREC:
         _aa.ILUFactorize(*_a);
         break;

      case SSOR_PREC:
         break;
   }
}


template<class T_>
void Prec<T_>::solve(Vect<T_>& x) const
{
   if (Verbosity>10)
      cout << "Solving linear system using preconditioning matrix." << endl;
   Vect<T_> b(x);
   switch (_type) {

      case DIAG_PREC:
         for (size_t i=0; i<_size; i++)
            x[i] *= _pivot[i];
         break;

      case DILU_PREC:
         _a->DILUSolve(_id,_pivot,b,x);
         break;

      case ILU_PREC:
         _aa.ILUSolve(b,x);
         break;

      case SSOR_PREC:
         _a->SSORSolve(b,x);
         break;
   }
}


template<class T_>
void Prec<T_>::solve(const Vect<T_>& b,
                     Vect<T_>&       x) const
{
   if (Verbosity>10)
      cout << "Solving linear system using preconditioning matrix." << endl;
   switch (_type) {

      case IDENT_PREC:
         x = b;
         break;

      case DIAG_PREC:
         for (size_t i=0; i<_size; i++)
            x[i] = b[i]*_pivot[i];
         break;

      case DILU_PREC:
         _a->DILUSolve(_id,_pivot,b,x);
         break;

      case ILU_PREC:
         _aa.ILUSolve(b,x);
         break;

      case SSOR_PREC:
         _a->SSORSolve(b,x);
         break;
   }
}


template<class T_>
void Prec<T_>::TransSolve(Vect<T_>& x) const
{
   if (Verbosity>10)
      cout << "Solving linear system using transpose of preconditioning matrix." << endl;
   switch (_type) {

      case DIAG_PREC:
         solve(x);
         break;

      case ILU_PREC:
         {
            Vect<T_> z(_size), temp(x);
            T_ tmp;
            for (size_t i=0; i<_size; i++) {
               z[i] = temp[i];
               tmp = _pivot[i] * z[i];
               for (size_t j=_id[i]+1; j<_row_ptr[i+1]; j++)
                  temp(_col_ind[j-1]) -= tmp * (*_a)(j);
            }
            for (int i=int(_size)-1; i>=0; i--) {
               x[i] = _pivot[i]*z[i];
               for (size_t j=_row_ptr[i]; j<_id[i]; j++)
                  z(_col_ind[j-1]) -= (*_a)(j) * x[i];
            }
         }
         break;
   }
}


template<class T_>
void Prec<T_>::TransSolve(const Vect<T_>& b,
                          Vect<T_>&       x) const
{
   x = b;
   TransSolve(x);
}


template<class T_>
T_& Prec<T_>::getPivot(size_t i) const
{
   return _pivot[i-1];
}


template<class T_>
int Prec<T_>::inv_diag()
{
   for (size_t i=1; i<=_a->getNbRows(); i++) {
      if ((*_a)(i,i) == static_cast<T_>(0))
         throw OFELIException("In Prec::inv_diag(): null diagonal term: " + itos(i));
      if ((_pivot[i-1]=static_cast<T_>(1)/(*_a)(i,i)) == static_cast<T_>(0))
         return -i;
   }
   return 0;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */

#endif
