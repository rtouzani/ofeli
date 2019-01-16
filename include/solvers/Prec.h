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

                             Preconditioner classes

  ==============================================================================*/

#ifndef __PRECOND_H
#define __PRECOND_H

#include <iostream>
using std::ostream;
using std::endl;

#include <iomanip>
using std::setw;

#include "OFELI_Config.h"
#include "linear_algebra/Vect.h"
#include "util/util.h"
#include "linear_algebra/Matrix.h"
#include "OFELIException.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Prec.h
 *  \brief Definition file for preconditioning classes.
 */

template<class T_> class SpMatrix;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_,class M_>
size_t inv_diag(size_t      n,
                const M_&   A,
                vector<T_>& diag)
{
   for (int i=1; i<=int(n); i++) {
      if (A(i,i)==static_cast<T_>(0))
         return i;
      if ((diag[i-1]=static_cast<T_>(1)/A(i,i)) == static_cast<T_>(0))
         return -i;
   }
   return 0;
}

template<class T_>
size_t inv_diag(size_t            n,
                const Matrix<T_>* A,
                vector<T_>&       diag)
{
   for (size_t i=1; i<=n; i++) {
      T_ d=(*A)(i,i);
      if (d==static_cast<T_>(0))
         return i;
      diag[i-1] = static_cast<T_>(1)/d;
      if (diag[i-1] == static_cast<T_>(0))
         return -i;
   }
   return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! \class Prec
 *  \ingroup Solver
 *  \brief To set a preconditioner.
 *  \details The preconditioner type is chosen in the constructor
 *
 * \tparam <T_> Data type (real_t, float, complex<real_t>, ...)
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

template<class T_> class Prec
{

 public:

/// \brief Default constructor.
    Prec() : _type(-1) { }

/** \brief Constructor that chooses preconditioner.
 *  @param [in] type Preconditioner type:
 *  <ul>
 *    <li><tt>IDENT_PREC</tt>: Identity preconditioner (No preconditioning)
 *    <li><tt>DIAG_PREC</tt>: Diagonal preconditioner
 *    <li><tt>DILU_PREC</tt>: Diagonal Incomplete factorization preconditioner
 *    <li><tt>ILU_PREC</tt>: Incomplete factorization preconditioner
 *    <li><tt>SSOR_PREC</tt>: SSOR (Symmetric Successive Over Relaxation) preconditioner
 *  </ul>
 */
    Prec(int type) : _type(type) { }

/** \brief Constructor using matrix of the linear system to precondition
 *  @param [in] A %Matrix to precondition
 *  @param [in] type Preconditioner type:
 *  <ul>
 *    <li><tt>IDENT_PREC</tt>: Identity preconditioner (No preconditioning)
 *    <li><tt>DIAG_PREC</tt>: Diagonal preconditioner
 *    <li><tt>DILU_PREC</tt>: Diagonal Incomplete factorization preconditioner
 *    <li><tt>ILU_PREC</tt>: Incomplete factorization preconditioner
 *    <li><tt>SSOR_PREC</tt>: SSOR (Symmetric Successive Over Relaxation) preconditioner
 *  </ul>
 */
    Prec(const SpMatrix<T_>& A,
         int                 type=DIAG_PREC) : _type(type) { setMatrix(A); }

/** \brief Constructor using matrix of the linear system to precondition
 *  @param [in] A Pointer to abstract Matrix class to precondition
 *  @param [in] type Preconditioner type:
 *  <ul>
 *    <li><tt>IDENT_PREC</tt>: Identity preconditioner (No preconditioning)
 *    <li><tt>DIAG_PREC</tt>: Diagonal preconditioner
 *    <li><tt>DILU_PREC</tt>: Diagonal Incomplete factorization preconditioner
 *    <li><tt>ILU_PREC</tt>: Incomplete factorization preconditioner
 *    <li><tt>SSOR_PREC</tt>: SSOR (Symmetric Successive Over Relaxation) preconditioner
 *  </ul>
 */
    Prec(const Matrix<T_>* A,
         int               type=DIAG_PREC) : _type(type) { setMatrix(A); }

/// \brief Destructor
    ~Prec() { }

/** \brief Define preconditioner type
 *  @param [in] type Preconditioner type:
    <ul>
       <li> <tt>IDENT_PREC</tt>: Identity preconditioner (No preconditioning)
       <li> <tt>DIAG_PREC</tt>: Diagonal preconditioner
       <li> <tt>DILU_PREC</tt>: Diagonal Incomplete factorization preconditioner
       <li> <tt>ILU_PREC</tt>: Incomplete factorization preconditioner
       <li> <tt>SSOR_PREC</tt>: SSOR (Symmetric Successive Over Relaxation) preconditioner
    </ul>
 */
    void setType(int type) { _type = type; }

/// \brief Define pointer to matrix for preconditioning (if this one is abstract)
/// @param [in] A %Matrix to precondition
    void setMatrix(const Matrix<T_>* A)
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

/// \brief Define the matrix for preconditioning
/// @param [in] A %Matrix to precondition (instance of class SpMatrix)
    void setMatrix(const SpMatrix<T_>& A)
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

/// \brief Solve a linear system with preconditioning matrix.
/// @param [in,out] x Right-hand side on input and solution on output.
    void solve(Vect<T_>& x) const
    {
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

/** \brief Solve a linear system with preconditioning matrix.
 *  @param [in] b Right-hand side
 *  @param [out] x Solution vector
 */
    void solve(const Vect<T_>& b,
               Vect<T_>&       x) const
    {

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

/// \brief Solve a linear system with transposed preconditioning matrix.
/// @param [in,out] x Right-hand side in input and solution in output.
    void TransSolve(Vect<T_>& x) const
    {
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


/** \brief Solve a linear system with transposed preconditioning matrix.
 *  @param [in] b Right-hand side vector
 *  @param [out] x Solution vector
 */
    void TransSolve(const Vect<T_>& b,
                    Vect<T_>&       x) const { x = b; TransSolve(x); }

/// Return i-th pivot of preconditioning matrix
    T_ & getPivot(size_t i) const { return _pivot[i-1]; }

 private:
   vector<T_>         _pivot;
   const SpMatrix<T_> *_a;
   SpMatrix<T_>       _aa;
   size_t             _size, _length;
   vector<size_t>     _id, _row_ptr, _col_ind;
   int                _type;

   size_t inv_diag()
   {
      for (size_t i=1; i<=_a->getNbRows(); i++) {
         if ((*_a)(i,i) == static_cast<T_>(0))
            throw OFELIException("In Prec::inv_diag(): null diagonal term: " + itos(i));
         if ((_pivot[i-1]=static_cast<T_>(1)/(*_a)(i,i)) == static_cast<T_>(0))
            return -i;
      }
      return 0;
   }
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
