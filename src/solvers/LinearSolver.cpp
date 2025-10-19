/*==============================================================================

                                    O  F  E  L  I

                          Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2025 Rachid Touzani

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

                       Implementation of class LinearSolver

  ==============================================================================*/

#include <string>
#include <iostream>
using std::ostream;
using std::endl;

#include <iomanip>
using std::setw;

#include "solvers/LinearSolver.h"

#if defined (USE_PETSC)
#include "linear_algebra/petsc/PETScMatrix.h"
#else
#include "util/util.h"
#include "io/output.h"
#include "solvers/CG.h"
#include "solvers/CGS.h"
#include "solvers/BiCG.h"
#include "solvers/BiCGStab.h"
#include "solvers/GMRes.h"
#endif
#include "linear_algebra/Matrix_impl.h"
#include "OFELIException.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */


LinearSolver::LinearSolver() : _fact(0), _max_it(1000), _matrix_set(0),
                               _s(DIRECT_SOLVER), _p(DIAG_PREC), _toler(sqrt(OFELI_EPSMCH)),
                               _x(nullptr), _b(nullptr), _A(nullptr)
{ }


LinearSolver::LinearSolver(int    max_it,
                           real_t tolerance) : _fact(0), _max_it(max_it),
                                               _matrix_set(0), _s(DIRECT_SOLVER), _p(DIAG_PREC),
                                               _toler(tolerance), _x(nullptr), _b(nullptr), _A(nullptr)
{ }


LinearSolver::LinearSolver(SpMatrix<real_t>&   A,
                           const Vect<real_t>& b,
                           Vect<real_t>&       x) : _fact(0), _max_it(1000), _matrix_set(1),
                                                    _s(CG_SOLVER), _p(DIAG_PREC), _toler(sqrt(OFELI_EPSMCH)),
                                                    _x(&x), _b(&b), _A(&A)
{ }


LinearSolver::LinearSolver(SkMatrix<real_t>&   A,
                           const Vect<real_t>& b,
                           Vect<real_t>&       x) : _fact(0), _max_it(1000), _matrix_set(1),
                                                    _s(DIRECT_SOLVER), _p(DIAG_PREC), _toler(sqrt(OFELI_EPSMCH)),
                                                    _x(&x), _b(&b), _A(&A)
{ }


LinearSolver::LinearSolver(TrMatrix<real_t>&   A,
                           const Vect<real_t>& b,
                           Vect<real_t>&       x) : _fact(0), _max_it(1000), _matrix_set(1),
                                                    _s(DIRECT_SOLVER), _p(DIAG_PREC), _toler(sqrt(OFELI_EPSMCH)),
                                                    _x(&x), _b(&b), _A(&A)
{ }


LinearSolver::LinearSolver(BMatrix<real_t>&    A,
                           const Vect<real_t>& b,
                           Vect<real_t>&       x) : _fact(0), _max_it(1000), _matrix_set(1),
                                                    _s(DIRECT_SOLVER), _p(DIAG_PREC), _toler(sqrt(OFELI_EPSMCH)),
                                                    _x(&x), _b(&b), _A(&A)
{ }


LinearSolver::LinearSolver(DMatrix<real_t>&    A,
                           const Vect<real_t>& b,
                           Vect<real_t>&       x) : _fact(0), _max_it(1000), _matrix_set(1),
                                                    _s(DIRECT_SOLVER), _p(DIAG_PREC), _toler(sqrt(OFELI_EPSMCH)),
                                                    _x(&x), _b(&b), _A(&A)
{ }


LinearSolver::LinearSolver(DSMatrix<real_t>&   A,
                           const Vect<real_t>& b,
                           Vect<real_t>&       x) : _fact(0), _max_it(1000), _matrix_set(1),
                                                    _s(DIRECT_SOLVER), _p(DIAG_PREC), _toler(sqrt(OFELI_EPSMCH)),
                                                    _x(&x), _b(&b), _A(&A)
{ }


LinearSolver::LinearSolver(SkSMatrix<real_t>&  A,
                           const Vect<real_t>& b,
                           Vect<real_t>&       x) : _fact(0), _max_it(1000), _matrix_set(1),
                                                    _s(DIRECT_SOLVER), _p(DIAG_PREC), _toler(sqrt(OFELI_EPSMCH)),
                                                    _x(&x), _b(&b), _A(&A)
{ }


LinearSolver::LinearSolver(SkMatrix<real_t>& A,
                           Vect<real_t>&     b,
                           Vect<real_t>&     x) : _fact(0), _max_it(1000), _matrix_set(1),
                                                  _s(DIRECT_SOLVER), _p(DIAG_PREC), _toler(sqrt(OFELI_EPSMCH)),
                                                  _x(&x), _b(&b), _A(&A)
{ }


void LinearSolver::setMatrix(OFELI::Matrix<real_t>* A)
{
   _A = A;
   _matrix_set = 1;
   _fact = 0;
}


void LinearSolver::setMatrix(SpMatrix<real_t>& A)
{
   _A = &A;
   _matrix_set = 1;
}


void LinearSolver::setMatrix(SkMatrix<real_t>& A)
{
   _A = &A;
   _matrix_set = 1;
   _fact = 0;
}


void LinearSolver::set(SpMatrix<real_t>&   A,
                       const Vect<real_t>& b,
                       Vect<real_t>&       x)
{
   _A = &A;
   _matrix_set = 1;
   _b = &b;
   _x = &x;
}


int LinearSolver::solve(SpMatrix<real_t>&   A,
                        const Vect<real_t>& b,
                        Vect<real_t>&       x,
                        Iteration           s,
                        Preconditioner      p)
{
   set(A,b,x);
   setSolver(s,p);
   return solve();
}


int LinearSolver::solve(Iteration      s,
                        Preconditioner p)
{
   setSolver(s,p);
   return solve();
}


int LinearSolver::solve()
{
   int ret=0;
   _nb_it = 0;
   if (!_matrix_set)
      throw OFELIException("In LinearSolver::solve(): No matrix has been defined.");
   if (_A->isDiagonal()) {
      if (_x) {
         for (size_t i=1; i<=_A->getNbRows(); i++)
            (*_x)(i) = (*_b)(i)/(*_A)(i,i);
      }
      else {
         for (size_t i=1; i<=_A->getNbRows(); i++)
            (*_x)(i) /= (*_A)(i,i);
      }
      return ret;
   }

   switch (_s) {

      case DIRECT_SOLVER:
         if (_fact)
            _A->Factor();
         _fact = 1;
         ret = _A->solve(*_b,*_x);
         break;

      case CG_SOLVER:
         _nb_it = CG(_A,_p,*_b,*_x,_max_it,_toler);
         break; 

      case BICG_SOLVER:
         _nb_it = BiCG(_A,_p,*_b,*_x,_max_it,_toler);
         break;

      case CGS_SOLVER:
         _nb_it = CGS(_A,_p,*_b,*_x,_max_it,_toler);
         break;

      case BICG_STAB_SOLVER:
         _nb_it = BiCGStab(_A,_p,*_b,*_x,_max_it,_toler);
         break;

      case GMRES_SOLVER:
         _nb_it = GMRes(_A,_p,*_b,*_x,_b->size()/5,_max_it,_toler);
         break;
   }
   if (_nb_it<0)
      ret = -_nb_it;
   return ret;
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */