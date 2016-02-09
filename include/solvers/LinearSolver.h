/*==============================================================================

                                    O  F  E  L  I

                          Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzani

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

                       Definition of class LinearSolver
       to select iterative or direct solver for a linear system of equations

  ==============================================================================*/

#ifndef __LINEAR_SOLVER_H
#define __LINEAR_SOLVER_H

#include <string>
#include <iostream>
using std::ostream;
using std::endl;

#include <iomanip>
using std::setw;

#include "OFELI_Config.h"

#if defined (USE_EIGEN)
#include "linear_algebra/Matrix.h"
#include <Eigen/IterativeLinearSolvers>
using namespace Eigen;
#else
#include "util/util.h"
#include "io/output.h"
#include "solvers/CG.h"
#include "solvers/CGS.h"
#include "solvers/BiCG.h"
#include "solvers/BiCGStab.h"
#include "solvers/GMRes.h"
#include "linear_algebra/SpMatrix.h"
#include "linear_algebra/SkMatrix.h"
#include "linear_algebra/SkSMatrix.h"
#include "linear_algebra/TrMatrix.h"
#include "linear_algebra/BMatrix.h"
#include "linear_algebra/DMatrix.h"
#include "linear_algebra/DSMatrix.h"
#endif

namespace OFELI {

template<class T_> class Matrix;
template<class T_> class SpMatrix;
template<class T_> class SkMatrix;
template<class T_> class SkSMatrix;
template<class T_> class TrMatrix;
template<class T_> class BMatrix;
template<class T_> class DMatrix;
template<class T_> class DSMatrix;

template<class T_> class LinearSolver
{

 public:

#if defined (USE_EIGEN)
    typedef Eigen::Matrix<T_,Eigen::Dynamic,1> VectorX;
#endif

/// \brief Default Constructor.
/// \details Initializes default parameters and pointers to 0.
    LinearSolver() : _fact(0), _verbose(0), _max_it(1000), _matrix_set(0),
                     _s(DIRECT_SOLVER), _p(DIAG_PREC), _toler(sqrt(OFELI_EPSMCH)),
                     _x(NULL), _b(NULL), _A(NULL)
    {  }

/** \brief Constructor with iteration parameters
 *  @param [in] max_it Maximal number of iterations
 *  @param [in] toler Tolerance for convergence (measured in relative weighted 2-Norm) in input,
 *  effective discrepancy in output.
 *  @param [in] verbose Information output parameter
 *  <ul>
 *    <li><tt>0</tt>: No output
 *    <li><tt>1</tt>: Output iteration information,
 *    <li><tt>2</tt> and greater: Output iteration information and solution at each iteration.
 *  </ul>
 */
    LinearSolver(int    max_it,
                 real_t tolerance,
                 int    verbose) : _fact(0), _verbose(verbose), _max_it(max_it),
                                   _matrix_set(0), _s(DIRECT_SOLVER), _p(DIAG_PREC),
                                   _toler(tolerance), _x(NULL), _b(NULL), _A(NULL)
    {  }

/** \brief Constructor using matrix, right-hand side and solution vector
 *  @param [in] A Reference to instance of class SpMatrix
 *  @param [in] b Vect instance that contains the right-hand side
 *  @param [in,out] x Vect instance that contains initial guess on input and solution on output
 */
    LinearSolver(      SpMatrix<T_>& A,
                 const Vect<T_>&     b,
                       Vect<T_>&     x) : _fact(0), _verbose(0), _max_it(1000), _matrix_set(1),
                                          _s(CG_SOLVER), _p(DIAG_PREC),
                                          _toler(sqrt(OFELI_EPSMCH)), _x(&x), _b(&b), _A(&A)
    {  }

/** \brief Constructor using skyline-stored matrix, right-hand side and solution vector
 *  @param [in] A SkMatrix instance that contains matrix
 *  @param [in] b Vect instance that contains the right-hand side
 *  @param [in,out] x Vect instance that contains initial guess on input and solution on output
 */
    LinearSolver(      SkMatrix<T_>& A,
                 const Vect<T_>&     b,
                       Vect<T_>&     x) : _fact(0), _verbose(0), _max_it(1000), _matrix_set(1),
                                          _s(DIRECT_SOLVER), _p(DIAG_PREC),
                                          _toler(sqrt(OFELI_EPSMCH)), _x(&x), _b(&b), _A(&A)
    {  }

/** \brief Constructor using a tridiagonal matrix, right-hand side and solution vector
 *  @param [in] A TrMatrix instance that contains matrix
 *  @param [in] b Vect instance that contains the right-hand side
 *  @param [in,out] x Vect instance that contains initial guess on input and solution on output
 */
    LinearSolver(      TrMatrix<T_>& A,
                 const Vect<T_>&     b,
                       Vect<T_>&     x) : _fact(0), _verbose(0), _max_it(1000), _matrix_set(1),
                                          _s(DIRECT_SOLVER), _p(DIAG_PREC),
                                          _toler(sqrt(OFELI_EPSMCH)), _x(&x), _b(&b), _A(&A)
    {  }

/** \brief Constructor using a banded matrix, right-hand side and solution vector
 *  @param [in] A BMatrix instance that contains matrix
 *  @param [in] b Vect instance that contains the right-hand side
 *  @param [in,out] x Vect instance that contains initial guess on input and solution on output
 */
    LinearSolver(      BMatrix<T_>& A,
                 const Vect<T_>&    b,
                       Vect<T_>&    x) : _fact(0), _verbose(0), _max_it(1000), _matrix_set(1),
                                         _s(DIRECT_SOLVER), _p(DIAG_PREC),
                                         _toler(sqrt(OFELI_EPSMCH)), _x(&x), _b(&b), _A(&A)
    {  }

/** \brief Constructor using a dense matrix, right-hand side and solution vector
 *  @param [in] A DMatrix instance that contains matrix
 *  @param [in] b Vect instance that contains the right-hand side
 *  @param [in,out] x Vect instance that contains initial guess on input and solution on output
 */
    LinearSolver(      DMatrix<T_>& A,
                 const Vect<T_>&    b,
                       Vect<T_>&    x) : _fact(0), _verbose(0), _max_it(1000), _matrix_set(1),
                                         _s(DIRECT_SOLVER), _p(DIAG_PREC),
                                         _toler(sqrt(OFELI_EPSMCH)), _x(&x), _b(&b), _A(&A)
    {  }

/** \brief Constructor using a dense symmetric matrix, right-hand side and solution vector
 *  @param [in] A DSMatrix instance that contains matrix
 *  @param [in] b Vect instance that contains the right-hand side
 *  @param [in,out] x Vect instance that contains initial guess on input and solution on output
 */
    LinearSolver(      DSMatrix<T_>& A,
                 const Vect<T_>&     b,
                       Vect<T_>&     x) : _fact(0), _verbose(0), _max_it(1000), _matrix_set(1),
                                          _s(DIRECT_SOLVER), _p(DIAG_PREC),
                                          _toler(sqrt(OFELI_EPSMCH)), _x(&x), _b(&b), _A(&A)
    {  }

/** \brief Constructor using skyline-stored symmetric matrix, right-hand side and solution vector
 *  @param [in] A SkMatrix instance that contains matrix
 *  @param [in] b Vect instance that contains the right-hand side
 *  @param [in,out] x Vect instance that contains initial guess on input and solution on output
 */
    LinearSolver(      SkSMatrix<T_>& A,
                 const Vect<T_>&      b,
                       Vect<T_>&      x) : _fact(0), _verbose(0), _max_it(1000), _matrix_set(1),
                                          _s(DIRECT_SOLVER), _p(DIAG_PREC),
                                          _toler(sqrt(OFELI_EPSMCH)), _x(&x), _b(&b), _A(&A)
    {  }

/** \brief Constructor using matrix, right-hand side
 *  @param [in] A SkMatrix instance that contains matrix
 *  @param [in,out] b Vect instance that contains the right-hand side on input and solution on output.
 *  The initial guess is the 0-vector
 */
    LinearSolver(SkMatrix<T_>& A,
                 Vect<T_>&     b,
                 Vect<T_>&     x) : _fact(0), _verbose(0), _max_it(1000), _matrix_set(1),
                                    _s(DIRECT_SOLVER), _p(DIAG_PREC),
                                    _toler(sqrt(OFELI_EPSMCH)), _x(&x), _b(&b), _A(&A)
    {  }

/// Destructor
    virtual ~LinearSolver() { }

/// \brief Set message level
/// \details Default value is 0
    void setVerbose(int verb) { _verbose = verb; }

/// \brief Set Maximum number of iterations
/// \details Default value is 1000
    void setMaxIter(int m) { _max_it = m; }

/// \brief Set tolerance value
    void setTolerance(real_t tol) { _toler = tol; }

/// \brief Set solution vector
    void setSolution(Vect<T_>& x)
    {
       _x = &x;
    }

/// \brief Set right-hand side vector
    void setRHS(Vect<T_>& b)
    {
       _b = &b;
    }

/// \brief Set matrix in the case of a pointer to Matrix
/// @param [in] A Pointer to abstract Matrix class
    void setMatrix(OFELI::Matrix<T_>* A)
    {
       _A = A;
       _matrix_set = 1;
       _fact = 0;
    }

/// \brief Set matrix in the case of a pointer to matrix
/// @param [in] A Pointer to abstract Matrix class
    void setMatrix(SpMatrix<T_>& A)
    {
       _A = &A;
       _matrix_set = 1;
    }

/// \brief Set matrix in the case of a skyline matrix
/// @param [in] A %Matrix as instance of class SkMatrix
    void setMatrix(SkMatrix<T_>& A)
    {
       _A = &A;
       _matrix_set = 1;
       _fact = 0;
    }

/** \brief Set matrix, right-hand side and initial guess
 *  @param [in] A Reference to matrix as a SpMatrix instance
 *  @param [in] b Vector containing right-hand side
 *  @param [in,out] x Vector containing initial guess on input and solution on output
 */
    void set(      SpMatrix<T_>& A,
             const Vect<T_>&     b,
                   Vect<T_>&     x)
    {
       _A = &A;
       _matrix_set = 1;
       _fact = 0;
       _b = &b;
       _x = &x;
    }

/** \brief Set solver and preconditioner
 *  @param [in] s Solver identification parameter.
 *  To be chosen in the enumeration variable Iteration:\n 
 *  <tt>DIRECT_SOLVER</tt>, <tt>CG_SOLVER</tt>, <tt>CGS_SOLVER<tt>, <tt>BICG_SOLVER</tt>,
 *  <tt>BICG_STAB_SOLVER</tt>, <tt>GMRES_SOLVER</tt>, <tt>QMR_SOLVER</tt>
 *  @param [in] p Preconditioner identification parameter. By default, the diagonal
 *  preconditioner is used.
 *  To be chosen in the enumeration variable Preconditioner:\n
 *  <tt>IDENT_PREC</tt>, <tt>DIAG_PREC</tt>, <tt>SSOR_PREC</tt>, <tt>ILU_PREC</tt> [Default:
 *  <tt>ILU_PREC</tt>]
 *  @note The argument <tt>p</tt> has no effect if the solver is <tt>DIRECT_SOLVER</tt>
 */
    void setSolver(Iteration      s, 
                   Preconditioner p=DIAG_PREC)
    { _s = s; _p = p; }

/// \brief Return solver code
    int getSolver() const { return _s; }

/** \brief Solve equations using system data, prescribed solver and preconditioner
 *  @param [in] A Reference to matrix as a SpMatrix instance
 *  @param [in] b Vector containing right-hand side
 *  @param [in,out] x Vector containing initial guess on input and solution on output
 *  @param [in] s Solver identification parameter
 *  To be chosen in the enumeration variable Iteration:\n 
 *  <tt>DIRECT_SOLVER</tt>, <tt>CG_SOLVER</tt>, <tt>CGS_SOLVER<tt>, <tt>BICG_SOLVER</tt>,
 *  <tt>BICG_STAB_SOLVER</tt>, <tt>GMRES_SOLVER</tt>, <tt>QMR_SOLVER</tt> [Default: <tt>CGS_SOLVER<tt>]
 *  @param [in] p Preconditioner identification parameter.
 *  To be chosen in the enumeration variable Preconditioner:\n
 *  <tt>IDENT_PREC</tt>, <tt>DIAG_PREC</tt>, <tt>SSOR_PREC</tt>, <tt>ILU_PREC</tt>,
 *  <tt>DILU_PREC</tt> [Default: <tt>DIAG_PREC<tt>]
 *  @remark The argument <tt>p</tt> has no effect if the solver is <tt>DIRECT_SOLVER</tt>
 *  @warning If the library <it>eigen</it> is used, only the preconditioners
 *  <tt>IDENT_PREC</tt>, <tt>DIAG_PREC</tt> and <tt>ILU_PREC</tt> are available.
 */
    int solve(      SpMatrix<T_>&  A,
              const Vect<T_>&      b,
                    Vect<T_>&      x,
                    Iteration      s,
                    Preconditioner p=DIAG_PREC)
    {
       set(A,b,x);
       int ret = solve(s,p);
       return ret;
    }

/** Solve equations using prescribed solver and preconditioner
 *  @param [in] s Solver identification parameter
 *  To be chosen in the enumeration variable Iteration:\n 
 *  <tt>DIRECT_SOLVER</tt>, <tt>CG_SOLVER</tt>, <tt>CGS_SOLVER<tt>, <tt>BICG_SOLVER</tt>,
 *  <tt>BICG_STAB_SOLVER</tt>, <tt>GMRES_SOLVER</tt>, <tt>QMR_SOLVER</tt> [Default: <tt>CGS_SOLVER<tt>]
 *  @param [in] p Preconditioner identification parameter.
 *  To be chosen in the enumeration variable Preconditioner:\n
 *  <tt>IDENT_PREC</tt>, <tt>DIAG_PREC</tt>, <tt>SSOR_PREC</tt>, <tt>DILU_PREC</tt>,
 *  <tt>ILU_PREC</tt> [Default: <tt>DIAG_PREC<tt>]
 *  @note The argument <tt>p</tt> has no effect if the solver is <tt>DIRECT_SOLVER</tt>
 */
    int solve(Iteration      s,
              Preconditioner p=DIAG_PREC);

/** \brief Solve equations all arguments must have given by other member functions
 *  \details Solver and preconditioner parameters must have been set by function setSolver.
 *  Otherwise, default values are set.
 */
    int solve() { return solve(_s,_p); }

/// Factorize matrix
    void setFact() { _fact = 1; }

/// Do not factorize matrix
    void setNoFact() { _fact = 0; }

 private:

   int            _fact, _verbose, _max_it, _matrix_set;
   Iteration      _s;
   Preconditioner _p;
   real_t         _toler;
   Vect<T_>       *_x;
   const Vect<T_> *_b;
   Matrix<T_>     *_A;
};


///////////////////////////////////////////////////////////////////////////////
//                         I M P L E M E N T A T I O N                       //
///////////////////////////////////////////////////////////////////////////////

template<class T_>
int LinearSolver<T_>::solve(Iteration      s,
                            Preconditioner p)
{
   _s = s, _p = p;
   int ret=0;
   try {
      if (!_matrix_set)
         THROW_RT("solve(): No matrix has been defined.");
   }
   CATCH("LinearSolver");

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
   SpMatrix<T_> &A=(SpMatrix<T_> &)(*_A);

#if defined (USE_EIGEN)
   VectorX x;
   try {
      if (_s==CG_SOLVER) {
         try {
            if (_p==IDENT_PREC) {
               ConjugateGradient<SparseMatrix<T_>,Lower,IdentityPreconditioner> im;
               im.setTolerance(_toler);
               im.compute(A.getEigenMatrix());
               x = im.solve(VectorX(*_b));
               ret = im.iterations();
            }
            else if (_p==DIAG_PREC) {
               ConjugateGradient<SparseMatrix<T_>,Lower,DiagonalPreconditioner<T_> > im;
               im.setTolerance(_toler);
               im.compute(A.getEigenMatrix());
               x = im.solve(VectorX(*_b));
               ret = im.iterations();
            }
            else if (_p==ILU_PREC) {
               ConjugateGradient<SparseMatrix<T_>,Lower,IncompleteLUT<T_> > im;
               im.setTolerance(_toler);
               im.compute(A.getEigenMatrix());
               x = im.solve(VectorX(*_b));
               ret = im.iterations();
            }
            else
               THROW_RT("solve(Iteration,Preconditioner): This preconditioner is not available in the eigen library.");
         }
         CATCH("LinearSolver");
      }
      else if (_s==BICG_STAB_SOLVER) {
         try {
            if (_p==IDENT_PREC) {
               Eigen::BiCGSTAB<SparseMatrix<T_>,IdentityPreconditioner> im;
               im.setTolerance(_toler);
               im.compute(A.getEigenMatrix());
               x = im.solve(VectorX(*_b));
               ret = im.iterations();
            }
            else if (_p==DIAG_PREC) {
               Eigen::BiCGSTAB<SparseMatrix<T_>,DiagonalPreconditioner<T_> > im;
               im.setTolerance(_toler);
               im.compute(A.getEigenMatrix());
               x = im.solve(VectorX(*_b));
               ret = im.iterations();
            }
            else if (_p==ILU_PREC) {
               Eigen::BiCGSTAB<SparseMatrix<T_>,IncompleteLUT<T_> > im;
               im.setTolerance(_toler);
               im.compute(A.getEigenMatrix());
               x = im.solve(VectorX(*_b));
               ret = im.iterations();
            }
            else
               THROW_RT("solve(Iteration,Preconditioner): This preconditioner is not available in the eigen library.");
         }
         CATCH("LinearSolver");
      }
      else
         THROW_RT("solve(Iteration,Preconditioner): This iterative solver is not available in the eigen library.");
   }
   CATCH("LinearSolver");
   _x->setSize(x.size(),1,1);
   *_x = x;
#else
   try {
      if (_s==DIRECT_SOLVER) {
         if (_fact)
            _A->Factor();
         _A->solve(*_b,*_x);
      }
      else if (_s==CG_SOLVER)
         ret = CG(A,_p,*_b,*_x,_max_it,_toler,_verbose);
      else if (_s==GMRES_SOLVER)
         ret = GMRes(A,_p,*_b,*_x,10,_max_it,_toler,_verbose);
      else if (_s==CGS_SOLVER)
         ret = CGS(A,_p,*_b,*_x,_max_it,_toler,_verbose);
      else if (_s==BICG_SOLVER)
         ret = BiCG(A,_p,*_b,*_x,_max_it,_toler,_verbose);
      else if (_s==BICG_STAB_SOLVER)
         ret = BiCGStab(A,_p,*_b,*_x,_max_it,_toler,_verbose);
      else
         THROW_RT("solve(Iteration,Preconditioner): This solver is not available.");
   }
   CATCH("LinearSolver");
#endif
   return ret;
}

} /* namespace OFELI */

#endif
