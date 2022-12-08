/*==============================================================================

                                    O  F  E  L  I

                          Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani

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

#include "OFELI_Config.h"
#include "linear_algebra/SpMatrix.h"
#include "linear_algebra/SkMatrix.h"
#include "linear_algebra/SkSMatrix.h"
#include "linear_algebra/TrMatrix.h"
#include "linear_algebra/BMatrix.h"
#include "linear_algebra/DMatrix.h"
#include "linear_algebra/DSMatrix.h"
#include "linear_algebra/Vect.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

#if defined(USE_PETSC)
class PETScMatrix<real_t>;
#endif

//-----------------------------------------------------------------------------
// Class LinearSolver
//-----------------------------------------------------------------------------

/*! \class LinearSolver
 *  \ingroup Solver
 *  \brief Class to solve systems of linear equations by direct or iterative methods
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class LinearSolver
{

 public:

#if defined (USE_EIGEN)
    typedef Eigen::Matrix<real_t,Eigen::Dynamic,1> VectorX;
#endif

/// \brief Default Constructor.
/// \details Initializes default parameters and pointers to 0.
    LinearSolver();

/** \brief Constructor with iteration parameters
 *  @param [in] max_it Maximal number of iterations
 *  @param [in] tolerance Tolerance for convergence (measured in relative weighted 2-Norm) in input,
 *  effective discrepancy in output.
 */
    LinearSolver(int    max_it,
                 real_t tolerance);

/** \brief Constructor using matrix, right-hand side and solution vector
 *  @param [in] A Reference to instance of class SpMatrix
 *  @param [in] b Vect instance that contains the right-hand side
 *  @param [in,out] x Vect instance that contains initial guess on input and solution on output
 */
    LinearSolver(SpMatrix<real_t>&   A,
                 const Vect<real_t>& b,
                 Vect<real_t>&       x);

/** \brief Constructor using skyline-stored matrix, right-hand side and solution vector
 *  @param [in] A SkMatrix instance that contains matrix
 *  @param [in] b Vect instance that contains the right-hand side
 *  @param [in,out] x Vect instance that contains initial guess on input and solution on output
 */
    LinearSolver(SkMatrix<real_t>&   A,
                 const Vect<real_t>& b,
                 Vect<real_t>&       x);

/** \brief Constructor using a tridiagonal matrix, right-hand side and solution vector
 *  @param [in] A TrMatrix instance that contains matrix
 *  @param [in] b Vect instance that contains the right-hand side
 *  @param [in,out] x Vect instance that contains initial guess on input and solution on output
 */
    LinearSolver(TrMatrix<real_t>&   A,
                 const Vect<real_t>& b,
                 Vect<real_t>&       x);

/** \brief Constructor using a banded matrix, right-hand side and solution vector
 *  @param [in] A BMatrix instance that contains matrix
 *  @param [in] b Vect instance that contains the right-hand side
 *  @param [in,out] x Vect instance that contains initial guess on input and solution on output
 */
    LinearSolver(BMatrix<real_t>&    A,
                 const Vect<real_t>& b,
                 Vect<real_t>&       x);

/** \brief Constructor using a dense matrix, right-hand side and solution vector
 *  @param [in] A DMatrix instance that contains matrix
 *  @param [in] b Vect instance that contains the right-hand side
 *  @param [in,out] x Vect instance that contains initial guess on input and solution on output
 */
    LinearSolver(DMatrix<real_t>&    A,
                 const Vect<real_t>& b,
                 Vect<real_t>&       x);

/** \brief Constructor using a dense symmetric matrix, right-hand side and solution vector
 *  @param [in] A DSMatrix instance that contains matrix
 *  @param [in] b Vect instance that contains the right-hand side
 *  @param [in,out] x Vect instance that contains initial guess on input and solution on output
 */
    LinearSolver(DSMatrix<real_t>&   A,
                 const Vect<real_t>& b,
                 Vect<real_t>&       x);

/** \brief Constructor using skyline-stored symmetric matrix, right-hand side and solution vector
 *  @param [in] A SkMatrix instance that contains matrix
 *  @param [in] b Vect instance that contains the right-hand side
 *  @param [in,out] x Vect instance that contains initial guess on input and solution on output
 */
    LinearSolver(SkSMatrix<real_t>&  A,
                 const Vect<real_t>& b,
                 Vect<real_t>&       x);

/** \brief Constructor using matrix, right-hand side
 *  @param [in] A SkMatrix instance that contains matrix
 *  @param [in] b Vect instance that contains the right-hand side
 *  @param [in,out] x Vect instance that contains the initial guess on input and solution on output
 */
    LinearSolver(SkMatrix<real_t>& A,
                 Vect<real_t>&     b,
                 Vect<real_t>&     x);

/// Destructor
    virtual ~LinearSolver() { }

/// \brief Set Maximum number of iterations
/// \details Default value is 1000
    void setMaxIter(int m) { _max_it = m; }

/// \brief Set tolerance value
    void setTolerance(real_t tol) { _toler = tol; }

/// \brief Set solution vector
    void setSolution(Vect<real_t>& x) { _x = &x; }

/// \brief Set right-hand side vector
    void setRHS(Vect<real_t>& b) { _b = &b; }

/// \brief Set matrix in the case of a pointer to Matrix
/// @param [in] A Pointer to abstract Matrix class
    void setMatrix(Matrix<real_t>* A);

/// \brief Set matrix in the case of a pointer to matrix
/// @param [in] A Pointer to abstract Matrix class
    void setMatrix(SpMatrix<real_t>& A);

/// \brief Set matrix in the case of a skyline matrix
/// @param [in] A %Matrix as instance of class SkMatrix
    void setMatrix(SkMatrix<real_t>& A);

/** \brief Set matrix, right-hand side and initial guess
 *  @param [in] A Reference to matrix as a SpMatrix instance
 *  @param [in] b Vector containing right-hand side
 *  @param [in,out] x Vector containing initial guess on input and solution on output
 */
    void set(SpMatrix<real_t>&   A,
             const Vect<real_t>& b,
             Vect<real_t>&       x);

/** \brief Set solver and preconditioner
 *  @param [in] s Solver identification parameter.
 *  To be chosen in the enumeration variable Iteration:\n 
 *  <tt>DIRECT_SOLVER</tt>, <tt>CG_SOLVER</tt>, <tt>CGS_SOLVER</tt>, <tt>BICG_SOLVER</tt>,
 *  <tt>BICG_STAB_SOLVER</tt>, <tt>GMRES_SOLVER</tt>, <tt>QMR_SOLVER</tt>
 *  @param [in] p Preconditioner identification parameter. By default, the diagonal
 *  preconditioner is used.
 *  To be chosen in the enumeration variable Preconditioner:\n
 *  <tt>IDENT_PREC</tt>, <tt>DIAG_PREC</tt>, <tt>SSOR_PREC</tt>, <tt>ILU_PREC</tt> [Default:
 *  <tt>ILU_PREC</tt>]
 *  @note The argument <tt>p</tt> has no effect if the solver is <tt>DIRECT_SOLVER</tt>
 */
    void setSolver(Iteration      s, 
                   Preconditioner p=DIAG_PREC) { _s = s; _p = p; }

/// \brief Return solver code
    Iteration getSolver() const { return _s; }

/// \brief Return solver preconditioner
    Preconditioner getPreconditioner() const { return _p; }

/** \brief Solve equations using system data, prescribed solver and preconditioner
 *  @param [in] A Reference to matrix as a SpMatrix instance
 *  @param [in] b Vector containing right-hand side
 *  @param [in,out] x Vector containing initial guess on input and solution on output
 *  @param [in] s Solver identification parameter
 *  To be chosen in the enumeration variable Iteration:\n 
 *  <tt>DIRECT_SOLVER</tt>, <tt>CG_SOLVER</tt>, <tt>CGS_SOLVER</tt>, <tt>BICG_SOLVER</tt>,
 *  <tt>BICG_STAB_SOLVER</tt>, <tt>GMRES_SOLVER</tt>, <tt>QMR_SOLVER</tt> [Default: <tt>CGS_SOLVER</tt>]
 *  @param [in] p Preconditioner identification parameter.
 *  To be chosen in the enumeration variable Preconditioner:\n
 *  <tt>IDENT_PREC</tt>, <tt>DIAG_PREC</tt>, <tt>SSOR_PREC</tt>, <tt>ILU_PREC</tt>,
 *  <tt>DILU_PREC</tt> [Default: <tt>DIAG_PREC</tt>]
 *  @remark The argument <tt>p</tt> has no effect if the solver is <tt>DIRECT_SOLVER</tt>
 *  @warning If the library <tt>eigen</tt> is used, only the preconditioners
 *  <tt>IDENT_PREC</tt>, <tt>DIAG_PREC</tt> and <tt>ILU_PREC</tt> are available.
 */
    int solve(SpMatrix<real_t>&   A,
              const Vect<real_t>& b,
              Vect<real_t>&       x,
              Iteration           s,
              Preconditioner      p=DIAG_PREC);

/** \brief Solve equations using prescribed solver and preconditioner
 *  @param [in] s Solver identification parameter
 *  To be chosen in the enumeration variable Iteration:\n 
 *  <tt>DIRECT_SOLVER</tt>, <tt>CG_SOLVER</tt>, <tt>CGS_SOLVER</tt>, <tt>BICG_SOLVER</tt>,
 *  <tt>BICG_STAB_SOLVER</tt>, <tt>GMRES_SOLVER</tt>, <tt>QMR_SOLVER</tt> [Default: <tt>CGS_SOLVER</tt>]
 *  @param [in] p Preconditioner identification parameter.
 *  To be chosen in the enumeration variable Preconditioner:\n
 *  <tt>IDENT_PREC</tt>, <tt>DIAG_PREC</tt>, <tt>SSOR_PREC</tt>, <tt>DILU_PREC</tt>,
 *  <tt>ILU_PREC</tt> [Default: <tt>DIAG_PREC</tt>]
 *  @note The argument <tt>p</tt> has no effect if the solver is <tt>DIRECT_SOLVER</tt>
 */
    int solve(Iteration      s,
              Preconditioner p=DIAG_PREC);

/** \brief Solve equations all arguments must have been given by other member functions
 *  \details Solver and preconditioner parameters must have been set by function setSolver.
 *  Otherwise, default values are set.
 */
    int solve();

/// Factorize matrix
    void setFact() { _fact = 1; }

/// Do not factorize matrix
    void setNoFact() { _fact = 0; }

/// Get number of performed iterations
    int getNbIter() const { return _nb_it; }

 private:

   int                _fact, _max_it, _nb_it, _matrix_set;
   Iteration          _s;
   Preconditioner     _p;
   real_t             _toler;
   Vect<real_t>       *_x;
   const Vect<real_t> *_b;
   Matrix<real_t>     *_A;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
