/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani

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

                   Definition of class 'EigenProblemSolver'

  ==============================================================================*/

#ifndef __EIGEN_SOLVER_H
#define __EIGEN_SOLVER_H

#include <stdlib.h>
#include <math.h>

#include <iostream>
using std::ostream;
using std::endl;
using std::cerr;

#include <iomanip>
using std::setw;

#include <string>
using std::string;

#include "OFELI_Config.h"
#include "linear_algebra/Assembly.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/DMatrix.h"
#include "linear_algebra/DSMatrix.h"
#include "linear_algebra/SkSMatrix.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */


/*! \file EigenProblemSolver.h
 *  \brief Definition file for class EigenProblemSolver.
 */

class Equa;

//-----------------------------------------------------------------------------
// Class EigenProblemSolver
//-----------------------------------------------------------------------------

/*! \class EigenProblemSolver
 *  \ingroup Solver
 *  \brief Class to find eigenvalues and corresponding eigenvectors of a given
 *  matrix in a generalized eigenproblem, <i>i.e.</i> Find scalars l and 
 *  non-null vectors v such that
 *       [K]{v} = l[M]{v}
 *  where [K] and [M] are symmetric matrices.
 *  The eigenproblem can be originated from a PDE. For this, we will refer
 *  to the matrices K and M as <i>Stiffness</i> and <i>Mass</i> matrices
 *  respectively. 
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

enum EigenMethod {
   SUBSPACE = 1, /*!< Subspace iteration method (for symmetric matrices)  */
   QR       = 2, /*!< QR reduction method (Martin, Parlett, Peters, Reinsch and Wilkinson.)                         */
};

class EigenProblemSolver
{

 public:

//----------------------------   BASIC OPERATIONS   ----------------------------

/// \brief Default constructor
    EigenProblemSolver();

/** \brief Constructor for a dense unsymmetric matrix that computes the eigenvalues.
 *  \details This constructor solves in place the eigenvalues problem and stores
 *  them in a vector.
 *  The eigenvectors can be obtained by calling the member function getEigenVector.
 *  Eigenvalue extraction uses the reduction to Hessenberg form of the matrix.
 *  @param [in] A Matrix for which eigenmodes are sought.
 *  @param [in] eigv Switch to compute or not eigenvectors [Default: <tt>false</tt>]
 *  @remark All eigenvalues are extracted.
 */
    EigenProblemSolver(DMatrix<real_t>& A,
                       bool             eigv=false);

/** \brief Constructor for a dense symmetric matrix that computes the eigenvalues.
 *  \details This constructor solves in place the eigenvalues problem and stores
 *  them in a vector.
 *  The eigenvectors can be obtained by calling the member function getEigenVector.
 *  @param [in] K Matrix for which eigenmodes are sought.
 *  @param [in] n Number of eigenvalues to extract. By default all eigenvalues are 
 *  computed.
 */
    EigenProblemSolver(DSMatrix<real_t>& K,
                       int               n=0);

/** \brief Constructor for Symmetric Skyline Matrices.
 *  @param [in] K "Stiffness" matrix
 *  @param [in] M "Mass" matrix
 *  @param [in] n Number of eigenvalues to extract. By default all eigenvalues are 
 *  computed.
 *  \note The generalized eigenvalue problem is defined by <tt>Kx = aMx</tt>,
 *  where <tt>K</tt> and <tt>M</tt> are referred to as stiffness and mass matrix.
 */
    EigenProblemSolver(SkSMatrix<real_t>& K,
                       SkSMatrix<real_t>& M,
                       int                n=0);

/** \brief Constructor for Symmetric Skyline Matrices.
 *  @param [in] K "Stiffness" matrix
 *  @param [in] M Diagonal "Mass" matrix stored as a Vect instance
 *  @param [in] n Number of eigenvalues to extract. By default all eigenvalues are 
 *  computed.
 *  \note The generalized eigenvalue problem is defined by <tt>Kx = aMx</tt>,
 *  where <tt>K</tt> and <tt>M</tt> are referred to as stiffness and mass matrix.
 */
    EigenProblemSolver(SkSMatrix<real_t>& K,
                       Vect<real_t>&      M,
                       int                n=0);

/** \brief Constructor for a dense symmetric matrix that compute the eigenvalues
 *  \details This constructor solves in place the eigenvalues problem and stores
 *  them in a vector. Extraction of eigenvalues and eigenvectors uses the subspace iteration method.
 *  The eigenvectors can be obtained by calling the member function getEigenVector.
 *  @param [in] A Matrix for which eigenmodes are sought.
 *  @param [in] ev Vector containing all computed (real) eigenvalues.
 *  @param [in] n Number of eigenvalues to extract. By default all eigenvalues are 
 *  computed.
 *  @remark The vector ev does not need to be sized before.
 */
    EigenProblemSolver(DSMatrix<real_t>& A,
                       Vect<real_t>&     ev,
                       int               n=0);

/** \brief Constructor for a dense matrix that compute the eigenvalues
 *  \details This constructor solves in place the eigenvalues problem and stores
 *  them separately in real and imaginary part vectors.
 *  Extraction of eigenvalues and eigenvectors uses the QR method after transforming the matrix
 *  in Hessenberg form.
 *  The eigenvectors can be obtained by calling the member function getEigenVector.
 *  @param [in] A Matrix for which eigenmodes are sought.
 *  @param [in] evr Vector containing real parts of eigenvalues.
 *  @param [in] evr Vector containing imaginary parts of eigenvalues.
 *  @param [in] eigv Switch to say if eigenvectors are to be computed [Default: <tt>false</tt>].
 *  @remark The vectors \c evr and \c evi do not need to be sized before.
 */
    EigenProblemSolver(DMatrix<real_t>& A,
                       Vect<real_t>&    evr,
                       Vect<real_t>&    evi,
                       bool             eigv=false);

/** \brief Consrtuctor using partial differential equation
 *  \details The used equation class must have been constructed using the Mesh instance
 *  @param [in] eq Reference to equation instance
 *  @param [in] lumped Mass matrix is lumped (\a true) or not (\a false) [Default: <tt>true</tt>]
 */
    EigenProblemSolver(Equa& eq,
                       bool  lumped=true);

/// \brief Destructor
    ~EigenProblemSolver();

//-------------------------------   MODIFIERS  ---------------------------------

/**
 * @brief Set numerical method for eigenvalue computation
 * 
 * @param [in] m Numerical method to compute eigenvalues:
 *   <ul>
 *     <li> SUBSPACE: Subspace iteration method (Valid for symmetric matrices only)
 *     <li> QR:       QR method (Reduction to Hessenberg form)
 *   </ul> 
 * @param [in] sym Boolean variable to say if matrix is symmetric or not
 */
   void set(EigenMethod m,
            bool        sym);

/** \brief Set matrix instances (Symmetric matrices).
 *  \details This function is to be used when the default constructor is applied.
 *  Case where the mass matrix is consistent.
 *  @param [in] K Stiffness matrix instance
 *  @param [in] M Mass matrix instance
 */
    void setMatrix(SkSMatrix<real_t>& K,
                   SkSMatrix<real_t>& M);

/** \brief Set matrix instances (Symmetric matrices).
 *  \details This function is to be used when the default constructor is applied.
 *  Case where the mass matrix is (lumped) diagonal and stored in a vector.
 *  @param [in] K Stiffness matrix instance
 *  @param [in] M Mass matrix instance where diagonal terms are stored as a vector.
 */
    void setMatrix(SkSMatrix<real_t>& K,
                   Vect<real_t>&      M);

/** \brief Set matrix instance (Symmetric matrix).
 *  \details This function is to be used when the default constructor is applied.
 *  Case of a standard (not generalized) eigen problem is to be solved
 *  @param [in] K Matrix instance
 */
    void setMatrix(DSMatrix<real_t>& K);

/** \brief Set matrix instance
 *  \details This function is to be used when the default constructor is applied.
 *  Case of a standard (not generalized) eigen problem is to be solved
 *  @param [in] K Matrix instance
 */
    void setMatrix(DMatrix<real_t>& K);

/** \brief Define partial differential equation to solve
 *  \details The used equation class must have been constructed using the Mesh instance
 *  @param [in] eq Reference to equation instance
 *  @param [in] lumped Mass matrix is lumped (\a true) or not (\a false) [Default: <tt>true</tt>]
 */
    void setPDE(Equa& eq,
                bool  lumped=true);

/** \brief Run the eigenproblem solver
 *  @param [in] nb Number of eigenvalues to be computed. By default, all eigenvalues are computed.
 *  @param [in] opt Toggle to compute or not corresponding eigenvectors. By default, no eigenvectors extracted.
 */
    int run(int  nb=0,
            bool opt=false);

/** \brief Assemble element arrays into global matrices
 *  \details This member function is to be called from finite element equation
 *  classes
 *  @param [in] el Reference to Element class
 *  @param [in] eK Pointer to element stiffness (or assimilated) matrix
 *  @param [in] eM Pointer to element mass (or assimilated) matrix
 */
    void Assembly(const Element& el,
                  real_t*        eK,
                  real_t*        eM);

/** \brief Assemble side arrays into global matrix and right-hand side
 *  \details This member function is to be called from finite element equation
 *  classes
 *  @param [in] sd Reference to Side class
 *  @param [in] sK Pointer to side stiffness
 */
    void SAssembly(const Side& sd,
                   real_t*     sK);

/** \brief Run the subspace iteration solver.
 *  \details This function rune the Bathe subspace iteration method.
 *  @param [in] nb_eigv Number of eigenvalues to be extracted
 *  @param [in] ss_dim Subspace dimension. Must be at least equal to the number
 *  eigenvalues to seek. [Default: <tt>nb_eigv</tt>]
 *  @return   1: Normal execution. Convergence has been achieved. 
 *            2: Convergence for eigenvalues has not been attained. 
 */ 
    int runSubSpace(size_t nb_eigv,
                    size_t ss_dim=0);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int runJacobi(DMatrix<real_t>& wm);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Define the subspace dimension.
 *  @param [in] dim Subspace dimension. Must be larger or equal to the
 *  number of wanted eigenvalues.
 *  By default this value will be set to the number of wanted eigenvalues
 */
    void setSubspaceDimension(int dim) { _dim = dim; }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief set
    void setSubspaceOption(int opt) { _opt = opt; }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief set maximal number of iterations.
    void setMaxIter(int max_it) { _max_it = max_it; }

/// \brief set tolerance value
/// @param [in] eps Convergence tolerance for eigenvalues [Default: 1.e-8]
    void setTolerance(real_t eps) { _epsj = _epsv = eps; }

//-----------------------------   INSPECTORS  ----------------------------------

/** \brief Check how many eigenvalues have been found using Sturm sequence method
 *  @param [out] nb_found number of eigenvalues actually found
 *  @param [out] nb_lost  number of eigenvalues missing
 *  @return   - 0, Successful completion of subroutine. 
 *            - 1, No convergent eigenvalues found. 
 */
    int checkSturm(int& nb_found,
                   int& nb_lost);

/// \brief Return actual number of performed iterations
    int getNbIter() const { return int(_max_it); }

/** \brief Return the n-th eigenvalue
 *  @param [in] n Index of eigenvalue (must be > 0)
 *  @param [in] i i=1: Real part, i=2: Imaginary part. [Default: 1]
 *  @return Eigenvalue
 */
    real_t getEigenValue(int n,
			 int i=1) const;

/** \brief Return the n-th eigenvector
 *  @param [in] n Label of eigenvector (They are stored in ascending order of eigenvalues)
 *  @param [out] v Vect instance where the eigenvector is stored.
 */
    void getEigenVector(int           n,
                        Vect<real_t>& v) const;

/** \brief Return the n-th eigenvector in the case of complex eigenvectors
 *  @param [in] n Label of eigenvector (They are stored in the order of eigenvalues)
 *  @param [out] vr Vect instance where the real part of the eigenvector is stored.
 *  @param [out] vi Vect instance where the imaginary part of the eigenvector is stored.
 */
    void getEigenVector(int           n,
                        Vect<real_t>& vr,
                        Vect<real_t>& vi) const;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    friend ostream & operator<<(ostream&                  s,
                                const EigenProblemSolver& es);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 private:

   Equa             *_theEqua;
   Mesh             *_theMesh;
   Matrix<real_t>   *_K, *_M;
   Vect<real_t>     *_lM, _rc, _evr, _evi, _evec;
   DSMatrix<real_t> _pK, _pM;
   DMatrix<real_t>  _ev;
   bool             _K_alloc, _M_alloc, _lM_alloc, _diag;
   size_t           _max_it, _dim, _nb_eq, _nb_eigv;
   int              _opt, _low, _high;
   real_t           _epsv, _epsj;
   EigenMethod      _method;
   bool             _sym, _eigv;

   int init();
   void Mxv(const Vect<real_t>& b, Vect<real_t>& c);

/*  Balance the matrix so that the rows with zero entries
 *  off the diagonal are isolated and the remaining columns and rows
 *  are resized to have one norm close to 1.                          
 */
   void balance(Vect<real_t>& scal);
 
// Reverse the balancing of balance for the eigenvactors
   void balance_back(const Vect<real_t>& scal);

// Transform the matrix to upper Hessenberg form
   void elmhes(Vect<int>& perm);

// Copy the Hessenberg matrix stored in mat to h
   void elmtrans(const Vect<int>& perm);

/* Reduce matrix to upper Hessenberg form by Householder transformations. 
 * All details of the transformations are stored in the remaining triangle 
 * of the Hessenberg matrix and in vector d.
 */
   void orthes(Vect<real_t>& d);

/* Compute the matrix v of accumulated transformations from the
 * information left by the Householder reduction of matrix mat to upper
 * Hessenberg form below the Hessenberg matrix in A and in the
 * vector d. The contents of the latter are destroyed.
*/
   void orttrans(Vect<real_t>& d);

// Compute the eigenvectors for the eigenvalues found in hqr2
   int hqrvec();

// Compute the eigenvalues and (if opt==true) the eigenvectors of an
// upper Hessenberg matrix
   int hqr2(Vect<int>& cnt);

// norm_1 normalizes the one norm of the column vectors in v.
// (special attention to complex vectors in v  is given)
   void norm_2();

   int eigen(int ortho, int ev_norm);

};

/// \fn ostream & operator<<(ostream& s, const EigenProblemSolver &es)
/// \brief Output eigenproblem information
/// \ingroup Solver
    ostream & operator<<(ostream&                  s,
                         const EigenProblemSolver& es);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
