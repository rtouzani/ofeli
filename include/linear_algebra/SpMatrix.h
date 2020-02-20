/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

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

                  Definition of class 'SpMatrix' for Sparse matrix

  ==============================================================================*/


#ifndef __SPMATRIX_H
#define __SPMATRIX_H

#include "linear_algebra/Matrix.h"
#ifdef USE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Sparse>
using namespace Eigen;
#endif

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file SpMatrix.h
 *  \brief Definition file for class SpMatrix.
 */

using std::pair;
using std::unique;

/*! \class SpMatrix
 *  \ingroup VectMat
 *  \brief To handle matrices in sparse storage format.
 *
 * This template class enables storing and manipulating a sparse matrix,
 * i.e. only nonzero terms are stored. Internally, the matrix is stored as a
 * vector instance and uses for the definition of its graph a <tt>Vect<size_t></tt>
 * instance row_ptr and a Vect<size_t> instance <tt>col_ind</tt> that contains
 * respectively addresses of first element of each row and column indices.
 *
 * To illustrate this, consider the matrix
 * \verbatim
            1   2   0
            3   4   0
            0   5   0
   \endverbatim
 *
 * Such a matrix is stored in the vector<real_t> instance {1,2,3,4,5}.
 * The vectors <tt>row_ptr</tt> and <tt>col_ind</tt> are respectively:
 * <tt>{0,2,4,5}</tt>, <tt>{1,2,1,2,2}</tt>
 *
 * When the library <tt>eigen</tt> is used in conjunction with <tt>OFELI</tt>,
 * the class uses the sparse matrix class of <tt>eigen</tt> and enables then 
 * access to specific solvers (see class LinearSolver)
 *
 * \tparam T_ Data type (double, float, complex<double>, ...)
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */

template<class T_> class SpMatrix : public Matrix<T_>
{
    using Matrix<T_>::_nb_rows;
    using Matrix<T_>::_nb_cols;
    using Matrix<T_>::_size;
    using Matrix<T_>::_length;
    using Matrix<T_>::_zero;
    using Matrix<T_>::_temp;
    using Matrix<T_>::_a;
    using Matrix<T_>::_diag;
    using Matrix<T_>::_ch;
    using Matrix<T_>::_dof_type;
    using Matrix<T_>::_is_diagonal;
    using Matrix<T_>::_theMesh;
    using Matrix<T_>::_row_ptr;
    using Matrix<T_>::_col_ind;

#ifdef USE_EIGEN
    typedef Eigen::Matrix<T_,Eigen::Dynamic,1> VectorX;
    typedef SparseMatrix<T_>                   SpMat;
    typedef Triplet<real_t>                    Tr;
#endif

 public:

/// \brief Default constructor.
/// \details Initialize a zero-dimension matrix
    SpMatrix();

/** \brief Constructor that initializes current instance as a dense matrix.
 *  \details Normally, for a dense matrix this is not the right class.
 *  @param [in] nr Number of matrix rows.
 *  @param [in] nc Number of matrix columns.
 */
    SpMatrix(size_t nr,
             size_t nc);

/** \brief Constructor that initializes current instance as a dense matrix.
 *  \details Normally, for a dense matrix this is not the right class.
 *  @param [in] size Number of matrix rows (and columns).
 *  @param [in] is_diagonal Boolean argument to say is the matrix is actually a diagonal matrix or not.
 */
    SpMatrix(size_t size,
             int    is_diagonal=false);

/** \brief Constructor using a Mesh instance.
 *  @param [in] mesh Mesh instance from which matrix graph is extracted.
 *  @param [in] dof Option parameter, with default value <tt>0</tt>.\n
 *  <tt>dof=1</tt> means that only one degree of freedom for each node (or element or side)
 *  is taken to determine matrix structure. The value <tt>dof=0</tt> means that matrix 
 *  structure is determined using all DOFs.
 *  @param [in] is_diagonal Boolean argument to say is the matrix is actually a diagonal matrix or not.
 */
    SpMatrix(Mesh&  mesh,
             size_t dof=0,
             int    is_diagonal=false);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#ifndef USE_EIGEN
 /** \brief Constructor using a Mesh instance and selecting a DOF (For use with Eigen Library).
 *  @param [in] dof Option parameter.\n
 *  <tt>dof=1</tt> means that only one degree of freedom for each node (or element or side)
 *  is taken to determine matrix structure. The value <tt>dof=0</tt> means that matrix 
 *  structure is determined using all DOFs.
 *  @param [in] mesh Mesh instance from which matrix graph is extracted.
 *  @param [in] code.
 */
   SpMatrix(size_t dof,
             Mesh&  mesh,
             int    code=0);

/** \brief Constructor using a Mesh instance.
 *  @param [in] mesh Mesh instance from which matrix graph is extracted.
 *  @param [in] dof Option parameter, with default value <tt>0</tt>.\n
 *  <tt>dof=1</tt> means that only one degree of freedom for each node (or element or side)
 *  is taken to determine matrix structure. The value <tt>dof=0</tt> means that matrix 
 *  structure is determined using all DOFs.
 *  @param [in] is_diagonal Boolean argument to say is the matrix is actually a diagonal matrix or not.
 */
    SpMatrix(size_t dof,
             size_t nb_eq,
             Mesh&  mesh);
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef USE_EIGEN
/** \brief Constructor for a square matrix using non zero row and column indices
 *  @param [in] I Vector containing pairs of row and column indices
 *  @param [in] opt Flag indicating if vectors I is cleaned and ordered
 *  (opt=1) or not (opt=0). In the latter case, this vector can have
 *  the same contents more than once and are not necessarily ordered
 */
    SpMatrix(const Vect<RC>& I,
             int             opt=1);
#endif

#ifndef USE_EIGEN
/** \brief Constructor for a square matrix using non zero row and column indices
 *  @param [in] I Vector containing pairs of row and column indices
 *  @param [in] a Vector containing matrix entries in the same order than
 *  the one given by <tt>I</tt>
 *  @param [in] opt Flag indicating if vector <tt>I</tt> is cleaned and ordered
 *  (<tt>opt=1</tt>: default) or not (<tt>opt=0</tt>). In the latter case, this vector can have
 *  the same contents more than once and are not necessarily ordered
 */
    SpMatrix(const Vect<RC>& I,
             const Vect<T_>& a,
             int             opt=1);
#endif

#ifndef USE_EIGEN
/** \brief Constructor for a rectangle matrix.
 *  @param [in] nr Number of rows
 *  @param [in] nc Number of columns
 *  @param [in] row_ptr Vector of row pointers (See the above description of this class).
 *  @param [in] col_ind Vector of column indices (See the above description of this class).\n
 */
    SpMatrix(size_t                nr,
             size_t                nc,
             const vector<size_t>& row_ptr,
             const vector<size_t>& col_ind);
#endif

#ifndef USE_EIGEN
/** \brief Constructor for a rectangle matrix
 *  @param [in] nr Number of rows
 *  @param [in] nc Number of columns
 *  @param [in] row_ptr Vector of row pointers (See the above description of this class).
 *  @param [in] col_ind Vector of column indices (See the above description of this class).\n
 *  @param [in] a vector instance containing matrix entries stored columnwise
 */
    SpMatrix(size_t                nr,
             size_t                nc,
             const vector<size_t>& row_ptr,
             const vector<size_t>& col_ind,
             const vector<T_>&     a);
#endif

#ifndef USE_EIGEN
/** \brief Constructor for a rectangle matrix
 *  @param [in] row_ptr Vector of row pointers (See the above description of this class).
 *  @param [in] col_ind Vector of column indices (See the above description of this class).
 */
    SpMatrix(const vector<size_t>& row_ptr,
             const vector<size_t>& col_ind);
#endif

#ifndef USE_EIGEN
/** \brief Constructor for a rectangle matrix
 *  @param [in] row_ptr Vector of row pointers (See the above description of this class).
 *  @param [in] col_ind Vector of column indices (See the above description of this class).
 *  @param [in] a vector instance that contain matrix entries stored row by row.\n
 *  Number of rows is extracted from vector <tt>row_ptr</tt>.
 */
    SpMatrix(const vector<size_t>& row_ptr,
             const vector<size_t>& col_ind,
             const vector<T_>&     a);
#endif

/// \brief Copy constructor.
    SpMatrix(const SpMatrix& m);

/// \brief Destructor.
    ~SpMatrix();

/// \brief Define matrix as identity
    void Identity();

/// \brief Define matrix as a dense one
    void Dense();

/// \brief Define matrix as a diagonal one
    void Diagonal();

/// \brief Define matrix as a diagonal one
/// with diagonal entries equal to <tt>a</tt>
    void Diagonal(const T_& a);

/** \brief Sets the matrix as the one for the Laplace equation in 1-D
 *  \details The matrix is initialized as the one resulting from P<sub>1</sub>
 *  finite element discretization of the classical elliptic operator
 *     -u'' = f
 *  with homogeneous Dirichlet boundary conditions
 *  \remark This function is available for real valued matrices only.
 *  @param [in] n Size of matrix (Number of rows)
 *  @param [in] h %Mesh size (assumed constant)
 */
    void Laplace1D(size_t n,
                   real_t h);

/** \brief Sets the matrix as the one for the Laplace equation in 2-D
 *  \details The matrix is initialized as the one resulting from P<sub>1</sub>
 *  finite element discretization of the classical elliptic operator
 *     -Delta u = f
 *  with homogeneous Dirichlet boundary conditions
 *  \remark This function is available for real valued matrices only.
 *  @param [in] nx Number of unknowns in the <tt>x</tt>-direction
 *  @param [in] ny Number of unknowns in the <tt>y</tt>-direction
 *  @remark The number of rows is equal to <tt>nx*ny</tt>
 */
    void Laplace2D(size_t nx,
                   size_t ny);

/** \brief Determine mesh graph and initialize matrix.
 *  \details This member function is called by constructor with the same arguments
 *  @param [in] mesh Mesh instance for which matrix graph is determined.
 *  @param [in] dof Option parameter, with default value <tt>0</tt>.\n
 *  <tt>dof=1</tt> means that only one degree of freedom for each node (or element or side)
 *  is taken to determine matrix structure. The value <tt>dof=0</tt> means that matrix
 *  structure is determined using all DOFs.
 */
    void setMesh(Mesh&  mesh,
                 size_t dof=0);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setMesh(size_t dof,
                 Mesh&  mesh,
                 int    code=0);

    void setMesh(size_t dof, 
                 size_t nb_eq,
                 Mesh&  mesh);

    void setMesh(Mesh&  mesh,
                 size_t dof,
                 size_t type);

/** \brief Activate extended graph option.
 *  \details An extended graph is the one for which a node is linked not only to
 *  its neighbor nodes but also to the neighbors of the neighbors.
 */
    void setExtendedGraph();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Activate 1-DOF per node option.
    void setOneDOF();

/// \brief Activate Sides option.
    void setSides();

/// \brief Store diagonal entries in a separate internal vector.
    void setDiag();

/** \brief Impose by a diagonal method an essential boundary condition.
 *  \details This member function modifies diagonal terms in matrix and terms
 *  in vector that correspond to degrees of freedom with nonzero code
 *  in order to impose a boundary condition.
 *  The penalty parameter is defined by default equal to 1.e20.
 *  It can be modified by member function \b setPenal(..).
 *  @param [in] mesh Mesh instance from which information is extracted.
 *  @param [in,out] b Vect instance that contains right-hand side.
 *  @param [in] u Vect instance that conatins imposed valued at DOFs where they are to be imposed.
 */
    void DiagPrescribe(Mesh&           mesh,
                       Vect<T_>&       b,
                       const Vect<T_>& u);

/** \brief Impose by a diagonal method an essential boundary condition using the Mesh instance
 *  provided by the constructor
 *  \details This member function modifies diagonal terms in matrix and terms
 *  in vector that correspond to degrees of freedom with nonzero code
 *  in order to impose a boundary condition.
 *  The penalty parameter is defined by default equal to 1.e20.
 *  It can be modified by member function \b setPenal(..).
 *  @param [in,out] b Vect instance that contains right-hand side.
 *  @param [in] u Vect instance that conatins imposed valued at DOFs where they are to be imposed.
 */
    void DiagPrescribe(Vect<T_>&       b,
                       const Vect<T_>& u);

/// \brief Set size of matrix (case where it's a square matrix).
/// @param [in] size Number of rows and columns.
    void setSize(size_t size);

/** \brief Set size (number of rows) of matrix
 *  @param [in] nr Number of rows
 *  @param [in] nc Number of columns
 */
    void setSize(size_t nr,
                 size_t nc);

/** \brief Set graph of matrix by giving a vector of its nonzero entries
 *  @param [in] I Vector containing pairs of row and column indices
 *  @param [in] opt Flag indicating if vector <tt>I</tt> is cleaned and ordered
 *  (<tt>opt=1</tt>: default) or not (<tt>opt=0</tt>). In the latter case, this vector can have
 *  the same contents more than once and are not necessarily ordered
 */
    void setGraph(const Vect<RC>& I,
                  int             opt=1);

/// \brief Get <tt>i</tt>-th row vector.
    Vect<T_> getRow(size_t i) const;

/// \brief Get <tt>j</tt>-th column vector.
    Vect<T_> getColumn(size_t j) const;

/** \brief Operator () (Non constant version)
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    T_ & operator()(size_t i,
                    size_t j);

/** \brief Operator <tt>()</tt> (Constant version)
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    T_ operator()(size_t i,
                  size_t j) const;

/** \brief Operator <tt>()</tt> with one argument (Constant version)
 *  \details Returns <tt>i</tt>-th position in the array storing matrix entries.
 *  The first entry is at location <tt>1</tt>.
 *  Entries are stored row by row.
 */
    T_ operator()(size_t i) const;

/** \brief Operator <tt>[]</tt> (Constant version).
 *  \details Returns <tt>i</tt>-th position in the array storing matrix entries.
 *  The first entry is at location <tt>0</tt>.
 *  Entries are stored row by row.
 */
    T_ operator[](size_t i) const;

/** \brief Operator <tt>*</tt> to multiply matrix by a vector
 *  @param [in] x Vect instance to multiply by
 *  @return Vector product of matrix by <tt>x</tt>
 */
    Vect<T_> operator*(const Vect<T_>& x) const;

/** \brief Operator <tt>*=</tt> to premultiply matrix by a constant
 *  @param [in] a Constant to multiply matrix by
 *  @return Resulting matrix
 */
    SpMatrix<T_>& operator*=(const T_& a);

/// \brief Get mesh instance whose reference will be stored in current instance of SpMatrix.
    void getMesh(Mesh& mesh);

/** \brief Multiply matrix by vector and save in another one.
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector that contains on output the result.
 */
    void Mult(const Vect<T_>& x,
              Vect<T_>&       y) const;

/** \brief Multiply matrix by vector <tt>x</tt> and add to <tt>y</tt>.
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector to add to the result. <tt>y</tt> contains on output the result.
 */
    void MultAdd(const Vect<T_>& x,
                 Vect<T_>&       y) const;

/** \brief Multiply matrix by vector <tt>a*x</tt> and add to <tt>y</tt>.
 *  @param [in] a Constant to multiply by matrix
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector to add to the result. <tt>y</tt> contains on output the result.
 */
    void MultAdd(T_              a,
                 const Vect<T_>& x,
                 Vect<T_>&       y) const;

/** \brief Multiply transpose of matrix by vector <tt>x</tt> and save in <tt>y</tt>.
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector that contains on output the result.
 */
    void TMult(const Vect<T_>& x,
               Vect<T_>&       y) const;

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                  a,
              const SpMatrix<T_>& m);

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m Pointer to Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                a,
              const Matrix<T_>* m);

/** \brief Assign a value to an entry of the matrix
 *  @param [in] i Row index
 *  @param [in] j Column index
 *  @param [in] val Value to assign to <tt>a(i,j)</tt>
 */
    void set(size_t    i,
             size_t    j,
             const T_& val);

/** \brief Add a value to an entry of the matrix
 *  @param [in] i Row index
 *  @param [in] j Column index
 *  @param [in] val Constant value to add to <tt>a(i,j)</tt>
 */
    void add(size_t    i,
             size_t    j,
             const T_& val);

/// \brief Operator =.
/// \details Assign constant value <tt>x</tt> to all matrix entries.
    void operator=(const T_& x);

/// \brief Return storage information.
/// @return Column index of the <tt>i</tt>-th stored element in matrix
    size_t getColInd(size_t i) const;

/// \brief Return Row pointer at position <tt>i</tt>.
    size_t getRowPtr(size_t i) const;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Factorize matrix
    int Factor();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** \brief Perform a diagonal incomplete LU factorization of matrix
 *  @param [out] id
 *  @param [out] pivot
 *  \return Return <tt>0</tt> if the factorization was normally achieved,
 *  <tt>n</tt> if the <tt>n</tt>-th pivot is null.
 */
#ifdef USE_EIGEN
    int DILUFactorize(vector<size_t>& id,
                      vector<T_>&     pivot) const;
#else
    int DILUFactorize(vector<size_t>& id,
                      vector<T_>&     pivot) const;
#endif

/** \brief Perform an Incomplete LU factorization of matrix
 *  \return Return <tt>0</tt> if the factorization was normally achieved,
 *  <tt>n</tt> if the <tt>n</tt>-th pivot is null.
 */
#ifdef USE_EIGEN
    int ILUFactorize(const SpMatrix<T_>& A);
#else
    int ILUFactorize(const SpMatrix<T_>& A);
#endif

/** \brief Solve a linear system with an diagonal incompletely factorized matrix
 *  @param [in] b Vect instance containing the right-hand side
 *  @param [out] x Vect instance containing on output the solution
 */
#ifndef USE_EIGEN
    void DILUSolve(const vector<size_t>& id,
                   const vector<T_>&     pivot,
                   const Vect<T_>&       b,
                   Vect<T_>&             x) const;
#endif

/** \brief Solve a linear system with an incompletely factorized matrix
 *  @param [in] b Vect instance containing the right-hand side
 *  @param [out] x Vect instance containing on output the solution
 */
#ifndef USE_EIGEN
    void ILUSolve(const Vect<T_>& b,
                  Vect<T_>&       x) const;
#endif

/** \brief Solve a linear system with an incompletely factorized matrix
 *  @param [in] b Vect instance containing the right-hand side
 *  @param [out] x Vect instance containing on output the solution
 */
#ifndef USE_EIGEN
    void SSORSolve(const Vect<T_>& b,
                   Vect<T_>&       x) const;
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int solve(Vect<T_>& b,
              bool      fact=false);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/** \brief Solve the linear system of equations.
 *  \details The default parameters are:
 *  <ul>
 *     <li><tt>CG_SOLVER</tt> for solver
 *     <li><tt>DIAG_PREC</tt> for preconditioner
 *     <li>Max. Number of iterations is <tt>1000</tt>
 *     <li>Tolerance is <tt>1.e-8</tt>
 *  </ul>
 *  To change these values, call function setSolver before this function
 *  @param [in] b Vector that contains right-hand side
 *  @param [out] x Vector that contains the obtained solution
 *  @param [in] fact Unused argument
 *  @return Number of actual performed iterations 
 */
    int solve(const Vect<T_>& b,
              Vect<T_>&       x,
              bool            fact=false);

/** \brief Choose solver and preconditioner for an iterative procedure
 *  @param [in] solver Option to choose iterative solver in an enumerated variable
 *  <ul>
 *     <li><tt>CG_SOLVER</tt>: Conjugate Gradient [default]
 *     <li><tt>CGS_SOLVER</tt>: Squared conjugate gradient
 *     <li><tt>BICG_SOLVER</tt>: Biconjugate gradient
 *     <li><tt>BICG_STAB_SOLVER</tt>: Biconjugate gradient stabilized
 *     <li><tt>GMRES_SOLVER</tt>: Generalized Minimal Residual
 *  </ul>
 *  Default value is <tt>CG_SOLVER</tt>
 *  @param [in] prec Option to choose preconditioner in an enumerated variable
 *  <ul>
 *     <li><tt>IDENT_PREC</tt>: Identity preconditioner (no preconditioning)
 *     <li><tt>DIAG_PREC</tt>: Diagonal preconditioner [default]
 *     <li><tt>SSOR_PREC</tt>: SSOR (Symmetric Successive Over Relaxation) preconditioner
 *     <li><tt>DILU_PREC</tt>: ILU (Diagonal Incomplete factorization) preconditioner
 *     <li><tt>ILU_PREC</tt>: ILU (Incomplete factorization) preconditioner
 *  </ul>
 *  Default value is <tt>DIAG_PREC</tt>
 *  @param [in] max_it Maximum number of allowed iterations. Default value is <tt>1000</tt>.
 *  @param [in] toler Tolerance for convergence. Default value is <tt>1.e-8</tt>
 */
    void setSolver(Iteration      solver=CG_SOLVER,
                   Preconditioner prec=DIAG_PREC,
                   int            max_it=1000,
                   real_t         toler=1.e-8);

/// brief Set all matrix entries to zero
    void clear();

/// \brief Return C-Array.
/// \details Non zero terms of matrix is stored row by row.
    T_ *get() const;

/** \brief  Return entry <tt>(i,j)</tt> of matrix if this one is stored, <tt>0</tt> otherwise
 *  @param [in] i Row index (Starting from 1)
 *  @param [in] j Column index (Starting from 1)
 */
    T_ get(size_t i,
           size_t j) const;

#ifdef USE_EIGEN
/// \brief Return reference to the matrix instance in Eigen library
    SpMat& getEigenMatrix();
#endif

#ifdef USE_EIGEN
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
#endif

/** \fn ostream & operator<<(ostream& s, const SpMatrix<T_> &A)
 *  \ingroup VectMat
 *  Output matrix in output stream
 */
    template<class TT_>
    friend ostream& operator<<(ostream& s, const SpMatrix<TT_>& A);

 private:

   size_t _dof;
   int _is_dense, _type, _max_it;
   int _one_dof, _sides, _extended;
   real_t _toler;
   vector<RC> _IJ;
   vector<size_t> _nbc;
   Iteration _solver;
   Preconditioner _prec;
#ifdef USE_EIGEN
   SpMat _A;
#else
   vector<T_> _pivot, _aL, _aU;
   vector<size_t> _id, _l_row_ptr, _u_row_ptr, _l_col_ind, _u_col_ind;
   size_t _lnnz, _unnz;
   int _col_index(size_t i, size_t j) const;
#endif
};

/** \fn Vect<T_> operator*(const SpMatrix<T_> &A, const Vect<T_> &b)
 *  \brief Operator * (Multiply vector by matrix and return resulting vector
 *  \ingroup VectMat
 *  @param [in] A SpMatrix instance to multiply by vector
 *  @param [in] b Vect instance 
 *  \return Vect instance containing <tt>A*b</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_>
Vect<T_> operator*(const SpMatrix<T_>& A,
                   const Vect<T_>&     b);


/** \fn ostream & operator<<(ostream& s, const SpMatrix<T_> &A)
 *  \ingroup VectMat
 *  \brief Output matrix in output stream
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_>
ostream& operator<<(ostream&            s,
                    const SpMatrix<T_>& A);

} /* namespace OFELI */

#endif

