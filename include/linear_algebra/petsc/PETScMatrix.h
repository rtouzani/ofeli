/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2026 Rachid Touzani

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

          Definition and implementation of template class 'PETScMatrix' 
                   for matrices using the PETSc library

  ==============================================================================*/


#ifndef __PETSC_MATRIX_H
#define __PETSC_MATRIX_H

#if defined(USE_PETSC)

#include <petscmat.h>
#include <petsc.h>

#include <vector>
using std::vector;

#include "mesh/Partition.h"
#include "linear_algebra/GraphOfMatrix.h"
#include "linear_algebra/petsc/PETScVect.h"
#include "util/util.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file PETScMatrix.h
 *  \brief Definition file for class PETScMatrix.
 */

/*! \class PETScMatrix
 *  \ingroup VectMat
 *  \brief To handle matrices in sparse storage format using the Petsc library.
 *
 * @warning This class is available only when OFELI has been installed with Petsc.
 *
 * \tparam T_ Data type (double, float, complex<double>, ...)
 */

template<class T_> class PETScVect;
class Mesh;
class Element;
class Side;
class Partition;
 
template<class T_> class PETScMatrix
{

 public:

/// \brief Default constructor.
/// \details Initialize a zero-dimension matrix
    PETScMatrix();

/** \brief Constructor that initializes current instance as a dense matrix.
 *  \details Normally, for a dense matrix this is not the right class.
 *  @param [in] nr Number of matrix rows.
 *  @param [in] nc Number of matrix columns.
 */
    PETScMatrix(size_t nr,
                size_t nc);

/** \brief Constructor that initializes current instance as a dense matrix.
 *  \details Normally, for a dense matrix this is not the right class.
 *  @param [in] size Number of matrix rows (and columns).
 */
    PETScMatrix(size_t size);

/** \brief Constructor using a Mesh instance
 *  @param [in] mesh Mesh instance from which matrix graph is extracted.
 *  @param [in] dof Option parameter, with default value <tt>0</tt>.\n
 *  <tt>dof=1</tt> means that only one degree of freedom for each node (or element or side)
 *  is taken to determine matrix structure. The value <tt>dof=0</tt> means that matrix 
 *  structure is determined using all DOFs.
 */
    PETScMatrix(Mesh&  mesh,
                size_t dof=0);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    PETScMatrix(size_t dof,
                Mesh&  mesh,
                int    code=0);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    PETScMatrix(size_t dof,
                size_t nb_eq,
                Mesh&  mesh);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Constructor for a square matrix using non zero row and column indices
 *  @param [in] I Vector containing pairs of row and column indices
 *  @param [in] opt Flag indicating if vectors I is cleaned and ordered
 *  (opt=1) or not (opt=0). In the latter case, this vector can have
 *  the same contents more than once and are not necessarily ordered
 */
    PETScMatrix(const vector<std::pair<size_t,size_t> >& I,
                      int                                opt=1);

/// \brief Copy constructor
    PETScMatrix(const PETScMatrix& m);

/// \brief Destructor
    ~PETScMatrix(void);

/// \brief Define matrix as identity matrix
    void Identity();

/// \brief Define matrix as a diagonal one
    void Diagonal();

/// \brief Define matrix as a diagonal one
/// with diagonal entries equal to <tt>a</tt>
    void Diagonal(const T_& a);

/** \brief
 *  @param [in] nnz
 */
    void setAIJ(const vector<int>& nnz);

/** \brief
 *  @param [in] diag_nnz
 *  @param [in] off_nnz
 */
    void setAIJ_MPI(const vector<int>& diag_nnz,
                    const vector<int>& off_nnz);

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
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setMesh(size_t dof, 
                 size_t nb_eq,
                 Mesh&  mesh);
    void setMesh(Mesh&  mesh,
                 size_t dof,
                 size_t type);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** \brief Activate extended graph option.
 *  \details An extended graph is the one
 *  for which a node is linked not only to its neighbor nodes but also to the
 *  neighbors of the neighbors.
 */
    void ExtendedGraph();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Set a Partition instance in the class
 *  \details This member function is to be used when parallel computing is considered.
 *  @param [in] p Reference to Partition instance
 */
    void setPartition(Partition& p);

/** \brief Set number of processors and processor rank
 *  @param [in] np Total number of processors.
 *  @param [in] r  Rank of current processor [Default: <tt>0</tt>
 *  @warning If this member function is not called, only one processor is used and then
 *  sequential computing is involved.
 */
    void setRank(int np,
                 int r=0) { _np = np; _rank = r; }

/// \brief Activate 1-DOF per node option.
    void setOneDOF();

/// \brief Activate Sides option.
    void setSides();

/// \brief Set matrix as a symmetric one
    void setSymmetric() { MatSetOption(_A,MAT_SYMMETRIC,PETSC_TRUE); }

/** \brief Impose by a diagonal method an essential boundary condition using the Mesh instance
 *  provided by the constructor
 *  \details This member function modifies diagonal terms in matrix and terms
 *  in vector that correspond to degrees of freedom with nonzero code
 *  in order to impose a boundary condition.
 *  The penalty parameter is defined by default equal to 1.e20.
 *  It can be modified by member function \b setPenal(..).
 *  @param [in,out] b PETScVect instance that contains right-hand side.
 *  @param [in] u PETScVect instance that conatins imposed valued at DOFs where they are to be imposed.
 */
    void DiagPrescribe(PETScVect<T_>&       b,
                       const PETScVect<T_>& u);

/// \brief Set size of matrix (case where it's a square matrix).
/// @param [in] size Number of rows and columns.
    void setSize(size_t size);

/** \brief Set size (number of rows) of matrix
 *  @param [in] nr Number of rows
 *  @param [in] nc Number of columns
 */
    void setSize(size_t nr,
                 size_t nc);

/** \brief Return the range of matrix rows owned by this processor
 *  @param [out] istart Index of the first local row
 *  @param [out] iend Index of the last local row
 */
    void getRange(int istart,
                  int iend)
    {
       MatGetOwnershipRange(_A,&_istart,&_iend);
       istart = _istart;
       iend = _iend;
    }
 
/** \brief Set graph of matrix by giving a vector of its nonzero entries
 *  @param [in] I Vector containing pairs of row and column indices
 *  @param [in] opt Flag indicating if vector <tt>I</tt> is cleaned and ordered
 *  (<tt>opt=1</tt>: default) or not (<tt>opt=0</tt>). In the latter case, this vector can have
 *  the same contents more than once and are not necessarily ordered
 */
    void setGraph(const vector<std::pair<size_t,size_t> >& I,
                        int                                opt=1);

/** \brief Operator <tt>()</tt>
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    T_ operator()(size_t i,
                  size_t j) const;

/// \brief Return number of matrix rows
    size_t getNbRows() const { return _nb_rows; }

/// \brief Return number of matrix columns
    size_t getNbColumns() const { return _nb_cols; }

/// \brief Return length of matrix
/// \details The length is the total number of stored elements in the matrix
    size_t getLength() const { return _length; }

/// \brief Get Mesh instance whose reference will be stored in current instance of PETScMatrix.
    void getMesh(Mesh& mesh);

/** \brief Multiply matrix by vector and save in another one.
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector that contains on output the result.
 */
    void Mult(const PETScVect<T_>& x,
                    PETScVect<T_>& y) const;

/** \brief Multiply matrix by vector <tt>x</tt> and add to <tt>y</tt>.
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector to add to the result. <tt>y</tt> contains on output the result.
 */
    void MultAdd(const PETScVect<T_>& x,
                       PETScVect<T_>& y) const;

/** \brief Multiply matrix by vector <tt>a*x</tt> and add to <tt>y</tt>.
 *  @param [in] a Constant to multiply by matrix
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector to add to the result. <tt>y</tt> contains on output the result.
 */
    void MultAdd(T_                   a,
                 const PETScVect<T_>& x,
                 PETScVect<T_>&       y) const;

/** \brief Assign a value to an entry of the matrix
 *  @param [in] i Row index
 *  @param [in] j Column index
 *  @param [in] a Value to assign to <tt>a(i,j)</tt>
 */
    void set(size_t    i,
             size_t    j,
             const T_& a);

/** \brief Add a value to an entry of the matrix
 *  @param [in] i Row index
 *  @param [in] j Column index
 *  @param [in] a Constant value to add to <tt>a(i,j)</tt>
 */
    void add(size_t    i,
             size_t    j,
             const T_& a);

/** \brief Assign values to a portion of the matrix
 *  @param [in] ir Vector of row indexes to assign (instance of class <tt>vector</tt>)
 *  @param [in] ic Vector of column indexes to assign (instance of class <tt>vector</tt>)
 *  @param [in] val Vector of values to assign (instance of class <tt>vector</tt>)
 */
    void set(vector<int>& ir,
             vector<int>& ic,
             vector<T_>&  val);

/// \brief Operator <tt>=</tt>
/// \details Assign constant value <tt>a</tt> to matrix diagonal entries
    void operator=(const T_& a);

/// \brief Set all matrix entries to zero
    void clear();

/** \brief Sets the matrix as the one for the Laplace equation in 1-D
 *  \details The matrix is initialized as the one resulting from P<sub>1</sub>
 *  finite element discretization of the classical elliptic operator
 *     -u'' = f
 *  with homogeneous Dirichlet boundary conditions
 *  \remark This function is available for real valued matrices only.
 *  @param [in] h %Mesh size (assumed constant)
 *  @param [in] mpi true if MPI is used for parallel computing, false if not (sequential),
 *              [Default: <tt>false</tt>]
 */
    void Laplace1D(real_t h,
                   bool   mpi=false);

/** \brief Sets the matrix as the one for the Laplace equation in 2-D
 *  \details The matrix is initialized as the one resulting from P<sub>1</sub>
 *  finite element discretization of the classical elliptic operator
 *     -Delta u = f
 *  with homogeneous Dirichlet boundary conditions
 *  \remark This function is available for real valued matrices only.
 *  @param [in] nx Number of unknowns in the <tt>x</tt>-direction
 *  @param [in] ny Number of unknowns in the <tt>y</tt>-direction
 *  @param [in] mpi true if MPI is used for parallel computing, false if not (sequential),
 *              [Default: <tt>false</tt>]
 */
    void Laplace2D(size_t nx,
                   size_t ny,
                   bool   mpi=false);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Factor matrix
    int Factor();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Solve the linear system of equations.
 *  \details The default parameters are: the Conjugate Gradient method and the Jacobi
 *  method for preconditioner.
 *  To change these values, call function setSolver before this function
 *  @param [in,out] b Vector that contains right-hand side on input and solution on output
 *
 *  @return Number of actual performed iterations 
 */
    int solve(PETScVect<T_>& b);

/** \brief Solve the linear system of equations.
 *  \details The default parameters are: the Conjugate Gradient method and the Jacobi
 *  method for preconditioner.
 *  To change these values, call function setSolver before this function
 *  @param [in] b Vector that contains right-hand side
 *  @param [out] x Vector that contains the obtained solution
 *  @return Number of actual performed iterations 
 */
    int solve(const PETScVect<T_>& b,
              PETScVect<T_>&       x);

/** \brief Choose solver and preconditioner for an iterative procedure
 *  @param [in] solver Option to choose iterative solver among the macros (see PETSc documentation for
 *  more details):
 *  <ul>
 *     <li><tt>KSPRICHARDSON</tt>: The Richardson iterative method (Default damping parameter is
 *                                 <tt>1.0</tt>)
 *     <li><tt>KSPCHEBYSHEV</tt>: The Chebyshev iterative method
 *     <li><tt>KSPCG</tt>: The conjugate gradient method [Default]
 *     <li><tt>KSPCGNE</tt>: The CG method for normal equations (without explicitly forming the product
 *         <tt>A^TA</tt>
 *     <li><tt>KSPGMRES</tt>: The GMRES iterative method (see A Generalized Minimal Residual
 *         Algorithm for Solving Nonsymmetric Linear Systems. Y. Saad and M. H. Schultz, SIAM J. Sci.
 *         Stat. Comput. Vo|. 7, No. 3, July 1986, pp. 856-869)
 *     <li><tt>KSPFGMRES</tt>: The Flexible GMRES method (with restart)
 *     <li><tt>KSPLGMRES</tt>: The 'augmented' standard GMRES method where the subspace uses
 *         approximations to the error from previous restart cycles
 *     <li><tt>KSPTCQMR</tt>: A variant of QMR (quasi minimal residual) developed by Tony Chan
 *     <li><tt>KSPBCGS</tt>: The BiCGStab (Stabilized version of BiConjugate Gradient Squared) method
 *     <li><tt>KSPIBCGS</tt>: The IBiCGStab (Improved Stabilized version of BiConjugate Gradient
 *         Squared) method in an alternative form to have only a single global reduction operation
 *         instead of the usual 3 (or 4)
 *     <li><tt>KSPFBCGS</tt>: The flexible BiCGStab method.
 *     <li><tt>KSPCGS</tt>: The CGS (Conjugate Gradient Squared) method
 *     <li><tt>KSPTFQMR</tt>: A transpose free QMR (quasi minimal residual)
 *     <li><tt>KSPCR</tt>: The conjugate residuals method
 *     <li><tt>KSPLSQR</tt>: The LSQR method
 *     <li><tt>KSPBICG</tt>: The Biconjugate gradient method (similar to running the conjugate
 *         gradient on the normal equations)
 *     <li><tt>KSPMINRES</tt>: The MINRES (Minimum Residual) method
 *     <li><tt>KSPSYMMLQ</tt>: The SYMMLQ method
 *     <li><tt>KSPGCR</tt>: The Generalized Conjugate Residual method
 *  </ul>
 *  @param [in] prec Option to choose preconditioner in an enumerated variable
 *  <ul>
 *     <li><tt>PCJACOBI</tt>: [Default] Jacobi (<i>i.e.</i> diagonal scaling) preconditioning
 *     <li><tt>PCBJACOBI</tt>: Block Jacobi preconditioning, each block is (approximately) solved
 *         with its own KSP object
 *     <li><tt>PCSOR</tt>: (S)SOR (successive over relaxation, Gauss-Seidel) preconditioning
 *     <li><tt>PCEISENSTAT</tt>: An implementation of SSOR (symmetric successive over relaxation,
 *         symmetric Gauss-Seidel) preconditioning that incorporates Eisenstat's trick to reduce
 *         the amount of computation needed
 *     <li><tt>PCICC</tt>: Incomplete Cholesky factorization preconditioners
 *     <li><tt>PCILU</tt>: Incomplete factorization preconditioners
 *     <li><tt>PCASM</tt>: Use the (restricted) additive Schwarz method, each 
 *         block is (approximately) solved with its own KSP object
 *     <li><tt>PCLU</tt>: Uses a direct solver, based on LU factorization, as a preconditioner
 *     <li><tt>PCCHOLESKY</tt>: Uses a direct solver, based on Cholesky factorization, as a preconditioner
 *  </ul>
 *  @param [in] toler Tolerance for convergence [Default: <tt>1.e-12</tt>]
 *  @param [in] max_it Maximum number of allowed iterations [Default: <tt>1000</tt>]
 */
    void setSolver(string solver,
                   string prec,
                   real_t toler=1.e-12,
                   int    max_it=1000);

/// \brief Return C-Array.
/// \details Non zero terms of matrix is stored row by row.
    T_ *get() const;

/** \brief  Return entry <tt>(i,j)</tt> of matrix if this one is stored, <tt>0</tt> otherwise
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    T_ get(size_t i,
            size_t j) const;

/** \brief Casting operator
 *  \details This member functions enables casting an instance of class PETScMatrix
 *  into the <tt>Petsc</tt> matrix type <tt>Mat</tt>. This is useful when one wants
 *  to usr any Petsc function that is not available in the wrapper (class PETScWrapper)
 *  or PETScMatrix.
 */
    operator Mat() const { return _A; }

/// \brief Get 1-norm of matrix
    PetscReal getNorm1() const;
 
/// \brief Get Frobenius norm of matrix
    PetscReal getFrobeniusNorm() const;

/// \brief Get infinity norm of matrix
    PetscReal getNormMax() const;

/** \brief Assembly of element matrix into global matrix.
 *  @param [in] el Reference to element instance
 *  @param [in] a Element matrix as a C-array
 */
    void Assembly(const Element& el,
                  T_*            a);

/** \brief Assembly of side matrix into global matrix.
 *  @param [in] sd Reference to side instance
 *  @param [in] a Side matrix as a C-array
 */
    void Assembly(const Side& sd,
                  T_*         a);

/// \brief Matrix assembly.
/// \details This function assembles matrix (begins and ends)
    void setAssembly();

/// \brief Activate MPI option
    void setMPI() { _mpi = true; }
    
 private:

    Mesh *_theMesh, *_theSubMesh;
    string _solver, _prec;
    int _type, _max_it, _dof_type, _is_diagonal, _istart, _iend;
    bool _extended, _one_dof, _sides, _mpi;
    bool _created, _set_nodes, _set_elements, _set_sides;
    size_t _nb_rows, _nb_cols, _size, _length, _dof;
    size_t _local_nb_rows, _local_nb_cols, _local_size;
    real_t _temp, _toler;
    PetscInt _np, _rank;
    PetscErrorCode _err;
    vector<int> _nnz, _cnnz, _col_ind, _row_ptr, _ir;
    vector<int> _diag_nnz, _off_diag_nnz;
    vector<T_> _v;
    Mat _A;

    void NodeGraph();

    void getRange()
    {
       MatGetOwnershipRange(_A,&_istart,&_iend);
    }

    static Mat create(size_t   nr,
                      size_t   nc,
                      MPI_Comm comm=PETSC_COMM_WORLD)
    {
       Mat m;
       MatCreate(comm,&m);
       MatSetSizes(m,PETSC_DECIDE,PETSC_DECIDE,nr,nc);
       MatSetFromOptions(m);
       MatSetUp(m);
       return m;
    }

    int col_index(size_t i, size_t j) const
    {
       for (int ii=0; ii<_nnz[i-1]; ii++) {
          size_t k=_col_ind[_cnnz[i-1]+ii];
          if (j==k+1)
             return k;
       }
       return -1;
    }
};


///////////////////////////////////////////////////////////////////////////////
//                         I M P L E M E N T A T I O N                       //
///////////////////////////////////////////////////////////////////////////////

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template<class T_>
PETScMatrix<T_>::PETScMatrix()
                : _solver(KSPCG), _prec(PCJACOBI), _max_it(1000), _mpi(false),
                  _created(false), _nb_rows(0), _nb_cols(0), _size(0),
                  _local_nb_rows(0), _local_nb_cols(0), _local_size(0), _toler(1.e-8),
                  _np(1), _rank(0)
{
}


template<class T_>
PETScMatrix<T_>::PETScMatrix(size_t nr,
                             size_t nc)
                : _solver(KSPCG), _prec(PCJACOBI), _max_it(1000), _mpi(false),
                  _nb_rows(nr), _nb_cols(nc), _local_nb_rows(nr), _local_nb_cols(nc),
                  _toler(1.e-8), _np(1), _rank(0)
{
   _size = 0;
   if (_nb_rows==_nb_cols)
      _size = _nb_rows;
   _length = _nb_rows * _nb_cols;
   _istart = 0; _iend = _nb_rows;
   _A = create(_nb_rows,_nb_cols);
   _created = true;
}


template<class T_>
PETScMatrix<T_>::PETScMatrix(size_t size)
                : _solver(KSPCG), _prec(PCJACOBI), _max_it(1000), _mpi(false),
                  _nb_rows(size), _nb_cols(size), _local_nb_rows(size), _local_nb_cols(size),
                  _toler(1.e-8),_np(1), _rank(0)

{
   _size = _nb_rows;
   _length = _size * _size;
   _istart = 0; _iend = _size;
   _A = create(_nb_rows,_nb_cols);
   _created = true;
}


template<class T_>
PETScMatrix<T_>::PETScMatrix(Mesh&  mesh,
                             size_t dof)
                : _solver(KSPCG), _prec(PCJACOBI), _max_it(1000), _toler(1.e-8),
                  _np(1), _rank(0)
{
   _created = false;
   setMesh(mesh,dof);
   _istart = 0; _iend = _nb_rows;
}


template<class T_>
PETScMatrix<T_>::PETScMatrix(size_t dof,
                             Mesh&  mesh,
                             int    code)
                : _solver(KSPCG), _prec(PCJACOBI), _max_it(1000), _toler(1.e-8),
                  _np(1), _rank(0)
{
   _created = false;
   setMesh(dof,mesh,mesh.getDOFSupport());
   _istart = 0; _iend = _nb_rows;
   _created = true;
}


template<class T_>
PETScMatrix<T_>::PETScMatrix(size_t dof,
                             size_t nb_eq,
                             Mesh&  mesh)
                : _solver(KSPCG), _prec(PCJACOBI), _max_it(1000), _toler(1.e-8),
                  _np(1), _rank(0)
{
   _created = false;
   setMesh(dof,nb_eq,mesh);
   _istart = 0; _iend = _nb_rows;
   _created = true;
}


template<class T_>
PETScMatrix<T_>::PETScMatrix(const vector<std::pair<size_t,size_t> >& I,
                                   int                                opt)
                : _solver(KSPCG), _prec(PCJACOBI), _max_it(1000), _created(true),
                  _toler(1.e-8), _np(1), _rank(0)
{
   setGraph(I,opt);
   _istart = 0; _iend = _nb_rows;
}


template<class T_>
PETScMatrix<T_>::PETScMatrix(const PETScMatrix& m)
{
   _dof = m._dof;
   _size = m._size;
   _length = m._length;
   _nnz.resize(_size+1);
   _cnnz.resize(_size+1);
   _col_ind.resize(_length);
   _col_ind = m._col_ind;
   _nnz = m._nnz;
   _cnnz = m._cnnz;
   _solver = m._solver;
   _prec = m._prec;
   _max_it = m._max_it;
   _toler = m._toler;
   _np = m._np;
   _rank = m._rank;
   setAssembly();
   _istart = 0; _iend = _nb_rows;
   _created = true;
}


template<class T_>
PETScMatrix<T_>::~PETScMatrix()
{
   if (_created)
      MatDestroy(&_A);
}


template<class T_>
void PETScMatrix<T_>::setAssembly()
{
   MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);
}


template<class T_>
void PETScMatrix<T_>::setSize(size_t size)
{
   if (_created)
      MatDestroy(&_A);
   _nb_rows = _nb_cols = _size = size;
   _local_nb_rows = _nb_rows, _local_nb_cols = _nb_cols;
   _length = _nb_rows*_nb_cols;
   _istart = 0; _iend = _nb_rows;
}


template<class T_>
void PETScMatrix<T_>::setSize(size_t nr,
                              size_t nc)
{
   if (_created)
      MatDestroy(&_A);
   _nb_rows = _local_nb_rows = nr;
   _nb_cols = _local_nb_cols = nc;
   _size = 0;
   if (_nb_rows==_nb_cols)
      _size = _nb_rows;
   _length = lsize_t(_nb_rows*_nb_cols);
   static int n=0;
   _istart = 0; _iend = _nb_rows;
   if (_np==1)
      MatCreateSeqAIJ(PETSC_COMM_SELF,_nb_rows,_nb_cols,0,&n,&_A);
   else
      MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,_nb_rows,_nb_cols,
                      0,PETSC_NULLPTR,0,PETSC_NULLPTR,&_A);
   MatSetFromOptions(_A);
}


template<class T_>
void PETScMatrix<T_>::setAIJ(const vector<int>& nnz)
{
   MatSetFromOptions(_A);
   MatSeqAIJSetPreallocation(_A,nnz.size(),&nnz[0]);
}


template<class T_>
void PETScMatrix<T_>::setAIJ_MPI(const vector<int>& diag_nnz,
                                 const vector<int>& off_nnz)
{
   MatSetFromOptions(_A);
   MatMPIAIJSetPreallocation(_A,diag_nnz.size(),&diag_nnz[0],off_nnz.size(),&off_nnz[0]);
}


template<class T_>
void PETScMatrix<T_>::setPartition(Partition& p)
{
   _theMesh = p.getMesh();
   _theSubMesh = &(p.getSubMesh(_rank));
   _dof_type = _theMesh->getDOFSupport();
   _size = _nb_rows = _nb_cols = _theMesh->getNbEq();
   if (_dof_type==ELEMENT_DOF)
      _size = _nb_rows = _nb_cols = _theMesh->getNbElements();
   NodeGraph();
   if (_created)
      MatDestroy(&_A);
   _diag_nnz.resize(_nb_rows);
   _off_diag_nnz.resize(_nb_rows);
   MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,_nb_rows,_nb_cols,
                   _diag_nnz.size(),&_diag_nnz[0],_off_diag_nnz.size(),&_off_diag_nnz[0],&_A);
   MatSetFromOptions(_A);
   MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);
}


template<class T_>
void PETScMatrix<T_>::setMesh(Mesh&  mesh,
                              size_t dof)
{
   _dof_type = mesh.getDOFSupport();
   _theMesh = &mesh;
   _dof = dof;
   _size = _nb_rows = _nb_cols = _theMesh->getNbEq();
   if (_dof_type==NODE_DOF) {
      if (_dof)
         _size = _nb_rows = _nb_cols = _theMesh->getNbNodes();
      else
         _size = _nb_rows = _nb_cols = _theMesh->getNbEq();
   }
   else if (_dof_type==SIDE_DOF) {
      if (_dof)
         _size = _nb_rows = _nb_cols = _theMesh->getNbSides();
      else
         _size = _nb_rows = _nb_cols = _theMesh->getNbEq();
   }
   else if (_dof_type==ELEMENT_DOF)
      _size = _nb_rows = _nb_cols = _theMesh->getNbElements();
   NodeGraph();
   if (_created)
      MatDestroy(&_A);
   MatCreateSeqAIJ(PETSC_COMM_SELF,_nb_rows,_nb_cols,PETSC_DEFAULT,&_nnz[0],&_A);
   MatSeqAIJSetColumnIndices(_A,&_col_ind[0]);
   MatSetFromOptions(_A);
   MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);
}


template<class T_>
void PETScMatrix<T_>::NodeGraph()
{
   std::vector<RC> pp;
   mesh_elements(*_theMesh) {
      for (size_t in=1; in<=The_element.getNbNodes(); in++) {
         Node *nd1=The_element(in);
         for (size_t k=1; k<=nd1->getNbDOF(); k++) {
            if (nd1->getDOF(k)!=0) {
               for (size_t jn=1; jn<=The_element.getNbNodes(); jn++) {
                  Node *nd2=The_element(jn);
                  for (size_t l=1; l<=nd2->getNbDOF(); l++) {
                     if (nd2->getDOF(l)!=0) {
                        pp.push_back(RC(nd1->getDOF(k)-1,nd2->getDOF(l)-1));
                        if (nd1->getDOF(k)!=nd2->getDOF(l))
                           pp.push_back(RC(nd2->getDOF(l)-1,nd1->getDOF(k)-1));
                     }
                  }
               }
            }
         }
      }
   }
   sort(pp.begin(),pp.end());
   vector<RC>::iterator new_end=std::unique(pp.begin(),pp.end());
   pp.erase(new_end,pp.end());
   size_t n=_theMesh->getNbEq();
   _length = pp.size();
   _col_ind.resize(_length);
   _nnz.resize(n+1);
   Clear(_nnz);
   for (size_t j=0; j<_length; j++) {
      _nnz[pp[j].first]++;
      _col_ind[j] = pp[j].second;
   }
   _cnnz.resize(n+1);
   _cnnz[0] = 0;
   for (size_t i=1; i<=n; i++)
      _cnnz[i] = _cnnz[i-1] + _nnz[i-1];
}


template<class T_>
void PETScMatrix<T_>::setGraph(const vector<std::pair<size_t,size_t> >& I,
                                     int                                opt)
{
   _length = I.size();
   _nb_rows = _nb_cols = 0;
   vector<RC> pp(I.size());
   for (size_t i=0; i<_length; i++) {
      _nb_rows = max(_nb_rows,I[i].first);
      _nb_cols = max(_nb_cols,I[i].second);
      pp[i] = RC(I[i].first-1,I[i].second-1);
   }
   _size = _nb_rows;
   if (opt==0) {
      sort(pp.begin(),pp.end());
      vector<RC>::iterator new_end = std::unique(pp.begin(),pp.end());
      pp.erase(new_end,pp.end());
   }
   _length = pp.size();
   _col_ind.resize(_length);
   _nnz.resize(_nb_rows);
   Clear(_nnz);
   for (size_t j=0; j<_length; j++)
      _nnz[pp[j].first]++;
   for (size_t k=1; k<_size; k++)
      _nnz[k] += _nnz[k-1];
   for (int ii=_size; ii>0; ii--)
      _nnz[ii] = _nnz[ii-1] + 1;
   MatCreateSeqAIJ(PETSC_COMM_SELF,_nb_rows,_nb_cols,PETSC_DEFAULT,&_nnz[0],&_A);
   MatSeqAIJSetColumnIndices(_A,&_col_ind[0]);
   MatSetFromOptions(_A);
}


template<class T_>
void PETScMatrix<T_>::setOneDOF()
{
   _one_dof = true;
}


template<class T_>
void PETScMatrix<T_>::setSides()
{
   _sides = true; 
}


template<class T_>
void PETScMatrix<T_>::Identity()
{
   Diagonal(1);
}


template<class T_>
void PETScMatrix<T_>::Diagonal()
{
   _is_diagonal = 1;
}


template<class T_>
void PETScMatrix<T_>::Diagonal(const T_& a)
{
   MatSeqAIJSetPreallocation(_A,1,NULL);
   MatSetFromOptions(_A);
   _is_diagonal = 1;
   _length = _nb_rows;
   for (size_t i=0; i<_nb_rows; ++i) {
      int ii=i, jj=i;
      MatSetValues(_A,1,&ii,1,&jj,1.,INSERT_VALUES);
   }
   MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);
}


template<>
inline void PETScMatrix<real_t>::Laplace1D(real_t h,
                                           bool   mpi)
{
   if (mpi)
      MatMPIAIJSetPreallocation(_A,3,NULL,3,NULL);
   MatSeqAIJSetPreallocation(_A,3,NULL);
   MatSetFromOptions(_A);
   _is_diagonal = 0;
   _length = _nb_rows;
   real_t v=2./h;
   for (size_t i=_istart; i<_iend; ++i) {
      int ii=i, jj=i;
      MatSetValues(_A,1,&ii,1,&jj,&v,INSERT_VALUES);
   }
   v = -1./h;
   for (size_t i=_istart+1; i<_iend; ++i) {
      int ii=i, jj=i-1;
      MatSetValues(_A,1,&ii,1,&jj,&v,INSERT_VALUES);
   }
   for (size_t i=_istart; i<_iend-1; ++i) {
      int ii=i, jj=i+1;
       MatSetValues(_A,1,&ii,1,&jj,&v,INSERT_VALUES);
   }
   MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);
}


template<>
inline void PETScMatrix<real_t>::Laplace2D(size_t nx,
                                           size_t ny,
                                           bool   mpi)
{
   _is_diagonal = 0;
   _size = nx*ny;
   _length = 5*_size;
   if (mpi)
      MatMPIAIJSetPreallocation(_A,5,NULL,5,NULL);
   MatSeqAIJSetPreallocation(_A,5,NULL);
   MatSetFromOptions(_A);
   if (mpi)
      getRange();
   for (size_t ii=_istart; ii<_iend; ii++) {
      real_t v=-1.0;
      int I=ii, J;
      size_t i=ii/ny, j=ii-i*ny;
      if (i>0) {
         J = ii - ny;
         MatSetValues(_A,1,&I,1,&J,&v,INSERT_VALUES);
      }
      if (i<nx-1) {
         J = I + ny;
         MatSetValues(_A,1,&I,1,&J,&v,INSERT_VALUES);
      }
      if (j>0) {
         J = I - 1;
         MatSetValues(_A,1,&I,1,&J,&v,INSERT_VALUES);
      }
      if (j<ny-1) {
         J = I + 1;
         MatSetValues(_A,1,&I,1,&J,&v,INSERT_VALUES);
      }
      v = 4.0;
      MatSetValues(_A,1,&I,1,&I,&v,INSERT_VALUES);
   }
   MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);
      //   MatSetUp(_A);
   if (mpi)
      MatGetOwnershipRange(_A,&_istart,&_iend);
   std::cout<<"$$ "<<mpi<<"  - "<<_istart<<"  "<<_iend<<endl;
}


template<class T_>
void PETScMatrix<T_>::DiagPrescribe(      PETScVect<T_>& b,
                                    const PETScVect<T_>& u)
{
   real_t p = 0;
   for (size_t j=_istart+1; j<=_iend; j++)
      p = max(p,std::abs(get(j,j)));
   size_t k=0;
   MESH_ND {
      for (size_t i=1; i<=TheNode.getNbDOF(); ++i) {
         size_t ii=TheNode.getDOF(i)-1;
         for (size_t j=0; j<_row_ptr[ii+1]-_row_ptr[ii]; j++,k++) {
            if (TheNode.getCode(i)>0) {
               b.set(ii+1,p*u[ii]);
               MatSetValue(_A,k,col_index(k,k),0,INSERT_VALUES);
	       if (ii+1==_col_ind[k])
                  MatSetValue(_A,k,col_index(k,k),p,INSERT_VALUES);
            }
	 }
      }
   }
   MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);
}


template<class T_>
T_ PETScMatrix<T_>::operator()(size_t i,
                               size_t j) const
{
   PetscInt ii=i-1, jj=col_index(i,j);
   PetscScalar a;
   MatGetValues(_A,1,&ii,1,&jj,&a);
   return a;
}


template<class T_>
void PETScMatrix<T_>::getMesh(Mesh& mesh)
{
      /*   if (_sides)
           SideGraph(mesh,_row_ptr,_col_ind);
           else {
           if (_dof)
           _size = _nb_rows = _nb_cols = mesh.getNbNodes();
           else
           _size = _nb_rows = _nb_cols = mesh.getNbEq();
           if (_type)
           XGraph(mesh,_row_ptr,_col_ind);
           //             XGraphScal(mesh,_row_ptr,_col_ind);
           else*/
   NodeGraph(mesh,_nnz,_col_ind);
      //   }
}


template<class T_>
void PETScMatrix<T_>::Mult(const PETScVect<T_>& x,
                                 PETScVect<T_>& y) const
{
   MatMult(_A,x,y);
   MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);
}


template<class T_>
void PETScMatrix<T_>::MultAdd(const PETScVect<T_>& x,
                                    PETScVect<T_>& y) const
{
   PETScVect<T_> z(y);
   MatMult(_A,x,y);
   y += z;
   MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);
}


template<class T_>
void PETScMatrix<T_>::MultAdd(      T_             a,
                              const PETScVect<T_>& x,
                                    PETScVect<T_>& y) const
{
   PETScVect<T_> z(y);
   MatMult(_A,x,y);
   y += a*z;
   MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);
}


template<class T_>
void PETScMatrix<T_>::set(      size_t i,
                                size_t j,
                          const T_&    a)
{
   MatSetValue(_A,i-1,col_index(i,j),a,INSERT_VALUES);
}


template<class T_>
void PETScMatrix<T_>::add(      size_t i,
                                size_t j,
                          const T_&    a)
{
   _err = MatSetValue(_A,i-1,col_index(i,j),a,ADD_VALUES);
}


template<class T_>
void PETScMatrix<T_>::clear()
{
   _err = MatZeroEntries(_A);
   MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);
}


template<class T_>
void PETScMatrix<T_>::operator=(const T_& a)
{
   _err = MatZeroEntries(_A);
   for (size_t i=_istart; i<_iend; i++)
      MatSetValue(_A,i,col_index(i+1,i+1),a,INSERT_VALUES);
   MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);
}


template<class T_>
int PETScMatrix<T_>::Factor()
{
   return -1;
}


template<class T_>
int PETScMatrix<T_>::solve(PETScVect<T_>& b)
{
   PETScVect<T_> x(b.size());
   int nb_it = Solve(b,x);
   b = x;
   _err = VecAssemblyBegin(b);
   _err = VecAssemblyEnd(b);
   return nb_it;
}


template<class T_>
int PETScMatrix<T_>::solve(const PETScVect<T_>& b,
                                 PETScVect<T_>& x)
{
   return 0;
}


template<class T_>
void PETScMatrix<T_>::setSolver(string solver,
                                string prec,
                                real_t toler,
                                int    max_it)
{
   _solver = solver;
   _prec = prec;
   _max_it = max_it;
   _toler = toler;
}


template<class T_>
PetscReal PETScMatrix<T_>::getNorm1() const
{
   PetscReal n;
   MatNorm(_A,NORM_1,&n);
   return n;
}


template<class T_>
PetscReal PETScMatrix<T_>::getFrobeniusNorm() const
{
   PetscReal n;
   MatNorm(_A,NORM_FROBENIUS,&n);
   return n;
}


template<class T_>
PetscReal PETScMatrix<T_>::getNormMax() const
{
   PetscReal n;
   MatNorm(_A,NORM_MAX,&n);
   return n;
}


template<class T_>
void PETScMatrix<T_>::Assembly(const Element& el,
                               T_*            a)
{
   _ir.clear();
   _v.clear();
   size_t kk=0;
   for (size_t i=1; i<=el.getNbNodes(); ++i) {
      Node *nd1=el(i);
      for (size_t k=1; k<=nd1->getNbDOF(); ++k) {
         _ir.push_back(nd1->getDOF(k)-1);
         for (size_t j=1; j<=el.getNbNodes(); ++j)
            for (size_t l=1; l<=el(j)->getNbDOF(); ++l)
               _v.push_back(a[kk++]);
      }
   }
   MatSetValues(_A,_ir.size(),&_ir[0],_ir.size(),&_ir[0],&_v[0],ADD_VALUES);
}


template<class T_>
void PETScMatrix<T_>::Assembly(const Side& sd,
                               T_*         a)
{
   _ir.clear();
   _v.clear();
   size_t kk=0;
   for (size_t i=1; i<=sd.getNbNodes(); ++i) {
      Node *nd1=sd(i);
      for (size_t k=1; k<=nd1->getNbDOF(); ++k) {
         _ir.push_back(nd1->getDOF(k)-1);
         for (size_t j=1; j<=sd.getNbNodes(); ++j)
            for (size_t l=1; l<=sd(j)->getNbDOF(); ++l)
               _v.push_back(a[kk++]);
      }
   }
   MatSetValues(_A,_ir.size(),&_ir[0],_ir.size(),&_ir[0],&_v[0],ADD_VALUES);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

///////////////////////////////////////////////////////////////////////////////
//                           ASSOCIATED  FUNCTIONS                           //
///////////////////////////////////////////////////////////////////////////////


/** \fn PETScVect<T_> operator*(const PETScMatrix<T_>& A, const PETScVect<T_>& x)
 *  \ingroup VectMat
 *  Multiply a matrix by a vector
 *  @param [in] A %Matrix to multiply by (Instance of class PETScMatrix)
 *  @param [in] x Vector to multiply by (Instance of class PETScVect)
 *  @return Vector product <tt>y = Ax</tt>
 */
template <class T_>
PETScVect<T_> operator*(const PETScMatrix<T_>& A,
                        const PETScVect<T_>&   x)
{
   VecAssemblyBegin(x);
   VecAssemblyEnd(x);
   PETScVect<T_> y(x.size());
   A.Mult(x,y);
   return y;
}


/** \fn ostream & operator<<(ostream& s, PETScMatrix<T_> &a)
 *  \ingroup VectMat
 *  \brief Output matrix in output stream
 */
template<class T_>
ostream& operator<<(ostream&         s,
                    PETScMatrix<T_>& A)
{
   A.setAssembly();
   PetscViewer viewer;
   PetscViewerASCIIOpen(PETSC_COMM_WORLD,PETSC_NULLPTR,&viewer);
   MatView(A,viewer);
   PetscViewerDestroy(&viewer);
   return s;
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif

#endif
