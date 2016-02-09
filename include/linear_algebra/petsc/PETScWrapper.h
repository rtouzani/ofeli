/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzaniset

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

                    Definition of Template class PETScWrapper
                         A wrapper for the PETSC library

  ==============================================================================*/

#ifndef __PETSC_WRAPPER_H
#define __PETSC_WRAPPER_H

#if defined (USE_PETSC)

#include "linear_algebra/petsc/PETScVect.h"
#include "linear_algebra/petsc/PETScMatrix.h"
#include "mesh/Mesh.h"
#include "mesh/Partition.h"

#include "petscksp.h"
#include "petscsnes.h"
#include "petscts.h"
#include "petsctao.h"

namespace OFELI {

/*! 
 *  \class PETScWrapper
 *  \ingroup VectMat
 *  \brief This class is a wrapper to be used when the library Petsc is installed and
 *  used with OFELI.
 *  \details When Petsc is used, an instance of class PETScWrapper is to be declared. It
 *  initializes the use of Petsc and enables calling solver functions in Petsc.
 *
 * \tparam T_ Data type (double, int, complex<double>, ...)
 *
 * When a linear system is invoked, the choice of iterative solvers can be made among the following methods 
 * (see Petsc documentation for more details):
 *  <ul>
 *     <li><tt>KSPRICHARDSON</tt>: The Richardson iterative method (Default damping parameter is
 *                                 <tt>1.0</tt>)
 *     <li><tt>KSPCHEBYSHEV</tt>: The Chebyshev iterative method
 *     <li><tt>KSPCG</tt>: The conjugate gradient method [Default]
 *     <li><tt>KSPCGNE</tt>: The CG method for normal equations (without explicitly forming the product
 *         <tt>A^TA</tt>
 *     <li><tt>KSPGMRES</tt>: [Default] The GMRES iterative method (see A Generalized Minimal Residual
 *         Algorithm for Solving Nonsymmetric Linear Systems. Y. Saad and M. H. Schultz, SIAM J. Sci.
 *         Stat. Comput., Vo|. 7, No. 3, July 1986, pp. 856-869)
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
 *
 * When a linear system is invoked, the choice of a preconditioner can be made among the following methods 
 * (see Petsc documentation for more details):
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
 */

template<class T_>
class PETScWrapper
{

 public:

/** \brief Constructor with program arguments
 *  @param [in] argc Count of number of command line arguments
 *  @param [in] args The command line arguments. Here is the list of arguments:
 *  <ul>
 *     <li> <tt>-start_in_debugger [noxterm,dbx,xdb,gdb,...]</tt>
 *          - Starts program in debugger
 *     <li> <tt>-on_error_attach_debugger [noxterm,dbx,xdb,gdb,...]</tt>
 *          - Starts debugger when error detected
 *     <li> <tt>-on_error_emacs &lt;machinename&gt;</tt> causes emacsclient to jump to error file
 *          - . -on_error_abort calls abort() when error detected (no traceback)
 *     <li> <tt>-on_error_mpiabort</tt> calls MPI_abort() when error detected
 *          - . -error_output_stderr prints error messages to stderr instead of the default stdout
 *     <li> <tt>-error_output_none</tt> does not print the error messages
 *          (but handles errors in the same way as if this was not called)
 *          - . -debugger_nodes [node1,node2,...] - Indicates nodes to start in debugger
 *     <li> <tt>-debugger_pause [sleeptime]</tt> (in seconds)
 *          - Pauses debugger
 *     <li> <tt>-stop_for_debugger</tt>
 *          - Print message on how to attach debugger manually to process and wait
 *            (<tt>-debugger_pause</tt>) seconds for attachment
 *     <li> <tt>-malloc</tt>
 *          - Indicates use of PETSc error-checking malloc (on by default for debug version of libraries)
 *     <li> <tt>-malloc no</tt>
 *          - Indicates not to use error-checking malloc
 *     <li> <tt>-malloc_debug</tt>
 *          - check for memory corruption at EVERY malloc or free
 *     <li> <tt>-malloc_dump</tt>
 *          - prints a list of all unfreed memory at the end of the run
 *     <li> <tt>-malloc_test</tt>
 *          - like <tt>-malloc_dump</tt> <tt>-malloc_debug</tt>, but only active for debugging builds
 *     <li> <tt>-fp_trap</tt>
 *          - Stops on floating point exceptions (Note that on the IBM RS6000 this slows code by at 
 *          least a factor of 10.)
 *     <li> <tt>-no_signal_handler</tt>
 *          - Indicates not to trap error signals
 *     <li> <tt>-shared_tmp</tt>
 *          - indicates <tt>/tmp</tt> directory is shared by all processors
 *     <li> <tt>-not_shared_tmp</tt>
 *          - each processor has own <tt>/tmp</tt>
 *     <li> <tt>-tmp</tt>
 *          - alternative name of <tt>/tmp</tt> directory
 *     <li> <tt>-get_total_flops</tt>
 *          - returns total flops done by all processors
 *     <li> <tt>-memory_info</tt>
 *          - Print memory usage at end of run
 *  </ul>
 *  @param [in] help String that contains message to display when argument
 *              <tt>-v</tt> is used
 *
 * @warning This class is available only when OFELI has been installed with Petsc
 *
 */
    PETScWrapper(int    argc,
                 char** args,
                 string help="")
                : _set_ksp(0), _set_snes(0), _set_ts(0), _set_tao(0), _tol(1.e-14)
    {
       PetscInt np;
       PetscInitialize(&argc,&args,(char *)0,help.c_str());
       MPI_Comm_size(PETSC_COMM_WORLD,&_size);
       MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);
       PetscOptionsGetInt(PETSC_NULL,"-np",&np,PETSC_NULL);
    }

/// \brief Destructor
/// \details Destroy the KSP context and release memory allocated by petsc
    ~PETScWrapper()
    {
       if (_set_ksp)
          KSPDestroy(&_ksp);
       _err = PetscFinalize();
    }

/** \brief Get an option as an integer number
 *  @param [in] s String to preprend the name of the option
 *  @param [out] n Obtained integer value  
 *  @param [out] set <tt>true</tt> if found, <tt>false</tt> if not.
 */
    PetscErrorCode getIntOption(string     s,
                                PetscInt&  n,
                                PetscBool& set) const
    {
       return PetscOptionsGetInt(NULL,s.c_str(),&n,&set);
    }

/** \brief Get an option as a bool variable
 *  @param [in] s String to preprend the name of the option
 *  @param [out] b Obtained boolean value  
 *  @param [out] set <tt>true</tt> if found, <tt>false</tt> if not.
 */
    PetscErrorCode getBoolOption(string     s,
                                 PetscBool& b,
                                 PetscBool& set) const
    {
       return PetscOptionsGetBool(NULL,s.c_str(),&b,&set);
    }

/// \brief Return wrapper size, <i>i.e.</i> number of processors
    PetscMPIInt size() const { return _size; }

/// \brief Set mesh
/// @param [in] ms Mesh instance
    void setMesh(Mesh& ms) {
       _theMesh = &ms;
       _A->setMesh(*_theMesh);
       _A->setRank(_rank);
    }

/// \brief Set mesh partition
/// \details This function is to be used for parallel computing
/// @param [in] p Partition instance
    void setPartition(Partition& p) {
       _part = &p;
       _A->setPartition(*_part);
    }

/// \brief Define problem matrix
/// @param [in] A PETScMatrix instance that contains matrix
    void setMatrix(PETScMatrix<T_>& A)
    {
       _set_ksp = 1;
       _A = &A;
       KSPCreate(PETSC_COMM_WORLD,&_ksp);
       KSPSetOperators(_ksp,*_A,*_A);

/*     Set linear solver defaults for this problem (optional).
       - By extracting the KSP and PC contexts from the KSP context,
         we can then directly call any KSP and PC routines to set
         various options.
       - The following four statements are optional; all of these
         parameters could alternatively be specified at runtime via
         KSPSetFromOptions();
 */
       KSPGetPC(_ksp,&_pc);
       KSPSetTolerances(_ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
       KSPSetInitialGuessNonzero(_ksp,PETSC_TRUE);
    }

/** \brief Set linear system features
 *  @param [in] A PETScMatrix instance that contains matrix
 *  @param [in] b Vector containing the right-hand side
 *  @param [in] s Option to choose iterative solver. See the definition of the class for 
 *              iterative methods
 *  @param [in] p Option to choose preconditioner. See the definition of the class for available
 *  preconditioners.
 *  @param [in] tol Tolerance for convergence of iteration process
 *              [Default: <tt>1.e-12</tt>]
 *  @param [in] max_it Maximum number of linear solver iterations [Default: <tt>1000</tt>]
 */
    void setLinearSystem(PETScMatrix<T_>& A,
                         PETScVect<T_>&   b,
                         string           s=KSPCG,
                         string           p=PCJACOBI,
                         real_t           tol=1.e-12,
                         size_t           max_it=1000)
    {
       setMatrix(A);
       _b = &b;
       setIterationMethod(s);
       setPreconditioner(p);
       setIterationParameters(tol,max_it);
    }

/** \brief Choose preconditioner for the iterative procedure
 *  @param [in] p Option to choose preconditioner. See the definition of the class for available
 *  preconditioners.
 */
    void setPreconditioner(string p)
    {
       PCSetType(_pc,p.c_str());
    }

/** \brief Choose iteration parameters
 *  @param [in] tol Tolerance for convergence of iteration process
 *  @param [in] max_it Maximum number of linear solver iterations
 */
    void setIterationParameters(real_t tol,
                                size_t max_it)
    {
       _tol = tol;
       KSPSetTolerances(_ksp,_tol,PETSC_DEFAULT,PETSC_DEFAULT,max_it);
    }

/** \brief Choose the iterative method.
 *  @param [in] m Option to choose iterative solver. See the definition of the class for available
 *  iterative solvers.
 */
    void setIterationMethod(string m)
    {
       KSPSetType(_ksp,m.c_str());
    }

/** \brief Solve the linear system
 *  \details If the member functions setIterationMethod and setPreconditioner have not been used,
 *  default methods are used
 *  @param [in,out] x Vector containing the initial guess on input and, if convergence is achieved,
 *  the solution on output
 */
    void solve(PETScVect<T_>& x)
    {
       MatAssemblyBegin(*_A,MAT_FINAL_ASSEMBLY);
       MatAssemblyEnd(*_A,MAT_FINAL_ASSEMBLY);
       _x = &x;
       KSPSolve(_ksp,*_b,*_x);
    }

/** \brief Solve the linear system
 *  \details If the member functions setIterationMethod and setPreconditioner have not been used,
 *  default methods are used
 *  @param [in] b Vector containing the right-hand side
 *  @param [in,out] x Vector containing the initial guess on input and, if convergence is achieved,
 *  the solution on output
 */
    void solve(const PETScVect<T_>& b,
                     PETScVect<T_>& x)
    {
       _b = &b;
       solve(x);
    }

/** \brief Check residual error
 *  \details This function computes the residual A*x - b and outputs the
 *  number of iterations
 *  @param [out] u Residual vector
 */
    void checkError(PETScVect<T_>& u) const
    {
       VecAXPY(*_x,-1.0,u);
       PetscReal norm;
       VecNorm(*_x,NORM_2,&norm);
       PetscInt it;
       KSPGetIterationNumber(_ksp,&it);
       cout << "Number of iterations: " << it << endl;
       cout << "Norm of error: " << norm << endl;
    }

/// \brief Return the number of iterations
    int getIterationNumber() const
    {
       int it;
       KSPGetIterationNumber(_ksp,&it);
       return it;
    }

/** \brief Set tolerance parameters for a linear system
 *  @param [in] rel_tol Relative convergence tolerance, relative decrease in the preconditioned residual norm
 *  @param [in] abs_tol Absolute convergence tolerance of the preconditioned residual norm
 *  @param [in] div_tol Divergence tolerance: Amount preconditioned residual norm 
 *  @param [in] max_it Maximum number of iterations
 */
    void setLSTolerances(real_t rel_tol,
                         real_t abs_tol,
                         real_t div_tol=PETSC_DEFAULT,
                         int    max_it=PETSC_DEFAULT) const
    {
       KSPSetTolerances(_ksp,rel_tol,abs_tol,div_tol,max_it);
    }

/// \brief Return the rank of the current processor
    PetscMPIInt getRank() const { return _rank; }

/// \brief
    void setPartition()
    {
      //       _par.set(*_theMesh,_size);
    }

/// \brief Output wrapper information
    template<class S_>
    friend ostream &operator<<(      ostream&          s,
                               const PETScWrapper<S_>& w);

 private:
   Mesh            *_theMesh;
   Partition       *_part;
   int             _init_petsc, _set_ksp, _set_snes, _set_ts, _set_tao;
   PetscErrorCode  _err;
   PetscInt        _size;
   PetscMPIInt     _rank;
   KSP             _ksp;
   SNES            _snes;
   TS              _ts;
   PC              _pc;
   real_t          _tol;
   PETScVect<T_>   *_x, *_b;
   PETScMatrix<T_> *_A;
   //   Partition       _par;
};


/** \fn ostream &operator<<(ostream &s, const PETScWrapper<T_> &w)
 *  \brief Output Petsc Wrapper data in output stream
 *  \ingroup VectMat
 */
template<class T_>
ostream &operator<<(      ostream&          s, 
                    const PETScWrapper<T_>& w)
{
   if (w._set_ksp) {
      s << "Linear Solver data:" << endl;
      KSPView(w._ksp,PETSC_VIEWER_STDOUT_SELF);
   }
   if (w._set_snes) {
      s << "Nonlinear Solver data:" << endl;
      SNESView(w._snes,PETSC_VIEWER_STDOUT_SELF);
   }
   if (w._set_ts) {
      s << "Time Stepping data:" << endl;
      TSView(w._ts,PETSC_VIEWER_STDOUT_SELF);
   }
   return s;
}

} /* namespace OFELI */

/*! @} */

#endif

#endif

/*

Available preconditioners:
PCJACOBI, PCBJACOBI, PCSOR, PCEISENSTAT, PCICC, PCILU, PCASM, PCKSP, PCLU, PCCHOLESKY

Available iterative solvers:
KSPRICHARDSON, KSPCHEBYSHEV, KSPCG, KSPGROPPCG, KSPPIPECG, KSPCGNE, KSPNASH, KSPSTCG
KSPGLTR, KSPGMRES, KSPFGMRES, KSPLGMRES, KSPDGMRES, KSPPGMRES, KSPTCQMR, KSPBCGS
KSPIBCGS, KSPFBCGS, KSPFBCGSR, KSPBCGSL, KSPCGS, KSPTFQMR, KSPCR, KSPPIPECR, KSPLSQR
KSPPREONLY, KSPQCG, KSPBICG, KSPBICGSL, KSPMINRES, KSPSYMMLQ, KSPLCD, KSPPYTHON, KSPGCR, KSPSPECEST

*/
