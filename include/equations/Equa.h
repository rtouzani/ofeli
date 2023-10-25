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

        Definition of abstract class 'Equa' for Finite Element Equations

  ==============================================================================*/

#ifndef __EQUA_H
#define __EQUA_H

#if defined (USE_PETSC)
#include "linear_algebra/petsc/PETScMatrix.h"
#else
#include "linear_algebra/DMatrix.h"
#include "linear_algebra/DSMatrix.h"
#include "linear_algebra/SpMatrix.h"
#include "linear_algebra/SkMatrix.h"
#include "linear_algebra/SkSMatrix.h"
#include "linear_algebra/BMatrix.h"
#include "linear_algebra/TrMatrix.h"
#endif

#include "mesh/Mesh.h"
#include "mesh/Grid.h"
#include "io/Prescription.h"
#include "linear_algebra/LocalMatrix.h"
#include "linear_algebra/LocalVect.h"
#include "solvers/LinearSolver.h"
#include "solvers/EigenProblemSolver.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \defgroup Equation General Purpose Equations
 *  \brief Gathers equation related classes
 */

/*! \file Equa.h
 *  \brief Definition file for abstract class Equa.
 */

/*! \enum PDE_Terms
 * Enumerate variable that selects various terms in partial differential equations
 */
enum PDE_Terms {
   CONSISTENT_MASS     =  0x00001000,     /*!< Consistent mass term                    */
   LUMPED_MASS         =  0x00002000,     /*!< Lumped mass term                        */
   MASS                =  0x00002000,     /*!< Consistent mass term                    */
   CAPACITY            =  0x00004000,     /*!< Consistent capacity term                */
   CONSISTENT_CAPACITY =  0x00004000,     /*!< Consistent capacity term                */
   LUMPED_CAPACITY     =  0x00008000,     /*!< Lumped capacity term                    */
   VISCOSITY           =  0x00010000,     /*!< Viscosity term                          */
   STIFFNESS           =  0x00020000,     /*!< Stiffness term                          */
   DIFFUSION           =  0x00040000,     /*!< Diffusion term                          */
   MOBILITY            =  0x00040000,     /*!< Mobility term                           */
   CONVECTION          =  0x00080000,     /*!< Convection term                         */
   DEVIATORIC          =  0x00100000,     /*!< Deviatoric term                         */
   DILATATION          =  0x00200000,     /*!< Dilatational term                       */
   ELECTRIC            =  0x00400000,     /*!< Electric term                           */
   MAGNETIC            =  0x00800000,     /*!< Magnetic term                           */
   LOAD                =  0x01000000,     /*!< Body load term                          */
   HEAT_SOURCE         =  0x02000000,     /*!< Body heat source term                   */
   BOUNDARY_TRACTION   =  0x04000000,     /*!< Boundary traction (pressure) term       */
   HEAT_FLUX           =  0x08000000,     /*!< Boundary heat flux term                 */
   CONTACT             =  0x10000000,     /*!< Signorini contact                       */
   BUOYANCY            =  0x20000000,     /*!< Buoyancy force term                     */
   LORENTZ_FORCE       =  0x40000000,     /*!< Lorentz force term                      */
   DAMPING             =  0x80000000      /*!< Damping term                            */
};


/*! \enum GENERIC_PDE_Terms
 * Enumerate variable that selects various terms in a generic partial differential equation
 */
enum GENERIC_PDE_Terms {
   NOTERM       =  0x00000000,     /*!< No term (empty equation)                      */
   L00          =  0x00001000,     /*!< 0th order in time and space, to LHS           */
   L10          =  0x00002000,     /*!< 1st order in time, 0th order in space, to LHS */
   L20          =  0x00004000,     /*!< 2nd order in time, 0th order in space, to LHS */
   L01          =  0x00008000,     /*!< 0th order in time, 1st order in space, to LHS */
   L02          =  0x00010000,     /*!< 0th order in time, 2nd order in space, to LHS */
   L11          =  0x00020000,     /*!< 1st order in time, 1st order in space, to LHS */
   BODY_RHS     =  0x00040000,     /*!< Given right-hand side on domain               */
   BOUNDARY_RHS =  0x00080000,     /*!< Given right-hand side on boundary             */
   NEUMANN      =  0x00008000,     /*!< Neumann boundary condition                    */
};


/*! \enum Analysis
 * Selects Analysis type
 */
enum Analysis {
   STATIONARY         =  0,    /*!< Steady State analysis                         */
   STEADY_STATE       =  0,    /*!< Steady state analysis                         */
   TRANSIENT          =  1,    /*!< Transient problem                             */
   TRANSIENT_ONE_STEP =  2,    /*!< Transient problem, perform only one time step */
   OPTIMIZATION       =  3,    /*!< Optimization problem                          */
   EIGEN              =  4     /*!< Eigenvalue problem                            */
};


/*! \enum TimeScheme
 * Selects Time integration scheme
 */
enum TimeScheme {
   NONE               =  0,    /*!< No time integration scheme                    */ 
   FORWARD_EULER      =  1,    /*!< Forward Euler scheme (Explicit)               */
   BACKWARD_EULER     =  2,    /*!< Backward Euler scheme (Implicit)              */
   CRANK_NICOLSON     =  3,    /*!< Crank-Nicolson scheme                         */
   HEUN               =  4,    /*!< Heun scheme                                   */
   NEWMARK            =  5,    /*!< Newmark scheme                                */
   LEAP_FROG          =  6,    /*!< Leap Frog scheme                              */
   ADAMS_BASHFORTH    =  7,    /*!< Adams-Bashforth scheme (2nd Order)            */
   AB2                =  7,    /*!< Adams-Bashforth scheme (2nd Order)            */
   RUNGE_KUTTA        =  8,    /*!< 4-th Order Runge-Kutta scheme (4th Order)     */
   RK4                =  8,    /*!< 4-th Order Runge-Kutta scheme                 */
   RK3_TVD            =  9,    /*!< 3-rd Order Runge-Kutta TVD scheme             */
   BDF2               = 10,    /*!< Backward Difference Formula (2nd Order)       */
   BUILTIN            = 11     /*!< Builtin scheme, implemented in equation class */
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/*! \enum PDE_Name
 * Choose partial differential equation to `
 */
enum PDE_Name {
   LAPLACE,                      /*!< Laplace equation                        */
   DIFFUSION_CONVECTION,         /*!< Diffusion Convection equation           */
   THERMAL_PHASE_CHANGE,         /*!< Thermal phase change problem (Stefan)   */
   INCOMPRESSIBLE_NAVIER_STOKES, /*!< Incompressible Navier-Stokes equations  */
   LINEARIZED_ELASTICITY,        /*!< Linearized elasticity equations         */
   PLANAR_TRUSS,                 /*!< 2-D truss equation                      */
   SPATIAL_BEAM                  /*!< 3-D beam equations                      */
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! \enum FEType
 * Choose Finite Element Type
 */
enum FEType {
   FE_2D_3N,                     /*!< 2-D elements, 3-Nodes (P1)              */
   FE_2D_6N,                     /*!< 2-D elements, 6-Nodes (P2)              */
   FE_2D_4N,                     /*!< 2-D elements, 4-Nodes (Q1)              */
   FE_3D_AXI_3N,                 /*!< 3-D Axisymmetric elements, 3-Nodes (P1) */
   FE_3D_4N,                     /*!< 3-D elements, 4-Nodes (P1)              */
   FE_3D_8N                      /*!< 3-D elements, 8-Nodes (Q1)              */
};

/*! \enum PDECoefType
 * Choose PDE Coefficient
 */
enum PDECoefType {
    C00                    =  0,   /*!    */
    C10                    =  1,   /*!    */
    C01                    =  2,   /*!    */
    C20                    =  3,   /*!    */
    C02                    =  4,   /*!    */
    RHO                    =  5,   /*!    */
    DENSITY                =  5,   /*!    */
    CP                     =  6,   /*!    */
    KAPPA                  =  7,   /*!    */
    THERMAL_CONDUCTIVITY   =  7,   /*!    */
    MU                     =  8,   /*!    */
    DYNAMIC_VISCOSITY      =  8,   /*!    */
    SIGMA                  =  9,   /*!    */
    ELECTRIC_CONDUCTIVITY  =  9,   /*!    */
    MMU                    = 10,   /*!    */
    MAGNETIC_PERMEABILITY  = 10,   /*!    */
    EPSILON                = 11,   /*!    */
    ELECTRIC_PERMITTIVITY  = 11,   /*!    */
    OMEGA                  = 12,   /*!    */
    ANGULAR_FREQUENCY      = 12,   /*!    */
    BETA                   = 13,   /*!    */
    VV                     = 14,   /*!    */
    YOUNG                  = 15,   /*!    */
    POISSON                = 16    /*!    */
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/*! \struct TimeIntegration
 *  Structure for time sequencing
 */
struct TimeIntegration {
TimeIntegration() : step(0), init(0.), final(1.), delta(0.1), scheme(FORWARD_EULER) { }
   int step;
   real_t time, init, final, delta;
   real_t theta, alpha, beta, time_parameter1, time_parameter2, time_parameter3;
   TimeScheme scheme;
};

/*! \struct ElementGeom
 *  Structure Element geometry data
 */
struct ElementGeom {
   real_t volume, area, length, size, det;
   Point<real_t> center, circumcenter;
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! \class Equa
 * \ingroup Equation
 * \brief Mother abstract class to describe equation.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class TimeStepping;


class Equa
{

 public:

/// \brief Default constructor
    Equa();

/// \brief Destructor
    virtual ~Equa();

/// \brief Define mesh and renumber DOFs after removing imposed ones
    void setMesh(Mesh &m);

/// \brief Return reference to Mesh instance
/// @return Reference to Mesh instance
    Mesh &getMesh() const;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** \brief Build equation.
 *  This member function is to be used if one wants to use a class that inherits
 *  from Equa that handles the whole solution process (including assembly).
 */
    virtual void build() = 0;
    virtual void build(EigenProblemSolver& e);
    virtual void build(TimeStepping& s);
    virtual int runOneTimeStep();
    virtual int runSteadyState();
    virtual int runTransient();

 /** \brief Solve the equation
 *  \details This function solves the thermal problem according to the argument
 *  @param [in] pa Problem type: must be chosen among enumerated values:
 *  <li>STATIONARY: For a stationary problem (Default value)
 *  <li>TRANSIENT: For a transient (time dependent) problem.
 *  <li>TRANSIENT_ONE_STEP: For a one time step in the transient problem
 */
    int run(Analysis   a=STATIONARY,
            TimeScheme s=NONE);

/** \brief 
 *  \details 
 *  @param [in] Df
 */
    void getTangent(Matrix<real_t>* Df);

/** \brief Set transient analysis settings
 *  \details Define a set of parameters for time integration
 *  @param [in] s Time inegration scheme
 *  @param [in] p1 First parameter for the scheme
 *  @param [in] p2 Second parameter for the scheme
 *  @param [in] p3 Third parameter for the scheme
 */
    void setTransient(TimeScheme s,
                      real_t     p1=1,
                      real_t     p2=1,
                      real_t     p3=1);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Return reference to linear solver instance
    LinearSolver &getLinearSolver();

/// \brief Return pointer to matrix
#if defined (USE_PETSC)
    PETScMatrix<real_t> *getMatrix() const;
#else
    Matrix<real_t> *getMatrix() const;
#endif

/** \brief Choose solver for the linear system
 *  @param [in] ls Solver of the linear system.
 *  To choose among the enumerated values: <tt>DIRECT_SOLVER</tt>, <tt>CG_SOLVER</tt>,
 *  <tt>GMRES_SOLVER</tt>
 *  <ul>
 *    <li> <tt>DIRECT_SOLVER</tt>, Use a facorization solver [default]
 *    <li> <tt>CG_SOLVER</tt>, Conjugate Gradient iterative solver
 *    <li> <tt>CGS_SOLVER</tt>, Squared Conjugate Gradient iterative solver
 *    <li> <tt>BICG_SOLVER</tt>, BiConjugate Gradient iterative solver
 *    <li> <tt>BICG_STAB_SOLVER</tt>, BiConjugate Gradient Stabilized iterative solver
 *    <li> <tt>GMRES_SOLVER</tt>, GMRES iterative solver
 *    <li> <tt>QMR_SOLVER</tt>, QMR iterative solver
 *  </ul>
 *  @param [in] pc Preconditioner to associate to the iterative solver.
 *  If the direct solver was chosen for the first argument this argument is not used.
 *  Otherwise choose among the enumerated values:
 *  <ul> 
 *    <li> <tt>IDENT_PREC</tt>, Identity preconditioner (no preconditioning [default])
 *    <li> <tt>DIAG_PREC</tt>, Diagonal preconditioner
 *    <li> <tt>ILU_PREC</tt>, Incomplete LU factorization preconditioner
 *  </ul>
 */
    void setSolver(Iteration      ls,
                   Preconditioner pc=IDENT_PREC);

/** \brief Choose type of matrix
 *  @param [in] t Type of the used matrix.
 *  To choose among the enumerated values: <tt>SKYLINE</tt>, <tt>SPARSE</tt>, <tt>DIAGONAL</tt>
 *  <tt>TRIDIAGONAL</tt>, <tt>SYMMETRIC</tt>, <tt>UNSYMMETRIC</tt>, <tt>IDENTITY</tt>
 */
    void setMatrixType(int t);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    virtual void setWithConvection(int f);
    void setTerms(PDE_Terms t);
    void set_rho(string exp);
    void set_Cp(string exp);
    void set_kappa(string exp);
    void set_mu(string exp);
    void set_sigma(string exp);
    void set_Mu(string exp);
    void set_epsilon(string exp);
    void set_omega(string exp);
    void set_beta(string exp);
    void set_v(string exp);
    void set_young(string exp);
    void set_poisson(string exp);
    virtual void set_00(real_t a=1.0) { }
    virtual void set_00(Fct& f) { }
    virtual void set_00(const string& f) { }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Solve the linear system with given matrix and right-hand side
 *  @param [in] A Pointer to matrix of the system
 *  @param [in] b Vector containing right-hand side
 *  @param [in,out] x Vector containing initial guess of solution on input, actual solution on output
 */
#if defined(USE_PETSC)
    int solveLinearSystem(PETScMatrix<real_t>* A,
                          PETScVect<real_t>&   b,
                          PETScVect<real_t>&   x);
#else
    int solveLinearSystem(Matrix<real_t>* A,
                          Vect<real_t>&   b,
                          Vect<real_t>&   x);
#endif

/** \brief Solve the linear system with given right-hand side
 *  @param [in] b Vector containing right-hand side
 *  @param [in,out] x Vector containing initial guess of solution on input, actual solution on output
 */
#if defined(USE_PETSC)
    int solveLinearSystem(PETScVect<real_t>& b,
                          PETScVect<real_t>& x);
#else
    int solveLinearSystem(Vect<real_t>& b,
                          Vect<real_t>& x);
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Set Analysis Type
    void setAnalysis(Analysis a);

/// \brief Set time iteration index
    void setTimeIndex(size_t step);

/// \brief Set initial time
    void setInitTime(real_t t);

/// \brief Set time step value
    void setTimeStep(real_t t);

/// \brief Return time step
    real_t getTimeStep() const;

/// \brief Set current time
    void setTime(real_t t);

/// \brief Set final time
    void setFinalTime(real_t t);

/// \brief Return number of fields
    size_t getNbFields() const;

/// \brief Set time integration scheme
    void setTimeIntegration(TimeScheme s);

/// \brief Set time integration parameters
    void setTimeIntegrationParam();

/// \brief Return time integration scheme
    TimeScheme getTimeIntegration() const;
    
/// \brief Return equation name
    string getEquationName() const;

/// \brief Return Finite Element Type
    string getFiniteElementType() const;

/** \brief Set equation input data
 *  @param [in] opt Parameter that selects data type for input. This parameter
 *  is to be chosen in the enumerated variable EqDataType
 *  @param [in] u Vect instance that contains input vector data
 *  List of data types contains <tt>INITIAL_FIELD</tt>, <tt>BOUNDARY_CONDITION_DATA</tt>, 
 *  <tt>SOURCE_DATA</tt> or <tt>FLUX</tt> with obvious meaning
 */
#if defined(USE_PETSC)
    virtual void setInput(EqDataType     opt,
                          PETScVect<real_t>& u);
#else
    virtual void setInput(EqDataType opt,
                          Vect<real_t>&  u);
#endif

    void set(Prescription& p);
    void setTolerance(real_t toler);
    void setPDECoef(PDECoefType t,
                    real_t      a);
    void setPDECoef(PDECoefType   t,
                    const string& s);
    void setPDECoef(PDECoefType t,
                    Fct&        f);
    real_t getPDECoef(PDECoefType      c,
                      const SpaceTime& p);
    Vect<real_t> getPDECoefV(PDECoefType      c,
                             const SpaceTime& p);
    real_t getPDECoef(PDECoefType c,
                      real_t      x,
                      real_t      y,
                      real_t      z,
                      real_t      t);
    Vect<real_t> getPDECoefV(PDECoefType c,
                             real_t      x,
                             real_t      y,
                             real_t      z,
                             real_t      t);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Print info on linear system solver
    void LinearSystemInfo();

 protected:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
   Mesh                   *_theMesh;
   Grid                   *_theGrid;
   size_t                 _nb_nodes, _nb_sides, _nb_boundary_sides, _nb_el, _nb_eq, _nb_dof, _nb_dof_total;
   size_t                 _ngx, _ngy, _ngz;
   real_t                 _hx, _hy, _hz;
   int                    _field_type, _terms;
   GENERIC_PDE_Terms      _gterms;
   int                    _matrix_type, _solver, _max_it;
   size_t                 _nb_fields, _nb_eigv;
   EigenProblemSolver     _ev;
   Iteration              _its;
   Preconditioner         _prec;
   bool                   _eigen;
   bool                   _constant_matrix, _constant_mesh, _set_matrix, _set_solver;
   int                    _sol_type, _init_type, _bc_type, _bf_type, _sf_type;
   string                 _equation_name, _finite_element;
   bool                   _rho_set, _Cp_set, _kappa_set, _mu_set, _sigma_set, _Mu_set;
   bool                   _epsilon_set, _omega_set, _beta_set, _v_set, _young_set, _poisson_set;
   Fct                    _rho_fct, _Cp_fct, _kappa_fct, _mu_fct, _sigma_fct, _Mu_fct;
   Fct                    _epsilon_fct, _omega_fct, _beta_fct, _v_fct, _young_fct, _poisson_fct;
   Fct                    _theFct;
   std::map<PDECoefType,real_t>        _coef_value;
   std::map<PDECoefType,Vect<real_t>*> _coef_vector;
   std::map<PDECoefType,string>        _coef_string;
   std::map<PDECoefType,Fct *>         _coef_fct;
   std::map<PDECoefType,int>           _set_coef;
   real_t                 _ex, _ey, _ez, _et;
   LinearSolver           _ls;
   real_t                 _toler;
   Prescription           *_prescription;
   TimeIntegration        _TimeInt;
   Analysis               _analysis;
   vector<Point<real_t> > _dSh;
   vector<real_t>         _sh, _wg;
   real_t                 _body_source, _bound_source;
#if defined(USE_PETSC)
   PETScMatrix<real_t>    *_A, *_CM, *_Df;
   PETScVect<real_t>      *_b, *_u, *_bc, *_bf, *_sf, *_pf, *_v, *_w, *_LM, _uu;
#else
   Matrix<real_t>         *_A, *_CM, *_Df;
   Vect<real_t>           *_b, *_u, *_bc, *_bf, *_sf, *_pf, *_v, *_w, *_LM, _uu;
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#if defined(USE_PETSC)
   void setMatrix(PETScMatrix<real_t> &A);
#else
   void setMatrix(SkSMatrix<real_t> &A);
   void setMatrix(SkMatrix<real_t> &A);
   void setMatrix(SpMatrix<real_t> &A);
#endif

   bool SolverIsSet() const;
   bool isConstantMatrix() const;
   bool isConstantMesh() const;
   void setConstantMatrix();
   void setConstantMesh();
   virtual void setTerms(int opt);
   void set_exprtk();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
