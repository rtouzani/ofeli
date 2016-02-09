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

       Definition of abstract class 'AbsEqua' for Finite Element Equations

  ==============================================================================*/


#ifndef __ABS_EQUA_H
#define __ABS_EQUA_H

#include <float.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#ifdef USE_PETSC
#include "linear_algebra/petsc/PETScMatrix.h"
#endif

#include "OFELI_Config.h"
#include "mesh/Mesh.h"
#include "mesh/Element.h"
#include "mesh/Side.h"

#if defined (USE_PETSC)
#include "linear_algebra/petsc/PETScMatrix.h"
#endif

#include "linear_algebra/DMatrix.h"
#include "linear_algebra/DSMatrix.h"
#include "linear_algebra/SpMatrix.h"
#include "linear_algebra/SkMatrix.h"
#include "linear_algebra/SkSMatrix.h"
#include "linear_algebra/BMatrix.h"
#include "linear_algebra/TrMatrix.h"

#include "linear_algebra/LocalMatrix.h"
#include "linear_algebra/LocalVect.h"
#include "mesh/Material.h"
#include "io/UserData.h"
#include "solvers/LinearSolver.h"
#include "util/Gauss.h"
#include "io/Prescription.h"
#include "solvers/EigenProblemSolver.h"

namespace OFELI {

/*! \defgroup Equation General Purpose Equations
 *  \brief Gathers equation related classes
 */

extern Material theMaterial;
class TimeStepping;
class Mesh;
class EigenProblemSolver;
class Prescription;

/*! \file AbsEqua.h
 *  \brief Definition file for abstract class AbsEqua.
 */

/*! \enum PDE_Terms
 * Enumerate variable that selects various terms in partial differential equations
 */
enum PDE_Terms {
   MASS                =  0,              /*!< Consistent mass term                    */
   CONSISTENT_MASS     =  0,              /*!< Consistent mass term                    */
   LUMPED_MASS         =  0x00001000,     /*!< Lumped mass term                        */
   CAPACITY            =  0x00002000,     /*!< Consistent capacity term                */
   CONSISTENT_CAPACITY =  0x00002000,     /*!< Consistent capacity term                */
   LUMPED_CAPACITY     =  0x00004000,     /*!< Lumped capacity term                    */
   VISCOSITY           =  0x00008000,     /*!< Viscosity term                          */
   STIFFNESS           =  0x00010000,     /*!< Stiffness term                          */
   DIFFUSION           =  0x00020000,     /*!< Diffusion term                          */
   CONVECTION          =  0x00040000,     /*!< Convection term                         */
   DEVIATORIC          =  0x00080000,     /*!< Deviatoric term                         */
   DILATATION          =  0x00100000,     /*!< Dilatational term                       */
   ELECTRIC            =  0x00200000,     /*!< Electric term                           */
   MAGNETIC            =  0x00400000,     /*!< Magnetic term                           */
   LOAD                =  0x00800000,     /*!< Body load term                          */
   HEAT_SOURCE         =  0x01000000,     /*!< Body heat source term                   */
   BOUNDARY_TRACTION   =  0x02000000,     /*!< Boundary traction (pressure) term       */
   HEAT_FLUX           =  0x04000000,     /*!< Boundary heat flux term                 */
   CONTACT             =  0x08000000,     /*!< Signorini contact                       */
   BUOYANCY            =  0x10000000,     /*!< Buoyancy force term                     */
   LORENTZ_FORCE       =  0x20000000      /*!< Lorentz force term                      */
};


/*! \enum EqDataType
 * Enumerate variable that selects equation data type
 */
enum EqDataType {
   INITIAL_FIELD                =    1,    /*!< Initial condition                      */
   SOLUTION                     =    1,    /*!< Solution vector (same as Initial)      */
   INITIAL_AUX_1                =    2,    /*!< Initial auxiliary field                */
   INITIAL_AUX_2                =    3,    /*!< Initial auxiliary field                */
   INITIAL_AUX_3                =    4,    /*!< Initial auxiliary field                */
   INITIAL_AUX_4                =    5,    /*!< Initial auxiliary field                */
   BOUNDARY_CONDITION           =    6,    /*!< Boundary condition data                */
   BODY_FORCE                   =    7,    /*!< Body force data                        */
   SOURCE                       =    7,    /*!< Source data (same as Body force)       */
   POINT_FORCE                  =    8,    /*!< Localized (at point) force             */
   BOUNDARY_FORCE               =    9,    /*!< Boundary force data                    */
   FLUX                         =    9,    /*!< Flux data (same as Boundary force)     */
   TRACTION                     =    9,    /*!< Traction data (same as Boundary force) */
   AUX_INPUT_FIELD_1            =   10,    /*!< Auxiliary input field 1                */
   AUX_INPUT_FIELD_2            =   11,    /*!< Auxiliary input field 2                */
   AUX_INPUT_FIELD_3            =   12,    /*!< Auxiliary input field 3                */
   AUX_INPUT_FIELD_4            =   13,    /*!< Auxiliary input field 4                */
   DISPLACEMENT_FIELD           =   14,    /*!< A displacement field                   */
   VELOCITY_FIELD               =   15,    /*!< A velocity field                       */
   TEMPERATURE_FIELD            =   16     /*!< A temperature field                    */
};


/*! \enum ArrayType
 * Selects local or global option for array as argument.
 */
enum ArrayType {
   LOCAL_ARRAY  = 0,   /*!< For a local array labeled with local numbering */
   GLOBAL_ARRAY = 1    /*!< For a local array labeled with global numbering */
};


/*! \enum TimeScheme
 * Selects time integration scheme
 */
enum TimeScheme {
   STATIONARY      =  0,    /*!< No time scheme: stationary                */
   FORWARD_EULER   =  1,    /*!< Forward Euler scheme (Explicit)           */
   BACKWARD_EULER  =  2,    /*!< Backward Euler scheme (Implicit)          */
   CRANK_NICOLSON  =  3,    /*!< Crank-Nicolson scheme                     */
   HEUN            =  4,    /*!< Heun scheme                               */
   NEWMARK         =  5,    /*!< Newmark scheme                            */
   LEAP_FROG       =  6,    /*!< Leap Frog scheme                          */
   ADAMS_BASHFORTH =  7,    /*!< Adams-Bashforth scheme (2nd Order)        */
   AB2             =  7,    /*!< Adams-Bashforth scheme (2nd Order)        */
   RUNGE_KUTTA     =  8,    /*!< 4-th Order Runge-Kutta scheme (4th Order) */
   RK4             =  8,    /*!< 4-th Order Runge-Kutta scheme             */
   RK3_TVD         =  9,    /*!< 3-rd Order Runge-Kutta TVD scheme         */
   BDF2            = 10     /*!< Backward Difference Formula (2nd Order)   */
};

/*! \enum PDE
 * Choose partial differential equation to `
 */
enum PDE {
   LAPLACE,                      /*!< Laplace equation                        */
   DIFFUSION_CONVECTION,         /*!< Diffusion Convection equation           */
   THERMAL_PHASE_CHANGE,         /*!< Thermal phase change problem (Stefan)   */
   INCOMPRESSIBLE_NAVIER_STOKES, /*!< Incompressible Navier-Stokes equations  */
   LINEARIZED_ELASTICITY,        /*!< Linearized elasticity equations         */
   PLANAR_TRUSS,                 /*!< 2-D truss equation                      */
   SPATIAL_BEAM                  /*!< 3-D beam equations                      */
};

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

/*! \enum AnalysisType
 * Choose analysis type
 */
enum AnalysisType {
   STEADY_STATE,                 /*!< Steady state analysis                   */
   TRANSIENT,                    /*!< Transient analysis                      */
   OPTIMIZATION                  /*!< Optimization analysis                   */
};


/*! \class AbsEqua
 * \ingroup Equation
 * \brief Mother abstract class to describe equation.
 *
 * \tparam <T_> Data type (real_t, float, complex<real_t>, ...)
 */

template<class T_> class SkMatrix;
template<class T_> class SkSMatrix;
template<class T_> class SpMatrix;

template<class T_>
class AbsEqua
{

 public:

/// \brief Default constructor
    AbsEqua() : _theMesh(NULL), _time_scheme(STATIONARY), _analysis(STEADY_STATE),
                _solver(-1), _step(0), _nb_fields(1), _final_time(0.),
                _A(NULL), _bc(NULL), _bf(NULL), _sf(NULL), _u(NULL), _v(NULL),
                _eigen(false), _sol_given(false), _bc_given(false),
                _init_given(false), _bf_given(false), _sf_given(false), _set_matrix(false)
    {
       _ls.setSolver(DIRECT_SOLVER);
    }

/// \brief Constructor with mesh instance
    AbsEqua(Mesh &mesh) : _theMesh(&mesh), _time_scheme(STATIONARY), _analysis(STEADY_STATE),
                          _solver(-1), _step(0), _nb_fields(1), _final_time(0.),
                          _A(NULL), _bc(NULL), _bf(NULL), _sf(NULL), _u(NULL), _v(NULL),
                          _eigen(false), _sol_given(false), _bc_given(false), _init_given(false),
                          _bf_given(false), _sf_given(false), _set_matrix(false)
    {  }

/// \brief Destructor
    virtual ~AbsEqua()
    {
       if (_A)
          delete _A;
    }

/// \brief Define mesh and renumber DOFs after removing imposed ones
    void setMesh(Mesh &m) { _theMesh = &m; _theMesh->removeImposedDOF(); }

/// \brief Return reference to Mesh instance
/// @return Reference to Mesh instance
    Mesh & getMesh() const { return *_theMesh; }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** \brief Build equation.
 *  This member function is to be used if one wants to use a class that inherits
 *  from AbsEqua that handles the whole solution process (including assembly).
 */
    virtual void build() { }
    virtual void build(EigenProblemSolver& e) { }
    virtual void build(TimeStepping& s) { }
    virtual int runOneTimeStep() { return 0; }
    virtual int run() { return 0; }
    virtual int runSteadyState() { return run(); }
    virtual int runTransient() { return 0; }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** \brief Set transient analysis settings
 *  \details Define a set of parameters for time integration
 *  @param [in] s Scheme
 *  @param [in] p1
 *  @param [in] p2
 *  @param [in] p3
 */
    void setTransient(int    s,
                      real_t p1=1,
                      real_t p2=1,
                      real_t p3=1)
    {
       _time_scheme = s;
       _time_parameter1 = p1;
       _time_parameter2 = p2;
       _time_parameter3 = p3;
       if (s==FORWARD_EULER)
          _theta = 0;
       else if (s==BACKWARD_EULER)
          _theta = 1;
       else if (s==CRANK_NICOLSON)
          _theta = 0.5;
   }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

// \brief Return number of (element or side) nodes.
//    virtual size_t getNbNodes() const;

// \brief Return number of (element or side) equations.
//    virtual size_t getNbEq() const;

/// \brief Return reference to linear solver instance
    LinearSolver<T_> &getLinearSolver() { return _ls; }

/** \brief Choose solver for the linear system
 *  @param [in] ls Solver of the linear system.
 *  To choose among the enumerated values: <tt>DIRECT_SOLVER</tt>, <tt>CG_SOLVER</tt>,
    <tt>GMRES_SOLVER</tt>
    <ul>
      <li> <tt>DIRECT_SOLVER</tt>, Use a facorization solver [default]
      <li> <tt>CG_SOLVER</tt>, Conjugate Gradient iterative solver
      <li> <tt>CGS_SOLVER</tt>, Squared Conjugate Gradient iterative solver
      <li> <tt>BICG_SOLVER</tt>, BiConjugate Gradient iterative solver
      <li> <tt>BICG_STAB_SOLVER</tt>, BiConjugate Gradient Stabilized iterative solver
      <li> <tt>GMRES_SOLVER</tt>, GMRES iterative solver
      <li> <tt>QMR_SOLVER</tt>, QMR iterative solver
    </ul>
 *  @param [in] pc Preconditioner to associate to the iterative solver.
 *  If the direct solver was chosen for the first argument this argument is not used.
 *  Otherwise choose among the enumerated values:
    <ul> 
      <li> <tt>IDENT_PREC</tt>, Identity preconditioner (no preconditioning [default])
      <li> <tt>DIAG_PREC</tt>, Diagonal preconditioner
      <li> <tt>ILU_PREC</tt>, Incomplete LU factorization preconditioner
    </ul>
 */
    void setSolver(Iteration      ls,
                   Preconditioner pc=IDENT_PREC)
    {
       if (_set_matrix == false) {
          _matrix_type = SPARSE;
          _A = new SpMatrix<T_>;
          _A->setMesh(*_theMesh,0);
       }
       _ls.setSolver(ls,pc);
    }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Set terms in equations
    virtual void setWithConvection(int f) { f = 0; }
    void setTerms(PDE_Terms t) { _terms = t; }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/** \brief Solve the linear system
 *  @param [in] A Pointer to matrix of the system (Instance of class SpMatrix)
 *  @param [in] b Vector containing right-hand side
 *  @param [in,out] x Vector containing initial guess of solution on input, actual solution on output
 */
    int SolveLinearSystem(Matrix<T_>* A,
                          Vect<T_>&   b,
                          Vect<T_>&   x)
    {
       _set_matrix = true;
       _ls.setMatrix(A);
       _ls.setRHS(b);
       _ls.setSolution(x);
       if (_ls.getSolver()<0) {
          if (_matrix_type==SKYLINE)
             _ls.setSolver(DIRECT_SOLVER);
          else if (_matrix_type==(SPARSE|SYMMETRIC))
             _ls.setSolver(CG_SOLVER);
           else if (_matrix_type==SPARSE)
             _ls.setSolver(GMRES_SOLVER);
      }
      return _ls.solve();
    }


#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Set Analysis Type
    void setAnalysis(int a) { _analysis = a; }

/// \brief Set time iteration index
    void setTimeIndex(size_t step) { _step = step; }

/// \brief Set time step duration
    void setTimeStep(real_t t) { _time_step = t; }

/// \brief Return time step
    real_t getTimeStep() const { return _time_step; }

/// \brief Set current time
    void setTime(real_t t) { _time = t; }

/// \brief Set final time
    void setFinalTime(real_t t) { _final_time = t; }

/// \brief Set verbose parameter
    void setVerbose(int v) { _verbose = v; }

/// \brief Return number of fields
    size_t getNbFields() const { return _nb_fields; }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setTimeIntegration(int scheme) { _time_scheme = scheme; }
    int getTimeIntegration() const { return _time_scheme; }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
    
#ifndef DOXYGEN_SHOULD_SKIP_THIS

/// \brief Return equation name
    string getEquationName() const { return _equation_name; }

/// \brief Return Finite Element Type
    string getFiniteElementType() const { return _finite_element; }

/** \brief Set equation input data
 *  @param [in] opt Parameter that selects data type for input. This parameter
 *  is to be chosen in the enumerated variable EqDataType
 *  @param [in] u Vect instance that contains input vector data
 *  List of data types contains <tt>INITIAL_FIELD</tt>, <tt>BOUNDARY_CONDITION_DATA</tt>, 
 *  <tt>SOURCE_DATA</tt> or <tt>FLUX</tt> with obvious meaning
 */
    virtual void setInput(EqDataType opt,
                          Vect<T_>&  u)
    {
       if (opt==INITIAL_FIELD) {
          _u = &u;
          _init_given = true;
       }
       else if (opt==SOLUTION) {
          _u = &u;
          _sol_given = true;
       }
       else if (opt==BOUNDARY_CONDITION) {
          _bc = &u;
          _bc_given = true;
       }
       else if (opt==SOURCE || opt==BODY_FORCE) {
          _bf = &u;
          _bf_given = true;
       }
       else if (opt==FLUX || opt==TRACTION) {
          _sf = &u;
          _sf_given = true;
       }
       else
          ;
    }

/// \brief Set prescription
/// @param [in] p Prescription instance
    void set(Prescription& p) { _prescription = &p; }

/// \brief Choose tolerance value
    void setTolerance(real_t toler) { _toler = toler; }

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
   
protected :
#ifndef DOXYGEN_SHOULD_SKIP_THIS
   Mesh               *_theMesh;
   int                _field_type, _time_scheme, _analysis, _verbose, _terms;
   int                _matrix_type, _solver, _max_it;
   size_t             _step, _nb_fields, _nb_eigv;
   real_t             _theta, _alpha, _beta, _time_step, _time, _final_time;
   real_t             _time_parameter1, _time_parameter2, _time_parameter3;
   EigenProblemSolver _ev;
   Matrix<T_>         *_A, *_CM;
   Vect<T_>           *_bc, *_bf, *_sf, *_u, *_v, *_b, *_LM, _uu;
   bool               _eigen, _sol_given, _bc_given, _init_given, _bf_given, _sf_given; 
   bool               _constant_matrix, _constant_mesh, _set_matrix;
   int                _sol_type, _init_type, _bc_type, _bf_type, _sf_type;
   string             _equation_name, _finite_element;
   LinearSolver<T_>   _ls;
   real_t             _toler;
   Prescription       *_prescription;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   void setMatrix(SkSMatrix<T_> &A)
   {
      _A = &A;
      _matrix_type = SKYLINE|SYMMETRIC;
      _A->setMesh(*_theMesh);
   }

   void setMatrix(SkMatrix<T_> &A)
   {
      _A = &A;
      _matrix_type = SKYLINE;
      _A->setMesh(*_theMesh);
   }

   void setMatrix(SpMatrix<T_> &A)
   {
      _A = &A;
      _matrix_type = SPARSE;
      _A->setMesh(*_theMesh);
   }
   bool isConstantMatrix() const { return _constant_matrix; }
   bool isConstantMesh() const { return _constant_mesh; }
   void setConstantMatrix() { _constant_matrix = true; }
   void setConstantMesh() { _constant_mesh = true; }
   virtual void setTerms(int opt) { _terms = opt; }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

} /* namespace OFELI */

#endif
