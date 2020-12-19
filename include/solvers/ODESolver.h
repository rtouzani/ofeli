/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

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

                          Definition of class 'ODESolver'

  ==============================================================================*/

#ifndef __ODE_SOLVER_H
#define __ODE_SOLVER_H

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
#include "linear_algebra/Vect.h"
#include "solvers/LinearSolver.h"
#include "solvers/Iter.h"
#include "equations/AbsEqua.h"
#include "io/Fct.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file ODESolver.h
 *  \brief Definition file for class ODESolver.
 */


//-----------------------------------------------------------------------------
// Class ODESolver
//-----------------------------------------------------------------------------

/*! \class ODESolver
 * \ingroup Solver
 * \brief To solve a system of ordinary differential equations
 * \details The class ODESolver enables solving by a numerical scheme a system
 * or ordinary differential equations taking one of the forms:
 * <ul>
 *    <li>A linear system of differential equations of the first-order:\n
 *          <b>  A<sub>1</sub>(t)u'(t) + A<sub>0</sub>(t)u(t) = f(t) </b>
 *    <li>A linear system of differential equations of the second-order:\n
 *          <b>  A<sub>2</sub>(t)u''(t) + A<sub>1</sub>(t)u'(t) + A<sub>0</sub>(t)u(t) = f(t)</b>
 *    <li>A system of ordinary differential equations of the form:\n
 *          <b> u'(t) = f(t,u(t))</b>
 * </ul>
 * 
 * The following time integration schemes can be used:
 * <ul>
 *    <li>Forward Euler scheme (value: \a FORWARD_EULER) for first-order systems
 *    <li>Backward Euler scheme (value: \a BACKWARD_EULER) for first-order linear systems
 *    <li>Crank-Nicolson (value: \a CRANK_NICOLSON) for first-order linear systems
 *    <li>Heun (value: \a HEUN) for first-order systems
 *    <li>2nd Order Adams-Bashforth (value: \a AB2) for first-order systems
 *    <li>4-th order Runge-Kutta (value: \a RK4) for first-order systems
 *    <li>2nd order Backward Differentiation Formula (value: \a BDF2) for linear first-order systems
 *    <li>Newmark (value: \a NEWMARK) for linear second-order systems with constant matrices
 * </ul>
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class ODESolver
{

 public:

//----------------------------   BASIC OPERATIONS   ----------------------------

/// \brief Default constructor
    ODESolver();

/// \brief Constructor providing the number of equations
    ODESolver(size_t nb_eq);

/** \brief Constructor using time discretization data
 *  @param [in] s Choice of the scheme: To be chosen in the enumerated variable
 *  \a Scheme (see the presentation of the class)
 *  @param [in] time_step Value of the time step. This value will be modified
 *  if an adaptive method is used. The default value for this parameter
 *  if the value given by the global variable \c theTimeStep
 *  @param [in] final_time Value of the final time (time starts at 0). The default
 *  value for this parameter is the value given by the global variable \c theFinalTime
 *  @param [in] nb_eq Number of differential equations (size of the system)
 *  [Default: <tt>1</tt>]
 */
    ODESolver(TimeScheme s,
              real_t     time_step=theTimeStep,
              real_t     final_time=theFinalTime,
              size_t     nb_eq=1);

/// \brief Destructor
    ~ODESolver();

//-------------------------------   MODIFIERS  ---------------------------------

/** \brief Define data of the differential equation or system
 *  @param [in] s Choice of the scheme: To be chosen in the enumerated variable
 *  \a Scheme (see the presentation of the class)
 *  @param [in] time_step Value of the time step. This value will be modified
 *  if an adaptive method is used. The default value for this parameter
 *  if the value given by the global variable \c theTimeStep
 *  @param [in] final_time Value of the final time (time starts at 0). The default
 *  value for this parameter is the value given by the global variable \c theFinalTime
 */
    void set(TimeScheme s,
             real_t     time_step=theTimeStep,
             real_t     final_time=theFinalTime);

/** \brief Set the number of equations [Default: <tt>1</tt>]
 *  \details This function is to be used if the default constructor was used
 */
    void setNbEq(size_t nb_eq) { _nb_eq = 1; }

/** \brief Define coefficients in the case of a scalar differential equation
 *  \details This function enables giving coefficients of the differential
 *  equation as an algebraic expression of time \a t (see the function fparse)
 *  @param [in] a0 Coefficient of the 0-th order term
 *  @param [in] a1 Coefficient of the 1-st order term
 *  @param [in] a2 Coefficient of the 2-nd order term
 *  @param [in] f Value of the right-hand side
 *  @note Naturally, the equation is of the first order if \a a2=0
 */
    void setCoef(real_t a0,
                 real_t a1,
                 real_t a2,
                 real_t f);

/** \brief Define coefficients in the case of a scalar differential equation
 *  @param [in] a0 Coefficient of the 0-th order term
 *  @param [in] a1 Coefficient of the 1-st order term
 *  @param [in] a2 Coefficient of the 2-nd order term
 *  @param [in] f Value of the right-hand side
 *  @note Naturally, the equation if of the first order if \a a2=0
 */
    void setCoef(string a0,
                 string a1,
                 string a2,
                 string f);

/// \brief Claim that ODE is linear
/// \details Claim that the defined ODE (or system of ODEs) is linear
    void setLinear();

/** \brief Set time derivative, given as an algebraic expression, for a nonlinear ODE
 *  \details This function enables prescribing the value of the 1-st derivative
 *  for a 1st order ODE or the 2nd one for a 2nd-order ODE. It is to be
 *  used for nonlinear ODEs of the form 
 *  y'(t) = f(t,y(t)) or y''(t) = f(t,y(t),y'(t))\n
 *  In the case of a system of ODEs, this function can be called once for each equation, 
 *  given in the order of the unknowns
 *  @param [in] f Expression of the function
 */
    void setF(string f);

/** \brief Set time derivative, given as an algebraic expression, for a nonlinear ODE
 *  \details This function enables prescribing the value of the 1-st derivative
 *  for a 1st order ODE or the 2nd one for a 2nd-order ODE. It is to be
 *  used for nonlinear ODEs of the form 
 *  y'(t) = f(t,y(t)) or y''(t) = f(t,y(t),y'(t))\n
 *  This function is to be used for the <tt>i</tt>-th equation of a system of ODEs
 *  @param [in] f Expression of the function
 *  @param [in] i Index of equation. Must be not larger than the number of equations
 */
    void setF(string f,
              int    i);

/** \brief Set time derivative of the function defining the ODE
 *  \details This function enables prescribing the value of the 1-st derivative
 *  for a 1st order ODE or the 2nd one for a 2nd-order ODE. It is to be
 *  used for nonlinear ODEs of the form 
 *  y'(t) = f(t,y(t)) or y''(t) = f(t,y(t),y'(t))\n
 *  In the case of a system of ODEs, this function can be called once for each equation, 
 *  given in the order of the unknowns
 */
    void setDF(string df,
               int    i,
               int    j);

/** \brief Set time derivative of the function defining the ODE
 *  \details This function enables prescribing the value of the 1-st derivative
 *  for a 1st order ODE or the 2nd one for a 2nd-order ODE. It is to be
 *  used for nonlinear ODEs of the form 
 *  y'(t) = f(t,y(t)) or y''(t) = f(t,y(t),y'(t))\n
 *  In the case of a system of ODEs, this function can be called once for each equation, 
 *  given in the order of the unknowns
 */
    void setdFdt(string df,
                 int    i);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setF(Fct& f);
    void setDF(Fct& df, int i=1, int j=1);
    void setdFdt(Fct& df, int i=1);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Set intermediate right-hand side vector for the Runge-Kutta method
 *  @param [in] f Value of right-hand side
 */
    void setRK4RHS(real_t f) { _d01 = f; }

/** \brief Set intermediate right-hand side vector for the Runge-Kutta method
 *  @param [in] f right-hand side vector
 */
    void setRK4RHS(Vect<real_t> &f);

/** \brief Set initial condition for a first-oder system of differential equations
 *  @param [in] u Vector containing initial condition for the unknown
 */
    void setInitial(Vect<real_t>& u);

/** \brief Set initial condition for a first-oder system of differential equations
 *  @param [in] u Initial condition for an unknown
 *  @param [in] i Index of the unknown
 */
    void setInitial(real_t u,
                    int    i);

/** \brief Set initial condition for a second-order system of differential equations
 *  \details Giving the right-hand side at initial time is somtimes required
 *  for high order methods like Runge-Kutta
 *  @param [in] u Vector containing initial condition for the unknown
 *  @param [in] v Vector containing initial condition for the time derivative of the unknown
 */
    void setInitial(Vect<real_t>& u,
                    Vect<real_t>& v);

/** \brief Set initial RHS for a system of differential equations
 *  \details Giving the right-hand side at initial time is somtimes required
 *  for high order methods like Runge-Kutta
 *  @param [in] f Vector containing right-hand side at initial time. This 
 *  vector is helpful for high order methods
 */
    void setInitialRHS(Vect<real_t>& f);

/** \brief Set initial condition for a second-order ordinary differential equation
 *  @param [in] u Initial condition (unknown) value
 *  @param [in] v Initial condition (time derivative of the unknown) value
 */
    void setInitial(real_t u,
                    real_t v);

/** \brief Set initial condition for a first-order ordinary differential equation
 *  @param [in] u Initial condition (unknown) value
 */
    void setInitial(real_t u);

/** \brief Set initial right-hand side for a single differential equation
 *  @param [in] f Value of right-hand side at initial time. This 
 *  value is helpful for high order methods
 */
    void setInitialRHS(real_t f);

/** \brief Define matrices for a system of first-order ODEs
 *  \details Matrices are given as references to class DMatrix.
 *  @param [in] A0 Reference to matrix in front of the 0-th order term (no time derivative)
 *  @param [in] A1 Reference to matrix in front of the 1-st order term (first time derivative)
 *  @remark This function has to be called at each time step
 */
    void setMatrices(DMatrix<real_t>& A0,
                     DMatrix<real_t>& A1);

/** \brief Define matrices for a system of second-order ODEs
 *  \details Matrices are given as references to class DMatrix.
 *  @param [in] A0 Reference to matrix in front of the 0-th order term (no time derivative)
 *  @param [in] A1 Reference to matrix in front of the 1-st order term (first time derivative)
 *  @param [in] A2 Reference to matrix in front of the 2-nd order term (second time derivative)
 *  @remark This function has to be called at each time step
 */
    void setMatrices(DMatrix<real_t>& A0,
                     DMatrix<real_t>& A1,
                     DMatrix<real_t>& A2);

/** \brief Define matrices for an implicit nonlinear system of first-order ODEs
 *  \details The system has the nonlinear implicit form 
 *       a1(u)' + a0(u) = 0 
 *  Vectors a0, a1 are given as references to class Vect.
 *  @param [in] a0 Reference to vector in front of the 0-th order term (no time derivative)
 *  @param [in] a1 Reference to vector in front of the 1-st order term (first time derivative)
 *  @remark This function has to be called at each time step
 */
    void seODEVectors(Vect<real_t>& a0,
                      Vect<real_t>& a1);

/** \brief Define matrices for an implicit nonlinear system of second-order ODEs
 *  \details The system has the nonlinear implicit form 
 *       a2(u)'' + a1(u)' + a0(u) = 0 
 *  Vectors a0, a1, a2 are given as references to class Vect.
 *  @param [in] a0 Reference to vector in front of the 0-th order term (no time derivative)
 *  @param [in] a1 Reference to vector in front of the 1-st order term (first time derivative)
 *  @param [in] a2 Reference to vector in front of the 2-nd order term (second time derivative)
 *  @remark This function has to be called at each time step
 */
    void seODEVectors(Vect<real_t>& a0,
                      Vect<real_t>& a1,
                      Vect<real_t>& a2);

/// \brief Set right-hand side vector for a system of ODE
/// @param [in] b Vect instance containing right-hand side for a linear system
/// of ordinary differential equations
    void setRHS(Vect<real_t>& b);

/// \brief Set right-hand side for a linear ODE
/// @param [in] f Value of the right-hand side for a linear ordinary differential equation
    void setRHS(real_t f);

/// \brief Set right-hand side value for a linear ODE
    void setRHS(string f);

/** \brief Define parameters for the Newmarxk scheme
 *  @param [in] beta Parameter beta [Default: <tt>0.25</tt>]
 *  @param [in] gamma Parameter gamma [Default: <tt>0.5</tt>]
 */
    void setNewmarkParameters(real_t beta,
                              real_t gamma);

/** \brief Say that matrix problem is constant
 *  \details This is useful if the linear system is solved by
 *  a factorization method but has no effect otherwise
 */
    void setConstantMatrix() { _constant_matrix = true; }

/** \brief Say that matrix problem is variable
 *  \details This is useful if the linear system is solved by
 *  a factorization method but has no effect otherwise
 */
    void setNonConstantMatrix() { _constant_matrix = false; }

/** \brief Set linear solver data
 *  @param [in] s Solver identification parameter.
 *  To be chosen in the enumeration variable Iteration:\n 
 *  DIRECT_SOLVER, CG_SOLVER, CGS_SOLVER, BICG_SOLVER, BICG_STAB_SOLVER,
 *  GMRES_SOLVER, QMR_SOLVER  [Default: <tt>DIRECT_SOLVER</tt>]
 *  @param [in] p Preconditioner identification parameter.
 *  To be chosen in the enumeration variable Preconditioner:\n
 *  IDENT_PREC, DIAG_PREC, ILU_PREC  [Default: <tt>DIAG_PREC</tt>]
 *  @note The argument \a p has no effect if the solver is DIRECT_SOLVER
 */
    void setLinearSolver(Iteration      s=DIRECT_SOLVER,
                         Preconditioner p=DIAG_PREC)
   { _s = s; _p = p; }

/** \brief Set maximal number of iterations
 *  \details This function is useful for a non linear ODE (or system of ODEs)
 *  if an implicit scheme is used
 *  @param [in] max_it Maximal number of iterations [Default: <tt>100</tt>]
 */
    void setMaxIter(int max_it) { _iter.setMaxIter(max_it); }

/** \brief Set tolerance value for convergence
 *  \details This function is useful for a non linear ODE (or system of ODEs)
 *  if an implicit scheme is used
 *  @param [in] toler Tolerance value [Default: <tt>1.e-8</tt>]
 */
    void setTolerance(real_t toler) { _iter.setTolerance(toler); }

/// \brief Run one time step
/// @return Value of new time step if this one is updated
    real_t runOneTimeStep();

/** \brief Run the time stepping procedure
 *  @param [in] opt Flag to say if problem matrix is constant while time stepping 
 *  (true) or not (Default value is false)
 *  @note This argument is not used if the time stepping scheme is explicit
 */
    void run(bool opt=false);

//-----------------------------   INSPECTORS  ----------------------------------

/// \brief Return number of equations
    size_t getNbEq() const { return _nb_eq; }

/// \brief Return LinearSolver instance
    LinearSolver<real_t> &getLSolver() { return _ls; }

/** \brief Get time derivative of solution
 *  \details Return approximate time derivative of solution in the case
 *  of a single equation
 *  @param [in] i Index of component whose time derivative is sought
 *  @return Time derivative of the i-th component of the solution
 *  @remark If we are solving one equation, this parameter is not used.
 */
    real_t getTimeDerivative(int i=1) const;

/** \brief Get time derivative of solution (for a system)
 *  \details Get approximate time derivative of solution in the case
 *  of an ODE system
 *  @param [out] y Vector containing time derivative of solution
 */
    void getTimeDerivative(Vect<real_t>& y) const;

/// \brief Return solution in the case of a scalar equation
    real_t get() const { return _y2; }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    friend ostream & operator<<(ostream&         s,
                                const ODESolver& de);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 private:

   enum DEType {
      SCALAR_LINEAR  = 0,      /*!< Linear Scalar Differential Equation        */
      SCALAR_NL      = 1,      /*!< Nonlinear Scalar Differential Equation     */
      VECTOR_LINEAR  = 2,      /*!< Linear Differential system of equations    */
      VECTOR_NL      = 3,      /*!< Nonlinear Differential system of equations */
   };

   AbsEqua<real_t> *_theEqua;
   size_t _order, _nb_eq, _nb_ssteps, _step, _sstep;
   int _sc;
   Iteration _s;
   Preconditioner _p;
   bool _fct_allocated, _dF_computed, _setF_called, _RK4_rhs;
   bool _linear, _a0, _a1, _a2, _constant_matrix, _regex, _explicit, _init, _lhs, _rhs, _rhsF;
   Vect<real_t> _x, _u, _v, *_w, _f0, _f1, _f2, _b, _f01, _f, *_bc, _bb, _vv, _dudt;
   Vect<real_t> *_du, _ddu, _dv, _ddv, _vF1, _vF2, _vF, _vDF1, _D, _k1, _k2, _k3, _k4;
   DMatrix<real_t> *_A0, *_A1, *_A2;
   real_t _time_step0, _time_step, _time, _final_time, _c0, _c1, _c2;
   real_t _y0, _y1, _dy1, _y2, _dy2, _ddy, _d0, _d1, _d2, _d01, _dydt;
   vector<string> _expA0, _expA1, _expA2, _var;
   vector<Fct *> _theF, _theDF, _thedFdt;
   vector<Fct> _theC;
   vector<real_t> _xv;
   real_t _beta, _gamma;
   LinearSolver<real_t> _ls;
   DEType _type;
   Iter<real_t> _iter;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   typedef void (ODESolver::* TSPtr)();
   static TSPtr TS[11];
   TSPtr _solve;

   typedef Vect<real_t>& (ODESolver::* FPtr)();
   static FPtr F[11];
   FPtr _set_f;

   std::map<TimeScheme,int> _sch;
   void setScheme();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

   void solveForwardEuler();
   void solveBackwardEuler();
   void solveCrankNicolson();
   void solveHeun();
   void solveAB2();
   void solveLeapFrog();
   void solveRK4();
   void solveRK3_TVD();
   void solveRKC();
   void solveNewmark();
   void solveBDF2();

// Functions to return the right-hand side for each time integration scheme
   Vect<real_t>& setF_ForwardEuler();
   Vect<real_t>& setF_BackwardEuler();
   Vect<real_t>& setF_CrankNicolson();
   Vect<real_t>& setF_Heun();
   Vect<real_t>& setF_AB2();
   Vect<real_t>& setF_LeapFrog();
   Vect<real_t>& setF_RK4();
   Vect<real_t>& setF_RK3_TVD();
   Vect<real_t>& setF_Newmark();
   Vect<real_t>& setF_BDF2();

   real_t eval(real_t t, real_t y, Fct* f);
   real_t eval(real_t t, const Vect<real_t>& y, Fct* f);

};

/// \fn ostream & operator<<(ostream& s, const ODESolver &de)
/// \brief Output differential system information
/// \ingroup Solver
    ostream & operator<<(ostream&         s,
                         const ODESolver& de);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
