/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2018 Rachid Touzani

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

                       Definition of class 'TimeStepping'

  ==============================================================================*/

#ifndef __TIME_STEPPING_H
#define __TIME_STEPPING_H

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
#include "equations/AbsEqua.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/Assembly.h"

#define MAX_NB_PDES  10

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

template class LinearSolver<real_t>;

/*! \file TimeStepping.h
 *  \brief Definition file for class TimeStepping.
 */

//-----------------------------------------------------------------------------
// Class TimeStepping
//-----------------------------------------------------------------------------

/*! \class TimeStepping
 * \ingroup Solver
 * \brief To solve time stepping problems, i.e. systems of linear ordinary 
 * differential equations of the form
 *       [A2]{y"} + [A1]{y'} + [A0]{y} = {b}
 * 
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 *
 * \details Features:
 * <ul>
 *    <li>The system may be first or second order (first and/or second order
 *        time derivatives
 *    <li>The following time integration schemes can be used:
 *      <ul>
 *         <li>For first order systems: The following schemes are implemented
 *             Forward Euler (value: \a FORWARD_EULER)\n
 *             Backward Euler (value: \a BACKWARD_EULER)\n
 *             Crank-Nicolson (value: \a CRANK_NICOLSON)\n
 *             Heun (value: \a HEUN)\n
 *             2nd Order Adams-Bashforth (value: \a AB2)\n
 *             4-th order Runge-Kutta (value: \a RK4)\n
 *             2nd order Backward Differentiation Formula (value: \a BDF2)
 *         <li>For second order systems: The following schemes are implemented
 *             Newmark (value: \a NEWMARK)
 *      </ul>
 * </ul>
 */

class TimeStepping
{

 public:

//----------------------------   BASIC OPERATIONS   ----------------------------

/// \brief Default constructor
    TimeStepping();

/** \brief Constructor using time discretization data
 *  @param [in] s Choice of the scheme: To be chosen in the enumerated variable
 *  \a TimeScheme (see the presentation of the class)
 *  @param [in] time_step Value of the time step. This value will be modified
 *  if an adaptive method is used. The default value for this parameter
 *  if the value given by the global variable \c theTimeStep
 *  @param [in] final_time Value of the final time (time starts at 0). The default
 *  value for this parameter is the value given by the global variable \c theFinalTime
 */
    TimeStepping(TimeScheme s,
                 real_t     time_step=theTimeStep,
                 real_t     final_time=theFinalTime);

/// \brief Destructor
    ~TimeStepping();

//-------------------------------   MODIFIERS  ---------------------------------

/** \brief Define data of the differential equation or system
 *  @param [in] s Choice of the scheme: To be chosen in the enumerated variable
 *  \a TimeScheme (see the presentation of the class)
 *  @param [in] time_step Value of the time step. This value will be modified
 *  if an adaptive method is used. The default value for this parameter
 *  if the value given by the global variable \c theTimeStep
 *  @param [in] final_time Value of the final time (time starts at 0). The default
 *  value for this parameter is the value given by the global variable \c theFinalTime
 */
    void set(TimeScheme s,
             real_t     time_step=theTimeStep,
             real_t     final_time=theFinalTime);

/** \brief Set reference to LinearSolver instance
 *  @param [in] ls Reference to LinearSolver instance
 */
    void setLinearSolver(LinearSolver<real_t>& ls) { _ls = &ls; }

/** \brief Define partial differential equation to solve
 *  \details The used equation class must have been constructed using the Mesh instance
 *  @param [in] eq Reference to equation instance
 *  @param [in] nl Toggle to say if the considered equation is linear [Default: <tt>0</tt>] or not
 */
    void setPDE(AbsEqua<real_t>& eq,
                bool             nl=false);

/** \brief Set intermediate right-hand side vector for the Runge-Kutta method
 *  @param [in] f Vector containing the RHS
 */
    void setRK4RHS(Vect<real_t>& f);

/** \brief Set intermediate right-hand side vector for the TVD Runge-Kutta 3 method
 *  @param [in] f Vector containing the RHS
 */
    void setRK3_TVDRHS(Vect<real_t> &f);

/** \brief Set initial condition for the system of differential equations
 *  @param [in] u Vector containing initial condition for the unknown
 *  @remark If a second-order differential equation is to be solved, use the
 *  the same function with two initial vectors (one for the unknown, the second
 *  for its time derivative)
 */
    void setInitial(Vect<real_t>& u);

/** \brief Set initial condition for a system of differential equations
 *  @param [in] u Vector containing initial condition for the unknown
 *  @param [in] v Vector containing initial condition for the time derivative
 *  of the unknown
 *  @note This function can be used to provide solution at previous time step
 *  if a restarting procedure is used.
 *  @note This member function is to be used only in the case of a second
 *  order system
 */
    void setInitial(Vect<real_t>& u,
                    Vect<real_t>& v);

/** \brief Set initial RHS for a system of differential equations when the used
 *  scheme requires it
 *  \details Giving the right-hand side at initial time is somtimes required
 *  for high order methods like Runge-Kutta
 *  @param [in] f Vector containing right-hand side at initial time. This 
 *  vector is helpful for high order methods
 *  @note This function can be used to provide solution at previous time step
 *  if a restarting procedure is used.
 */
    void setInitialRHS(Vect<real_t>& f);

/// \brief Set right-hand side vector
    void setRHS(Vect<real_t>& b);

/// \brief Set vector containing boundary condition to enforce
    void setBC(Vect<real_t>& u);

/** \brief Define parameters for the Newmark scheme
 *  @param [in] beta Parameter beta [Default: <tt>0.25</tt>]
 *  @param [in] gamma Parameter gamma [Default: <tt>0.5</tt>]
 */
    void setNewmarkParameters(real_t beta,
                              real_t gamma);

/** \brief Say that matrix problem is constant
 *  \details This is useful if the linear system is solved by
 *  a factorization method but has no effect otherwise
 */
    void setConstantMatrix() { _constant_matrix[_nb_pdes-1] = true; }

/** \brief Say that matrix problem is variable
 *  \details This is useful if the linear system is solved by
 *  a factorization method but has no effect otherwise
 */
    void setNonConstantMatrix() { _constant_matrix[_nb_pdes-1] = false; }

/** \brief Set linear solver data
 *  @param [in] s Solver identification parameter.
 *  To be chosen in the enumeration variable Iteration:\n 
 *  DIRECT_SOLVER, CG_SOLVER, CGS_SOLVER, BICG_SOLVER, BICG_STAB_SOLVER,
 *  GMRES_SOLVER, QMR_SOLVER [Default: <tt>DIRECT_SOLVER</tt>]
 *  @param [in] p Preconditioner identification parameter.
 *  To be chosen in the enumeration variable Preconditioner:\n
 *  IDENT_PREC, DIAG_PREC, ILU_PREC [Default: <tt>DIAG_PREC</tt>]
 *  @note The argument \a p has no effect if the solver is DIRECT_SOLVER
 */
    void setLinearSolver(Iteration      s=DIRECT_SOLVER,
                         Preconditioner p=DIAG_PREC);

/** \brief Set vectors defining a nonlinear first order system of ODEs
 *  \details The ODE system has the form
 *     a1(u)' + a0(u) = 0
 *  @param [in] a0 Vect instance defining the 0-th order term
 *  @param [in] A0 Matrix instance
 */
    void setNLTerm0(Vect<real_t>&   a0,
                    Matrix<real_t>& A0);

/** \brief Set vectors defining a nonlinear second order system of ODEs
 *  \details The ODE system has the form
 *     a2(u)'' + a1(u)' + a0(u) = 0
 *  @param [in] a0 Vect instance defining the 0-th order term
 *  @param [in] a1 Vect instance defining the first order term
 *  @param [in] a2 Vect instance defining the second order term
 */
    void setNLTerm(Vect<real_t>& a0,
                   Vect<real_t>& a1,
                   Vect<real_t>& a2);

/** \brief Set verbosity parameter:
 *  \details
 *  <ul>
 *     <li> = 0, No output
 *     <li> = 1, Print step label and time value
 *     <li> = 2, Print step label, time value, time step and integration scheme
 *  </ul>
 */
    void setVerbose(int v=0) { _verb = v; }

/// \brief Run one time step
/// @return Value of new time step if this one is updated
    real_t runOneTimeStep();

/** \brief Run the time stepping procedure
 *  @param [in] opt Flag to say if problem matrix is constant while time stepping 
 *  (true) or not (Default value is false)
 *  @note This argument is not used if the time stepping scheme is explicit
 */
    void run(bool opt=false);

/** \brief Assemble element arrays into global matrix and right-hand side
 *  \details This member function is to be called from finite element equation
 *  classes
 *  @param [in] el Reference to Element class
 *  @param [in] b Pointer to element right-hand side
 *  @param [in] A0 Pointer to matrix of 0-th order term (involving no time derivative)
 *  @param [in] A1 Pointer to matrix of first order term (involving time first derivative)
 *  @param [in] A2 Pointer to matrix of second order term (involving time second derivative)
 *  [Default: <tt>NULL</tt>]

 */
    void Assembly(const Element& el,
                  real_t*        b,
                  real_t*        A0,
                  real_t*        A1,
                  real_t*        A2=NULL);

/** \brief Assemble side arrays into global matrix and right-hand side
 *  \details This member function is to be called from finite element equation
 *  classes
 *  @param [in] sd Reference to Side class
 *  @param [in] b Pointer to side right-hand side
 *  @param [in] A Pointer to matrix [Default: <tt>NULL</tt>]
 */
    void SAssembly(const Side& sd,
                   real_t*     b,
                   real_t*     A=NULL);

//-----------------------------   INSPECTORS  ----------------------------------

/// \brief Return LinearSolver instance
    LinearSolver<real_t> &getLSolver() { return *_ls; }

    friend ostream & operator<<(ostream&            s,
                                const TimeStepping& ts);

 private:

   AbsEqua<real_t> *_theEqua[MAX_NB_PDES];
   Mesh *_theMesh[MAX_NB_PDES];
   size_t _order, _nb_eq[MAX_NB_PDES], _nb_dof[MAX_NB_PDES], _nb_ssteps, _step, _sstep;
   int _nb_pdes, _ind, _verb, _sc, _non_linear, _max_it;
   Iteration _s[MAX_NB_PDES];
   Preconditioner _p[MAX_NB_PDES];
   bool _constant_matrix[MAX_NB_PDES], _regex, _explicit[MAX_NB_PDES], _set_bc[MAX_NB_PDES], _nl[MAX_NB_PDES];
   Vect<real_t> _u[MAX_NB_PDES], _v[MAX_NB_PDES], *_w[MAX_NB_PDES], _f0[MAX_NB_PDES], _f1[MAX_NB_PDES];
   Vect<real_t> *_f2[MAX_NB_PDES], _b[MAX_NB_PDES], *_f01[MAX_NB_PDES], _f[MAX_NB_PDES], *_bc[MAX_NB_PDES];
   Vect<real_t> _bb[MAX_NB_PDES], _vv[MAX_NB_PDES];
   Vect<real_t> *_du[MAX_NB_PDES], _ddu[MAX_NB_PDES], _dv[MAX_NB_PDES], _ddv[MAX_NB_PDES], _D[MAX_NB_PDES];
   Vect<real_t> _k1[MAX_NB_PDES], _k2[MAX_NB_PDES], _k3[MAX_NB_PDES], _k4[MAX_NB_PDES];
   Vect<real_t> *_a0, *_a1, *_a2;
   Matrix<real_t> *_A[MAX_NB_PDES];
   real_t _time_step0, _time_step, _time, _final_time, _c0, _c1, _c2, _toler;
   real_t _beta, _gamma, _nl_toler;
   int _max_nl_it;
   LinearSolver<real_t> *_ls;

   typedef void (TimeStepping::* TSPtr)();
   static TSPtr TS[11];
   TSPtr _solve;

   typedef void (TimeStepping::* ASPtr)(const Element&,real_t*,real_t*,real_t*,real_t*);
   static ASPtr AS[11];
   ASPtr _assemb;

   typedef void (TimeStepping::* ASSPtr)(const Side&,real_t*,real_t*);
   static ASSPtr ASS[11];
   ASSPtr _sassemb;

   typedef void (TimeStepping::* MassPtr)();
   static MassPtr MASS[2];
   MassPtr _mass_opt;

   typedef Vect<real_t>& (TimeStepping::* RHSPtr)();
   static RHSPtr RHS[11];
   RHSPtr _set_rhs;

   typedef void (TimeStepping::* PreSolvePtr)();
   static PreSolvePtr PS[11];
   PreSolvePtr _presolve;

// Functions to solve the linear system issued from time integration schemes
   void solveStationary();
   void solveForwardEuler();
   void solveBackwardEuler();
   void solveCrankNicolson();
   void solveHeun();
   void solveAB2();
   void solveLeapFrog();
   void solveRK4();
   void solveRK3_TVD();
   void solveNewmark();
   void solveBDF2();

// Functions to assemble the linear system for time integration schemes
   void AssembleStationary(const Element& el,
                           real_t*        eb,
                           real_t*        eA0,
                           real_t*        eA1,
                           real_t*        eA2=NULL);
   void SAssembleStationary(const Side& sd,
                            real_t*     sb,
                            real_t*     sA=NULL);
   void AssembleForwardEuler(const Element& el,
                             real_t*        eb,
                             real_t*        eA0,
                             real_t*        eA1,
                             real_t*        eA2=NULL);
   void SAssembleForwardEuler(const Side& sd,
                              real_t*     sb,
                              real_t*     sA=NULL);
   void AssembleBackwardEuler(const Element& el,
                              real_t*        eb,
                              real_t*        eA0,
                              real_t*        eA1,
                              real_t*        eA2=NULL);
   void SAssembleBackwardEuler(const Side& sd,
                               real_t*     sb,
                               real_t*     sA=NULL);
   void AssembleCrankNicolson(const Element& el,
                              real_t*        eb,
                              real_t*        eA0,
                              real_t*        eA1,
                              real_t*        eA2=NULL);
   void SAssembleCrankNicolson(const Side& sd,
                               real_t*     sb,
                               real_t*     sA=NULL);
   void AssembleHeun(const Element& el,
                     real_t*        eb,
                     real_t*        eA0,
                     real_t*        eA1,
                     real_t*        eA2=NULL);
   void SAssembleHeun(const Side& sd,
                      real_t*     sb,
                      real_t*     sA=NULL);
   void AssembleAB2(const Element& el,
                    real_t*        eb,
                    real_t*        eA0,
                    real_t*        eA1,
                    real_t*        eA2=NULL);
   void SAssembleAB2(const Side& sd,
                     real_t*     sb,
                     real_t*     sA=NULL);
   void AssembleLeapFrog(const Element& el,
                         real_t*        eb,
                         real_t*        eA0,
                         real_t*        eA1,
                         real_t*        eA2=NULL);
   void SAssembleLeapFrog(const Side& d,
                          real_t*     sb,
                          real_t*     sA=NULL);
   void AssembleRK4(const Element& el,
                    real_t*        eb,
                    real_t*        eA0,
                    real_t*        eA1,
                    real_t*        eA2=NULL);
   void AssembleRK3_TVD(const Element& el,
                        real_t*        eb,
                        real_t*        eA0,
                        real_t*        eA1,
                        real_t*        eA2=NULL);
   void SAssembleRK4(const Side& sd,
                     real_t*     sb,
                     real_t*     sA=NULL);
   void SAssembleRK3_TVD(const Side& sd,
                         real_t*     sb,
                         real_t*     sA=NULL);
   void AssembleNewmark(const Element& el,
                        real_t*        eb,
                        real_t*        eA0,
                        real_t*        eA1,
                        real_t*        eA2=NULL);
   void SAssembleNewmark(const Side& sd,
                         real_t*     sb,
                         real_t*     eA=NULL);
   void AssembleBDF2(const Element& el,
                     real_t*        eb,
                     real_t*        eA0,
                     real_t*        eA1,
                     real_t*        eA2=NULL);
   void SAssembleBDF2(const Side& sd,
                      real_t*     sb,
                      real_t*     sA=NULL);

// Functions to return the right-hand side for each time integration scheme
   Vect<real_t>& setRHS_Stationary();
   Vect<real_t>& setRHS_ForwardEuler();
   Vect<real_t>& setRHS_BackwardEuler();
   Vect<real_t>& setRHS_CrankNicolson();
   Vect<real_t>& setRHS_Heun();
   Vect<real_t>& setRHS_AB2();
   Vect<real_t>& setRHS_LeapFrog();
   Vect<real_t>& setRHS_RK4();
   Vect<real_t>& setRHS_RK3_TVD();
   Vect<real_t>& setRHS_Newmark();
   Vect<real_t>& setRHS_BDF2();

// Functions to compute a predictor
   void PreSolve_Stationary() { }
   void PreSolve_ForwardEuler() { }
   void PreSolve_BackwardEuler() { }
   void PreSolve_CrankNicolson() { }
   void PreSolve_Heun() { }
   void PreSolve_AB2() { }
   void PreSolve_LeapFrog() { }
   void PreSolve_RK4() { }
   void PreSolve_RK3_TVD() { }
   void PreSolve_Newmark();
   void PreSolve_BDF2() { }

   void eval(real_t t);
   void insertBC(const Vect<real_t>& b, Vect<real_t>& v);
   void insertBC0(const Vect<real_t>& b, Vect<real_t>& v);
};

/** \fn ostream & operator<<(ostream& s, const TimeStepping &ts)
 * \brief Output differential system information
 * \ingroup Solver
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
    ostream & operator<<(ostream&            s,
                         const TimeStepping& ts);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
