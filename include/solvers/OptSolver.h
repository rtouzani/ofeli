/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2019 Rachid Touzani

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

                          Definition of class 'OptSolver'

  ==============================================================================*/

#ifndef __OPT_SOLVER_H
#define __OPT_SOLVER_H

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
#include "MyOpt.h"

namespace OFELI {

/*! \file OptSolver.h
 *  \brief Definition file for class OptSolver.
 */

//-----------------------------------------------------------------------------
// Class OptSolver
//-----------------------------------------------------------------------------

/*! \class OptSolver
 * \ingroup Solver
 * \brief To solve an optimization problem with bound constraints
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
class OptSolver
{

 public:

/*! \enum OptMethod
 * \brief Choose optimization algorithm
 */
   enum OptMethod {
      GRADIENT             = 0,   /*!< Gradient method                                  */
      TRUNCATED_NEWTON     = 1,   /*!< Truncated Newton method                          */
      SIMULATED_ANNEALING  = 2,   /*!< Simulated annealing global optimization method   */
      NELDER_MEAD          = 3,   /*!< Nelder-Mead global optimization method           */
   };

//----------------------------   BASIC OPERATIONS   ----------------------------

/// \brief Default constructor
    OptSolver();

/** \brief Constructor using vector of optimization variables
 *  @param [in] x Vector having as size the number of optimization variables.
 *  It contains the initial guess for the optimization algorithm.
 *  @remark After using the member function run, the vector \a x contains the
 *  obtained solution if the optimization procedure was successful
 */
    OptSolver(Vect<real_t>& x);

/** \brief Constructor using vector of optimization variables
 *  @param [in] opt Reference to instance of user defined optimization class.
 *  This class inherits from abstract class MyOpt. It must contain the member function
 *      \c double \c Objective(Vect<double> &x)
 *  which returns the value of the objective for a given solution vector \c x. The user
 *  defined class must contain, if the optimization algorithm requires it the member function
 *      \c Gradient(Vect<double> &x, Vect<double> &g)
 *  which stores the gradient of the objective in the vector \c g for a given optimization 
 *  vector \c x.
 *  The user defined class must also contain, if the optimization algorithm requires it the
 *  member function
 *  @param [in] x Vector having as size the number of optimization variables.
 *  It contains the initial guess for the optimization algorithm.
 *  @remark After using the member function run, the vector \a x contains the
 *  obtained solution if the optimization procedure was successful
 */
    OptSolver(MyOpt&        opt,
              Vect<real_t>& x);

/// \brief Destructor
    ~OptSolver();

//-------------------------------   MODIFIERS  ---------------------------------

/// \brief Return the total number of function evaluations.
    int getNbFctEval() const { return _nb_obj_eval; }

/** \brief Choose optimization method
 *  @param [in] m Enumerated value to choose the optimization algorithm to use.
 *                Must be chosen among the enumerated values:
 *                <ul>
 *                <li>GRADIENT: Gradient steepest descent method with projection 
 *                              for bounded constrained problems
 *                <li>TRUNCATED_NEWTON: The Nash's Truncated Newton Algorithm, due to
 *                              S.G. Nash (Newton-type Minimization via the 
 *                              Lanczos method, SIAM J. Numer. Anal. 21 (1984) 770-778).
 *                <li>SIMULATED_ANNEALING: Global optimization simulated annealing 
 *                    method. See Corana et al.'s article: "Minimizing Multimodal Functions 
 *                    of Continuous Variables with the Simulated Annealing Algorithm" in the
 *                    September 1987 (vol. 13, no. 3, pp. 262-280)
 *                    issue of the ACM Transactions on Mathematical Software.
 *                <li>NELDER_MEAD: Global optimization Nelder-Mead method due to
 *                    John Nelder, Roger Mead (A simplex method for function minimization,
 *                    Computer Journal, Volume 7, 1965, pages 308-313). As implemented by
 *                    R. ONeill (Algorithm AS 47: Function Minimization Using a Simplex Procedure,
 *                    Applied Statistics, Volume 20, Number 3, 1971, pages 338-345).
 *                 </ul>
 */
    void setOptMethod(OptMethod m);

/** \brief Prescribe boundary conditions as constraints
 *  \details This member function is useful in the case of optimization problems
 *  where the optimization variable vector is the solution of a partial differential
 *  equation. For this case, Dirichlet boundary conditions can be prescribed as
 *  constraints for the optimization problem
 *  @param [in] bc Vector containing the values to impose on degrees of freedom.
 *  This vector must have been constructed using the Mesh instance.
 *  @remark Only degrees of freedom with positive code are taken into account as
 *  prescribed
 */
    void setBC(const Vect<real_t>& bc);

/** \brief Define the objective function to minimize by an algebraic expression
 *  @param [in] exp Regular expression defining the objective function
 */
    void setObjective(string exp);

/** \brief Define a component of the gradient of the objective function to minimize
 *  by an algebraic expression
 *  @param [in] exp Regular expression defining the objective function
 *  @param [in] i Component of gradient [Default: <tt>1</tt>]
 */
    void setGradient(string exp,
                     int    i=1);

/** \brief Choose user defined optimization class
 *  @param [in] opt Reference to inherited user specified optimization class
 */
    void setOptClass(MyOpt& opt) { _opt = &opt; }

/** \brief Define upper bound for optimization variable
 *  \details Case of a one-variable problem
 *  @param [in] ub Upper bound
 */
    void setUpperBound(real_t ub);

/** \brief Define upper bounds for optimization variables
 *  @param [in] ub Vector containing upper values for variables
 */
    void setUpperBounds(Vect<real_t>& ub);

/** \brief Define lower bound for optimization variable
 *  \details Case of a one-variable problem
 *  @param [in] lb Lower value
 */
    void setLowerBound(real_t lb);

/** \brief Set verbosity parameter
 * @param [in] verb Verbosity parameter
 */
    void setVerbosity(int verb) { _verb = verb; }

/** \brief Define lower bounds for optimization variables
 * @param [in] lb Vector containing lower values for variables
 */
    void setLowerBounds(Vect<real_t>& lb);

/** \brief Set Simulated annealing options
 *  @remark This member function is useful only if simulated annealing is used.
 *  @param [in] rt The temperature reduction factor. The value suggested by
 *  Corana et al. is \a .85. See Goffe et al. for more advice.
 *  @param maxevl [in] The maximum number of function evaluations. If it is
 *  exceeded, the return \a code=1.
 *  @param [in] ns Number of cycles. After \a ns*nb_var function evaluations, each
 *  element of \a vm is adjusted so that approximately half of all function evaluations
 *  are accepted. The suggested value is 20.
 *  @param [in] nt Number of iterations before temperature reduction. After
 *  \a nt*ns*n function evaluations, temperature \a (t) is changed
 *  by the factor \a rt. Value suggested by Corana et al. is
 *  \a max(100,5*nb_var). See Goffe et al. for further advice.
 *  @param [in] neps Number of final function values used to decide upon termination.
 *  See \c eps. Suggested value is \a 4
 *  @param [in] maxevl
 *  @param [in] t The initial temperature. See Goffe et al. for advice.
 *  @param [in] vm The step length vector. On input it should encompass the
 *  region of interest given the starting value \a x. For point
 *  x[i], the next trial point is selected is from \a x[i]-vm[i]
 *  to \a x[i]+vm[i]. Since \a vm is adjusted so that about half
 *  of all points are accepted, the input value is not very
 *  important (i.e. is the value is off, \a OptimSA adjusts \a vm to the
 *  correct value).
 *  @param [in] xopt
 *  @param [in] fopt
 */
    void setSAOpt(real_t        rt,
                  int           ns,
                  int           nt,
                  int&          neps,
                  int           maxevl,
                  real_t        t,
                  Vect<real_t>& vm,
                  Vect<real_t>& xopt,
                  real_t&       fopt);

/** \brief Set error tolerance
 *  @param [in] toler Error tolerance for termination. If the final function
 *  values from the last neps temperatures differ from the
 *  corresponding value at the current temperature by less than
 *  eps and the final function value at the current temperature
 *  differs from the current optimal function value by less than
 *  toler, execution terminates and the value \a 0 is returned.
 */
    void setTolerance(real_t toler) { _toler = toler; }

/// \brief Set maximal number of iterations
    void setMaxIterations(int n) { _max_it = n; }

/// \brief Return number of objective function evaluations
    int getNbObjEval() const { return _nb_obj_eval; }

/// \brief Return the final temperature
/// \details This function is meaningful only if the Simulated Annealing algorithm is used
    real_t getTemperature() const { return _t; }

/// @brief Return the number of accepted objective function evaluations
/// \details This function is meaningful only if the Simulated Annealing algorithm is used
    int getNbAcc() const { return _nacc; }

/// \brief Return the total number of trial function evaluations that
/// would have been out of bounds
/// \details This function is meaningful only if the Simulated Annealing algorithm is used
    int getNbOutOfBounds() const { return _nobds; }

/// \brief Return Optimal value of the objective
    real_t getOptObj() const { return _fopt; }

/** \brief Run the optimization algorithm
 *  \details This function runs the optimization procedure using default values for parameters.
 *  To modify these values, user the function run with arguments
 */
    int run();

/** \brief Run the optimization algorithm
 *  @param [in] toler Tolerance value for convergence testing
 *  @param [in] max_it Maximal number of iterations to achieve convergence
 */
    int run(real_t toler,
            int    max_it);

/** \brief Return solution in the case of a one variable optimization
 *  \details In the case of a one variable problem, the solution value is 
 *  returned, if the optimization procedure was successful
 */
    real_t getSolution() const { return (*_x)[0]; }

/** \brief Get solution vector
 *  \details The vector \a x contains the solution of the optimization problem. Note that if
 *  the constructor using an initial vector was used, the vector will contain the solution 
 *  once the member function run has beed used (If the optimization procedure was successful)
 *  @param [out] x solution vector 
 */
    void getSolution(Vect<real_t>& x) const { x = *_x; }

/// Output class information
    friend ostream & operator<<(ostream&         s,
                                const OptSolver& os);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    friend int OptimTN(OptSolver&, Vect<real_t>&, const Vect<real_t>&, const Vect<real_t>&,
                       int&, int&, int&, real_t, int);
    friend int lmqnbc(OptSolver&, Vect<real_t>&, real_t&, Vect<real_t>&, const Vect<real_t>&,
                      const Vect<real_t>&, vector<int>&, int, int&, size_t, real_t&, real_t&,
                      real_t&, real_t&, int&, int&);
    friend int modlnp(OptSolver&, int, Vect<real_t>&, Vect<real_t>&, Vect<real_t>&, Vect<real_t>&,
                      Vect<real_t>&, Vect<real_t>&, Vect<real_t>&, Vect<real_t>& g, Vect<real_t>&,
                      int, int&, int&, int&, int&, real_t&, real_t&, real_t&, int&, int,
                      vector<int>&, real_t, real_t&, real_t&, real_t&, int, Vect<real_t>*);
    friend void gtims(OptSolver&, const Vect<real_t>&, Vect<real_t>&, const Vect<real_t>&,
                      const Vect<real_t>&, Vect<real_t>*, int&, real_t&, real_t, real_t);
    friend int linder(OptSolver&, real_t, real_t&, real_t&, real_t&, real_t&, real_t&,
                      const vector<real_t>&, real_t, Vect<real_t>&, real_t&, real_t&,
                      Vect<real_t>&, int&, Vect<real_t>*);
    friend int OptimPG(OptSolver&, Vect<real_t>&, const Vect<real_t>&, const Vect<real_t>&, int&,
                       int&, int&, real_t, int);
    friend int OptimNM(OptSolver&, Vect<real_t>&, real_t&, real_t, Vect<real_t>&, int, int, int&, int&);
    friend int OptimSA(OptSolver&, Vect<real_t>&, real_t&, real_t& , int&, int&, int&, int&,
                       const Vect<real_t>&, const Vect<real_t>&, const Vect<real_t>&, int&, real_t&,
                       Vect<real_t>&, real_t&, int&, int&, int&);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 private:

   size_t _size;
   int _verb, _eval, _nacc, _ns, _nt, _nobds, _nb_obj_eval, _nb_grad_eval;
   int _conv, _max_eval, _neps, _nb_restart, _max_it;
   vector<int> _pivot;
   Vect<real_t> *_x, _c, _vm, _lb, _ub, _g, _step;
   OptMethod _opt_method;
   real_t _accrcy, _toler, _reqmin;
   real_t _f, _rt, _fopt, _t;
   MyOpt *_opt;
   bool _exp, _sa_opt, _tn_opt, _obj_type, _x_set, _method_set;
   string _exp_obj, _var;
   Vect<string> _exp_grad, _exp_hess;

   void setTNDefaults();
   void setSADefaults();
   void setNMDefaults();
   void setPGDefaults();
   real_t *obj_function(const Vect<real_t>& x);
   void *obj_gradient(const Vect<real_t>& x, Vect<real_t>& g);

   real_t Objective(Vect<real_t>& x);
   void Gradient(Vect<real_t>& x, Vect<real_t>& g);

};

/// \fn ostream & operator<<(ostream& s, const OptSolver &de)
/// \brief Output differential system information
/// \ingroup Solver
    ostream & operator<<(ostream&         s,
                         const OptSolver& os);

} /* namespace OFELI */

#endif
