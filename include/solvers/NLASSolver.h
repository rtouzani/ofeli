/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2025 Rachid Touzani

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

                      Definition of for class 'NLASSolver'

  ==============================================================================*/

#ifndef __NLAS_SOLVER_H
#define __NLAS_SOLVER_H

#include <stdlib.h>
#include <math.h>
#include <vector>
#include <functional>
using std::vector;
using std::function;

#include <iostream>
using std::ostream;
using std::endl;
using std::cerr;

#include <iomanip>
using std::setw;

#include <string>
using std::string;

#include "OFELI_Config.h"
#include "equations/Equa.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/Matrix.h"
#include "solvers/MyNLAS.h"
#include "io/Fct.h"


namespace OFELI {

/*! \file NLASSolver.h
 *  \brief Definition file for class NLASSolver.
 */


/*! \class NLASSolver
 * \ingroup Solver
 * \brief To solve a system of nonlinear algebraic equations of the form
 *       f(u) = 0
 * 
 * \details Features:
 * <ul>
 *    <li>The nonlinear problem is solved by the Newton's method in the general case,
 *        and in the one variable case, either by the bisection or the Regula Falsi method
 *    <li>The function and its gradient are given:
 *        <ul>
 *          <li>Either by regular expressions
 *          <li>Or by user defined functions
 *          <li>Or by a user defined class. This feature enables defining the function 
 *              and its gradient through a PDE class for instance
 *        </ul>
 * </ul>
 */

class NLASSolver
{

 public:

//----------------------------   BASIC OPERATIONS   ----------------------------

/// \brief Default constructor.
    NLASSolver();

/** \brief Constructor defining the iterative method to solve the equation
 *  @param [in] nl Choose an iterative procedure to solve the nonlinear system of equations:
 *  To be chosen among the enumerated values: <tt>BISECTION</tt>, <tt>REGULA_FALSI</tt> or 
 *  <tt>NEWTON</tt>.
 *  @param [in] nb_eq Number of equations [Default: <tt>1</tt>]
 */
    NLASSolver(NonLinearIter nl,
               size_t        nb_eq=1);

/** \brief Constructor defining a one-variable problem
 *  @param [in] x Variable containing on input initial guess and on output
 *  solution, if convergence is achieved
 *  @param [in] nl Iterative procedure to solve the nonlinear system of equations:
 *  To be chosen among the enumerated values:  <tt>BISECTION</tt>, <tt>REGULA_FALSI</tt> or 
 *  <tt>NEWTON</tt>. 
 */
    NLASSolver(real_t&       x,
               NonLinearIter nl=NEWTON);

/** \brief Constructor defining a multi-variable problem
 *  @param [in] x Variable containing on input initial guess and on output
 *  solution, if convergence is achieved
 *  @param [in] nl Iterative procedure to solve the nonlinear system of equations:
 *  The only possible value (default one) in the current version is <tt>NEWTON</tt>. 
 */
    NLASSolver(Vect<real_t>& x,
               NonLinearIter nl=NEWTON);

/** \brief Constructor using a user defined class
 *  @param [in] my_nlas Reference to instance of user defined class.
 *  This class inherits from abstract class MyNLAS. It must contain the member function
 *      \c Vect<double> Function(const Vect<double>& x)
 *  which returns the value of the nonlinear function, as a vector, for a given solution vector 
 * \c x. The user defined class must contain, if the iterative scheme requires it the member function
 *      \c Vect<double> Gradient(const Vect<real_t>& x)
 *  which returns the gradient as a \c n*n vector, each index \c (i,j) containing the j-th partial 
 * derivative of the i-th function.
 *  @param [in] nl Iterative procedure to solve the nonlinear system of equations:
 *  To be chosen among the enumerated values:  <tt>BISECTION</tt>, <tt>REGULA_FALSI</tt> or 
 *  <tt>NEWTON</tt>. 
 */
    NLASSolver(MyNLAS&       my_nlas,
               NonLinearIter nl=NEWTON);

/// \brief Destructor
    ~NLASSolver();

//-------------------------------   MODIFIERS  ---------------------------------

/// \brief Set Maximal number of iterations
/// \details Default value of this parameter is 100
    void setMaxIter(int max_it) { _max_it = max_it; }

/// \brief Set tolerance value for convergence
/// \details Default value of this parameter is 1.e-8
    void setTolerance(real_t toler) { _toler = toler; }

/** \brief Define an iterative procedure
 *  To be chosen among the enumerated values: <tt>BISECTION</tt>, <tt>REGULA_FALSI</tt> or 
 *  <tt>NEWTON</tt>. 
 */
    void set(NonLinearIter nl);

/** \brief Define number of equations
 */
    void setNbEq(size_t nb_eq);

/** \brief Define the function associated to the equation to solve
 *  \details This function can be used in the case where a user defined function is
 *  to be given. To be used in the one-variable case.
 *  @param [in] f Function given as a function of one real variable and returning
 *  a real number. This function can be defined by the calling program as a
 *  C-function and then cast to an instance of class function
 */
    void setFunction(function<real_t(real_t)> f);

/** \brief Define the function associated to the equation to solve
 *  \details This function can be used in the case where a user defined function is
 *  to be given.
 *  @param [in] f Function given as a function of many variables, stored in an input vector,
 *  and returns a vector. This function can be defined by the calling program as a
 *  C-function and then cast to an instance of class function
 */
    void setFunction(function<Vect<real_t>(Vect<real_t>)> f);

/** \brief Define the function associated to the derivative of the equation to solve 
 *  @param [in] g Function given as a function of one real variable and returning
 *  a real number. This function can be defined by the calling program as a
 *  C-function and then cast to an instance of class function
 */
    void setGradient(function<real_t(real_t)> g);

/** \brief Define the function associated to the gradient of the equation to solve 
 *  @param [in] g Function given as a function of many variables, stored in an input vector.
 *  and returns a \c n*n vector (\n n is the number of variables). This function can be defined
 *  by the calling program as a C-function and then cast to an instance of class function
 */
    void setGradient(function<Vect<real_t>(Vect<real_t>)> g);

/** \brief Set function for which zero is sought (case of one equation)
 *  @param [in] exp Regular expression defining the function using the symbol \c x as
 *  a variable
 */
    void setf(string exp);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setf(Fct& f);
    void setDf(Fct& df, size_t i=1, size_t j=1);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Set pzrtial derivative of function for which zero is sought (case of many equations)
 *  @param [in] exp Regular expression defining the partial derivative. In this expression, the
 *  variables are \c x1, \c x2, ... \c x10 (up to 10 variables)
 *  @param [in] i Component of function [Default: <tt>=1</tt>]
 *  @param [in] j Index of the partial derivative [Default: <tt>=1</tt>]
 */
    void setDf(string exp,
               size_t i=1,
               size_t j=1);

/** \brief Define a PDE
 *  \details The solver can be used to solve a nonlinear PDE. In this case, the PDE is
 *  defined as an instance of a class inheriting of Equa.
 *  @param [in] eq Pointer to equation instance
 */
    void setPDE(Equa& eq);

/** \brief Set initial guess for the iterations
 *  @param [in] u Vector containing initial guess for the unknown
 */
    void setInitial(Vect<real_t> &u);

/** \brief Set initial guess for a unique unknown
 *  @param [in] x Rference to value of initial guess
 */
    void setInitial(real_t& x);

/** \brief Set initial guesses bisection or Regula falsi algorithms
 *  @param [in] a Value of first initial guess
 *  @param [in] b Value of second initial guess
 *  @note The function has to have opposite signs at these values
 *  <tt>i.e.</tt> f(a)f(b)<0.
 *  @warning This function makes sense only in the case of a unique function
 *  of one variable
 */
    void setInitial(real_t a,
                    real_t b);

/// \brief Run the solution procedure
    int run();

//-----------------------------   INSPECTORS  ----------------------------------

/// \brief Return solution (Case of a scalar equation)
    real_t get() const { return _y; }

/// \brief Return solution (case of a nonlinear system of equations)
/// @param [out] u Vector that contains on output the solution
    void get(Vect<real_t> &u) const;

/// \brief Return number of iterations
    int getNbIter() const { return _nb_it; }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    friend ostream & operator<<(ostream& s, const NLASSolver &nl);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

private:

   enum {
     FUNCTION   = 0,
     EXPRESSION = 1,
     CLASS      = 2
   };

   bool _fct_allocated, _df_computed, _cv, _f_given, _df_given, _u_set, _ab_given;
   Equa *_theEqua;
   int _nl, _max_it, _nl_it, _it, _nb_it, _fct_type;
   Vect<real_t> *_u, _v, _f, _w;
   real_t _toler, *_x, _y, _a, _b, _g;
   Matrix<real_t> *_Df;
   Mesh *_theMesh;
   MyNLAS *_my_nlas;
   size_t _nb_eq, _nb_fct_def;
   function<real_t(real_t)> _fct1, _grad1; 
   function<Vect<real_t>(Vect<real_t>)> _fct, _grad;
   vector<Fct *> _theFct, _theDFct;
   typedef void (NLASSolver::* NLPtr)();
   static NLPtr NL[7];
   NLPtr _nlp;

   real_t Function(const Vect<real_t>& u, int i=1);
   real_t Function(real_t x);
   void Gradient(const Vect<real_t>& u, size_t i=1, size_t j=1);
   real_t Gradient(real_t x);
   void solveBisection();
   void solveRegulaFalsi();
   void solvePicard();
   void solveSecant();
   void solveNewton();
};

//-----------------------------------------------------------------------------
// Associated functions
//-----------------------------------------------------------------------------

/// \fn ostream & operator<<(ostream& s, const NLASSolver &nl)
/// \brief Output nonlinear system information
/// \ingroup Solver
    ostream & operator<<(ostream&          s,
                         const NLASSolver& nl);

} /* namespace OFELI */

#endif
