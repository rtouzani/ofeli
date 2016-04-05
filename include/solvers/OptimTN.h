/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2014 Rachid Touzani

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
        Definition of function 'OptimTN' for truncared Newton optimization
  ==============================================================================*/

#ifndef __OPTIM_TN_H
#define __OPTIM_TN_H

#include "OFELI_Config.h"
#include <float.h>
#include <iostream>
using std::ostream;
using std::endl;

#include <iomanip>
using std::setw;
using std::ios;
using std::setprecision;

#include "linear_algebra/Vect.h"
using std::abs;

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/* ************************************************************************* */
/*                                                                           */
/*                           TRUNCATED NEWTON METHOD                         */
/*                                                                           */
/*   Written by:  Stephen G. Nash                                            */
/*                Operations Research and Applied Statistics Dept.           */
/*                George Mason University                                    */
/*                Fairfax, VA 22030                                          */
/*                                                                           */
/* ************************************************************************* */

/*! \file OptimTN.h
 *  \ingroup Solver
 *  \brief Function to solve an optimization problem using the Truncated Newton method.
 *
 */

/** \fn int OptimTN(OPT_ &theOpt, Vect<real_t> &x, Vect<real_t> &low, Vect<real_t> &up,
                    Vect<int> &pivot, int max_it=-1, real_t toler=100*OFELI_EPSMCH, int msg_lvl=5)
    \ingroup Solver
   \brief Truncated Newton optimization solver.
  
   \details Solves a bounds-constrained optimization problem using the Nash's Truncated 
   Newton Algorithm (See paper by S.G. Nash, Newton-type Minimization via the 
   Lanczos method, SIAM J. Numer. Anal. 21 (1984) 770-778). All vector variables 
   are instances of class Vect<real_t>.

   @param [in] theOpt Instance of class \b OPT_ that is implemented by the user and that provides the objective function.
   @param [in,out] x Vector that contains an initial guess of the solution and as output the final optimization variables if optimization has succeeded
   @param [in] low Vector of the same size as <tt>x</tt> that contains for each variable 
   the lower bound to impose. Note that Dirichlet boundary conditions are treated as equality conditions (i.e. lower and upper bounds) and that these ones can be imposed via an auxiliary optimization file (BCAsConstraint)
   @param [in] up Vector of the same size as <tt>x</tt> that contains for each variable 
   the upper bound to impose. Note that Dirichlet boundary conditions are treated as equality conditions (i.e. lower and upper bounds) and that these ones can be imposed via an auxiliary optimization function BCAsConstraint)
   @param [in] pivot Vector of the same size as <tt>x</tt> that contains on return for each variable an integer value that says if the corresponding constraint was reached (different from <tt>0</tt>) or not (<tt>= 0</tt>). Note that Dirichlet boundary conditions 
   are treated as equality constraints
   @param [in] max_it Maximum number of iterations for convergence
   @param [in] toler Tolerance for convergence (measured in relative weighted 2-Norm of projected gradient)
   @param [in] msg_lvl Output message level. Must be between <tt>0</tt> and <tt>10</tt>
   @return  Number of performed iterations

   \tparam <OPT_> Class that provides the objective function. This class is defined by the user.
 
   The \b OPT_ class:
 
   This class is defined by the user. It must have the member function :
  
   \b void Objective(Vect<real_t> &x, real_t &f, Vect<real_t> &g) 

   Here above, <tt>x</tt> is the optimization variable vector, <tt>f</tt> is the value of the objective to 
   calculate for the given <tt>x</tt> and <tt>g</tt> is the gradient vector for <tt>x</tt>.
  
   The function \b BCAsConstraint:
  
   This function is defined by the user:

   \b void BCAsConstraint(const Mesh &m, const Vect<real_t> &bc, Vect<real_t> &up, Vect<real_t> &low)

   This function imposes Dirichlet boundary conditions in an optimization problem 
   as optimization constraints. If such conditions are to be present, this function 
   has to be invoked by giving on input <tt>bc(i)</tt> as the value to impose for the 
   <tt>i</tt>-th optimization variable.
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class OPT_>
int lmqnbc(      OPT_&   theOpt,
                 size_t  dim,
                 real_t* x,
                 real_t& f,
                 real_t* g,
                 real_t* low,
                 real_t* up,
                 int*    pivot,
           const int&    msg_lvl,
           const int&    max_it,
                 size_t& max_fun,
                 real_t& eta,
                 real_t& stepmx,
                 real_t& accrcy,
                 real_t& xtol);

void monit(      size_t  dim,
                 real_t* x,
           const real_t& f,
                 real_t* g,
                 int     niter,
                 int     nftotl,
                 int     nfeval,
                 int*    pivot);

void ztime(size_t  dim,
           real_t* x,
           int*    pivot);

real_t stpmax(const real_t& stepmx,
              const real_t& pe,
                    size_t  dim,
                    real_t* x,
                    real_t* p,
                    int*    pivot,
                    real_t* low,
                    real_t* up);

void modz(size_t  dim,
          real_t* x,
          real_t* p,
          int*    pivot,
          real_t* low,
          real_t* up,
          real_t& flast,
          real_t& fnew);

int cnvtst(const real_t& alpha,
           const real_t& pnorm,
           const real_t& toleps,
           const real_t& xnorm,
           const real_t& difnew,
           const real_t& rtleps,
           const real_t& ftest,
           const real_t& gtg,
           const real_t& peps,
           const real_t& gtpnew,
                 real_t& fnew,
                 real_t& flast,
                 real_t* g,
                 int*    pivot,
                 size_t  dim,
           const real_t& accrcy);

int crash(size_t  dim,
	  real_t* x,
	  int*    pivot,
	  real_t* low,
	  real_t* up);

template<class OPT_>
int modlnp (OPT_ &theOpt, int modet, real_t *zsol, real_t *gv, real_t *r,
            real_t *v, real_t *diagb, real_t *emat, real_t *x, real_t *g,
            real_t *zk, size_t dim, const int &max_it, int &nfeval, int &nmodif,
            int &nlincg, int &upd1, real_t &yksk, real_t &gsk, real_t &yrsr,
            int &lreset, const int &bounds, int *pivot, const real_t &accrcy,
            real_t &gtp, real_t &gnorm, real_t &xnorm, const int &lhyr, real_t *w);

void ndia3 (size_t dim, real_t *e, real_t *v, real_t *gv, real_t *r, real_t vgv, int modet);

void negvec(size_t  dim,
            real_t* v);

void lsout (int loc, int test, real_t xmin, real_t fmin, real_t gmin,
            real_t xw, real_t fw, real_t gw, real_t u, real_t a, real_t b,
            real_t tol, real_t eps, real_t scxbd, real_t xlamda);

int chkucp (const int &lwtest, const int &max_fun, size_t dim, real_t &alpha,
            const real_t &eta, real_t &peps, real_t &rteps, real_t &rtol,
            real_t &rtolsq, const real_t &stepmx, real_t &tst, const real_t &xtol,
            real_t &xnorm, real_t *x, real_t &sm, real_t &tiny, const real_t &accrcy);

real_t step1(real_t fnew,
             real_t fm,
             real_t gtp,
             real_t smax);

template<class OPT_>
void gtims (OPT_ &theOpt, real_t *v, real_t *gv, size_t dim, real_t *x, real_t *g,
            real_t *w, int &first, real_t &delta, const real_t &accrcy, const real_t &xnorm);

void ssbfgs (size_t dim, real_t *sj, real_t *hjv, real_t *hjyj, const real_t &yjsj,
             const real_t &yjhyj, const real_t &vsj, const real_t &vhyj, real_t *hjp1v);

void mslv (real_t *g, real_t *y, size_t dim, real_t *sk, real_t *yk, real_t *diagb,
           real_t *sr, real_t *yr, real_t *hyr, real_t *hg, real_t *hyk, int &upd1,
           real_t &yksk, real_t &gsk, real_t &yrsr, const int &lreset, const int &first);

void msolve (real_t *g, real_t *y, size_t dim, real_t *w, int &upd1, real_t &yksk, real_t &gsk,
             real_t &yrsr, const int &lreset, int first, int lhyr);

int initp3 (real_t *diagb, real_t *emat, size_t dim, int &lreset, real_t &yksk, const real_t &yrsr,
            real_t *bsk, real_t *sk, real_t *yk, real_t *sr, real_t *yr, int &modet, int &upd1);

void initpc (real_t *diagb, real_t *emat, size_t dim, real_t *w, int &modet,
             int &upd1, real_t &yksk, real_t &gsk, real_t &yrsr, int &lreset);

template<class OPT_>
int linder (OPT_ &theOpt, size_t dim, const real_t &sm, real_t &reltol, real_t &abstol,
            real_t &tnytol, real_t &eta, real_t &xbnd, real_t *p, const real_t &gtp, real_t *x,
            real_t &f, real_t &alpha, real_t *g, int &nftotl, real_t *w);

int getptc(real_t& big,
	   real_t& rtsmll,
	   real_t& reltol,
           real_t& abstol,
	   real_t& tnytol,
	   real_t& fpresn,
	   real_t& eta,
           real_t& rmu,
	   real_t& xbnd,
	   real_t& u,
	   real_t& fu,
	   real_t& gu,
           real_t& xmin,
	   real_t& fmin,
	   real_t& gmin,
	   real_t& xw,
	   real_t& fw,
           real_t& gw,
	   real_t& a,
	   real_t& b,
	   real_t& oldf,
	   real_t& b1,
           real_t& scxbnd,
	   real_t& e,
	   real_t& step,
	   real_t& factor,
           int&    braktd,
	   real_t& gtest1,
	   real_t& gtest2,
	   real_t& tol,
           int&    ientry,
	   int&    itest);

template<class OPT_>
void objective(OPT_&   theOpt,
	       size_t  dim,
	       real_t* px,
	       real_t& f,
	       real_t* pg);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

template<class OPT_>
inline int OptimTN(OPT_&         theOpt,
                   Vect<real_t>& x,
                   Vect<real_t>& low,
                   Vect<real_t>& up,
                   Vect<int>&    pivot,
                   int           max_it,
                   real_t        toler,
                   int           msg_lvl)
/* ------------------------------------------------------------------------- */
/*                                                                           */
/* This template function solves the optimization problem                    */
/*                                                                           */
/*   minimize     f(x)                                                       */
/*      x                                                                    */
/*   subject to   low <= x <= up                                             */
/*                                                                           */
/* where x is a vector of n real_t variables. The method used is             */
/* a truncated-newton algorithm (see "Newton-type Minimization via           */
/* the Lanczos method" by S.G. Nash (SIAM J. Numer. Anal. 21 (1984),         */
/* pp. 770-778). This algorithm finds a local minimum of f(x).  it does      */
/* not assume that the function f is convex (and so cannot guarantee a       */
/* global solution), but does assume that the function is bounded below.     */
/* it can solve problems having any number of variables, but it is           */
/* especially useful when the number of variables (n) is large.              */
/*                                                                           */
/* Function Arguments   :                                                    */
/*                                                                           */
/* x       - (real) vector of length at least n; on input, an initial        */
/*           estimate of the solution; on output, the computed solution.     */
/* g       - (real) vector of length at least n; on output, the final        */
/*           value of the gradient                                           */
/* f       - (real) on input, a rough estimate of the value of the           */
/*           objective function at the solution; on output, the value        */
/*           of the objective function at the solution                       */
/* w       - (real) work vector of length at least 14*n                      */
/*                                                                           */
/* low, up - (real) vectors of length at least n containing                  */
/*           the lower and upper bounds on the variables.  if                */
/*           there are no bounds on a particular variable, set               */
/*           the bounds to -1.d38 and 1.d38, respectively.                   */
/* pivot  - (int) work vector of length at least n, used                     */
/*           to record which variables are at their bounds.                  */
/*                                                                           */
/* Return value (int)                                                        */
/*        ( 0 => normal return                                               */
/*        ( 2 => more than max_fun evaluations                               */
/*        ( 3 => line search failed to find lower point (may not be serious) */
/*        (-1 => error in input parameters                                   */
/*                                                                           */
/* This is an easy-to-use driver for the main optimization routine           */
/* lmqnbc.  more experienced users who wish to customize performance         */
/* of this algorithm should call lmqbc directly.                             */
/*                                                                           */
/* ------------------------------------------------------------------------- */
/* This template function sets up all the parameters for the Truncated-Newton*/
/* algorithm.  the parameters are:                                           */
/*                                                                           */
/* eta     - severity of the linesearch                                      */
/* max_fun - maximum allowable number of function evaluations                */
/* xtol    - desired accuracy for the solution x                             */
/* stepmx  - maximum allowable step in the linesearch                        */
/* accrcy  - accuracy of computed function values                            */
/* msg_lvl - controls quantity of printed output                             */
/*           0 = none, 1 = one line per major iteration.                     */
/* max_it  - maximum number of inner iterations per step                     */
/*---------------------------------------------------------------------------*/
{
    int ret;
    size_t max_fun;
    real_t xtol, accrcy, stepmx, eta;

//  Set parameters for the optimization routine

    size_t dim=x.size();
    if (max_it==-1)
      max_it = dim / 2;
    if (max_it>50)
      max_it = 50;
    if (max_it<=0)
      max_it = 1;
    max_fun = dim * 150;
    eta = 0.25;
    stepmx = 10.;
    accrcy = toler;
    xtol = sqrt(accrcy);
    real_t f=0;

//  Define pointers
    real_t *px = new real_t [dim];
    real_t *g = new real_t [dim];
    real_t *plow = new real_t [dim];
    real_t *pup = new real_t [dim];
    int *ppivot = new int [dim];
    for (size_t k=0; k<dim; k++) {
       px[k] = x[k];
       plow[k] = low[k];
       pup[k] = up[k];
       ppivot[k] = pivot[k];
    }

//  Minimize function
    ret = lmqnbc<OPT_>(theOpt, dim, px, f, g, plow, pup, ppivot, msg_lvl, max_it, max_fun,
                       eta, stepmx, accrcy, xtol);

//  Print results
    if (ret)
      cout << "\n\nError Code = " << ret << endl;
    cout << "\n\nOptimal Function Value = " << f << endl;
    for (size_t j=0; j<dim; j++) {
       x[j] = px[j];
       pivot[j] = ppivot[j];
    }
    delete [] px; delete [] g; delete [] plow; delete [] pup; delete [] ppivot;
    return ret;
}


#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class OPT_>
inline int lmqnbc(      OPT_ &theOpt,
                        size_t  dim,
                        real_t* x,
                        real_t& f,
                        real_t* g,
                        real_t* low,
                        real_t* up,
                        int*    pivot,
                  const int&    msg_lvl,
                  const int&    max_it,
                        size_t& max_fun,
                        real_t& eta,
                        real_t& stepmx,
                        real_t& accrcy,
                        real_t& xtol)
/*------------------------------------------------------------------------------*/
/*                                                                              */
/* This function is a bounds-constrained Truncated-Newton method.               */
/* the truncated-Newton method is preconditioned by a limited-memory            */
/* quasi-newton method (this preconditioning strategy is developed              */
/* in this routine) with a further diagonal scaling (see function ndia3).       */
/* for further details on the parameters, see function tnbc.                    */
/*                                                                              */
/*------------------------------------------------------------------------------*/
{
    size_t i;
    int numf=0, nwhy, ioldg, modet, niter=0, icycle, nlincg, idiagb;
    int nfeval, nmodif, ireset, nftotl, nm1, ipk, isk, iyk;
    int conv, newcon, lreset, upd1;
    real_t oldf, fnew, peps, rtol, yksk, tiny, yrsr, alpha, fkeep;
    real_t sm, flast, gnorm, ftest, pnorm, rteps, xnorm;
    real_t fm=0., pe, difold, difnew, reltol, gtpnew, toleps, gtg;
    real_t epsred, abstol, oldgtp, rtleps, rtolsq, tnytol, gsk, spe;

//  Check that initial x is feasible and that the bounds are consistent
    real_t *w = new real_t [14*dim];
    int ier = crash(dim, x, pivot, low, up);
    if (ier) {
      cout << "There is no feasible point; Terminating algorithm." << endl;
      delete [] w;
      return ier;
    }
    if (msg_lvl >= 1)
      cout << "\n\n   NIT  NF   CG           F                  GTG\n\n" << endl;

//  Initialize variables
    upd1 = 1;
    ireset = nfeval = nmodif = nlincg = 0;
    conv = 0;
    nm1 = dim - 1;

//  Within this routine the array w(loldg) is shared by w(lhyr)
    int lhyr = 9*dim, lwtest = 14*dim;

//  Check parameters and set constants
    nwhy = chkucp(lwtest, max_fun, dim, alpha, eta, peps, rteps, rtol, rtolsq,
                  stepmx, ftest, xtol, xnorm, x, sm, tiny, accrcy);
    if (nwhy < 0)
      goto L160;

    objective<OPT_>(theOpt,dim,x,fnew,g);
    nftotl = 1;
    oldf = fnew;
    gtg = Dot(dim, g, g);
    flast = fnew;

// Test the Lagrange multipliers to see if they are non-negative.
// because the constraints are only lower bounds, the components
// of the gradient corresponding to the active constraints are the
// lagrange multipliers. Afterwards, the projected gradient is formed.

   for (i=0; i<dim; i++)
      if (pivot[i] != 2)
         if (-pivot[i] * g[i] < 0.)
            pivot[i] = 0;
   ztime(dim, g, pivot);
   gtg = Dot(dim, g, g);
   if (msg_lvl >= 1)
      monit(dim, x, fnew, g, niter, nftotl, nfeval, pivot);

// Check if the initial point is a local minimum
   ftest = fabs(fnew) + 1.;
   if (gtg < OFELI_EPSMCH * 1e-4 * ftest * ftest)
      goto L130;

// Set initial values to other parameters
   icycle = nm1;
   toleps = rtol + rteps;
   rtleps = rtolsq + OFELI_EPSMCH;
   gnorm = sqrt(gtg);
   difnew = 0.;
   epsred = 0.05;
   fkeep = fnew;

// Set the diagonal of the approximate hessian to unity
   idiagb = 6*dim;
   for (i=0; i<dim; i++)
      w[idiagb++] = 1.;

/* ..................start of main iterative loop.......... */

// Compute the new search direction
   modet = msg_lvl - 3;
   modlnp<OPT_>(theOpt, modet, &w[12*dim], &w[0], &w[dim], &w[3*dim], &w[6*dim],
                &w[13*dim], x, g, &w[2*dim], dim, max_it, nfeval, nmodif, nlincg,
                upd1, yksk, gsk, yrsr, lreset, 1, pivot, accrcy, gtpnew, gnorm,
                xnorm, lhyr, w);
L20:
   Copy(dim, g, &w[9*dim]);
   pnorm = Nrm2(dim, &w[12*dim]);
   oldf = fnew;
   oldgtp = gtpnew;

// Prepare to compute the step length
   pe = pnorm + OFELI_EPSMCH;

// Compute the absolute and relative tolers for the linear search
   reltol = rteps * (xnorm + 1.) / pe;
   abstol = -OFELI_EPSMCH * ftest / (oldgtp - OFELI_EPSMCH);

// Compute the smallest allowable spacing between points in the linear search
   tnytol = OFELI_EPSMCH * (xnorm + 1.) / pe;
   spe = stpmax(stepmx, pe, dim, x, &w[12*dim], pivot, low, up);

// Set the initial step length
   alpha = step1 (fnew, fm, oldgtp, spe);

// Perform the linear search
   nwhy = linder<OPT_>(theOpt, dim, sm, reltol, abstol, tnytol, eta, spe,
                       &w[12*dim], oldgtp, x, fnew, alpha, g, numf, w);
   newcon = 0;
   if (fabs(alpha-spe) <= OFELI_EPSMCH * 10.) {
      newcon = 1;
      nwhy = 0;
      modz(dim, x, &w[12*dim], pivot, low, up, flast, fnew);
      flast = fnew;
   }

   if (msg_lvl >= 3)
      cout << "        Linesearch results:  alpha,pnorm : " << alpha << " "
           << pnorm << endl;
   ++niter;
   nftotl += numf;

// If required, print the details of this iteration
   if (msg_lvl >= 1)
      monit(dim, x, fnew, g, niter, nftotl, nfeval, pivot);
   if (nwhy < 0)
      goto L160;
   if (nwhy == 0 || nwhy == 2)
      goto L40;

// The linear search has failed to find a lower point
   nwhy = 3;
   goto L140;

L40:
   if (nwhy <= 1)
      goto L50;
   objective<OPT_>(theOpt,dim,x,fnew,g);
   ++nftotl;

//  Terminate if more than max_fun evaluations have been made
L50:
   nwhy = 2;
   if (nftotl > int(max_fun))
      goto L150;
   nwhy = 0;

// Set up parameters used in convergence and resetting tests

   difold = difnew;
   difnew = oldf - fnew;

// If this is the first iteration of a new cycle, compute the
// percentage reduction factor for the resetting test

   if (icycle == 1) {
      if (difnew > difold * 2.)
         epsred += epsred;
      if (difnew < difold * .5)
         epsred *= .5;
   }
   Copy(dim, g, &w[0]);
   ztime(dim, &w[0], pivot);
   gtg = Dot(dim, &w[0], &w[0]);
   gnorm = sqrt(gtg);
   ftest = fabs(fnew) + 1.;
   xnorm = Nrm2(dim, x);

// Test for convergence
   conv = cnvtst(alpha, pnorm, toleps, xnorm, difnew, rtleps, ftest,
                 gtg, peps, gtpnew, fnew, flast, g, pivot, dim, accrcy);
   if (conv)
      goto L130;
   ztime(dim, g, pivot);

// Compute the change in the iterates and the corresponding change in the gradients

   if (newcon)
      goto L90;
   isk = 4*dim; ipk = 12*dim;
   iyk = 5*dim; ioldg = 9*dim;
   for (i=0; i<dim; ++i) {
      w[iyk] = g[i] - w[ioldg];
      w[isk] = alpha * w[ipk];
      ++ipk; ++isk; ++iyk; ++ioldg;
   }

// Set up parameters used in updating the preconditioning strategy

   yksk = Dot(dim, &w[5*dim], &w[4*dim]);
   lreset = 0;
   if (icycle == nm1 || difnew < epsred * (fkeep - fnew))
      lreset = 1;
   if (lreset)
      goto L80;
   yrsr = Dot(dim, &w[8*dim], &w[7*dim]);
   if (yrsr <= 0.)
      lreset = 1;
L80:
   upd1 = 0;

// Compute the new search direction

L90:
   if (upd1 && msg_lvl >= 3)
      cout << "upd1 is true - Trivial Preconditioning" << endl;
   if (newcon && msg_lvl >= 3)
      cout << "newcon is true - Constraint added in linesearch" << endl;
   modet = msg_lvl - 3;
   modlnp<OPT_>(theOpt, modet, &w[12*dim], &w[0], &w[dim], &w[3*dim], &w[6*dim],
                &w[13*dim], x, g, &w[2*dim], dim, max_it, nfeval, nmodif,
                 nlincg, upd1, yksk, gsk, yrsr, lreset, 1, pivot, accrcy, gtpnew,
                 gnorm, xnorm, lhyr,w);
   if (newcon)
      goto L20;
    if (lreset)
      goto L110;

//  Compute the accumulated step and its corresponding gradient difference

    Xpy(dim, &w[4*dim], &w[7*dim]);
    Xpy(dim, &w[5*dim], &w[8*dim]);
    ++icycle;
    goto L20;

//  Reset

L110:
    ++ireset;

//  Initialize the sum of all the changes in x

    Copy(dim, &w[4*dim], &w[7*dim]);
    Copy(dim, &w[5*dim], &w[8*dim]);
    fkeep = fnew;
    icycle = 1;
    goto L20;

/* ...............end of main iteration....................... */

L130:
    f = fnew;
    delete [] w;
    return 0;
L140:
    oldf = fnew;

//  Local search could be installed here

L150:
    f = oldf;
    if (msg_lvl >= 1)
      monit(dim, x, f, g, niter, nftotl, nfeval, pivot);

L160:
    delete [] w;
    return nwhy;
}


inline void monit(size_t dim, real_t *x, const real_t &f, real_t *g, int niter,
                  int nftotl, int nfeval, int *pivot)
/*------------------------------------------------------------------------------*/
/*                    Print results of current iteration                        */
/*------------------------------------------------------------------------------*/
{
    if (x) { }
    real_t gtg = 0.;
    cout.setf(ios::scientific);
    for (size_t i=0; i<dim; i++)
       if (pivot[i] == 0)
         gtg += g[i] * g[i];
    cout << setw(5) << niter << setw(5) << nftotl << setw(5) << nfeval
         << "  " << setprecision(8) << setw(18) << f << "  "
         << setprecision(8) << setw(18) << gtg << endl;
}


inline void ztime(size_t  dim,
                  real_t* x,
                  int*    pivot)
//------------------------------------------------------------------------------
//      This function multiplies the vector x by the constraint matrix z
//------------------------------------------------------------------------------
{
   for (size_t i=0; i<dim; ++i)
      if (pivot[i])
         x[i] = 0.0;
}


inline real_t stpmax(const real_t &stepmx, const real_t &pe, size_t dim, real_t *x,
                     real_t *p, int *pivot, real_t *low, real_t *up)
/*------------------------------------------------------------------------------*/
/*                     Compute the maximum allowable step length                */
/*                                                                              */
/*            spe is the standard (unconstrained) max step                      */
/*------------------------------------------------------------------------------*/
{
   real_t t, spe = stepmx / pe;
   for (size_t i=0; i<dim; ++i) {
      if (pivot[i] == 0)
         if (p[i] != 0.) {
            if (p[i] <= 0.) {
               t = low[i] - x[i];
               if (t > spe * p[i])
                  spe = t / p[i];
            }
            else {
               t = up[i] - x[i];
               if (t < spe * p[i])
                  spe = t / p[i];
            }
         }
   }
   return (spe);
}


inline void modz(size_t dim, real_t *x, real_t *p, int *pivot, real_t *low, real_t *up,
                 real_t &flast, real_t &fnew)
/*------------------------------------------------------------------------------*/
/*      Update the constraint matrix if a new constraint is encountered         */
/*------------------------------------------------------------------------------*/
{
   real_t tol;
   for (size_t i=0; i<dim; i++) {
      if (!pivot[i]) {
         if (p[i] != 0.) {
            if (p[i] <= 0.) {
               tol = OFELI_EPSMCH * 10. * (fabs(low[i]) + 1.);
               if (x[i] - low[i] <= tol) {
                  flast = fnew;
                  pivot[i] = -1;
                  x[i] = low[i];
               }
            }
            else {
               tol = OFELI_EPSMCH * 10. * (fabs(up[i]) + 10.);
               if (up[i] - x[i] <= tol) {
                  flast = fnew;
                  pivot[i] = 1;
                  x[i] = up[i];
               }
            }
         }
      }
   }
}


inline int cnvtst (const real_t &alpha, const real_t &pnorm, const real_t &toleps,
                   const real_t &xnorm, const real_t &difnew, const real_t &rtleps,
                   const real_t &ftest, const real_t &gtg, const real_t &peps,
                   const real_t &gtpnew, real_t &fnew, real_t &flast, real_t *g,
                   int *pivot, size_t dim, const real_t &accrcy)
/*------------------------------------------------------------------------------*/
/*                                                                              */
/*                         Test for convergence                                 */
/*                                                                              */
/* For details, See Gill, Murray, and Wright (1981, p. 308) and                 */
/* Fletcher (1981, p. 116). The multiplier tests (here, testing                 */
/* the sign of the components of the gradient) may still need to                */
/* modified to incorporate tolers for zero.                                 */
/*                                                                              */
/*------------------------------------------------------------------------------*/
{
   real_t t, cmax=0.;
   int conv=0, imax=0;
   int ltest = flast - fnew <= gtpnew * -0.5;
   for (size_t i=0; i<=dim; i++) {
      if (pivot[i] != 0 && pivot[i] != 2) {
         t = -pivot[i] * g[i];
         if (t < 0.) {
            conv = 0;
            if (!ltest)
               if (cmax > t) {
                  cmax = t;
                  imax = int(i);
               }
         }
      }
   }
   if (imax!=0) {
      pivot[imax-1] = 0;
      flast = fnew;
      return conv;
   }
   conv = 0;
   if ((alpha * pnorm >= toleps * (xnorm + 1.) || fabs(difnew) >=
       rtleps * ftest || gtg >= peps * ftest * ftest) && gtg >=
       accrcy * 1e-4 * ftest * ftest)
      return 0;
   return 1;
}


inline int crash(size_t  dim,
                 real_t* x,
                 int*    pivot,
                 real_t* low,
                 real_t* up)
/*------------------------------------------------------------------------------*/
/*                                                                              */
/* This function initializes the constraint information, and ensures that the   */
/* initial point satisfies  low <= x <= up.                                     */
/* The constraints are checked for consistency.                                 */
/*                                                                              */
/*------------------------------------------------------------------------------*/
{
   int er = 0;
   for (size_t i=0; i<dim; i++) {
      if (x[i] < low[i])
         x[i] = low[i];
      if (x[i] > up[i])
         x[i] = up[i];
      pivot[i] = 0;
      if (x[i] == low[i])
         pivot[i] = -1;
      if (x[i] == up[i])
         pivot[i] = 1;
      if (up[i] == low[i])
         pivot[i] = 2;
      if (low[i] > up[i])
         er = -int(i) - 1;
   }
   return er;
}


template<class OPT_>
inline int modlnp (OPT_ &theOpt, int modet, real_t *zsol, real_t *gv, real_t *r,
                   real_t *v, real_t *diagb, real_t *emat, real_t *x, real_t *g,
                   real_t *zk, size_t dim, const int &max_it, int &nfeval, int &nmodif,
                   int &nlincg, int &upd1, real_t &yksk, real_t &gsk, real_t &yrsr,
                   int &lreset, const int &bounds, int *pivot, const real_t &accrcy,
                   real_t &gtp, real_t &gnorm, real_t &xnorm, const int &lhyr,
                   real_t *w)
/*------------------------------------------------------------------------------*/
/*                                                                              */
/* This function performs a preconditioned conjugate-gradient                   */
/* iteration in order to solve the newton equations for a search                */
/* direction for a Truncated-Newton algorithm. When the value of the            */
/* quadratic model is sufficiently reduced, the iteration is terminated.        */
/*                                                                              */
/* parameters                                                                   */
/*                                                                              */
/* modet       - integer which controls amount of output                        */
/* zsol        - computed search direction                                      */
/* g           - current gradient                                               */
/* gv,gz1,v    - scratch vectors                                                */
/* r           - residual                                                       */
/* diagb,emat  - diagonal preconditoning matrix                                 */
/* niter       - nonlinear iteration #                                          */
/* feval       - value of quadratic function                                    */
/*                                                                              */
/*------------------------------------------------------------------------------*/
{
   real_t beta=0, qnew=0, alpha=0, delta=0, rzold=0, rnorm=0, qtest=0, pr;
   nmodif = 0;

// Initialization

// General initialization

   if (modet > 0)
      cout << "\n\nEntering modlnop" << endl;
   if (max_it == 0)
      return 0;

   int first = 1;
   real_t rhsnrm = gnorm, tol = 1e-12, qold = 0.;

// Initialization for preconditioned conjugate-gradient algorithm
   initpc (diagb, emat, dim, w, modet, upd1, yksk, gsk, yrsr, lreset);
   for (size_t i=0; i<dim; i++) {
      r[i] = -g[i];
      v[i] = zsol[i] = 0.;
   }

// Main Iteration

   int k;
   for (k=1; k<=max_it; k++) {
      ++nlincg;
      if (modet > 1)
         cout << "\n\n### Iteration " << k << " ###" << endl;

//    CG iteration to solve system of equations

      if (bounds)
         ztime(dim, r, pivot);
      msolve (r, zk, dim, w, upd1, yksk, gsk, yrsr, lreset, first, lhyr);
      if (bounds)
         ztime(dim, zk, pivot);
      real_t rz = Dot(dim, r, zk);
      if (rz / rhsnrm < tol)
         goto L80;
      if (k==1)
         beta = 0.;
      else if (k > 1)
         beta = rz / rzold;
      for (size_t j=0; j<dim; j++)
         v[j] = zk[j] + beta * v[j];
      if (bounds)
         ztime(dim, v, pivot);
      gtims<OPT_>(theOpt, v, gv, dim, x, g, w, first, delta, accrcy, xnorm);
      if (bounds)
         ztime(dim, gv, pivot);
      ++nfeval;
      real_t vgv = Dot(dim, v, gv);
      if (vgv / rhsnrm < tol)
         goto L50;
      ndia3(dim, emat, v, gv, r, vgv, modet);

//    Compute linear step length
      alpha = rz / vgv;
      if (modet >= 1)
         cout << "alpha = " << alpha << endl;

//    Compute current solution and related vectors
      Axpy(dim,  alpha, v, zsol);
      Axpy(dim, -alpha, gv, r);

//    Test for convergence
      gtp = Dot(dim, zsol, g);
      pr =  Dot(dim, r, zsol);
      qnew = (gtp + pr) * 0.5;
      qtest = k * (1. - qold / qnew);
      if (qtest < 0.)
         goto L70;
      qold = qnew;
      if (qtest <= 0.5)
         goto L70;

//    Perform cautionary test
      if (gtp > 0.)
         goto L40;
      rzold = rz;
   }

// Terminate algorithm
   --k;
   goto L70;

//  Truncate algorithm in case of an emergency
L40:
    if (modet>=-1)
       cout << "g(t) positive at iteration " << k << " - truncating method" << endl;
    Axpy(dim, -alpha, v, zsol);
    gtp = Dot(dim, zsol, g);
    goto L90;
L50:
    if (modet>-2)
      cout << "         *** Hessian not positive-definite ***" << endl;

    if (k > 1)
      goto L70;

    msolve(g, zsol, dim, w, upd1, yksk, gsk, yrsr, lreset, first, lhyr);
    negvec(dim, zsol);
    if (bounds)
       ztime(dim, zsol, pivot);
    gtp = Dot(dim, zsol, g);
L70:
    if (modet >= -1)
       cout << "\n        modlan truncated after " << k << " iterations.  rnorm = "
            << rnorm << endl;
    goto L90;
L80:
    if (modet >= -1)
      cout << "Preconditioning not positive-definite." << endl;
    if (k > 1)
       goto L70;
    Copy(dim, g, zsol);
    negvec(dim, zsol);
    if (bounds)
       ztime(dim, zsol, pivot);
    gtp = Dot(dim, zsol, g);
    goto L70;

//  Store (or restore) diagonal preconditioning
L90:
    Copy(dim, emat, diagb);
    return 0;
}


inline void ndia3 (size_t  dim, 
                   real_t* e,
                   real_t* v,
                   real_t* gv,
                   real_t* r,
                   real_t  vgv,
                   int     modet)
/*------------------------------------------------------------------------------*/
/*                                                                              */
/* Update the preconditioing matrix based on a diagonal version of the BFGS     */
/* quasi-newton update.                                                         */
/*                                                                              */
/*------------------------------------------------------------------------------*/
{
   real_t vr = Dot(dim, v, r);
   for (size_t i=0; i<dim; i++) {
      e[i] = e[i] - r[i] * r[i] / vr + gv[i] * gv[i] / vgv;
      if (e[i] <= 1e-6) {
         if (modet > -2)
            cout << " *** emat negative :  " << e[i] << endl;
         e[i] = 1.;
      }
   }
}


/*==============================================================================*/
/*                      Service routines for optimization                       */
/*==============================================================================*/


inline void negvec (size_t  dim,
                    real_t* v)
{
   for (size_t i=0; i<dim; ++i)
      v[i] = -v[i];
}


inline void lsout (int    loc,
                   int    test,
                   real_t xmin,
                   real_t fmin,
                   real_t gmin,
                   real_t xw,
                   real_t fw,
                   real_t gw,
                   real_t u,
                   real_t a,
                   real_t b,
                   real_t tol,
                   real_t eps,
                   real_t scxbd,
                   real_t xlamda)
/*------------------------------------------------------------------------------*/
/*                            Error printouts for getptc                        */
/*------------------------------------------------------------------------------*/
{
   xlamda = 0;
   real_t yu = xmin + u, ya = a + xmin, yb = b + xmin, yw = xw + xmin, ybnd = scxbd + xmin;
   cout << "\n\n\nOutput from linear search" << endl;
   cout << "  Tol and eps : " << tol << "  " << eps << endl;
   cout << "  Current upper and lower bounds : " << ya << "  " << yb << endl;
   cout << "  Strict upper bound : " << ybnd << endl;
   cout << "  xw, fw, gw : " << yw << ", " << fw << ", " << gw << endl;
   cout << "  xmin, fmin, gmin : " << xmin << ", " << fmin << ", " << gmin << endl;
   cout << "  New estimate : " << yu << endl;
   cout << "  iloc, itest : " << loc << ", " << test << endl;
}


inline real_t step1(real_t fnew,
                    real_t fm,
                    real_t gtp,
                    real_t smax)
/*------------------------------------------------------------------------------*/
/* step1 returns the length of the initial step to be taken along the           */
/* vector p in the next linear search.                                          */
/*------------------------------------------------------------------------------*/
{
   static real_t d, alpha;
   d = fabs(fnew-fm);
   alpha = 1.;
   if (d * 2. <= -gtp && d >= OFELI_EPSMCH)
      alpha = d * -2. / gtp;
   if (alpha>=smax)
      alpha = smax;
   return alpha;
}


inline int chkucp(const int&    lwtest,
                  const int&    max_fun,
                        size_t  dim,
                        real_t& alpha,
                  const real_t& eta,
                        real_t& peps,
                        real_t& rteps,
                        real_t& rtol,
                        real_t& rtolsq,
                  const real_t& stepmx,
                        real_t& tst,
                  const real_t& xtol,
                        real_t& xnorm,
                        real_t* x,
                        real_t& sm,
                        real_t& tiny,
                  const real_t &accrcy)
/*------------------------------------------------------------------------------*/
/* Checks parameters and sets constants which are common to both derivative and */
/* non-derivative algorithms                                                    */
/*------------------------------------------------------------------------------*/
{
   int lw = 14*int(dim);
   int nwhy;

   sm = OFELI_EPSMCH*OFELI_EPSMCH;
   tiny = sm;
   nwhy = -1;
   rteps = sqrt(OFELI_EPSMCH);
   rtol = xtol;
   if (fabs(rtol)<accrcy)
      rtol = rteps*10.;

// Check for errors in the input parameters
   if (lw<lwtest || dim<1 || rtol<0. || eta>=1. || eta<0. || stepmx<rtol || max_fun<1)
      return nwhy;
   nwhy = 0;

// Set constants for later

   rtolsq = rtol*rtol;
   peps = pow(accrcy,0.6666);
   xnorm = Nrm2(dim, x);
   alpha = tst = 0.;
   return nwhy;
}


template<class OPT_>
inline void gtims(      OPT_&   theOpt,
                        real_t* v,
                        real_t* gv,
                        size_t  dim,
                        real_t* x,
                        real_t* g,
                        real_t* w,
                        int&    first,
                        real_t& delta,
                  const real_t& accrcy,
                  const real_t &xnorm)
/*------------------------------------------------------------------------------*/
/* This function computes the product of the matrix g times the vector          */
/* v and stores the result in the vector gv (finite-difference version)         */
/*------------------------------------------------------------------------------*/
{
   static real_t dinv, f;
   int ihg;
   size_t i;

   if (first) {
     delta = sqrt(accrcy) * (xnorm + 1.);
     first = 0;
   }
   dinv = 1. / delta;
   ihg = 10*dim;
   for (i=0; i<dim; i++)
      w[ihg++] = x[i] + delta * v[i];
   objective<OPT_>(theOpt,dim,&w[10*dim],f,gv);
   for (i=0; i<dim; i++)
      gv[i] = (gv[i] - g[i]) * dinv;
}


inline void ssbfgs (size_t dim, real_t *sj, real_t *hjv, real_t *hjyj, const real_t &yjsj,
                    const real_t &yjhyj, const real_t &vsj, const real_t &vhyj,
                    real_t *hjp1v)
/*------------------------------------------------------------------------------*/
/*                              Self-Scaled BFGS                                */
/*------------------------------------------------------------------------------*/
{
   real_t gamma = 1.;
   real_t delta = (gamma * yjhyj / yjsj + 1.) * vsj / yjsj - gamma * vhyj / yjsj;
   real_t beta = -gamma * vsj / yjsj;
   for (size_t i=0; i<dim; ++i)
      hjp1v[i] = gamma * hjv[i] + delta * sj[i] + beta * hjyj[i];
}


inline void mslv (real_t *g, real_t *y, size_t dim, real_t *sk, real_t *yk,
                  real_t *diagb, real_t *sr, real_t *yr, real_t *hyr, real_t *hg,
                  real_t *hyk, int &upd1, real_t &yksk, real_t &gsk, real_t &yrsr,
                  const int &lreset, const int &first)
/*------------------------------------------------------------------------------*/
/* This function acts as a preconditioning step for the linear                  */
/* conjugate-gradient routine. It is also the method of computing the search    */
/* direction from the gradient for the non-linear conjugate-gradient code.      */
/* it represents a two-step self-scaled bfgs formula.                           */
/*------------------------------------------------------------------------------*/
{
    real_t ghyk=0, ghyr=0, yksr=0, ykhyk=0, ykhyr=0, yrhyr=0, rdiagb=0, gsr=0;
    size_t i;

    if (!upd1) {
       gsk = Dot(dim, g, sk);
       if (!lreset) {

//     Compute hg and hy where h is the inverse of the diagonals

          for (i=0; i<dim; i++) {
             rdiagb = 1. / diagb[i];
             hg[i] = g[i] * rdiagb;
             if (first) {
                hyk[i] = yk[i] * rdiagb;
                hyr[i] = yr[i] * rdiagb;
             }
          }
          if (first) {
             yksr = Dot(dim, yk, sr);
             ykhyr = Dot(dim, yk, hyr);
          }
          gsr = Dot(dim, g, sr);
          ghyr = Dot(dim, g, hyr);
          if (first)
             yrhyr = Dot(dim, yr, hyr);
          ssbfgs(dim, sr, hg, hyr, yrsr, yrhyr, gsr, ghyr, hg);
          if (first)
             ssbfgs(dim, sr, hyk, hyr, yrsr, yrhyr, yksr, ykhyr, hyk);
          ykhyk = Dot(dim, hyk, yk);
          ghyk = Dot(dim, hyk, g);
          ssbfgs(dim, sk, hg, hyk, yksk, ykhyk, gsk, ghyk, y);
       }

//     Compute gh and hy where h is the inverse of the diagonals

       for (i=0; i<dim; i++) {
          rdiagb = 1. / diagb[i];
          hg[i] = g[i] * rdiagb;
          if (first)
             hyk[i] = yk[i] * rdiagb;
       }
       if (first)
          ykhyk = Dot(dim, yk, hyk);
       ghyk = Dot(dim, g, hyk);
       ssbfgs(dim, sk, hg, hyk, yksk, ykhyk, gsk, ghyk, y);
    }
    for (i=0; i<dim; i++)
       y[i] = g[i] / diagb[i];
}


inline void msolve(      real_t* g,
                         real_t* y,
                         size_t  dim,
                         real_t* w,
		         int&    upd1,
                         real_t& yksk,
		         real_t& gsk,
		         real_t& yrsr,
		   const int&    lreset,
                         int     first,
	                 int     lhyr)
/*------------------------------------------------------------------------------*/
/*   This function sets upt the arrays for mslv                                 */
/*------------------------------------------------------------------------------*/
{
   mslv(g, y, dim, &w[4*dim], &w[5*dim], &w[6*dim], &w[7*dim], &w[8*dim],
        &w[lhyr], &w[10*dim], &w[11*dim], upd1, yksk, gsk, yrsr, lreset, first);
}


/*==============================================================================*/
/*                   Functions to initialize preconditioner                     */
/*==============================================================================*/

inline int initp3(      real_t* diagb,
		        real_t* emat,
		        size_t  dim,
		        int&    lreset,
		        real_t& yksk,
                  const real_t& yrsr,
		        real_t* bsk,
		        real_t* sk,
		        real_t* yk,
                        real_t* sr,
		        real_t* yr,
		        int&    modet,
		        int&    upd1)
{
    size_t i;
    real_t srds, yrsk, d1, dn, td, sds;

    if (upd1)
       goto L90;
    if (lreset)
       goto L60;
    for (i=0; i<dim; i++)
       bsk[i] = diagb[i] * sr[i];
    sds = Dot(dim, sr, bsk);
    srds = Dot(dim, sk, bsk);
    yrsk = Dot(dim, yr, sk);
    for (i=0; i<dim; i++) {
       td = diagb[i];
       bsk[i] = td * sk[i] - bsk[i] * srds / sds + yr[i] * yrsk / yrsr;
       emat[i] = td - td * td * sr[i] * sr[i] / sds + yr[i] * yr[i] / yrsr;
    }
    sds = Dot(dim, sk, bsk);
    for (i=0; i<dim; i++)
       emat[i] = emat[i] - bsk[i] * bsk[i] / sds + yk[i] * yk[i] / yksk;
    goto L110;
L60:
    for (i=0; i<dim; i++)
       bsk[i] = diagb[i] * sk[i];
    sds = Dot(dim, &sk[1], &bsk[1]);
    for (i=0; i<dim; i++) {
       td = diagb[i];
       emat[i] = td - td * td * sk[i] * sk[i] / sds + yk[i] * yk[i] / yksk;
    }
    goto L110;
L90:
    Copy(dim, diagb, emat);
L110:
    if (modet < 1)
      return 0;
    d1 = dn = emat[0];
    for (i=0; i<dim; i++) {
       if (emat[i] < d1)
          d1 = emat[i];
       if (emat[i] > dn)
          dn = emat[i];
    }
    cout << "\n\n        dmin = " << d1 << ",  dmax = " << dn
         << ",  cond = " << dn/d1 << endl;
    return 0;
}


inline void initpc(real_t *diagb, real_t *emat, size_t dim, real_t *w, int &modet,
                   int &upd1, real_t &yksk, real_t &gsk, real_t &yrsr, int &lreset)
{
   gsk = 0;
   initp3(diagb, emat, dim, lreset, yksk, yrsr, &w[11*dim], &w[4*dim], &w[5*dim],
          &w[7*dim], &w[8*dim], modet, upd1);
}



/*==============================================================================*/
/*               Line Search algorithms of Gill and Murray                      */
/*==============================================================================*/


template<class OPT_>
inline int linder(OPT_ &theOpt, size_t dim, const real_t &sm, real_t &reltol,
                  real_t &abstol, real_t &tnytol, real_t &eta, real_t &xbnd,
                  real_t *p, const real_t &gtp, real_t *x, real_t &f, real_t &alpha,
                  real_t *g, int &nftotl, real_t *w)
{
    static int numf, l, itcnt, nprnt, lg, lx, ientry, lsprnt, tst, flag;
    static real_t oldf, fmin, gmin, step, xmin, a, b, e, u, b1, gtest1, gtest2;
    static real_t fu, gu, fw, gw, factor, scxbnd, xw, fpresn, rtsmll;
    static real_t big, tol, rmu, ualpha;
    static int braktd;

    lx = 0; lg = lx + dim;
    lsprnt = 0; nprnt = 10000;
    rtsmll = sqrt(sm);
    big = 1./sm;

//  Set the estimated relative precision in f(x)
    fpresn = OFELI_EPSMCH * 10.;
    numf = 0;
    u = alpha;
    fu = f;
    fmin = f;
    gu = gtp;
    rmu = 1e-4;

//  First entry sets up the initial interval of uncertainty
    ientry = 1;

//  Test for too many iterations
    for (itcnt=1; itcnt<=20; itcnt++) {
      flag = 0;
      getptc(big, rtsmll, reltol, abstol, tnytol, fpresn, eta, rmu,
             xbnd, u, fu, gu, xmin, fmin, gmin, xw, fw, gw, a, b,
             oldf, b1, scxbnd, e, step, factor, braktd, gtest1, gtest2,
             tol, ientry, tst);
      if (lsprnt >= nprnt)
        lsout (ientry, tst, xmin, fmin, gmin, xw, fw, gw, u, a, b, tol, reltol, scxbnd, xbnd);

//    If test=1, the algorithm requires the function value to be calculated

      if (tst == 1) {
        ualpha = xmin + u;
        l = lx;
        for (size_t i=0; i<dim; i++)
           w[l++] = x[i] + ualpha * p[i];
        objective<OPT_>(theOpt,dim,&w[lx],fu,&w[lg]);
        ++numf;
        gu = Dot(dim, &w[lg], p);

//      The gradient vector corresponding to the best point is overwritten if fu
//      is less than fmin and fu is sufficiently lower than f at the origin.

        if (fu <= fmin && fu <= oldf - ualpha * gtest1)
          Copy(dim, &w[lg], g);
        flag = 1;
      }

//    If test=2 or 3 a lower point could not be found

      else {
        nftotl = numf;
        flag = 1;
        if (tst==0) {

//        A successful search has been made
          flag = 0;
          f = fmin;
          alpha = xmin;
          for (size_t j=0; j<dim; ++j)
             x[j] += alpha*p[j];
        }
        return flag;
      }
    }
    return flag;
}


inline int getptc(real_t& big,
                  real_t& rtsmll,
                  real_t& reltol,
                  real_t& abstol,
                  real_t& tnytol,
                  real_t& fpresn,
                  real_t& eta,
                  real_t& rmu,
                  real_t& xbnd,
                  real_t& u,
                  real_t& fu,
                  real_t& gu,
                  real_t& xmin,
                  real_t& fmin,
                  real_t& gmin,
                  real_t& xw,
                  real_t& fw,
                  real_t& gw,
                  real_t& a,
                  real_t& b,
                  real_t& oldf,
                  real_t& b1,
                  real_t& scxbnd,
                  real_t& e,
                  real_t& step,
                  real_t& factor,
                  int&    braktd,
                  real_t& gtest1,
                  real_t& gtest2,
                  real_t& tol,
                  int& ientry,
                  int& itest)
/*------------------------------------------------------------------------------*/
/*                                                                              */
/* An algorithm for finding a steplength, called repeatedly by functions which  */
/* require a step length to be computed using cubic interpolation. The          */
/* parameters contain information about the interval in which a lower point is  */
/* to be found and from this getptc computes a point at which the function can  */
/* be evaluated by the calling program. The value of the integer parameters     */
/* ientry determines the path taken through the code.                           */
/*                                                                              */
/*------------------------------------------------------------------------------*/
{
   real_t d;
   static real_t abgw, absr, p, q, r, s, sc, denom, a1, d1, d2;
   static int convrg;
   static real_t xmidpt, twotol, sumsq, abgmin, chordm, chordu;

// Branch to appropriate section of code depending on the value of ientry

   switch (ientry) {
      case 1:  goto L10;
      case 2:  goto L20;
   }

// ientry=1  : Check input parameters
L10:
   itest = 2;
   if (u<=0. || xbnd<=tnytol || gu>0.)
      return 0;
   itest = 1;
   if (xbnd<abstol)
      abstol = xbnd;
   tol = abstol;
   twotol = 2 * tol;

// a and b define the interval of uncertainty, x and xw are points
// with lowest and second lowest function values so far obtained.
// initialize a,smin,xw at origin and corresponding values of
// function and projection of the gradient along direction of search
// at values for latest estimate at minimum.

   a = xw = xmin = 0.;
   oldf = fmin = fu;
   fw = fu;
   gw = gu;
   gmin = gu;
   step = u;
   factor = 5.;

// The minimum has not yet been bracketed
   braktd = 0;

// Set up xbnd as a bound on the step to be taken. (xbnd is not computed
// explicitly but scxbnd is its scaled value.)  set the upper bound
// on the interval of uncertainty initially to xbnd + tol(xbnd).

   scxbnd = xbnd;
   b = scxbnd + reltol*fabs(scxbnd) + abstol;
   e = b + b;
   b1 = b;

// Compute the constants required for the two convergence criteria
   gtest1 = -rmu * gu;
   gtest2 = -eta * gu;

// Set ientry to indicate that this is the first iteration

   ientry = 2;
   goto L210;

// ientry = 2

// update a,b,xw, and xmin
L20:
   if (fu>fmin)
      goto L60;

// if function value not increased, new point becomes next
// origin and other points are scaled accordingly

   chordu = oldf - (xmin + u) * gtest1;
   if (fu <= chordu)
      goto L30;

/* The new function value does not satisfy the sufficient decrease                */
/* criterion. prepare to move the upper bound to this point and                   */
/* force the interpolation scheme to either bisect the interval of                */
/* uncertainty or take the linear interpolation step which estimates              */
/* the root of f(alpha)=chord(alpha).                                             */

   chordm = oldf - xmin * gtest1;
   gu = -gmin;
   denom = chordm - fmin;
   if (fabs(denom) >= OFELI_EPSMCH)
      goto L25;

   denom = 1e-15;
   if (chordm - fmin < 0.)
      denom = -denom;

L25:
   if (xmin != 0.)
      gu = gmin * (chordu - fu) / denom;
   fu = u * 0.5 * (gmin + gu) + fmin;
   if (fu < fmin)
      fu = fmin;
   goto L60;

L30:
   fw = fmin; 
   fmin = fu;
   gw = gmin;
   gmin = gu;
   xmin += u;
   a -= u;
   b -= u;
   xw = -u;
   scxbnd -= u;
   if (gu <= 0.)
      goto L40;
   b = 0.;
   braktd = 1;
   goto L50;

L40:
   a = 0.;

L50:
   tol = fabs(xmin)*reltol + abstol;
   goto L90;

// If function value increased, origin remains unchanged
//  but new point may now qualify as w.
L60:
   if (u < 0.)
      goto L70;
   b = u;
   braktd = 1;
   goto L80;
L70:
   a = u;
L80:
   xw = u;
   fw = fu;
   gw = gu;
L90:
   twotol = tol + tol;
   xmidpt = (a+b) * 0.5;

//  Check termination criteria
   convrg = ((fabs(xmidpt)<=twotol-(b-a)*0.5 || fabs(gmin)<=gtest2)) && fmin<oldf
             && ((fabs(xmin-xbnd) > tol || !braktd));
   if (!convrg)
      goto L100;
   itest = 0;
   if (xmin != 0.)
      return 0;

// if the function has not been reduced, check to see that the relative
// change in f(x) is consistent with the estimate of the delta-
// unimodality constant, tol.  if the change in f(x) is larger than
// expected, reduce the value of tol
   itest = 3;
   if (fabs(oldf-fw) <= fpresn * (fabs(oldf) + 1.))
      return 0;
   tol *= 0.1;
   if (tol < tnytol)
      return 0;
   reltol *= 0.1;
   abstol *= 0.1;
   twotol *= 0.1;

// Continue with the computation of a trial step length
L100:
   r = q = s = 0.;
   if (fabs(e) <= tol)
      goto L150;

// Fit cubic through xmin and xw
   r = (fmin-fw)*3.0/xw + gmin + gw;
   absr = fabs(r);
   q = absr;
   if (gw == 0. || gmin == 0.)
      goto L140;

// Compute the square root of (r*r - gmin*gw) in a way which avoids underflow and overflow
   abgw = fabs(gw);
   abgmin = fabs(gmin);
   s = sqrt(abgmin*abgw);
   if (gw/abgw*gmin > 0.)
      goto L130;

// Compute the square root of r*r + s*s
   sumsq = 1.0;
   p = 0.;
   if (absr >= s)
      goto L110;

//  There is a possibility of overflow
   if (s > rtsmll)
      p = s*rtsmll;
   if (absr>=p) {
      d = absr/s;
      sumsq = d*d + 1.;
   }
   sc = s;
   goto L120;

// There is a possibility of underflow
L110:
   if (absr>rtsmll)
      p = absr*rtsmll;
   if (s>=p) {
      d = s / absr;
      sumsq = d*d + 1.;
   }
   sc = absr;

L120:
   sumsq = sqrt(sumsq);
   q = big;
   if (sc<big/sumsq)
      q = sc * sumsq;
   goto L140;

// Compute the square root of r*r - s*s
L130:
   q = sqrt(fabs((r+s)*(r-s)));
   if (r>=s || r<=-s)
      goto L140;
   r = q = 0.;
   goto L150;

// Compute the minimum of fitted cubic
L140:
   if (xw < 0.)
      q = -q;
   s = xw * (gmin - r - q);
   q = gw - gmin + q + q;
   if (q > 0.)
      s = -s;
   if (q <= 0.)
      q = -q;
   r = e;
   if (b1 != step || braktd)
      e = step;

// Construct an artificial bound on the estimated steplength
L150:
   a1 = a; b1 = b;
   step = xmidpt;
   if (braktd)
      goto L160;
   step = -factor*xw;
   if (step > scxbnd)
      step = scxbnd;
   if (step != scxbnd)
      factor *= 5.;
   goto L170;

//  If the minimum is bracketed by 0 and xw the step must lie within (a,b)

L160:
   if ((a!=0. || xw>=0.) && (b!=0. || xw<=0.))
      goto L180;

// If the minimum is not bracketed by 0 and xw the step must lie within (a1,b1)
   d1 = xw;
   d2 = a;
   if (a == 0.)
      d2 = b;
// This line might be: if (a .eq. 0.0) d2 = e
   u = -d1/d2;
   step = d2*5*(1./u + 0.1)/11.;
   if (u < 1.)
      step = d2 * 0.5 * sqrt(u);

L170:
   if (step <= 0.)
      a1 = step;
   else
      b1 = step;

/* Reject the step obtained by interpolation if it lies outside the    */
/* required interval or it is greater than half the step obtained      */
/* during the last-but-one iteration.                                  */

L180:
   if ((fabs(s) <= fabs(0.5*q*r)) || s<=q*a1 || s>=q*b1)
      goto L200;

// A cubic interpolation step
   step = s/q;

// The function must not be evalutated too close to a or b
   if (step-a>=twotol && b-step>=twotol)
      goto L210;
   if (xmidpt > 0.)
      goto L190;
   step = -tol;
   goto L210;

L190:
   step = tol;
   goto L210;

L200:
   e = b - a;

// If the step is too large, replace by the scaled bound (so as to
// compute the new point on the boundary)

L210:
   if (step >= scxbnd) {
      step = scxbnd;
      scxbnd -= (reltol*fabs(xbnd) + abstol) / (reltol + 1.);
   }
   u = step;
   if (fabs(step)<tol && step<0.)
      u = -tol;
   if (fabs(step)<tol && step>=0.)
      u = tol;
   itest = 1;
   return 0;
}


template<class OPT_>
inline void objective(OPT_&   theOpt,
                      size_t  dim,
                      real_t* px,
                      real_t &f,
                      real_t* pg)
{
   Vect<real_t> x(dim,px), g(dim,pg);
   theOpt.Objective(x,f,g);
   for (size_t i=0; i<dim; i++) {
      px[i] = x[i];
      pg[i] = g[i];
   }
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
