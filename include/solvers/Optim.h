/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani

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

                        Prototypes for optimization functions

  ==============================================================================*/

#ifndef __OPTIM_H
#define __OPTIM_H

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

namespace OFELI {

class OptSolver;

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/*==============================================================================*/
/*                                    OptimTN                                   */
/*==============================================================================*/
int OptimTN(OptSolver&          opt,
            Vect<real_t>&       x,
            const Vect<real_t>& lb,
            const Vect<real_t>& ub,
            int&                nb_obj_eval,
            int&                nb_grad_eval,
            int&                max_it,
            real_t              toler);

int lmqnbc(OptSolver&          opt,
           Vect<real_t>&       x,
           real_t&             f,
           Vect<real_t>&       g,
           const Vect<real_t>& lb,
           const Vect<real_t>& ub,
           vector<int>&        pivot,
           int&                max_it,
           size_t              max_fun,
           real_t&             eta,
           real_t&             stepmx,
           real_t&             accrcy,
           real_t&             xtol,
           int&                nb_obj_eval,
           int&                nb_grad_eval);

void monit(real_t              f,
           const Vect<real_t>& g,
           const vector<int>&  pivot,
           int                 niter,
           int                 nftotl,
           int                 nb_obj_eval);

void ztime(vector<real_t>&    x,
           const vector<int>& pivot);

/// Compute the maximum allowable step length
/// spe is the standard (unconstrained) max step
real_t stpmax(real_t                stepmx,
              real_t                pe,
              const Vect<real_t>&   x,
              const vector<real_t>& p,
              const vector<int>&    pivot,
              const Vect<real_t>&   lb,
              const Vect<real_t>&   ub);

/// Update the constraint matrix if a new constraint is encountered
void modz(Vect<real_t>&         x,
          const vector<real_t>& p,
          vector<int>&          pivot,
          const Vect<real_t>&   lb,
          const Vect<real_t>&   ub,
          real_t&               flast,
          real_t&               fnew);

int cnvtst(real_t              alpha,
           real_t              pnorm,
           real_t              toleps,
           real_t              xnorm,
           real_t              difnew,
           real_t              rtleps,
           real_t              ftest,
           real_t              gtg,
           real_t              peps,
           real_t              gtpnew,
           real_t&             fnew,
           real_t&             flast,
           const Vect<real_t>& g,
           vector<int>&        pivot,
           real_t              accrcy);

/* This function performs a preconditioned conjugate-gradient
 * iteration in order to solve the newton equations for a search
 * direction for a Truncated-Newton algorithm. When the value of the
 * quadratic model is sufficiently reduced, the iteration is terminated.
 *
 * parameters:
 * modet       - integer which controls amount of output
 * zsol        - computed search direction
 * gv,gz1,v    - scratch vectors
 * r           - residual
 * diagb,emat  - diagonal preconditoning matrix
 * niter       - nonlinear iteration
 * feval       - value of quadratic function
 */
int modlnp(OptSolver&    opt,
           Vect<real_t>& zsol,
           Vect<real_t>& gv,
           Vect<real_t>& r,
           Vect<real_t>& v,
           Vect<real_t>& diagb,
           Vect<real_t>& emat,
           Vect<real_t>& x,
           Vect<real_t>& g,
           Vect<real_t>& zk,
           int           max_it,
           int&          nb_obj_eval,
           int&          nb_grad_eval,
           int&          nmodif,
           int&          nlincg,
           int&          upd1,
           real_t&       yksk,
           real_t&       gsk,
           real_t&       yrsr,
           int&          lreset,
           int           bounds,
           vector<int>&  pivot,
           real_t        accrcy,
           real_t&       gtp,
           real_t&       gnorm,
           real_t&       xnorm,
           int           lhyr,
           Vect<real_t>* w);

/// Update the preconditioing matrix based on a diagonal version of the BFGS
/// quasi-newton update.
void ndia3(vector<real_t>&       e,
           const vector<real_t>& v,
           const vector<real_t>& gv,
           const vector<real_t>& r,
           real_t                vgv,
           int                   modet);

/*  This function initializes the constraint information, and ensures that the
 *  initial point satisfies  low <= x <= up.
 *  The constraints are checked for consistency.
 */
int crash(Vect<real_t>&       x,
          vector<int>&        pivot,
          const Vect<real_t>& lb,
          const Vect<real_t>& ub);

void lsout(int    loc,
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
           real_t xlamda);

/// step1 returns the length of the initial step to be taken along the
/// vector p in the next linear search.
real_t step1(real_t fnew,
             real_t fm,
             real_t gtp,
             real_t smax);

/// Check parameters and set constants which are common to both derivative and
/// non-derivative algorithms
int chkucp(int                 lwtest,
           int                 max_fun,
           real_t&             alpha,
           real_t              eta,
           real_t&             peps,
           real_t&             rteps,
           real_t&             rtol,
           real_t&             rtolsq,
           real_t              stepmx,
           real_t&             tst,
           real_t              xtol,
           real_t&             xnorm,
           const Vect<real_t>& x,
           real_t&             sm,
           real_t&             tiny,
           real_t              accrcy);

/// Compute the product of the matrix g times the vector v and store the result
/// in the vector gv (finite-difference version)
void gtims(OptSolver&          opt,
           const Vect<real_t>& v,
           Vect<real_t>&       gv,
           const Vect<real_t>& x,
           const Vect<real_t>& g,
           Vect<real_t>*       w,
           int&                first,
           real_t&             delta,
           real_t              accrcy,
           real_t              xnorm);

void ssbfgs(const vector<real_t>& sj,
            const vector<real_t>& hjv,
            const vector<real_t>& hjyj,
            real_t                yjsj,
            real_t                yjhyj,
            real_t                vsj,
            real_t                vhyj,
            vector<real_t>&       hjp1v);

/*  This function acts as a preconditioning step for the linear
 *  conjugate-gradient routine. It is also the method of computing the search
 *  direction from the gradient for the non-linear conjugate-gradient code.
 *  it represents a two-step self-scaled bfgs formula.
 */
void mslv(const vector<real_t>& g,
          vector<real_t>&       y,
          vector<real_t>&       sk,
          vector<real_t>&       yk,
          const vector<real_t>& diagb,
          vector<real_t>&       sr,
          vector<real_t>&       yr,
          vector<real_t>&       hyr,
          vector<real_t>&       hg,
          vector<real_t>&       hyk,
          int&                  upd1,
          real_t&               yksk,
          real_t&               gsk,
          real_t&               yrsr,
          int                   lreset,
          int                   first);

void msolve(const vector<real_t>& g,
            vector<real_t>&       y,
            Vect<real_t>*         w,
            int&                  upd1,
            real_t&               yksk,
            real_t&               gsk,
            real_t&               yrsr,
            int                   lreset,
            int                   first,
            int                   lhyr);

int initp3(const vector<real_t>& diagb,
           vector<real_t>&       emat,
           int&                  lreset,
           real_t&               yksk,
           real_t                yrsr,
           vector<real_t>&       bsk,
           const vector<real_t>& sk,
           const vector<real_t>& yk,
           const vector<real_t>& sr,
           const vector<real_t>& yr,
           int&                  modet,
           int&                  upd1);

void initpc(const vector<real_t>& diagb,
            vector<real_t>&       emat,
            Vect<real_t>*         w,
            int&                  modet,
            int&                  upd1,
            real_t&               yksk,
            real_t&               gsk,
            real_t&               yrsr,
            int&                  lreset);

int linder(OptSolver&            opt,
           real_t                sm,
           real_t&               reltol,
           real_t&               abstol,
           real_t&               tnytol,
           real_t&               eta,
           real_t&               xbnd,
           const vector<real_t>& p,
           real_t                gtp,
           Vect<real_t>&         x,
           real_t&               f,
           real_t&               alpha,
           Vect<real_t>&         g,
           int&                  nftotl,
           Vect<real_t>*         w);

/* An algorithm for finding a steplength, called repeatedly by functions which
 * require a step length to be computed using cubic interpolation. The
 * parameters contain information about the interval in which a lower point is
 * to be found and from this getptc computes a point at which the function can
 * be evaluated by the calling program. The value of the integer parameters
 * ientry determines the path taken through the code.
 */
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


/*==============================================================================*/
/*                                    OptimSA                                   */
/*==============================================================================*/

int OptimSA(OptSolver&          opt,
            Vect<real_t>&       x,
            real_t&             rt,
            real_t&             toler,
            int&                ns,
            int&                nt,
            int&                neps,
            int&                max_eval,
            const Vect<real_t>& lb,
            const Vect<real_t>& ub,
            const Vect<real_t>& c,
            real_t&             t,
            Vect<real_t>&       vm,
            real_t&             fopt,
            int&                nacc,
            int&                nb_eval,
            int&                nobds);


/*==============================================================================*/
/*                                    OptimNM                                   */
/*==============================================================================*/

int OptimNM(OptSolver&    opt,
            Vect<real_t>& x,
            real_t&       ynewlo,
            real_t        reqmin,
            Vect<real_t>& step,
            int           conv,
            int           max_eval,
            int&          nb_eval,
            int&          nb_restart);

/*==============================================================================*/
/*                                    OptimPG                                   */
/*==============================================================================*/

int OptimPG(OptSolver&          opt,
            Vect<real_t>&       x,
            const Vect<real_t>& lb,
            const Vect<real_t>& ub,
            int&                nb_eval,
            int                 max_it,
            real_t              toler);

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */

#endif
