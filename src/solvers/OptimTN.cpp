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

                       Implementation of class 'OptimTN'

  ==============================================================================*/

#include "solvers/OptSolver.h"
#include "solvers/Optim.h"
#include "util/util.h"

#include <iostream>
using std::cout;

namespace OFELI {

int OptimTN(OptSolver&          opt,
            Vect<real_t>&       x,
            const Vect<real_t>& lb,
            const Vect<real_t>& ub,
            int&                nb_obj_eval,
            int&                nb_grad_eval,
            int&                max_it,
            real_t              toler,
            int                 verb)
{
   int ret;
   size_t n=x.size();

// Set parameters for the optimization routine
   if (max_it==-1) {
      max_it = n / 2;
      if (n<=4)
         max_it = 50;
   }
   if (max_it>50)
      max_it = 50;
   if (max_it<=0)
      max_it = 1;
   real_t f = 0;

// Minimize function
   Vect<real_t> g(n);
   vector<int> pivot(n);
   real_t accrcy = toler;
   size_t max_fun = 150*n;
   real_t eta=0.25, stepmx=10., xtol=sqrt(accrcy);
   ret = lmqnbc(opt,x,f,g,lb,ub,pivot,verb,max_it,max_fun,eta,stepmx,accrcy,xtol,
                nb_obj_eval,nb_grad_eval);

// Print results
   if (ret)
      cout << "\n\nError Code = " << ret << endl;
   cout << "\n\nOptimal Function Value = " << f << endl;
   return ret;
}


int lmqnbc(OptSolver&          opt,
           Vect<real_t>&       x,
           real_t&             f,
           Vect<real_t>&       g,
           const Vect<real_t>& lb,
           const Vect<real_t>& ub,
           vector<int>&        pivot,
           int                 verb,
           int&                max_it,
           size_t              max_fun,
           real_t&             eta,
           real_t&             stepmx,
           real_t&             accrcy,
           real_t&             xtol,
           int&                nb_obj_eval,
           int&                nb_grad_eval)
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
   int numf=0, niter=0, icycle, newcon, lreset;
   real_t fm=0., fnew, peps, rtol, yksk, tiny, yrsr, alpha, fkeep;
   real_t oldf, flast, sm, gnorm, ftest, pnorm, rteps, xnorm;
   real_t pe, difold, difnew, reltol, gtpnew, toleps;
   real_t epsred, abstol, oldgtp, rtleps, rtolsq, tnytol, gsk, spe;

// Check that initial x is feasible and that the bounds are consistent
   int ier = crash(x,pivot,lb,ub);
   if (ier) {
      cout << "There is no feasible point; Terminating algorithm." << endl;
      return ier;
   }
   if (verb>0)
      cout << "\n\n   NIT  NF   CG           F                  GTG\n" << endl;

// Initialize variables
   nb_obj_eval = 0;
   int upd1=1, ireset=0, nmodif=0, nlincg=0, conv=0;
   Vect<real_t> w[15];
   size_t n=x.size();
   for (size_t i=0; i<15; i++)
      w[i].setSize(n);

// Within this routine the array w(loldg) is shared by w(lhyr)
   int lhyr=9*n, lwtest=14*n;

// Check parameters and set constants
   int nwhy = chkucp(lwtest,max_fun,alpha,eta,peps,rteps,rtol,rtolsq,stepmx,ftest,
                     xtol,xnorm,x,sm,tiny,accrcy);
   if (nwhy<0)
      return nwhy;

   oldf = flast = fnew = opt.Objective(x);
   opt.Gradient(x,g);
   int nftotl = 1;

/* Test the Lagrange multipliers to see if they are non-negative.
   because the constraints are only lower bounds, the components
   of the gradient corresponding to the active constraints are the
   Lagrange multipliers. Afterwards, the projected gradient is formed. */

   for (size_t i=0; i<n; i++) {
      if (pivot[i] != 2)
         if (-pivot[i]*g[i] < 0.)
            pivot[i] = 0;
   }
   ztime(g,pivot);
   real_t gtg = (g,g);
   if (verb>0)
      monit(fnew,g,pivot,niter,nftotl,nb_obj_eval);

// Check if the initial point is a local minimum
   ftest = fabs(fnew) + 1.;
   if (gtg < 1.e-4*OFELI_EPSMCH*ftest*ftest)
      goto L130;

// Set initial values to other parameters
   icycle = n - 1;
   toleps = rtol + rteps;
   rtleps = rtolsq + OFELI_EPSMCH;
   gnorm = sqrt(gtg);
   difnew = 0.;
   epsred = 0.05;
   fkeep = fnew;

// Set the diagonal of the approximate hessian to unity
   w[6] = 1.;

/* ..................start of main iterative loop.......... */

// Compute the new search direction
   modlnp(opt,verb-3,w[12],w[0],w[1],w[3],w[6],w[13],x,g,w[2],max_it,nb_obj_eval,
          nb_grad_eval,nmodif,nlincg,upd1,yksk,gsk,yrsr,lreset,1,pivot,accrcy,gtpnew,
          gnorm,xnorm,lhyr,w);

L20:
   w[9] = g;
   pnorm = Nrm2(w[12]);
   oldf = fnew;
   oldgtp = gtpnew;

// Prepare to compute the step length
   pe = pnorm + OFELI_EPSMCH;

// Compute the absolute and relative tolers for the linear search
   reltol = rteps * (xnorm + 1.) / pe;
   abstol = -OFELI_EPSMCH * ftest / (oldgtp - OFELI_EPSMCH);

// Compute the smallest allowable spacing between points in the linear search
   tnytol = OFELI_EPSMCH * (xnorm + 1.) / pe;
   spe = stpmax(stepmx,pe,x,w[12],pivot,lb,ub);

// Set the initial step length
   alpha = step1(fnew,fm,oldgtp,spe);

// Perform the linear search
   nwhy = linder(opt,sm,reltol,abstol,tnytol,eta,spe,w[12],oldgtp,x,fnew,alpha,
                 g,numf,w);
   newcon = 0;
   if (fabs(alpha-spe) <= OFELI_EPSMCH*10.) {
      newcon = 1;
      nwhy = 0;
      modz(x,w[12],pivot,lb,ub,flast,fnew);
      flast = fnew;
   }

   if (verb>2)
      cout << "        Linesearch results: alpha, pnorm: " << alpha << " " << pnorm << endl;
   ++niter;
   nftotl += numf;

// If required, print the details of this iteration
   if (verb>0)
      monit(fnew,g,pivot,niter,nftotl,nb_obj_eval);
   if (nwhy<0)
      return nwhy;
   if (nwhy==0 || nwhy==2)
      goto L40;

// The linear search has failed to find a lower point
   nwhy = 3;
   goto L140;

L40:
   if (nwhy > 1) { 
      fnew = opt.Objective(x);
      opt.Gradient(x,g);
      ++nftotl;
   }

// Terminate if more than max_fun evaluations have been made
   nwhy = 2;
   if (nftotl <= int(max_fun)) {
      nwhy = 0;

//    Set up parameters used in convergence and resetting tests
      difold = difnew;
      difnew = oldf - fnew;

//    If this is the first iteration of a new cycle, compute the
//    percentage reduction factor for the resetting test
      if (icycle==1) {
         if (difnew > difold * 2.0)
            epsred += epsred;
         if (difnew < difold * 0.5)
            epsred *= 0.5;
      }
      w[0] = g;
      ztime(w[0],pivot);
      gtg = w[0]*w[0];
      gnorm = sqrt(gtg);
      ftest = fabs(fnew) + 1.;
      xnorm = Nrm2(x);

//    Test for convergence
      conv = cnvtst(alpha,pnorm,toleps,xnorm,difnew,rtleps,ftest,gtg,peps,gtpnew,
                    fnew,flast,g,pivot,accrcy);
      if (conv)
         goto L130;
      ztime(g,pivot);

//    Compute the change in the iterates and the corresponding change in the gradients
      if (!newcon) {
         size_t isk=4*n, ipk=12*n, iyk=5*n, ioldg=9*n;
         for (size_t i=0; i<n; ++i) {
            w[0][iyk++] = g[i] - w[0][ioldg++];
            w[0][isk++] = alpha * w[0][ipk++];
         }

//       Set up parameters used in updating the preconditioning strategy
         yksk = (w[5],w[4]);
         lreset = 0;
         if (icycle==int(n)-1 || difnew<epsred*(fkeep-fnew))
            lreset = 1;
         if (!lreset) {
            yrsr = w[8]*w[7];
            if (yrsr <= 0.)
               lreset = 1;
         }
         upd1 = 0;
      }

//    Compute new search direction
      if (upd1 && verb>2)
         cout << "upd1 is true - Trivial Preconditioning" << endl;
      if (newcon && verb>=3)
         cout << "newcon is true - Constraint added in linesearch" << endl;
      modlnp(opt,verb-3,w[12],w[0],w[1],w[3],w[6],w[13],x,g,w[2],max_it,
             nb_obj_eval,nb_grad_eval,nmodif,nlincg,upd1,yksk,gsk,yrsr,lreset,1,
             pivot,accrcy,gtpnew,gnorm,xnorm,lhyr,w);   
      if (newcon)
         goto L20;
      if (lreset)
         goto L110;

//    Compute the accumulated step and its corresponding gradient difference
      Xpy(w[4],w[7]);
      Xpy(w[5],w[8]);
      ++icycle;
      goto L20;

//    Reset
L110:
      ++ireset;

//    Initialize the sum of all the changes in x
      w[7] = w[4];
      w[8] = w[5];
      fkeep = fnew;
      icycle = 1;
      goto L20;

/* ...............end of main iteration....................... */

L130:
      f = fnew;
      return 0;

L140:
      oldf = fnew;
   }

//  Local search could be installed here
   f = oldf;
   if (verb>0)
      monit(f,g,pivot,niter,nftotl,nb_obj_eval);
   return nwhy;
}


void monit(real_t              f,
           const Vect<real_t>& g,
           const vector<int>&  pivot,
           int                 niter,
           int                 nftotl,
           int                 nb_obj_eval)
{
   cout.setf(ios::scientific);
   real_t gtg = 0.;
   for (size_t i=0; i<g.size(); i++)
      if (pivot[i] == 0)
         gtg += g[i]*g[i];
   cout << setw(5) << niter << setw(5) << nftotl << setw(5) << nb_obj_eval
        << "  " << setprecision(8) << setw(18) << f << "  "
        << setprecision(8) << setw(18) << gtg << endl;
}


void ztime(vector<real_t>&    x,
           const vector<int>& pivot)
{
   for (size_t i=0; i<x.size(); ++i)
      if (pivot[i])
         x[i] = 0.0;
}


real_t stpmax(real_t                stepmx,
              real_t                pe,
              const Vect<real_t>&   x,
              const vector<real_t>& p,
              const vector<int>&    pivot,
              const Vect<real_t>&   lb,
              const Vect<real_t>&   ub)
{
   real_t spe = stepmx/pe;
   for (size_t i=0; i<x.size(); ++i) {
      if (pivot[i] == 0) {
         if (p[i] != 0.) {
            if (p[i] <= 0.) {
               real_t t = lb[i] - x[i];
               if (t > spe*p[i])
                  spe = t/p[i];
            }
            else {
               real_t t = ub[i] - x[i];
               if (t < spe*p[i])
                  spe = t/p[i];
            }
         }
      }
   }
   return (spe);
}


void modz(Vect<real_t>&         x,
          const vector<real_t>& p,
          vector<int>&          pivot,
          const Vect<real_t>&   lb,
          const Vect<real_t>&   ub,
          real_t&               flast,
          real_t&               fnew)
{
   real_t tol;
   for (size_t i=0; i<x.size(); i++) {
      if (pivot[i]==0) {
         if (p[i] != 0.) {
            if (p[i] <= 0.) {
               tol = OFELI_EPSMCH * 10. * (fabs(lb[i])+1.);
               if (x[i]-lb[i] <= tol) {
                  flast = fnew;
                  pivot[i] = -1;
                  x[i] = lb[i];
               }
            }
            else {
               tol = OFELI_EPSMCH * 10. * (fabs(ub[i])+10.);
               if (ub[i]-x[i] <= tol) {
                  flast = fnew;
                  pivot[i] = 1;
                  x[i] = ub[i];
               }
            }
         }
      }
   }
}


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
           real_t              accrcy)
{
   real_t cmax=0.;
   int conv=0, imax=0;
   int ltest = flast - fnew <= gtpnew * -0.5;
   for (size_t i=0; i<=g.size(); i++) {
      if (pivot[i] != 0 && pivot[i] != 2) {
         real_t t = -pivot[i] * g[i];
         if (t < 0.) {
            conv = 0;
            if (!ltest) {
               if (cmax > t) {
                  cmax = t;
                  imax = i;
               }
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
   if ((alpha*pnorm >= toleps*(xnorm+1.) || fabs(difnew) >=
      rtleps*ftest || gtg >= peps*ftest*ftest) && gtg >= accrcy*1e-4*ftest*ftest)
      return 0;
   return 1;
}


int modlnp(OptSolver&    opt,
           int           modet,
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
           Vect<real_t>* w)
{
   real_t beta=0, qnew=0, alpha=0, delta=0, rzold=0, rnorm=0, qtest=0, pr;
   nmodif = 0;
   size_t n=x.size();

// Initialization

// General initialization
   if (modet>0)
      cout << "\n\nEntering modlnop" << endl;
   if (max_it == 0)
      return 0;

   int first = 1;
   real_t rhsnrm=gnorm, tol=1e-12, qold=0.;

// Initialization for preconditioned conjugate-gradient algorithm
   initpc(diagb,emat,w,modet,upd1,yksk,gsk,yrsr,lreset);
   for (size_t i=0; i<n; i++)
      r[i] = -g[i];
   zsol = 0.;
   v = 0.;

// Main Iteration
   int k;
   for (k=1; k<=max_it; k++) {
      ++nlincg;
      if (modet>1)
         cout << "\n\n### Iteration " << k << " ###" << endl;

//    CG iteration to solve system of equations
      if (bounds)
        ztime(r,pivot);
      msolve(r,zk,w,upd1,yksk,gsk,yrsr,lreset,first,lhyr);
      if (bounds)
        ztime(zk,pivot);
      real_t rz = r*zk;
      if (rz/rhsnrm < tol)
         goto L80;
      if (k==1)
         beta = 0.;
      else
         beta = rz/rzold;
      for (size_t j=0; j<x.size(); j++)
         v[j] = zk[j] + beta*v[j];
      if (bounds)
        ztime(v,pivot);
      gtims(opt,v,gv,x,g,w,first,delta,accrcy,xnorm);
      if (bounds)
        ztime(gv,pivot);
      ++nb_obj_eval, ++nb_grad_eval;
      real_t vgv = v*gv;
      if (vgv/rhsnrm < tol)
         goto L50;
      ndia3(emat,v,gv,r,vgv,modet);

//    Compute linear step length
      alpha = rz / vgv;
      if (modet >= 1)
         cout << "alpha = " << alpha << endl;

//    Compute current solution and related vectors
      Axpy(alpha,v,zsol);
      Axpy(-alpha,gv,r);

//    Test for convergence
      gtp = zsol*g;
      pr =  r*zsol;
      qnew = (gtp + pr)*0.5;
      qtest = k * (1. - qold/qnew);
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
   max_it = k;
   goto L70;

// Truncate algorithm in case of an emergency
L40:
   if (modet>=-1)
      cout << "g(t) positive at iteration " << k << " - truncating method" << endl;
   Axpy(-alpha,v,zsol);
   gtp = zsol*g;
   goto L90;

L50:
  if (modet>-2)
     cout << "         *** Hessian not positive-definite ***" << endl;

  if (k > 1)
     goto L70;

   msolve(g,zsol,w,upd1,yksk,gsk,yrsr,lreset,first,lhyr);
   for (size_t i=0; i<n; i++)
      zsol[i] = -zsol[i];
   if (bounds)
     ztime(zsol,pivot);
   gtp = zsol*g;

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
   zsol = g;
   for (size_t i=0; i<n; i++)
      zsol[i] = -zsol[i];
   if (bounds)
     ztime(zsol,pivot);
   gtp = zsol*g;
   goto L70;

// Store (or restore) diagonal preconditioning
L90:
   diagb = emat;
   return 0;
}


void ndia3(vector<real_t>&       e,
           const vector<real_t>& v,
           const vector<real_t>& gv,
           const vector<real_t>& r,
           real_t                vgv,
           int                   modet)
{
   real_t vr = v*r;
   for (size_t i=0; i<r.size(); i++) {
      e[i] -= r[i]*r[i]/vr - gv[i]*gv[i]/vgv;
      if (e[i] <= 1e-6) {
         if (modet > -2)
            cout << " *** emat negative:  " << e[i] << endl;
         e[i] = 1.;
      }
   }
}


int crash(Vect<real_t>&       x,
          vector<int>&        pivot,
          const Vect<real_t>& lb,
          const Vect<real_t>& ub)
{
   int er = 0;
   for (size_t i=0; i<x.size(); i++) {
      if (x[i] < lb[i])
         x[i] = lb[i];
      if (x[i] > ub[i])
         x[i] = ub[i];
      pivot[i] = 0;
      if (x[i] == lb[i])
         pivot[i] = -1;
      if (x[i] == ub[i])
         pivot[i] = 1;
      if (ub[i] == lb[i])
         pivot[i] = 2;
      if (lb[i] > ub[i])
         er = -int(i) - 1;
   }
   return er;
}


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
           real_t xlamda)
{
   xlamda = 0;
   real_t yu=xmin+u, ya=a+xmin, yb=b+xmin, yw=xw+xmin, ybnd=scxbd+xmin;
   cout << "\n\n\nOutput from linear search" << endl;
   cout << "  Tol and eps: " << tol << "  " << eps << endl;
   cout << "  Current upper and lower bounds: " << ya << "  " << yb << endl;
   cout << "  Strict upper bound: " << ybnd << endl;
   cout << "  xw, fw, gw: " << yw << ", " << fw << ", " << gw << endl;
   cout << "  xmin, fmin, gmin: " << xmin << ", " << fmin << ", " << gmin << endl;
   cout << "  New estimate: " << yu << endl;
   cout << "  iloc, itest: " << loc << ", " << test << endl;
}


real_t step1(real_t fnew,
             real_t fm,
             real_t gtp,
             real_t smax)
{
   real_t alpha=1., d=fabs(fnew-fm);
   if (2.*d <= -gtp && d >= OFELI_EPSMCH)
      alpha = -2.*d/gtp;
   if (alpha>=smax)
      alpha = smax;
   return alpha;
}


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
           real_t              accrcy)
{
   int lw=14*x.size(), nwhy=-1;
   sm = OFELI_EPSMCH*OFELI_EPSMCH;
   tiny = sm;
   rteps = sqrt(OFELI_EPSMCH);
   rtol = xtol;
   if (fabs(rtol)<accrcy)
      rtol = rteps*10.;

// Check for errors in the input parameters
   if (lw<lwtest || x.size()<1 || rtol<0. || eta>=1. || eta<0. || stepmx<rtol || max_fun<1)
      return nwhy;
   nwhy = 0;

// Set constants for later
   rtolsq = rtol*rtol;
   peps = pow(accrcy,0.6666);
   xnorm = Nrm2(x);
   alpha = tst = 0.;
   return nwhy;
}


void gtims(OptSolver&          opt,
           const Vect<real_t>& v,
           Vect<real_t>&       gv,
           const Vect<real_t>& x,
           const Vect<real_t>& g,
           Vect<real_t>*       w,
           int&                first,
           real_t&             delta,
           real_t              accrcy,
           real_t              xnorm)
{
   if (first) {
      delta = sqrt(accrcy)*(xnorm + 1.);
      first = 0;
   }
   real_t dinv = 1./delta;
   w[10] = x + delta*v;
   opt.Gradient(w[10],gv);
   gv = (gv - g)*dinv;
}


void ssbfgs(const vector<real_t>& sj,
            const vector<real_t>& hjv,
            const vector<real_t>& hjyj,
            real_t                yjsj,
            real_t                yjhyj,
            real_t                vsj,
            real_t                vhyj,
            vector<real_t>&       hjp1v)
{
   real_t gamma = 1.;
   real_t delta = (gamma*yjhyj/yjsj + 1.) * vsj/yjsj - gamma*vhyj/yjsj;
   real_t beta = -gamma*vsj/yjsj;
   for (size_t i=0; i<sj.size(); ++i)
      hjp1v[i] = gamma*hjv[i] + delta*sj[i] + beta*hjyj[i];
}


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
          int                   first)
{
    real_t ghyk=0, ghyr=0, yksr=0, ykhyk=0, ykhyr=0, yrhyr=0, rdiagb=0, gsr=0;
    size_t n=g.size();
    if (!upd1) {
       gsk = g*sk;
       if (!lreset) {

//        Compute hg and hy where h is the inverse of the diagonals
          for (size_t i=0; i<n; i++) {
             rdiagb = 1. / diagb[i];
             hg[i] = g[i]*rdiagb;
             if (first) {
                hyk[i] = yk[i] * rdiagb;
                hyr[i] = yr[i] * rdiagb;
             }
          }
          if (first) {
             yksr = yk*sr;
             ykhyr = yk*hyr;
          }
          gsr = g*sr; ghyr = g*hyr;
          if (first)
             yrhyr = yr*hyr;
          ssbfgs(sr, hg, hyr, yrsr, yrhyr, gsr, ghyr, hg);
          if (first)
             ssbfgs(sr, hyk, hyr, yrsr, yrhyr, yksr, ykhyr, hyk);
          ykhyk = hyk*yk; ghyk = hyk*g;
          ssbfgs(sk,hg,hyk,yksk,ykhyk,gsk,ghyk,y);
       }

//     Compute gh and hy where h is the inverse of the diagonals
       for (size_t i=0; i<n; i++) {
          rdiagb = 1./diagb[i];
          hg[i] = g[i]*rdiagb;
          if (first)
             hyk[i] = yk[i]*rdiagb;
       }
       if (first)
          ykhyk = yk*hyk;
       ghyk = g*hyk;
       ssbfgs(sk,hg,hyk,yksk,ykhyk,gsk,ghyk,y);
    }
    for (size_t i=0; i<n; i++)
       y[i] = g[i] / diagb[i];
}


void msolve(const vector<real_t>& g,
            vector<real_t>&       y,
            Vect<real_t>*         w,
            int&                  upd1,
            real_t&               yksk,
            real_t&               gsk,
            real_t&               yrsr,
            int                   lreset,
            int                   first,
            int                   lhyr)
{
   mslv(g,y,w[4],w[5],w[6],w[7],w[8],w[lhyr/g.size()],w[10],
        w[11],upd1,yksk,gsk,yrsr,lreset,first);
}


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
           int&                  upd1)
{
   size_t n=diagb.size();
   if (!upd1) {
      if (!lreset) {
         for (size_t i=0; i<n; i++)
            bsk[i] = diagb[i]*sr[i];
         real_t sds=sr*bsk, srds=sk*bsk, yrsk=yr*sk;
         for (size_t i=0; i<n; i++) {
            real_t td = diagb[i];
            bsk[i] = td*sk[i] - bsk[i]*srds/sds + yr[i]*yrsk/yrsr;
            emat[i] = td - td*td*sr[i]*sr[i]/sds + yr[i]*yr[i]/yrsr;
         }
         sds = sk*bsk;
         for (size_t i=0; i<n; i++)
            emat[i] -= bsk[i]*bsk[i]/sds - yk[i]*yk[i]/yksk;
         goto L110;
      }
      for (size_t i=0; i<n; i++)
         bsk[i] = diagb[i]*sk[i];
      real_t sds = 0.;
      for (size_t i=1; i<n; i++)
         sds += sk[i] * bsk[i];
      for (size_t i=0; i<n; i++) {
         real_t td = diagb[i];
         emat[i] = td - td*td*sk[i]*sk[i]/sds + yk[i]*yk[i]/yksk;
      }
      goto L110;
   }

   emat = diagb;

L110:
   if (modet<1)
      return 0;
   real_t d1=emat[0];
   real_t dn=d1;
   for (size_t i=0; i<n; i++) {
      if (emat[i] < d1)
         d1 = emat[i];
      if (emat[i] > dn)
         dn = emat[i];
   }
   cout << "\n\n        dmin = " << d1 << ",  dmax = " << dn
        << ",  cond = " << dn/d1 << endl;
   return 0;
}


void initpc(const vector<real_t>& diagb,
            vector<real_t>&       emat,
            Vect<real_t>*         w,
            int&                  modet,
            int&                  upd1,
            real_t&               yksk,
            real_t&               gsk,
            real_t&               yrsr,
            int&                  lreset)
{
   gsk = 0;
   initp3(diagb,emat,lreset,yksk,yrsr,w[11],w[4],w[5],w[7],w[8],modet,upd1);
}


/*==============================================================================*/
/*               Line Search algorithms of Gill and Murray                      */
/*==============================================================================*/
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
           Vect<real_t>*         w)
{
   size_t n=x.size();
   int tst, flag=0;
   real_t oldf, gmin, step, xmin=0., a, b, e, b1, gtest1, gtest2;
   real_t fw, gw, scxbnd, xw, tol;

   int lsprnt=0, nprnt=10000;
   real_t rtsmll=sqrt(sm), big=1./sm;

// Set the estimated relative precision in f(x)
   real_t fpresn = OFELI_EPSMCH * 10.;
   int numf=0;
   real_t u=alpha, fu=f, fmin=f, gu=gtp, rmu=1e-4;

// First entry sets up the initial interval of uncertainty
   int ientry = 1;

// Test for too many iterations
   for (int itcnt=1; itcnt<=20; itcnt++) {
      flag = 0;
      int braktd;
      real_t factor;
      getptc(big,rtsmll,reltol,abstol,tnytol,fpresn,eta,rmu,xbnd,u,fu,gu,xmin,
             fmin,gmin,xw,fw,gw,a,b,oldf,b1,scxbnd,e,step,factor,braktd,gtest1,
             gtest2,tol,ientry,tst);
      if (lsprnt>=nprnt)
         lsout(ientry,tst,xmin,fmin,gmin,xw,fw,gw,u,a,b,tol,reltol,scxbnd,xbnd);

//    If test=1, the algorithm requires the function value to be calculated
      if (tst==1) {
         real_t ualpha = xmin + u;
         for (size_t i=0; i<n; i++)
            w[0][i] = x[i] + ualpha*p[i];
         fu = opt.Objective(w[0]);
         opt.Gradient(w[0],w[1]);
         ++numf;
         gu = w[1]*p;

//       The gradient vector corresponding to the best point is overwritten if fu
//       is less than fmin and fu is sufficiently lower than f at the origin.
         if (fu<=fmin && fu<=oldf-ualpha*gtest1)
            g = w[1];
         flag = 1;
      }

//    If test=2 or 3 a lower point could not be found
      else {
         nftotl = numf;
         flag = 1;
         if (tst==0) {

//          A successful search has been made
            flag = 0;
            f = fmin;
            alpha = xmin;
            for (size_t j=0; j<x.size(); ++j)
               x[j] += alpha*p[j];
         }
         return flag;
      }
   }
   return flag;
}


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
           int&    itest)
{
   real_t d, abgw, absr, p, q, r, s, sc, denom, a1, d1, d2;
   real_t xmidpt, twotol, sumsq, abgmin, chordm, chordu;
   int convrg;

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
   twotol = 2*tol;

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

   denom = 1.e-15;
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
   fw = fmin; fmin = fu;
   gw = gmin; gmin = gu;
   xmin += u;
   a -= u; b -= u;
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
   reltol *= 0.1; abstol *= 0.1; twotol *= 0.1;

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
   d1 = xw, d2 = a;
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
   if (step>=scxbnd) {
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


} /* namespace OFELI */
