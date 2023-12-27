/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2024 Rachid Touzani

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

                       Implementation of class 'OptSolver'

  ==============================================================================*/

#include "solvers/OptSolver.h"
#include "util/util.h"

#include <iostream>
using std::cout;

namespace OFELI {


/** \brief Simulated annealing optimization solver.

   \details Simulated annealing is a global optimization method that distinguishes
   between different local optima. Starting from an initial point, the
   algorithm takes a step and the function is evaluated. When minimizing a
   function, any downhill step is accepted and the process repeats from this
   new point. An uphill step may be accepted. Thus, it can escape from local
   optima. This uphill decision is made by the Metropolis criteria. As the
   optimization process proceeds, the length of the stx decline and the
   algorithm closes in on the global optimum. Since the algorithm makes very
   few assumptions regarding the function to be optimized, it is quite
   robust with respect to non-quadratic surfaces. The degree of robustness
   can be adjusted by the user. In fact, simulated annealing can be used as
   a local optimizer for difficult functions.
  
   This implementation of simulated annealing was used in "Global Optimization
   of Statistical Functions with Simulated Annealing," Goffe, Ferrier and
   Rogers, Journal of Econometrics, vol. 60, no. 1/2, Jan./Feb. 1994, pp. 65-100.
   Briefly, we found it competitive, if not superior, to multiple
   restarts of conventional optimization routines for difficult optimization problems.
  
   For more information on this routine, contact its author:
        Bill Goffe, bgoffe@whale.st.usm.edu
  
  
    \b Synopsis:
    This function implements the continuous simulated annealing global
    optimization algorithm described in Corana et al.'s article
    "Minimizing Multimodal Functions of Continuous Variables with the
    "Simulated Annealing" Algorithm" in the September 1987 (vol. 13,
    no. 3, pp. 262-280) issue of the ACM Transactions on Mathematical
    Software.
  
    A very quick (perhaps too quick) overview of OptimSA:\n
    OptimSA tries to find the global optimum of an N dimensional function.
    It moves both up and downhill and as the optimization process
    proceeds, it focuses on the most promising area.\n
    To start, it randomly chooses a trial point within the step length
    vm of the user selected starting point. The function is evaluated at
    this trial point and its value is compared to its value at the initial point.\n
    In a maximization problem, all uphill moves are accepted and the
    algorithm continues from that trial point. Downhill moves may be
    accepted; the decision is made by the Metropolis criteria. It uses T
    (temperature) and the size of the downhill move in a probabilistic
    manner. The smaller \a t and the size of the downhill move are, the more
    likely that move will be accepted. If the trial is accepted, the
    algorithm moves on from that point. If it is rejected, another point
    is chosen instead for a trial evaluation.\n
    Each element of \a vm periodically adjusted so that half of all
    function evaluations in that direction are accepted.

    A fall in \a t is imposed upon the system with the rt variable by
    \a t(i+1) = rt*t(i) where \a i is the \a i-th iteration. Thus, as \a t declines,
    downhill moves are less likely to be accepted and the percentage of
    rejections rise. Given the scheme for the selection for \a vm, \a vm falls.
    Thus, as \a t declines, \a vm falls and OptimSA focuses upon the most promising
    area for optimization.
  
    \b The \b importance \b of \b the \b parameter \a t:

    The parameter \a t is crucial in using OptimSA successfully. It influences
    \a vm, the step length over which the algorithm searches for optima. For
    a small initial \a t, the step length may be too small; thus not enough
    of the function might be evaluated to find the global optima. The user
    should carefully examine vm in the intermediate output (set msg_lvl =
    1) to make sure that vm is appropriate. The relationship between the
    initial temperature and the resulting step length is function
    dependent.\n
    To determine the starting temperature that is consistent with
    optimizing a function, it is worthwhile to run a trial run first. Set
    rt = 1.5 and t = 1.0. With rt > 1.0, the temperature increases and vm
    rises as well. Then select the \a T that produces a large enough vm.
  
    For modifications to the algorithm and many details on its use,
    (particularly for econometric applications) see Goffe, Ferrier
    and Rogers, "Global Optimization of Statistical Functions with
    Simulated Annealing," Journal of Econometrics, vol. 60, no. 1/2,
    Jan./Feb. 1994, pp. 65-100.\n
    For more information, contact\n
                Bill Goffe\n
                Department of Economics and International Business\n
                University of Southern Mississippi\n
                Hattiesburg, MS  39506-5072\n
                (601) 266-4484 (office)\n
                (601) 266-4920 (fax)\n
                bgoffe@whale.st.usm.edu (Internet)
  
    As far as possible, the parameters here have the same name as in
    the description of the algorithm on pp. 266-8 of Corana et al.
  
    \b Note: The suggested values generally come from Corana et al. To
             drastically reduce runtime, see Goffe et al., pp. 90-1 for
             suggestions on choosing the appropriate \p rt and \p nt.
  
    @param [in] theOpt Instance of class \a OPT_ that is implemented by the user and that provides the objective function.

    @param [in] x The starting values for the variables of the function to be optimized.

    @param [in] rt The temperature reduction factor. The value suggested by
    Corana et al. is \a .85. See Goffe et al. for more advice.

    @param [in] toler Error tolerance for termination. If the final function
    values from the last neps temperatures differ from the
    corresponding value at the current temperature by less than
    eps and the final function value at the current temperature
    differs from the current optimal function value by less than
    eps, execution terminates and the value \a 0 is returned.

    @param [in] ns Number of cycles. After \a ns*n function evaluations, each
    element of vm is adjusted so that approximately half of
    all function evaluations are accepted. The suggested value is 20.

    @param [in] nt Number of iterations before temperature reduction. After
    \a nt*ns*n function evaluations, temperature \a (t) is changed
    by the factor \a rt. Value suggested by Corana et al. is
    \a max(100,5*n). See Goffe et al. for further advice.

    @param [in] neps Number of final function values used to decide upon termination.
    See \c eps. Suggested value is \a 4

    @param max_eval [in] The maximum number of function evaluations. If it is
    exceeded, the return \a code=1.

    @param [in] lb The lower bound for the allowable solution variables.

    @param [in] ub The upper bound for the allowable solution variables.
    If the algorithm chooses \a x(i) \a < \a lb(i) or \a x(i) \a > \a ub(i)
    \a i=1,...,n, a point is from inside is randomly selected.
    This focuses the algorithm on the region inside \a ub and \a lb.
    Unless the user wishes to concentrate the search to a particular region, 
    \a ub and \a lb should be set to very large positive
    and negative values, respectively. Note that the starting
    vector \a x should be inside this region. Also note that \a lb and
    \a ub are fixed in position, while \a vm is centered on the last
    accepted trial set of variables that optimizes the function.
    @param c [in] Vector that controls the step length adjustment. The suggested
    value for all elements is 2.

    @param [in] msg_lvl controls printing inside \a OptimSA.
    <ul>
       <li>0 - Nothing printed.
       <li>1 - Function value for the starting value and
               summary results before each temperature
               reduction. This includes the optimal
               function value found so far, the total
               number of moves (broken up into uphill,
               downhill, accepted and rejected), the
               number of out of bounds trials, the
               number of new optima found at this
               temperature, the current optimal x and
               the step length \a vm. Note that there are
               \a n*ns*nt function evalutations before each
               temperature reduction. Finally, notice is
               is also given upon achieveing the termination criteria.
       <li>2 - Each new step length \a (vm), the current optimal
               \a x and the current trial \a x (x). This
               gives the user some idea about how far \a x
               strays from \a x as well as how \a vm is adapting
               to the function.
       <li>3 - Each function evaluation, its acceptance or
               rejection and new optima. For many problems,
               this option will likely require a small tree
               if hard copy is used. This option is best
               used to learn about the algorithm. A small
               value for \a MAXEVL is thus recommended when
               using \a msg_lvl=3.
               Suggested value: 1
    </ul>
    Note: For a given value of \a msg_lvl, the lower valued
          options (other than 0) are utilized.

    @param [in] t On input, the initial temperature. See Goffe et al. for advice.
    On output, the final temperature.

    @param [in] vm The step length vector. On input it should encompass the
    region of interest given the starting value \a x. For point
    x[i], the next trial point is selected is from \a x[i]-vm[i]
    to \a x[i]+vm[i]. Since \a vm is adjusted so that about half
    of all points are accepted, the input value is not very
    important (i.e. is the value is off, \a OptimSA adjusts \a vm to the
    correct value).
  
    @param [out] fopt The optimal value of the function.

    @param [out] nacc The number of accepted function evaluations.

    @param [out] nfcnev The total number of function evaluations. In a minor
    point, note that the first evaluation is not used in the
    core of the algorithm; it simply initializes the algorithm.

    @param [out] nobds The total number of trial function evaluations that
    would have been out of bounds of \a lb and \a ub. Note that
    a trial point is randomly selected between \a lb and \a ub.

    @return
    <ul>
        <li>0 - Normal return; termination criteria achieved.\n
        <li>1 - Number of function evaluations \a (nfcnev) is
                greater than the maximum number \a (max_eval).\n
        <li>2 - The starting value \a (x) is not inside the
                bounds (\a lb and \a ub).\n
        <li>3 - The initial temperature is not positive.\n
                99 - Should not be seen; only used internally.\n
     </ul>

 */
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
            int&                nb_obj_eval,
            int&                nobds)
{
   bool quit=false;
   int nrej, nnew, nup, totmov;
   real_t p, fp, pp;
   size_t n=x.size();
   Vect<real_t> fstar(neps), xp(n), nacp(n);

   real_t ranmar(void);
   
   srand((unsigned)time(nullptr));
   nacc = nobds = nb_obj_eval = 0;
   nacp = 0;
   fstar = 1.e20;

// If the initial temperature is not positive, notify the user and
// return to the calling routine.
   if (t <= 0.) {
      cout << "The initial temperature is not positive, reset the variable t.\n";
      return 3;
   }

// If the initial value is out of bounds, notify the user and return to the calling routine.
   for (size_t i=0; i<n; ++i) {
      if (x[i]>ub[i] || x[i]<lb[i]) {
         cout << "The starting value (x) is outside the bounds (lb and ub)\n";
         cout << "Execution terminated without any optimization.\n";
         cout << "Respecify x, ub or lb so that lb(i) < x(i) < ub(i), i=1,n.\n";
         return 2;
      }
   }

// Evaluate the function with input x and return value as f.
   real_t f = -opt.Objective(x);

// If the function is to be minimized, switch the sign of the function.
// Note that all intermediate and final output switches the sign back
// to eliminate any possible confusion for the user.
   ++nb_obj_eval;
   fopt = f;
   fstar[0] = f;
   if (Verbosity >= 1) {
      cout << endl;
      cout << "Initial x\n" << x;
      cout << "Initial f: " << -f << endl;
   }

// Start the main loop. Note that it terminates if
//       (i) the algorithm successfully optimizes the function or
//      (ii) there are too many function evaluations (more than MAXEVL).

L100:
   int lnobds=0, ndown=0;
   nup = nrej = nnew = 0;
   for (int m=1; m<=nt; ++m) {
      for (int j=0; j<ns; ++j) {
         for (size_t k=0; k<n; ++k) {

//          Generate xp, the trial value of x. Note use of vm to choose xp.
            for (size_t i=0; i<n; ++i) {
               if (i==k)
//                  xp[i] = x[i] + (2*ranmar() - 1) * vm[i];
                  xp[i] += (2*rand() - 1) * vm[i];

//             If xp is out of bounds, select a point in bounds for the trial.
               if (x[i]<lb[i] || x[i]>ub[i]) {
//                x[i] = lb[i] + (ub[i] - lb[i]) * ranmar();
                  x[i] = lb[i] + (ub[i] - lb[i]) * rand();
                  ++lnobds, ++nobds;
                  if (Verbosity >= 3) {
                     cout << endl;
                     cout << "Current x\n" << x;
                     cout << "Current f : " << -f << endl;
                     cout << "Trial x\n" << xp;
                  }
               }
            }

//          Evaluate the function with the trial point xp and return as fp
            fp = -opt.Objective(xp);
            ++nb_obj_eval;
            if (Verbosity >= 3) {
               cout << endl << "Current x\n" << x;
               cout << "Current f : " << -f << endl;
               cout << "Trial x\n" << xp;
               cout << "Resulting f : " << -fp << endl;
            }

//          If too many function evaluations occur, terminate the algorithm
            if (nb_obj_eval >= max_eval) {
               cout << "Too many function evaluations. Consider increasing max_eval or toler,\n";
               cout << "or decreasing nt or rt. These results are likely to be poor\n";
               return 1;
            }

//          Accept the new point if the function value increases
            if (fp >= f) {
               if (Verbosity >= 3)
                  cout << "Point accepted\n";
               x = xp;
               f = fp;
               ++nacc, ++nacp[k], ++nup;

//             If greater than any other point, record as new optimum
               if (fp > fopt) {
                  if (Verbosity >= 3)
                     cout << "New Optimum\n";
                  x = xp;
                  fopt = fp;
                  ++nnew;
               }

//          If the point is lower, use the Metropolis criteria to decide on
//          acceptance or rejection.
            } else {
               p = exprep((fp-f)/t);
//             pp = ranmar_();
               pp = rand();
               if (pp < p) {
                  if (Verbosity >= 3)
                     cout << "Though higher, point accepted\n";
                  x = xp;
                  f = fp;
                  ++nacc, ++nacp[k], ++ndown;
               } else {
                  ++nrej;
                  if (Verbosity >= 3)
                     cout << "Higher point rejected\n";
               }
            }
         }
      }

//    Adjust vm so that approximately half of all evaluations are accepted.
      for (size_t i=0; i<n; ++i) {
         real_t ratio = real_t(nacp[i]) / real_t(ns);
         if (ratio > 0.6)
            vm[i] *= c[i] * (ratio - 0.6) / 0.4 + 1.;
         else if (ratio < 0.4)
            vm[i] /= c[i] * ((0.4 - ratio) / 0.4) + 1.;
         if (vm[i] > (ub[i]-lb[i]))
            vm[i] = ub[i] - lb[i];
      }
      if (Verbosity >= 2) {
         cout << "Intermediate results after step length adjustment\n";
         cout << "New step length (vm)\n" << vm;
         cout << "Current optimal x\n" << x;
         cout << "Current x\n" << x;
      }
      nacp = 0;
   }
   if (Verbosity >= 1) {
      totmov = nup + ndown + nrej;
      cout << "Intermediate results before next temperature reduction\n";
      cout << "Current Temperature:            " << t << endl;
      cout << "Min function value so far:      " << -fopt << endl;
      cout << "Total moves:                    " << totmov << endl;
      cout << "  Downhill:                     " << nup << endl;
      cout << "  Accepted Uphill:              " << ndown << endl;
      cout << "  Rejected Uphill:              " << nrej << endl;
      cout << "  Trials out of bounds:         " << lnobds << endl;
      cout << "  New minima this temperature:  " << nnew << endl;
      cout << "Current optimal x\n" << x;
      cout << "Step length (vm)\n" << vm;
   }

// Check termination criteria.
   quit = false;
   fstar[0] = f;
   if ((fopt - fstar[0]) <= toler)
      quit = true;
   for (int i=0; i<neps; ++i)
      if (fabs(f-fstar[i]) > toler)
         quit = false;

// Terminate SA if appropriate.
   if (quit) {
      fopt = -fopt;
      if (Verbosity >= 1)
         cout << "SA Achieved termination criteria.\n";
      return 0;
   }

// If termination criteria is not met, prepare for another loop.
   t *= rt;
   for (int i=neps-1; i>=1; --i)
      fstar[i] = fstar[i-1];
   f = fopt;
   goto L100;
}


/*			
    Required Functions (included):
      exprep - Replaces the function \a exp to avoid under- and overflows.\n
      ranmar - The actual random number generator. Note that
               rmarin must run first (OptimSA does this). It produces uniform
               random numbers on [0,1]. These routines are from
               Usenet's comp.lang.fortran. For a reference, see
               "Toward a Universal Random Number Generator"
               by George Marsaglia and Arif Zaman, Florida State
               University Report: FSU-SCRI-87-50 (1987).
               It was later modified by F. James and published in
               "A Review of Pseudo-random Number Generators." For
               further information, contact stuart@ads.com. These
               routines are designed to be portable on any machine
               with a 24-bit or more mantissa. I have found it produces
               identical results on a IBM 3081 and a Cray Y-MP.\n

    Machine Specific Features:
      1. \a exprep may have to be modified if used on non-IBM type main-
         frames. Watch for under- and overflows in exprep.
      2. Some format statements use G25.18; this may be excessive for
         some machines.
      3. \a rmarin and ranmar are designed to be protable; they should not
      cause any problems.
*/	  


real_t ranmar(int i1, int i2)
{
   static real_t c_=0.021602869033813477, cd_=0.45623308420181274, cm_=0.9999998211860657;
   static real_t u_[97];
   real_t uni = u_[i1-1] - u_[i2-1];
   if (uni < 0.)
      uni += 1.;
   u_[i1-1] = uni;
   --i1;
   if (i1 == 0)
      i1 = 97;
   --i2;
   if (i2 == 0)
      i2 = 97;
   c_ -= cd_;
   if (c_ < 0.)
      c_ += cm_;
   uni -= c_;
   if (uni < 0.)
      uni += 1.;
   return uni;
}

} /* namespace OFELI */
