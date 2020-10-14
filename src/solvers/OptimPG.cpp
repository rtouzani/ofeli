
/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

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

             Function to solve a bound constrained optimization problem
                         by the gradient projected method

         The implementation of this code is mainly inspired by the Matlab
         implementation of C.T. Kelley

  ==============================================================================*/

#include "OFELI_Config.h"
#include "solvers/OptSolver.h"
#include "util/util.h"

#include <iostream>
using std::cout;

namespace OFELI {


void ProjActiveSet(Vect<real_t>&       x,
                   const Vect<real_t>& lb,
                   const Vect<real_t>& ub,
                   Vect<real_t>&       px)
{
   for (size_t i=0; i<x.size(); ++i)
      px[i] = fmax(lb[i],fmin(ub[i],x[i]));
}


int OptimPG(OptSolver&          opt,
            Vect<real_t>&       x,
            const Vect<real_t>& lb,
            const Vect<real_t>& ub,
            int&                nb_obj_eval,
            int&                nb_grad_eval,
            int&                max_it,
            real_t              toler)
{
   size_t n = x.size();
   Vect<real_t> px(n), xc(n), xt(n), gc(n), pgc(n), pl(n), xcc(n);

// put initial iterate in feasible set
   xc = x;
   ProjActiveSet(xc,lb,ub,px);
   if ((xc-px).getNorm2()>0) {
      cout << "Warning in OptimPG(..): Initial iterate not feasible." << endl;
      ProjActiveSet(xc,lb,ub,xc);
   }
   real_t alp = 1.e-4;
   real_t fc = opt.Objective(xc);
   opt.Gradient(xc,gc);
   nb_obj_eval = nb_grad_eval = 1;
   xcc = xc - gc;
   ProjActiveSet(xcc,lb,ub,xt);
   ProjActiveSet(xcc,lb,ub,pgc);
   pgc = xc - pgc;

   int it=1;
   real_t lambda=0.;
   while ((pgc.getWNorm2()>toler) & (it<=max_it)) {
      lambda = 1;
      xcc = xc - lambda*gc;
      ProjActiveSet(xcc,lb,ub,xt);
      real_t ft = opt.Objective(xt);
      nb_obj_eval++, it++;
      size_t iarm=0;
      pl = xc - xt;
      real_t fgoal = fc - (pl,pl)*(alp/lambda);

//    Simple line search
      while (ft > fgoal) {
         lambda *= 0.1;
         iarm++;
         xcc = xc - lambda*gc;
         ProjActiveSet(xcc,lb,ub,xt);
         pl = xc - xt;
         ft = opt.Objective(xt);
         nb_obj_eval++;
         if (iarm>10) {
            cout << "OptimPG: Armijo error in gradient projection" << endl;
            max_it = it;
            return 1;
         }
         fgoal = fc - (pl,pl)*(alp/lambda);
      }
      xc = xt;
      fc = opt.Objective(xc);
      opt.Gradient(xc,gc);
      nb_obj_eval++, nb_grad_eval++;
      xcc = xc - gc;
      ProjActiveSet(xcc,lb,ub,pgc);
      pgc = xc - pgc;
      x = xc;
   }
   max_it = it;
   return 0;
}

} /* namespace OFELI */
