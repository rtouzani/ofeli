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

                       Implementation of class 'OptSolver'

  ==============================================================================*/

#include "solvers/OptSolver.h"
#include "solvers/Optim.h"
#include <limits>
#include <algorithm>
#include "OFELIException.h"

namespace OFELI {

OptSolver::OptSolver()
          : _size(1), _verb(1), _nb_obj_eval(0), _nb_grad_eval(0), _max_eval(100000),
            _max_it(1000), _toler(1.e-10), _exp(true), _sa_opt(false), _tn_opt(false),
            _obj_type(0), _x_set(true), _method_set(false), _var("")
{
   _x = new Vect<real_t>(1);
   _lb.setSize(_size);
   _ub.setSize(_size);
   _lb = -std::numeric_limits<real_t>::max();
   _ub =  std::numeric_limits<real_t>::max();
}


OptSolver::OptSolver(Vect<real_t>& x)
          : _size(x.size()), _verb(1), _nb_obj_eval(0), _nb_grad_eval(0), _max_eval(100000),
            _max_it(1000), _x(&x), _toler(1.e-10), _exp(true), _sa_opt(false), _tn_opt(false),
            _obj_type(0), _x_set(false), _method_set(false)
{
   _lb.setSize(_size);
   _ub.setSize(_size);
   _lb = -std::numeric_limits<real_t>::max();
   _ub =  std::numeric_limits<real_t>::max();
}


OptSolver::OptSolver(MyOpt&        opt,
                     Vect<real_t>& x)
          : _size(x.size()), _verb(1), _nb_obj_eval(0), _nb_grad_eval(0), _max_eval(100000),
            _max_it(1000), _x(&x), _toler(1.e-10), _opt(&opt), _exp(false), _sa_opt(false),
            _tn_opt(false), _obj_type(0), _x_set(false), _method_set(false)
{
   _lb.setSize(_size);
   _ub.setSize(_size);
   _lb = -std::numeric_limits<real_t>::max();
   _ub =  std::numeric_limits<real_t>::max();
}


OptSolver::~OptSolver()
{
   if (_x_set)
      delete _x;
}


void OptSolver::setOptMethod(OptMethod m)
{
   _opt_method = m;
   _method_set = true;
   if (_opt_method==TRUNCATED_NEWTON)
      setTNDefaults();
   else if (_opt_method==SIMULATED_ANNEALING)
      setSADefaults();
   else if (_opt_method==NELDER_MEAD)
      setNMDefaults();
   else if (_opt_method==GRADIENT)
      setPGDefaults();
}


void OptSolver::setUpperBound(real_t ub)
{
   _ub[0] = ub;
}


void OptSolver::setUpperBounds(Vect<real_t>& ub)
{
   _ub = ub;
}


void OptSolver::setLowerBound(real_t lb)
{
   _lb[0] = lb;
}


void OptSolver::setLowerBounds(Vect<real_t>& lb)
{
   _lb = lb;
}


real_t OptSolver::Objective(const Vect<real_t>& x)
{
   if (_exp) {
      int err;
      theParser.Parse(_exp_obj.c_str(),_var.c_str());
      real_t v = theParser.Eval(x);
      if ((err=theParser.EvalError()))
         throw OFELIException("In OptSolver::Objective(x): Illegal algebraic expression "+itos(err));
      return v;
   }
   else
      return _opt->Objective(x);
}


void OptSolver::Gradient(const Vect<real_t>& x,
                         Vect<real_t>&       g)
{
   if (_exp) {
      int err;
      for (size_t i=0; i<_size; i++) {
         theParser.Parse(_exp_grad[i].c_str(),_var.c_str());
         g[i] = theParser.Eval(x);
         if ((err=theParser.EvalError()))
            throw OFELIException("In OptSolver::Gradient(x,g): Illegal algebraic expression "+itos(err));
      }
   }
   else
      _opt->Gradient(x,g);
}


void OptSolver::setObjective(string exp)
{
   _exp_obj = exp;
   _exp_grad.setSize(_size);
   _exp_hess.setSize(_size,_size);
   _var = "";
   for (size_t i=1; i<_size; i++)
      _var += "x" + itos(i) + ",";
   _var += "x" + itos(_size);
}


void OptSolver::setGradient(string exp,
                            int    i)
{
   if (_opt_method==SIMULATED_ANNEALING)
      throw OFELIException("In OptSolver::setGradient(exp,i): Providing the gradient is useless for the simulated annealing method.");
   if (i>int(_size) || i<=0)
      throw OFELIException("In OptSolver::setGradient(exp,i): Index is out of bounds");
   _exp_grad[i-1] = exp;
}


void OptSolver::setTNDefaults()
{
   _verb = 1;
   _max_it = 200;
   _toler = OFELI_EPSMCH;
}


void OptSolver::setPGDefaults()
{
   _verb = 1;
   _max_it = 100;
   _toler = OFELI_EPSMCH;
}


void OptSolver::setNMDefaults()
{
   _verb = 1;
   _max_it = 100;
   _toler = OFELI_EPSMCH;
   _reqmin = 0.,
   _step.setSize(_size);
   _step = 0.1;
   _conv = 1;
   _max_eval = 100000;
   _nb_restart = 1;
}


void OptSolver::setSADefaults()
{
   _t = 1.;
   _rt = 1.5;
   _ns = 20;
   _nt = fmax(100,5*_size);
   _vm.setSize(_size);
   _vm = 1.;
   _neps = 4;
   _max_eval = 100000;
   _c.setSize(_size);
   _c = 2;
}


void OptSolver::setSAOpt(real_t        rt,
                         int           ns,
                         int           nt,
                         int&          neps,
                         int           maxevl,
                         real_t        t,
                         Vect<real_t>& vm,
                         Vect<real_t>& xopt,
                         real_t&       fopt)
{
   _rt = rt;
   _ns = ns;
   _nt = nt;
   _t = t;
   _vm = vm;
   _neps = neps;
   _max_eval = maxevl;
// *  @param c [in] Vector that controls the step length adjustment. The suggested
// *  value for all elements is 2.
   _c.setSize(_size);
   _c = 2;
}


void OptSolver::setBC(const Vect<real_t>& bc)
{
   _lb = -std::numeric_limits<real_t>::max();
   _ub =  std::numeric_limits<real_t>::max();
   mesh_nodes(bc.getMesh()) {
      for (size_t i=1; i<=The_node.getNbDOF(); i++) {
          size_t k = The_node.getDOF(i) - 1;
          if (The_node.getCode(i)>0)
             _lb[k] = _ub[k] = bc[k];
      }
   }
}


int OptSolver::run()
{
   if (_method_set==false)
      throw OFELIException("In OptSolver::run(): No optimization method has been chosen.");
   int ret = 0;
   if (_opt_method==TRUNCATED_NEWTON) {
      if (_verb>0)
         cout << "Solving the optimization problem by the Truncated Newton method ..." << endl;
      ret = OptimTN(*this,*_x,_lb,_ub,_nb_obj_eval,_nb_grad_eval,_max_it,_toler,_verb);
   }
   else if (_opt_method==SIMULATED_ANNEALING)
      ret = OptimSA(*this,*_x,_rt,_toler,_ns,_nt,_neps,_max_eval,_lb,_ub,_c,_verb,_t,
                    _vm,_fopt,_nacc,_nb_obj_eval,_nobds);
   else if (_opt_method==NELDER_MEAD)
      ret = OptimNM(*this,*_x,_fopt,_reqmin,_step,_conv,_max_eval,_nb_obj_eval,_nb_restart);
   else if (_opt_method==GRADIENT)
      ret = OptimPG(*this,*_x,_lb,_ub,_nb_obj_eval,_nb_grad_eval,_max_it,_toler,_verb);
   else
      ;
   return ret;
}


int OptSolver::run(real_t toler,
                   int    max_it,
                   int    verb)
{
   _toler = toler;
   _max_it = max_it;
   _verb = verb;
   return run();
}


ostream& operator<<(ostream&         s,
                    const OptSolver& os)
{
   string om = "Gradient";
   if (os._opt_method==OptSolver::TRUNCATED_NEWTON)
      om = "Truncated Newton";
   else if (os._opt_method==OptSolver::SIMULATED_ANNEALING)
      om = "Simulated Annealing";
   else if (os._opt_method==OptSolver::NELDER_MEAD)
      om = "Nelder Mead";
   else if (os._opt_method==OptSolver::GRADIENT)
      om = "Gradient";
   s << "\n\nSUMMARY OF OPTIMIZATION SOLVER:" << endl;
   s << "Optimization method:\t\t\t" << om << endl;
   s << "Problem size:\t\t\t\t" << setw(6) << os._size << endl;
   if (os._opt_method!=OptSolver::SIMULATED_ANNEALING &&
       os._opt_method!=OptSolver::NELDER_MEAD)
      s << "Number of performed iterations:\t\t" << setw(6) << os._max_it << endl;
   if (os._opt_method!=OptSolver::SIMULATED_ANNEALING &&
       os._opt_method!=OptSolver::NELDER_MEAD)
   s << "Tolerance for convergence testing:\t" << setw(6) << os._toler << endl;
   s << "Number of objective evaluations:\t" << setw(6) << os._nb_obj_eval << endl;
   s << "Number of gradient evaluations:\t\t" << setw(6) << os._nb_grad_eval << endl;
   return s;
}

} /* namespace OFELI */
