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

                       Implementation of class 'OptSolver'

  ==============================================================================*/

#include "solvers/OptSolver.h"
#include "solvers/Optim.h"
#include "linear_algebra/Vect_impl.h"
#include <limits>
#include <algorithm>
#include "OFELIException.h"

using std::to_string;

namespace OFELI {

OptSolver::OptSolver()
          : _size(1), _nb_in_const(0), _nb_eq_const(0), _nb_obj_eval(0), _nb_grad_eval(0),
            _max_eval(100000), _max_it(1000), _toler(1.e-10), _sa_opt(false), _tn_opt(false),
            _obj_type(0), _x_set(true), _method_set(false),
            _fct_allocated(false), _grad_allocated(false), _hessian_allocated(false)
{
   _x = new Vect<real_t>(1);
   _lb.setSize(_size);
   _ub.setSize(_size);
   _lb = -std::numeric_limits<real_t>::max();
   _ub =  std::numeric_limits<real_t>::max();
   _var.push_back("x");
   _theDFct.resize(1);
   _theDDFct.resize(1);
   _theFct = nullptr;
}


OptSolver::OptSolver(Vect<real_t>& x)
          : _size(x.size()), _nb_in_const(0), _nb_eq_const(0), _nb_obj_eval(0), _nb_grad_eval(0),
            _max_eval(100000), _max_it(1000), _x(&x), _toler(1.e-10), _sa_opt(false), _tn_opt(false),
            _obj_type(0), _x_set(false), _method_set(false),
            _fct_allocated(false), _grad_allocated(false), _hessian_allocated(false)
{
   set(x);
   if (_size==1)
      _var.push_back("x");
   else {
      for (size_t j=0; j<_size; ++j)
         _var.push_back("x"+to_string(j+1));
   }
   _theDFct.resize(_size);
   _theDDFct.resize(_size*_size);
   _theFct = nullptr;
}


OptSolver::OptSolver(MyOpt&        opt,
                     Vect<real_t>& x)
          : _size(x.size()), _nb_in_const(0), _nb_eq_const(0), _nb_obj_eval(0), _nb_grad_eval(0),
            _max_eval(100000), _max_it(1000), _x(&x), _toler(1.e-10), _opt(&opt), _sa_opt(false),
            _tn_opt(false), _obj_type(0), _x_set(false), _method_set(false),
            _fct_allocated(false), _grad_allocated(false), _hessian_allocated(false)
{
   set(x);
   if (_size==1)
      _var.push_back("x");
   else {
      for (size_t j=0; j<_size; ++j)
         _var.push_back("x"+to_string(j+1));
   }
   _theDFct.resize(_size);
   _theDDFct.resize(_size*_size);
   _theFct = nullptr;
}


OptSolver::~OptSolver()
{
   if (_x_set)
      delete _x;
   if (_fct_allocated)
      delete _theFct;
   if (_grad_allocated) {
      for (size_t i=0; i<_size; ++i)
         delete _theDFct[i];
   }
   if (_hessian_allocated) {
      for (size_t i=0; i<_size*_size; ++i)
         delete _theDDFct[i];
   }
}


void OptSolver::set(Vect<real_t>& x)
{
   if (_x_set)
      delete _x;
   _x = &x;
   _x_set = false;
   _size = _x->size();
   _lb.setSize(_size);
   _ub.setSize(_size);
   _lb = -std::numeric_limits<real_t>::max();
   _ub =  std::numeric_limits<real_t>::max();
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


real_t OptSolver::Objective(Vect<real_t>& x)
{
   if (_type==EXPRESSION || _type==FCT)
      return eval(x,_theFct);
   else
      return _opt->Objective(x);
}


void OptSolver::Gradient(Vect<real_t>& x,
                         Vect<real_t>& g)
{
   if ((_type==FCT || _type==EXPRESSION) && _grad_computed) {
      _xv = x;
      for (size_t i=0; i<_size; i++)
         g[i] = _theFct->D(_xv,i+1);
   }
   else if ((_type==FCT || _type==EXPRESSION) && !_grad_computed) {
      for (size_t i=0; i<_size; i++)
         g[i] = eval(x,_theDFct[i]);
   }
   else
      _opt->Gradient(x,g);
}


void OptSolver::setObjective(string exp)
{
   _theFct = new Fct(exp,_var);
   _fct.set(exp,_var);
   _type = EXPRESSION;
   _fct_allocated = true;
   _grad_computed = true;
}


void OptSolver::setObjective(Fct& f)
{
   _theFct = &f;
   _fct = f;
   _type = FCT;
   _grad_computed = true;
}


void OptSolver::setGradient(string exp,
                            int    i)
{
   if (_opt_method==SIMULATED_ANNEALING)
      throw OFELIException("In OptSolver::setGradient(exp,i): Providing the gradient "
			   "is useless for the simulated annealing method.");
   if (i>int(_size) || i<=0)
      throw OFELIException("In OptSolver::setGradient(exp,i): Index is out of bounds");
   _theDFct[i-1] = new Fct(exp,_var);
   _grad_allocated = true;
   _grad_computed = false;
}


void OptSolver::setHessian(string exp,
                           int    i,
                           int    j)
{
   if (_opt_method!=NEWTON)
      throw OFELIException("In OptSolver::setHessian(exp,i,j): Providing the hessian "
			   "is available for Newton's method only.");
   if (i>int(_size) || i<=0)
      throw OFELIException("In OptSolver::setHessian(exp,i,j): First index is out of bounds");
   if (j>int(_size) || j<=0)
      throw OFELIException("In OptSolver::setHessian(exp,i,j): Second index is out of bounds");
   _theDDFct[_size*(i-1)+j-1] = new Fct(exp,_var);
   _hessian_allocated = true;
   _hessian_computed = false;
}


void OptSolver::setGradient(Fct& f,
                            int  i)
{
   if (_opt_method==SIMULATED_ANNEALING)
      throw OFELIException("In OptSolver::setGradient(Fct,i): Providing the gradient is useless for the simulated annealing method.");
   if (i>int(_size) || i<=0)
      throw OFELIException("In OptSolver::setGradient(Fct,i): Index is out of bounds");
   _theDFct[i-1] = &f;
   _grad_computed = false;
}


void OptSolver::setHessian(Fct&   f,
                           int    i,
                           int    j)
{
   if (_opt_method!=NEWTON)
      throw OFELIException("In OptSolver::setHessian(exp,i,j): Providing the hessian "
			   "is available for Newton's method only.");
   if (i>int(_size) || i<=0)
      throw OFELIException("In OptSolver::setHessian(exp,i,j): First index is out of bounds");
   if (j>int(_size) || j<=0)
      throw OFELIException("In OptSolver::setHessian(exp,i,j): Second index is out of bounds");
   _theDDFct[_size*(i-1)+j-1] = &f;
   _hessian_computed = false;
}


void OptSolver::setIneqConstraint(Fct&   f,
				  real_t penal)
{
   if (_theFct==nullptr)
      throw OFELIException("In OptSolver::setIneqConstraint(f,penal): setObjective must be called first");
   string s = _theFct->expr + " + " + to_string(0.5*penal) + "*max(" + f.expr + ",0.)^2";
   vector<string> var = _theFct->var;
   _theFct->set(s,var);
   _nb_in_const++;
}


void OptSolver::setEqConstraint(Fct&   f,
                                real_t penal)
{
   if (_theFct==nullptr)
      throw OFELIException("In OptSolver::setEqConstraint(f,penal): setObjective must be called first");
   string s = _theFct->expr + " + " + to_string(0.5*penal) + "*(" + f.expr + ")^2";
   vector<string> var = _theFct->var;
   _theFct->set(s,var);
   _nb_eq_const++;
}


void OptSolver::setIneqConstraint(string exp,
				  real_t penal)
{
   if (_theFct==nullptr)
      throw OFELIException("In OptSolver::setIneqConstraint(exp,penal): setObjective must be called first");
   string s = _theFct->expr + " + " + to_string(0.5*penal) + "*max(" + exp + ",0.)^2";
   vector<string> var = _theFct->var;
   _theFct->set(s,var);
   _nb_in_const++;
}


void OptSolver::setEqConstraint(string exp,
                                real_t penal)
{
   if (_theFct==nullptr)
      throw OFELIException("In OptSolver::setEqConstraint(exp,penal): setObjective must be called first");
   string s = _theFct->expr + " + " + to_string(0.5*penal) + "*(" + exp + ")^2";
   vector<string> var = _theFct->var;
   _theFct->set(s,var);
   _nb_eq_const++;
}


void OptSolver::setTNDefaults()
{
   _max_it = 200;
   _toler = OFELI_EPSMCH;
}


void OptSolver::setPGDefaults()
{
   _max_it = 100;
   _toler = OFELI_EPSMCH;
}


void OptSolver::setNMDefaults()
{
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
// c Vector that controls the step length adjustment. The suggested value for all elements is 2.
   _c.setSize(_size);
   _c = 2;
}


void OptSolver::setBC(const Vect<real_t>& bc)
{
   _lb = -std::numeric_limits<real_t>::max();
   _ub =  std::numeric_limits<real_t>::max();
   size_t k = 0;
   node_loop(&(bc.getMesh())) {
     for (size_t i=1; i<=The_node.getNbDOF(); i++, k++) {
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
      if (Verbosity>0)
         cout << "Solving the optimization problem by the Truncated Newton method ..." << endl;
      ret = OptimTN(*this,*_x,_lb,_ub,_nb_obj_eval,_nb_grad_eval,_max_it,_toler);
   }
   else if (_opt_method==SIMULATED_ANNEALING)
      ret = OptimSA(*this,*_x,_rt,_toler,_ns,_nt,_neps,_max_eval,_lb,_ub,_c,_t,
                    _vm,_fopt,_nacc,_nb_obj_eval,_nobds);
   else if (_opt_method==NELDER_MEAD)
      ret = OptimNM(*this,*_x,_fopt,_reqmin,_step,_conv,_max_eval,_nb_obj_eval,_nb_restart);
   else if (_opt_method==GRADIENT)
      ret = OptimPG(*this,*_x,_lb,_ub,_nb_obj_eval,_nb_grad_eval,_max_it,_toler);
   else
      ;
   return ret;
}


int OptSolver::run(real_t toler,
                   int    max_it)
{
   _toler = toler;
   _max_it = max_it;
   return run();
}


real_t OptSolver::getObjective()
{
   return _fct(*_x);
}


real_t OptSolver::eval(const Vect<real_t>& x,
                       Fct*                f)
{
   _xv = x;
   return (*f)(_xv);
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
