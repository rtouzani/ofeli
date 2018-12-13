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

                    Implementation of class 'NLASSolver'

  ==============================================================================*/

#include "solvers/NLASSolver.h"

namespace OFELI {

NLASSolver::NLPtr NLASSolver::NL [] = {
   NULL,                             // Not used
   &NLASSolver::solveBisection,      // Bisection method
   &NLASSolver::solveRegulaFalsi,    // Regula Falsi method
   &NLASSolver::solvePicard,         // Picard iteration method
   &NLASSolver::solveSecant,         // Secant method
   &NLASSolver::solveNewton          // Newton's method
};


NLASSolver::NLASSolver()
           : _cv(false), _f_given(false), _grad_given(false), _u_set(true),
             _ab_given(false), _theEqua(NULL), _nl(int(NEWTON)), _verb(1), _max_it(100),
             _nb_eq(0), _fct_type(FUNCTION), _toler(1.e-8), _var(""), _Df(NULL),
             _theMesh(NULL), _my_nlas(NULL)
{
   _Df = NULL;
   _u = new Vect<real_t>(1);
}


NLASSolver::NLASSolver(NonLinearIter nl,
                       int           nb_eq)
           : _cv(false), _f_given(false), _grad_given(false), _u_set(false),
             _ab_given(false), _theEqua(NULL), _nl(int(nl)), _verb(1), _max_it(100),
             _nb_eq(nb_eq), _fct_type(FUNCTION), _toler(1.e-8), _Df(NULL),
             _theMesh(NULL), _my_nlas(NULL)
{
}


NLASSolver::NLASSolver(real_t&       x,
                       NonLinearIter nl)
           : _cv(false), _f_given(false), _grad_given(false), _u_set(false),
             _ab_given(false), _theEqua(NULL), _nl(int(nl)), _verb(1), _max_it(100),
             _nb_eq(1), _fct_type(FUNCTION), _toler(1.e-8), _x(&x), _Df(NULL),
             _theMesh(NULL), _my_nlas(NULL)
{
   set(nl);
}


NLASSolver::NLASSolver(Vect<real_t>& u,
                       NonLinearIter nl)
           : _cv(false), _f_given(false), _grad_given(false), _u_set(false),
             _ab_given(false), _theEqua(NULL), _nl(int(nl)), _verb(1), _max_it(100),
             _nb_eq(u.size()), _fct_type(FUNCTION), _u(&u), _toler(1.e-8), _var(""),
             _Df(NULL), _theMesh(NULL), _my_nlas(NULL)
{
   _v.setSize(_nb_eq);
   set(nl);
}


NLASSolver::NLASSolver(MyNLAS&       my_nlas,
                       NonLinearIter nl)
           : _cv(false), _f_given(false), _grad_given(false), _u_set(false),
             _ab_given(false), _theEqua(NULL), _nl(int(nl)), _verb(1), _max_it(100),
             _nb_eq(0), _fct_type(FUNCTION), _toler(1.e-8), _var(""), _Df(NULL),
             _theMesh(NULL), _my_nlas(NULL)
{
   set(nl);
}


NLASSolver::~NLASSolver()
{
   if (_u_set)
      delete _u;
}


void NLASSolver::setInitial(real_t a,
                            real_t b)
{
   if (_nb_eq>1)
      throw OFELIException("In NLASSolver::setInitial(double,double):\n"
                           "This function is available for one-variable problems only.");
   _a = a;
   _b = b;
   _ab_given = true;
}


void NLASSolver::setPDE(AbsEqua<real_t>& eq)
{
   _theEqua = &eq;
   _theMesh = &(_theEqua->getMesh());
   _Df = new DMatrix<real_t>(_nb_eq);
   _theEqua->getTangent(_Df);
   _fct_type = CLASS;
}


void NLASSolver::eval(real_t x)
{
   int err;
   theParser.Parse(_f_exp[0].c_str(),"x");
   _f_exp[0] = theParser.Eval(x);
   if ((err=theParser.EvalError()))
      throw OFELIException("In NLASSolver::eval(double):\n"
                           "Illegal algebraic expression.");
}


void NLASSolver::setFunction(function<real_t(real_t)> f)
{
   _fct_type = FUNCTION;
   _f_given = true;
   _fct1 = f;
   _nb_eq = 1;
}


void NLASSolver::setFunction(function<Vect<real_t>(Vect<real_t>)> f)
{
   _fct_type = FUNCTION;
   _f_given = true;
   _fct = f;
}


void NLASSolver::setGradient(function<real_t(real_t)> g)
{
   if (_f_given==false)
      throw OFELIException("In NLASSolver::setGradient(g):\n"
                           "Function must be given first.");
   _grad_given = true;
   _grad1 = g;
}


void NLASSolver::setGradient(function<Vect<real_t>(Vect<real_t>)> g)
{
   if (_f_given==false)
      throw OFELIException("In NLASSolver::setGradient(g):\n"
                           "Function must be given first.");
   _grad_given = true;
   _grad = g;
}


void NLASSolver::setf(string exp)
{
   static int i=0;
   if (i==0) {
      _fct_type = EXPRESSION;
      _f_given = true;
      _f_exp.resize(_nb_eq);
      _Df_exp.setSize(_nb_eq,_nb_eq);
      _var = "";
      for (int j=1; j<_nb_eq; ++j)
         _var += "x" + itos(j) + ",";
      _var += "x" + itos(_nb_eq);
      if (_nb_eq==1)
         _var = "x";
   }
   _f_exp[i++] = exp;
}


void NLASSolver::setDf(string exp,
                       int    i,
                       int    j)
{
   if (_f_given==false)
      throw OFELIException("In NLASSolver::setDf(exp,i,j):\n"
                           "Function must be given first.");
   if (_nl!=NEWTON && _nl!=SECANT)
      throw OFELIException("In NLASSolver::setDf(exp,i,j):\n"
                           "Providing the gradient is useless for the chosen algorithm.");
   if (i<=0 || j>int(_nb_eq) || j>int(_nb_eq) || j<=0)
      throw OFELIException("In NLASSolver::setDf(exp,i):\n"
                           "Index (" + itos(i) + "," + itos(j) + ") is out of bounds");
   _Df_exp(i,j) = exp;
   _grad_given = true;
}

  
void NLASSolver::set(NonLinearIter nl)
{
   _nl_it = int(nl);
   if (_nl_it<0 || _nl_it>6) {
      _nlp = NULL;
      throw OFELIException("In NLASSolver::set(NonLinearIter):\n"
                           "Iterative method not available.");
   }
   else
      _nlp = NL[int(_nl_it)+1];
}


void NLASSolver::setInitial(real_t& x)
{
   _x = &x;
}


void NLASSolver::setInitial(Vect<real_t>& u)
{
   _u = &u;
   _nb_eq = _u->size();
   _v.setSize(_nb_eq);
}


real_t NLASSolver::Function(const Vect<real_t>& u,
                            int                 i)
{
   if (_fct_type==EXPRESSION) {
      int err;
      theParser.Parse(_f_exp[i-1].c_str(),_var.c_str());
      real_t v = theParser.Eval(u);
      if ((err=theParser.EvalError()))
         throw OFELIException("In NLASSolver::Function(u,i):\n"
                              "Illegal algebraic expression "+itos(err));
      return v;
   }

   else if (_fct_type==FUNCTION) {
      return _fct(u)(i);
   }

   else
      return _my_nlas->Function(u,i);
}


void NLASSolver::Gradient(const Vect<real_t>& u,
                          int                 i,
                          int                 j)
{
   if (_fct_type==EXPRESSION) {
      int err;
      theParser.Parse(_Df_exp(i,j).c_str(),_var.c_str());
      (*_Df)(i,j) = theParser.Eval(&u[0]);
      if ((err=theParser.EvalError()))
         throw OFELIException("In NLASSolver::Gradient(u,i):\n"
                              "Illegal algebraic expression "+itos(err));
   }

   else if (_fct_type==FUNCTION) {
     (*_Df)(i,j) = _grad(u)(i,j);
   }

   else
      (*_Df)(i,j) = _my_nlas->Gradient(u,i,j);
   _grad_given = true;
}


real_t NLASSolver::Function(real_t x)
{
   if (_fct_type==EXPRESSION) {
      int err;
      theParser.Parse(_f_exp[0].c_str(),_var.c_str());
      real_t v = theParser.Eval(x);
      if ((err=theParser.EvalError()))
         throw OFELIException("In NLASSolver::Function(x):\n"
                              "Illegal algebraic expression "+itos(err));
      return v;
   }

   if (_fct_type==FUNCTION) {
      return _fct1(x);
   }

   else
      return 0.;
}


real_t NLASSolver::Gradient(real_t x)
{
   if (_fct_type==EXPRESSION) {
      int err;
      theParser.Parse(_Df_exp[0].c_str(),_var.c_str());
      real_t g = theParser.Eval(x);
      if ((err=theParser.EvalError()))
         throw OFELIException("In NLASSolver::Gradient(x,g):\n"
                              "Illegal algebraic expression "+itos(err));
      _grad_given = true;
      return g;
   }

   else if (_fct_type==FUNCTION) {
      _grad_given = true;
      return _grad1(x);
   }

   else
      return 0.;
}


void NLASSolver::run()
{
   if (_verb > 0)
      cout << "Running the nonlinear solver ... " << endl;
   if (_verb > 1)
      cout << "Running iterations ..." << endl;
   (this->*_nlp)();
}


void NLASSolver::solveBisection()
{
   if (_nb_eq>1)
      throw OFELIException("In NLASSolver::solveBisection():\n"
                           "This method is valid for one variable problems only.");
   if (_ab_given==false)
      throw OFELIException("In NLASSolver::solveBisection():\n"
                           "Values of a and b must be given before.");
   if (Function(_a)*Function(_b)>0.)
      throw OFELIException("In NLASSolver::solveBisection():\n"
                           "Values of a and b must give opposite sign for f.");
   _it = 0;
   *_x = 0.5*(_a+_b);
   while (++_it < _max_it) {
      (Function(*_x)*Function(_a)<0.) ? _b = *_x : _a = *_x;
      _y = 0.5*(_a+_b);
      if (_verb>1)
         cout << "Iteration " << _it+1 << ", Solution: " << _y << endl;
      if (fabs(*_x-_y)/fabs(*_x) < _toler) {
         _cv = true;
         _nb_it = _it, _it = _max_it;
      }
      *_x = _y;
   }
   if (_verb>0) {
      if (_cv)
         cout << "Convergence after " << _nb_it << " iterations." << endl;
      else
         cout << "No convergence after " << _nb_it << " iterations." << endl;
   }
}


void NLASSolver::solveRegulaFalsi()
{
   if (_nb_eq>1)
      throw OFELIException("In NLASSolver::solveRegulaFalsi():\n"
                           "This method is valid for one variable problems only.");
   if (_ab_given==false)
      throw OFELIException("In NLASSolver::solveRegulaFalsi():\n"
                           "Values of a and b must be given before");
   if (Function(_a)*Function(_b)>0.)
      throw OFELIException("In NLASSolver::solveRegulaFalsi():\n"
                           "Values of a and b must give opposite sign for f.");
   real_t fa=Function(_a), fb=Function(_b);
   *_x = (fa*_b-fb*_a)/(fa-fb);
   _it = 0;
   while (++_it < _max_it) {
      (Function(*_x)*fa<0.) ? _b = *_x : _a = *_x;
      fa = Function(_a), fb = Function(_b);
      _y = (fa*_b-fb*_a)/(fa-fb);
      if (_verb>1)
         cout << "Iteration " << _it+1 << ", Solution: " << _y << endl;
      if (fabs(*_x-_y)/fabs(*_x) < _toler) {
         _cv = true;
         _nb_it = _it, _it = _max_it;
      }
      *_x = _y;
   }
   if (_verb>0) {
      if (_cv)
         cout << "Convergence after " << _nb_it << " iterations." << endl;
      else
         cout << "No convergence after " << _nb_it << " iterations." << endl;
   }
}


void NLASSolver::solvePicard()
{
}


void NLASSolver::solveSecant()
{
   if (_grad_given==false)
      throw OFELIException("In NLASSolver::solveSecant():\n"
                           "The Secant method requires providing "
                           "the gradient of the given function.");
   if (_nb_eq==1) {
      if (_verb>1)
         cout << "Initial guess: " << *_x << endl;
      _g = Gradient(*_x);
      _it = 0;
      while (++_it < _max_it) {
         _y = *_x - Function(*_x)/_g;
         if (_verb>1)
            cout << "Iteration " << _it+1 << ", Solution: : " << _y << endl;
         if (fabs(*_x-_y)/fabs(*_x) < _toler) {
            _nb_it = _it, _it = _max_it;
            _cv = true;
         }
         *_x = _y;
      }
   }
   else {
   }
   if (_verb>0) {
      if (_cv)
         cout << "Convergence after " << _nb_it << " iterations." << endl;
      else
         cout << "No convergence after " << _nb_it << " iterations." << endl;
   }
}


void NLASSolver::solveNewton()
{
   if (_grad_given==false)
      throw OFELIException("In NLASSolver::solveNewton():\n"
                           "The Newton's method requires providing "
                           "the gradient of the given function.");
   _it = 0;
   _cv = false;
   if (_nb_eq==1) {
      if (_verb>1)
         cout << "Initial guess: " << *_x << endl;
      while (++_it < _max_it) {
         _y = *_x - Function(*_x)/Gradient(*_x);
         if (_verb>1)
            cout << "Iteration " << _it+1 << ", Solution: " << _y << endl;
         if (fabs(*_x-_y)/fabs(*_x) < _toler) {
            _nb_it = _it, _it = _max_it;
            _cv = true;
         }
         *_x = _y;
      }
   }
   else {
      _Df = new DMatrix<real_t>(_nb_eq);
      Vect<real_t> b(_nb_eq);
      while (++_it < _max_it) {
         _Df->reset();
         for (int i=1; i<=_nb_eq; ++i) {
            b(i) = -Function(*_u,i);
            for (int j=1; j<=_nb_eq; j++)
               Gradient(*_u,i,j);
         }
         _Df->solve(b,_v);
         if (_v.getWNorm2()/_u->getWNorm2() < _toler) {
            _nb_it = _it, _it = _max_it;
            _cv = true;
         }
         *_u += _v;
         if (_verb>1)
            cout << "Iteration " << _it+1 << ", Solution:\n" << *_u << endl;
      }
      delete _Df;
   }
   if (_verb>0) {
      if (_cv)
         cout << "Convergence after " << _nb_it << " iterations." << endl;
      else
         cout << "No convergence after " << _nb_it << " iterations." << endl;
   }
}


void NLASSolver::get(Vect<real_t> &u) const
{
   u = *_u;
}


ostream& operator<<(ostream&          s,
                    const NLASSolver& nl)
{
   s << "\nNUMERICAL SOLUTION OF A NONLINEAR ALGEBRAIC SYSTEM\n\n";
   s << "Number of equations:                            " << nl._nb_eq << endl;
   string Nl;
   if (nl._nl==0)
      Nl = "Bisection";
   else if (nl._nl==1)
      Nl = "Regula Falsi";
   else if (nl._nl==2)
      Nl = "Picard iteration";
   else if (nl._nl==3)
      Nl = "Secant method";
   else if (nl._nl==4)
      Nl = "Newton";
   else
      ;
   s << "Nonlinear iteration method:                     " << Nl << endl;
   if (nl._cv)
      s << "Number of performed iterations for convergence: " << nl._nb_it << endl;
   else
      s << "NO CONVERGENCE !" << endl;
   return s;
}

} /* namespace OFELI */
