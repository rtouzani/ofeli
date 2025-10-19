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

                    Implementation of class 'NLASSolver'

  ==============================================================================*/

#include "solvers/NLASSolver.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/DMatrix_impl.h"

namespace OFELI {

NLASSolver::NLPtr NLASSolver::NL [] = {
   nullptr,                          // Not used
   &NLASSolver::solveBisection,      // Bisection method
   &NLASSolver::solveRegulaFalsi,    // Regula Falsi method
   &NLASSolver::solvePicard,         // Picard iteration method
   &NLASSolver::solveSecant,         // Secant method
   &NLASSolver::solveNewton          // Newton's method
};


NLASSolver::NLASSolver()
           : _fct_allocated(false), _df_computed(false), _cv(false), _f_given(false), _df_given(false),
             _u_set(true), _ab_given(false), _theEqua(nullptr), _nl(int(NEWTON)), _max_it(100),
             _fct_type(FUNCTION), _toler(1.e-8), _Df(nullptr), _theMesh(nullptr), _my_nlas(nullptr),
             _nb_eq(0), _nb_fct_def(0)
{
   _Df = nullptr;
   _u = new Vect<real_t>(1);
}


NLASSolver::NLASSolver(NonLinearIter nl,
                       size_t        nb_eq)
           : _fct_allocated(false), _df_computed(false), _cv(false), _f_given(false), _df_given(false),
             _u_set(false), _ab_given(false), _theEqua(nullptr), _max_it(100), _fct_type(FUNCTION), 
             _toler(1.e-8), _Df(nullptr), _theMesh(nullptr), _my_nlas(nullptr), _nb_eq(nb_eq), _nb_fct_def(0)
{
   set(nl);
}


NLASSolver::NLASSolver(real_t&       x,
                       NonLinearIter nl)
           : _fct_allocated(false), _df_computed(false), _cv(false), _f_given(false), _df_given(false), 
             _u_set(false), _ab_given(false), _theEqua(nullptr), _max_it(100), _fct_type(FUNCTION),
             _toler(1.e-8), _x(&x), _Df(nullptr), _theMesh(nullptr), _my_nlas(nullptr), _nb_eq(1),
             _nb_fct_def(0)
{
   set(nl);
}


NLASSolver::NLASSolver(Vect<real_t>& u,
                       NonLinearIter nl)
           : _fct_allocated(false), _df_computed(false), _cv(false), _f_given(false), _df_given(false),
             _u_set(false),_ab_given(false), _theEqua(nullptr), _nl(int(nl)), _max_it(100),
             _fct_type(FUNCTION), _u(&u), _toler(1.e-8), _Df(nullptr), _theMesh(nullptr),  _nb_eq(u.size()),
             _nb_fct_def(0)
{
   _my_nlas = nullptr;
   _v.setSize(_nb_eq);
   _w.setSize(_nb_eq);
   set(nl);
}


NLASSolver::NLASSolver(MyNLAS&       my_nlas,
                       NonLinearIter nl)
           : _fct_allocated(false), _cv(false), _f_given(false), _df_given(false), _u_set(false),
             _ab_given(false), _theEqua(nullptr), _nl(int(nl)), _max_it(100), _fct_type(FUNCTION),
             _toler(1.e-8), _Df(nullptr), _theMesh(nullptr), _nb_eq(0), _nb_fct_def(0)
{
   _my_nlas = &my_nlas;
   set(nl);
}


NLASSolver::~NLASSolver()
{
   if (_u_set)
      delete _u;
   if (_fct_allocated) {
      for (size_t i=0; i<_nb_eq; ++i) {
         delete _theFct[i];
         delete _theDFct[i];
      }
   }
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


void NLASSolver::setPDE(Equa& eq)
{
   _theEqua = &eq;
   _theMesh = &(_theEqua->getMesh());
   _Df = new DMatrix<real_t>(_nb_eq);
   _theEqua->getTangent(_Df);
   _fct_type = CLASS;
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
      throw OFELIException("In NLASSolver::setGradient(g):\nFunction must be given first.");
   _df_given = true;
   _grad1 = g;
}


void NLASSolver::setGradient(function<Vect<real_t>(Vect<real_t>)> g)
{
   if (_f_given==false)
      throw OFELIException("In NLASSolver::setGradient(g):\nFunction must be given first.");
   _df_given = true;
   _grad = g;
}


void NLASSolver::setf(Fct& f)
{
   if (_nb_fct_def+1>_nb_eq)
      throw OFELIException("In NLASSolver::setf(Fct):\nToo many function definitions.");
   if (_nb_fct_def==0) {
      _fct_type = EXPRESSION;
      _f_given = true;
      _theFct.resize(_nb_eq);
      _theDFct.resize(_nb_eq*_nb_eq);
      _df_computed = _df_given = true;
   }
   _theFct[_nb_fct_def++] = &f;
}


void NLASSolver::setDf(Fct&   df,
                       size_t i,
                       size_t j)
{
   if (_f_given==false)
      throw OFELIException("In NLASSolver::setDf(Fct,i,j):\nFunction must be given first.");
   if (_nl!=NEWTON && _nl!=SECANT)
      throw OFELIException("In NLASSolver::setDf(Fct,i,j):\n"
                           "Providing the gradient is useless for the chosen algorithm.");
   if (i==0 || j>_nb_eq || j>_nb_eq || j==0)
      throw OFELIException("In NLASSolver::setDf(Fct,i,j):\n"
                           "Index (" + to_string(i) + "," + to_string(j) + ") is out of bounds");
   _theDFct[_nb_eq*(i-1)+j-1] = &df;
   _df_computed = false, _df_given = true;
}


void NLASSolver::setf(string exp)
{
   if (_nb_fct_def+1>_nb_eq)
      throw OFELIException("In NLASSolver::setf(string):\nToo many function definitions.");
   if (_nb_fct_def==0) {
      _fct_type = EXPRESSION;
      _f_given = true;
      _theFct.resize(_nb_eq);
      _theDFct.resize(_nb_eq*_nb_eq);
      _df_computed = _df_given = true;
   }
   vector<string> var(_nb_eq);
   if (_nb_eq==1)
      var[0] = "x";
   else {
      for (size_t j=0; j<_nb_eq; ++j)
         var[j] = "x" + to_string(j+1);
   }
   _theFct[_nb_fct_def++] = new Fct(exp,var);
   _fct_allocated = true;
}


void NLASSolver::setDf(string exp,
                       size_t i,
                       size_t j)
{
   if (_f_given==false)
      throw OFELIException("In NLASSolver::setDf(exp,i,j):\nFunction must be given first.");
   if (_nl!=NEWTON && _nl!=SECANT)
      throw OFELIException("In NLASSolver::setDf(exp,i,j):\n"
                           "Providing the gradient is useless for the chosen algorithm.");
   if (i==0 || j>_nb_eq || j>_nb_eq || j==0)
      throw OFELIException("In NLASSolver::setDf(exp,i,j):\n"
                           "Index (" + to_string(i) + "," + to_string(j) + ") is out of bounds");
   vector<string> var(_nb_eq);
   if (_nb_eq==1)
      var[0] = "x";
   else {
      for (size_t j=0; j<_nb_eq; ++j)
         var[j] = "x" + to_string(j+1);
   }
   _theDFct[_nb_eq*(i-1)+j-1] = new Fct(exp,var);
   _fct_allocated = true;
   _df_computed = false, _df_given = true;
}


void NLASSolver::set(NonLinearIter nl)
{
   _nl_it = int(nl);
   if (_nl_it<0 || _nl_it>6) {
      _nlp = nullptr;
      throw OFELIException("In NLASSolver::set(NonLinearIter):\n"
                           "Iterative method not available.");
   }
   else
      _nlp = NL[_nl_it+1];
}


void NLASSolver::setNbEq(size_t nb_eq)
{
   _nb_eq = nb_eq;
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
   _w.setSize(_nb_eq);
}


real_t NLASSolver::Function(real_t x)
{
   if (_fct_type==EXPRESSION)
      return (*_theFct[0])(x);
   else if (_fct_type==FUNCTION)
      return _fct1(x);
   else
      return 0.;
}


real_t NLASSolver::Function(const Vect<real_t>& u,
                            int                 i)
{
   if (_fct_type==EXPRESSION)
      return (*_theFct[i-1])(u);
   else if (_fct_type==FUNCTION)
      return _fct(u)(i);
   else
      return _my_nlas->Function(u,i);
}


real_t NLASSolver::Gradient(real_t x)
{
   if (_fct_type==EXPRESSION) {
      if (_df_computed)
         return _theFct[0]->D(x);
      else
         return (*_theDFct[0])(x);
   }
   else if (_fct_type==FUNCTION)
      return _grad1(x);
   else
      return 0.;
}


void NLASSolver::Gradient(const Vect<real_t>& x,
                          size_t              i,
                          size_t              j)
{
   if (_fct_type==EXPRESSION) {
      if (_df_computed)
         (*_Df)(i,j) = _theFct[i-1]->D(x,j);
      else
         (*_Df)(i,j) = (*_theDFct[_nb_eq*(i-1)+j-1])(x);
   }
   else if (_fct_type==FUNCTION)
      (*_Df)(i,j) = _grad(x)(i,j);
   else
      (*_Df)(i,j) = _my_nlas->Gradient(x,i,j);
   _df_given = true;
}


int NLASSolver::run()
{
   if (Verbosity > 0)
      cout << "Running the nonlinear solver ... " << endl;
   if (Verbosity > 1)
      cout << "Running iterations ..." << endl;
   (this->*_nlp)();
   return 0;
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
      if (Verbosity>1)
         cout << "Iteration " << _it+1 << ", Solution: " << _y << endl;
      real_t xx = fabs(*_x);
      if (xx<OFELI_EPSMCH)
         xx = 1.;
      if (fabs(*_x-_y)/xx < _toler) {
         _cv = true;
         _nb_it = _it, _it = _max_it;
      }
      *_x = _y;
   }
   if (Verbosity>0) {
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
      if (Verbosity>1)
         cout << "Iteration " << _it+1 << ", Solution: " << _y << endl;
      real_t xx = fabs(*_x);
      if (xx<OFELI_EPSMCH)
         xx = 1.;
      if (fabs(*_x-_y)/xx < _toler) {
         _cv = true;
         _nb_it = _it, _it = _max_it;
      }
      *_x = _y;
   }
   if (Verbosity>0) {
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
   if (!_df_given && !_df_computed)
      throw OFELIException("In NLASSolver::solveSecant():\n"
                           "The Secant method requires providing "
                           "the gradient of the given function.");
   _cv = false;
   _it = 0;
   if (_nb_eq==1) {
      if (Verbosity>2)
         cout << "Initial guess: " << *_x << endl;
      while (++_it < _max_it) {
         _g = Gradient(*_x);
         if (fabs(_g)<OFELI_EPSMCH) {
            if (_it==1)
               _g = 1.;
            else
               throw OFELIException("In NLASSolver::solveSecant(): Null derivative.");
         }
         _y = *_x - Function(*_x)/_g;
         if (Verbosity>1)
            cout << "Iteration " << _it+1 << ", Solution: : " << _y << endl;
         real_t xx = fabs(*_x);
         if (xx<OFELI_EPSMCH)
            xx = 1.;
         if (fabs(*_x-_y)/xx < _toler) {
            _nb_it = _it, _it = _max_it;
            _cv = true;
         }
         *_x = _y;
      }
   }
   else {
      Vect<real_t> b(_nb_eq);
      _Df = new DMatrix<real_t>(_nb_eq);
      for (size_t i=1; i<=_nb_eq; ++i) {
         b(i) = -Function(*_u,i);
         for (size_t j=1; j<=_nb_eq; j++)
            Gradient(*_u,i,j);
      }
      while (++_it < _max_it) {
         for (size_t i=1; i<=_nb_eq; ++i) {
            b(i) = -Function(*_u,i);
            for (size_t j=1; j<=_nb_eq; j++)
               Gradient(*_u,i,j);
         }
         _Df->solve(b,_v);
         if (_v.getWNorm2()/_u->getWNorm2() < _toler) {
            _nb_it = _it, _it = _max_it;
            _cv = true;
         }
         *_u += _v;
         if (Verbosity>1)
            cout << "Iteration " << _it+1 << ", Solution:\n" << *_u << endl;
      }
      delete _Df;
   }
   if (Verbosity>0) {
      if (_cv)
         cout << "Convergence after " << _nb_it << " iterations." << endl;
      else
         cout << "No convergence after " << _nb_it << " iterations." << endl;
   }
}


void NLASSolver::solveNewton()
{
   if (!_df_given && !_df_computed)
      throw OFELIException("In NLASSolver::solveNewton():\n"
                           "The Newton's method requires providing "
                           "the gradient of the given function.");
   _it = 0;
   _cv = false;
   if (_nb_eq==1) {
      if (Verbosity>2)
         cout << "Initial guess: " << *_x << endl;
      while (++_it < _max_it) {
         _g = Gradient(*_x);
         if (fabs(_g)<OFELI_EPSMCH) {
            if (_it==1)
               _g = 1.;
            else
               throw OFELIException("In NLASSolver::solveNewton(): Null derivative.");
         }
         _y = *_x - Function(*_x)/_g;
         if (Verbosity>1)
            cout << "Iteration " << _it+1 << ", Solution: " << _y << endl;
         real_t xx = fabs(*_x);
         if (xx<OFELI_EPSMCH)
            xx = 1.;
         if (fabs(*_x-_y)/xx < _toler) {
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
         for (size_t i=1; i<=_nb_eq; ++i) {
            b(i) = -Function(*_u,i);
            for (size_t j=1; j<=_nb_eq; j++)
               Gradient(*_u,i,j);
         }
         _Df->solve(b,_v);
         if (_v.getWNorm2()/_u->getWNorm2() < _toler) {
            _nb_it = _it, _it = _max_it;
            _cv = true;
         }
         *_u += _v;
         if (Verbosity>1)
            cout << "Iteration " << _it+1 << ", Solution:\n" << *_u << endl;
      }
      delete _Df;
   }
   if (Verbosity>0) {
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
