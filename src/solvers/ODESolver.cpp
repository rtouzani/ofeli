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

                       Implementation of class 'ODESolver'

  ==============================================================================*/

#include "solvers/ODESolver.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/DMatrix_impl.h"
#include "OFELIException.h"

using std::to_string;
using std::cout;

namespace OFELI {

ODESolver::TSPtr ODESolver::TS [] = {
   nullptr,                              // Stationary problem (not used)
   &ODESolver::solveForwardEuler,        // Forward Euler
   &ODESolver::solveBackwardEuler,       // Backward Euler
   &ODESolver::solveCrankNicolson,       // Crank-Nicolson
   &ODESolver::solveHeun,                // Heun
   &ODESolver::solveNewmark,             // Newmark
   &ODESolver::solveLeapFrog,            // Leap Frog
   &ODESolver::solveAB2,                 // Adams Bashforth 2
   &ODESolver::solveRK4,                 // Runge Kutta 4
   &ODESolver::solveRK3_TVD,             // Runge Kutta 3 TVD
   &ODESolver::solveBDF2                 // Backward Differentiation Formula
};


ODESolver::FPtr ODESolver::F [] = {
   nullptr,
   &ODESolver::setF_ForwardEuler,
   &ODESolver::setF_BackwardEuler,
   &ODESolver::setF_CrankNicolson,
   &ODESolver::setF_Heun,
   &ODESolver::setF_Newmark,
   &ODESolver::setF_LeapFrog,
   &ODESolver::setF_AB2,
   &ODESolver::setF_RK4,
   &ODESolver::setF_RK3_TVD,
   &ODESolver::setF_BDF2
};


ODESolver::ODESolver()
          : _order(0), _nb_eq(1), _step(0), _s(DIRECT_SOLVER), _p(DIAG_PREC), _my_ode(nullptr),
            _fct_allocated(false), _dF_computed(false), _setF_called(false), _RK4_rhs(false),
            _linear(false), _a0(false), _a1(false), _a2(false), _constant_matrix(false), _regex(false),
            _explicit(false), _init(false), _lhs(false), _rhs(false), _time(0.), _beta(0.25), _gamma(0.5),
            _nb_fct_def(0), _rhs_count(0)
{
}


ODESolver::ODESolver(size_t nb_eq)
          : _order(0), _nb_eq(nb_eq), _step(0), _s(DIRECT_SOLVER), _p(DIAG_PREC), _my_ode(nullptr),
            _fct_allocated(false), _dF_computed(false), _setF_called(false), _RK4_rhs(false), _linear(false),
            _a0(false), _a1(false), _a2(false), _constant_matrix(false), _regex(false), _explicit(false),
            _init(false), _lhs(false), _rhs(false), _time(0.), _beta(0.25), _gamma(0.5), _nb_fct_def(0),
	    _rhs_count(0)
{
   setNbEq(_nb_eq);
}


ODESolver::ODESolver(TimeScheme s,
                     real_t     time_step,
                     real_t     final_time,
                     size_t     nb_eq)
          : _order(0), _nb_eq(nb_eq), _step(0), _s(DIRECT_SOLVER), _p(DIAG_PREC), _my_ode(nullptr),
            _fct_allocated(false), _dF_computed(false), _setF_called(false), _RK4_rhs(false), _linear(false),
            _a0(false), _a1(false), _a2(false), _constant_matrix(false), _regex(false), _explicit(false),
            _init(false), _lhs(false), _rhs(false), _time_step(time_step), _time(0.), _final_time(final_time),
            _beta(0.25), _gamma(0.5), _nb_fct_def(0), _rhs_count(0)
{
   set(s,time_step,final_time);
   setNbEq(_nb_eq);
}


ODESolver::ODESolver(MyODE&     my_ode,
                     TimeScheme s,
                     real_t     time_step,
                     real_t     final_time,
                     size_t     nb_eq)
          : _order(0), _nb_eq(nb_eq), _step(0), _s(DIRECT_SOLVER), _p(DIAG_PREC), _my_ode(&my_ode),
            _fct_allocated(false), _dF_computed(false), _setF_called(false), _RK4_rhs(false), _linear(false),
            _a0(false), _a1(false), _a2(false), _constant_matrix(false), _regex(false), _explicit(false),
            _init(false), _lhs(false), _rhs(true), _time_step(time_step), _time(0.), _final_time(final_time),
            _beta(0.25), _gamma(0.5), _nb_fct_def(0), _rhs_count(0)
{
   set(s,time_step,final_time);
   setNbEq(_nb_eq);
   _type = SCALAR_NL;
   if (nb_eq>1)
      _type = VECTOR_NL;
}


ODESolver::~ODESolver()
{
   if (_fct_allocated) {
      for (size_t i=0; i<_nb_eq; ++i) {
         delete _theF[i];
         delete _theDF[i];
         delete _thedFdt[i];
      }
   }
}


void ODESolver::setNbEq(size_t n)
{
   _nb_eq = n;
   _x.setSize(_nb_eq);
   _u.setSize(_nb_eq);
   _v.setSize(_nb_eq);
   _b.setSize(_nb_eq);
   _f0.setSize(_nb_eq);
   _f1.setSize(_nb_eq);
   _f01.setSize(_nb_eq);
   _f2.setSize(_nb_eq);
   _vF1.setSize(_nb_eq);
   _vF2.setSize(_nb_eq);
   _ddu.setSize(_nb_eq);
   _dudt.setSize(_nb_eq);
   _theF.resize(_nb_eq);
   _theDF.resize(_nb_eq*_nb_eq);
   _thedFdt.resize(_nb_eq);
   _f1.clear();
   _f2.clear();
   _f01.clear();
   _type = SCALAR_NL;
   if (_nb_eq>1)
      _type = VECTOR_LINEAR;
   _xv.resize(_nb_eq+1);
   _var.push_back("t");
   if (_nb_eq==1)
      _var.push_back("y");
   else {
      for (size_t j=1; j<=_nb_eq; ++j)
         _var.push_back("y"+to_string(j));
   }
}


void ODESolver::setCoef(real_t a0,
                        real_t a1,
                        real_t a2,
                        real_t f)
{
   _regex = false;
   _type = SCALAR_LINEAR;
   _c0 = a0; _c1 = a1; _c2 = a2;
   _d2 = f;
   _order = 2;
   if (_c2==0)
      _order = 1;
   if (_c1==0 && _c2==0)
      _order = 0;
   _lhs = true;
}


void ODESolver::setCoef(string a0,
                        string a1,
                        string a2,
                        string f)
{
   if (_nb_eq>1)
      throw OFELIException("In ODESolver::setCoef(...): This function is available for scalar ode's only.");
   vector<string> var;
   var.push_back("t");
   _theC.resize(3);
   _theC[0].set(a0,var);
   _theC[1].set(a1,var);
   _theC[2].set(a2,var);
   _nb_eq = 1;
   _type = SCALAR_LINEAR;
   _regex = true;
   _lhs = true;
}


void ODESolver::setMatrices(DMatrix<real_t>& A0,
                            DMatrix<real_t>& A1)
{
   if (A0.getNbRows()!=_nb_eq || A1.getNbRows()!=_nb_eq)
      throw OFELIException("In ODESolver::setMatrices(...): Matrix size is "
                           "different from system size");
   _A0 = &A0;
   _A1 = &A1;
   _order = 1;
   _lhs = true;
}


void ODESolver::setMatrices(DMatrix<real_t>& A0,
                            DMatrix<real_t>& A1,
                            DMatrix<real_t>& A2)
{
   setMatrices(A0,A1);
   if (A2.getNbRows()!=_nb_eq)
      throw OFELIException("In ODESolver::setMatrices(...): Matrix size is "
                           "different from system size");
   _A2 = &A2;
   _order = 2;
   _lhs = true;
}


void ODESolver::setLinear()
{
   _linear = true;
   _type = SCALAR_LINEAR;
   if (_nb_eq>1)
      _type = VECTOR_LINEAR;
}


void ODESolver::setF(string f)
{
   if (_nb_fct_def+1>_nb_eq)
      throw OFELIException("In ODESolver::setF(string):\nToo many function definitions.");
   _dF_computed = true;
   _theF[_nb_fct_def++] = new Fct(f,_var);
   _fct_allocated = true;
   _regex = true;
   _lhs = _rhs = true;
   _type = SCALAR_NL;
   if (_nb_eq>1)
      _type = VECTOR_NL;
   _setF_called = true;
}


void ODESolver::setDF(string df,
                      int    i,
                      int    j)
{
   if (_setF_called==false)
      throw OFELIException("In ODESolver::setDF(string,int,int): Function setF must be called before this one.");
   _regex = true;
   _theDF[_nb_eq*(i-1)+j-1] = new Fct(df,_var);
   _lhs = _rhs = true;
}


void ODESolver::setdFdt(string df,
                        int    i)
{
   if (_setF_called==false)
      throw OFELIException("In ODESolver::setdFdt(string,int): Function setF must be called before this one.");
   _regex = true;
   _thedFdt[i-1] = new Fct(df,_var);
   _y0 = _u[0];
   _lhs = _rhs = true;
}


void ODESolver::setF(Fct& f)
{
   if (_nb_fct_def+1>_nb_eq)
      throw OFELIException("In ODESolver::setF(Fct):\nToo many function definitions.");
   _dF_computed = true;
   _theF[_nb_fct_def++] = &f;
   _fct_allocated = false;
   _regex = true;
   if (_nb_eq==1)
      _y0 = _u[0];
   _lhs = _rhs = true;
   _setF_called = true;
   _type = SCALAR_NL;
   if (_nb_eq>1)
      _type = VECTOR_NL;
}


void ODESolver::setDF(Fct& df,
                      int  i,
                      int  j)
{
   if (_setF_called==false)
      throw OFELIException("In ODESolver::setDF(Fct,int,int):\nFunction must be given first.");
   if (i<=0 || i>int(_nb_eq) || j>int(_nb_eq) || j<=0)
      throw OFELIException("In ODESolver::setDF(Fct,i,j):\n"
                           "Index (" + to_string(i) + "," + to_string(j) + ") is out of bounds");
   _theDF[_nb_eq*(i-1)+j-1] = &df;
   _dF_computed = false;
}


void ODESolver::setdFdt(Fct& df,
                        int  i)
{
   if (_setF_called==false)
      throw OFELIException("In ODESolver::setdFdt(Fct,int):\nFunction must be given first.");
   if (i<=0 || i>int(_nb_eq))
      throw OFELIException("In ODESolver::setdFdt(Fct,int,int):\nIndex " + to_string(i) + " is out of bounds");
   _thedFdt[i-1] = &df;
   _dF_computed = false;
}


void ODESolver::setRK4RHS(Vect<real_t>& f)
{
   _f01 = f;
   _RK4_rhs = true;
}


void ODESolver::set(TimeScheme s,
                    real_t     time_step,
                    real_t     final_time)
{
   _final_time = final_time;
   _time_step = _time_step0 = time_step;
   _nb_ssteps = 1;
   _explicit = false;
   if (s&FORWARD_EULER || s&HEUN || s&RK4 || s&LEAP_FROG || s&AB2)
      _explicit = true;
   if (s&HEUN)
      _nb_ssteps = 2;
   if (s&RK4)
      _nb_ssteps = 4;
   _sc = int(s);
   if (s<1 || s>10) {
      _solve = nullptr;
      throw OFELIException("In ODESolver::set(...): Time integration scheme not available.");
   }
   else
      _solve = TS[int(s)];
   _constant_matrix = false;
}


void ODESolver::setInitial(real_t u)
{
   _u[0] = _y0 = u;
   _nb_eq = 1;
   _init = true;
}


void ODESolver::setInitial(real_t u,
                           int    i)
{
   _u[i-1] = u;
   _order = 1;
   _init = true;
}


void ODESolver::setInitial(real_t u,
                           real_t v)
{
   _u[0] = _y0 = u;
   _dy1 = v;
   _nb_eq = 1;
   _init = true;
}


void ODESolver::setInitialRHS(real_t f)
{
   _d0 = _d1 = f;
}


void ODESolver::setInitial(Vect<real_t>& u)
{
   if (u.size()!=_nb_eq)
      throw OFELIException("In ODESolver::setInitial(Vect<real_t>): "
                           "Vector size is different from system size");
   _w = &u;
   _u = *_w;
   _order = 1;
   _init = true;
}


void ODESolver::setInitial(Vect<real_t>& u,
                           Vect<real_t>& v)
{
   setInitial(u);
   if (v.size()!=_nb_eq)
      throw OFELIException("In ODESolver::setInitial(Vect<real_t>,Vect<real_t>): "
                           "Vector size is different from system size");
   _du = &v;
   _order = 2;
   _init = true;
}


void ODESolver::setInitialRHS(Vect<real_t>& f)
{
   _f0 = _f1 = f;
}


void ODESolver::setRHS(real_t f)
{
   _d2 = f;
   _rhs = true;
}


void ODESolver::setRHS(string f)
{
   if (_time>0)
      _rhs_count = 0;
   if (_rhs_count+1>_nb_eq)
      throw OFELIException("In ODESolver::setRHS(string):\nToo many function definitions.");
   _theF[_rhs_count] = new Fct(f,"t");
   _fct_allocated = true;
   _f1(_rhs_count+1) = (*_theF[_rhs_count])(_time);
   _regex = true;
   _rhs = true;
   _type = SCALAR_LINEAR;
   if (_rhs_count>0)
      _type = VECTOR_LINEAR;
   _rhs_count++;
}


void ODESolver::setRHS(Vect<real_t>& f)
{
   if (f.size()!=_nb_eq)
      throw OFELIException("In ODESolver::setRHS(Vect<real_t>): Vector size is "
                           "different from system size");
   _f2 = f;
   _rhs = true;
}


void ODESolver::setNewmarkParameters(real_t beta,
                                     real_t gamma)
{
   _beta = beta;
   _gamma = gamma;
}


real_t ODESolver::getTimeDerivative(int i) const
{
   if (_nb_eq==1)
      return _dydt;
   else
      return _dudt[i-1];
}


void ODESolver::getTimeDerivative(Vect<real_t>& y) const
{
   y = _dudt;
}


real_t ODESolver::runOneTimeStep()
{
   if (_step==0) {
      if (_init==false)
         throw OFELIException("In ODESolver::runOneTimeStep(): No given initial condition.");
      if (_lhs==false && _rhs==false)
         throw OFELIException("In ODESolver::runOneTimeStep(): ODE not properly defined.");
   }
   _step++;
   _time = theTime;
   if (Verbosity>1)
      cout << "Running time step: " << _step << ", time: " << _time << " ..." << endl;
   if (_regex) {
      if (_type==SCALAR_LINEAR) {
         _c0 = _theC[0](_time);
         _c1 = _theC[1](_time);
         _c2 = _theC[2](_time);
         _d2 = (*_theF[0])(_time);
         _d01 = (*_theF[0])(_time-0.5*_time_step);
         _order = 2;
         if (_c2==0)
            _order = 1;
      }
      else if (_type==VECTOR_LINEAR || _type==VECTOR_NL) {
         if (_setF_called==false) {
            for (size_t i=0; i<_nb_eq; ++i) {
               _f1[i]  = (*_theF[i])(_time-_time_step);
               _f2[i]  = (*_theF[i])(_time);
               _f01[i] = (*_theF[i])(_time-0.5*_time_step);
            }
         }
      }
   }
   (this->*_solve)();
   theTimeStep = _time_step;
   return _time_step;
}


void ODESolver::run(bool opt)
{
   if (Verbosity>0)
      cout << "Running time integration ..." << endl;
   _constant_matrix = opt;
   TimeLoop
      runOneTimeStep();
}


real_t ODESolver::evalF(real_t t,
                        real_t y)
{
   real_t f = 0.;
   if (_regex) {
      _xv[0] = t, _xv[1] = y;
      f = (*_theF[0])(_xv);
   }
   else if (_my_ode!=nullptr)
      f = _my_ode->Function(t,y);
   return f;
}


real_t ODESolver::evalDF(real_t t,
                         real_t y)
{
   real_t f = 0.;
   if (_regex) {
      _xv[0] = t, _xv[1] = y;
      f = (*_theDF[0])(_xv);
   }
   else if (_my_ode!=nullptr)
      f = _my_ode->Jacobian(t,y);
   return f;
}


real_t ODESolver::evalF(real_t              t,
                        const Vect<real_t>& y,
                        size_t              i)
{
   real_t f = 0.;
   if (_regex) {
      _xv[0] = t;
      for (size_t j=0; j<_nb_eq; ++j)
         _xv[j+1] = y[j];
      f = (*_theF[i-1])(_xv);
   }
   else if (_my_ode!=nullptr)
      f = _my_ode->Function(t,y,i);
   return f;
}


real_t ODESolver::evalDF(real_t              t,
                         const Vect<real_t>& y,
                         size_t              i,
                         size_t              j)
{
   real_t f = 0.;
   if (_regex) {
      _xv[0] = t;
      for (size_t j=0; j<_nb_eq; ++j)
         _xv[j+1] = y[j];
      f = (*_theDF[i-1])(_xv);
   }
   else if (_my_ode!=nullptr)
      f = _my_ode->Jacobian(t,y,i);
   return f;
}


void ODESolver::getPhase(Vect<real_t>& x,
                         Vect<real_t>& v,
                         size_t        i)
{
   size_t n = _y.size();
   x.setSize(n);
   v.setSize(n);
   for (size_t j=0; j<n; ++j) {
      x[j] = _y[j][i-1];
      v[j] = _dy[j][i-1];
   }
}


void ODESolver::solveForwardEuler()
{
   if (_order==2)
      throw OFELIException("In ODESolver::solveForwardEuler(): Forward Euler scheme is "
                           "implemented for first order equations only.");
   if (_type==SCALAR_LINEAR) {
      _y1 = _y0 + _time_step*(_d1-_c0*_y0)/_c1;
      _y2 = _y0 = _y1;
      _d1 = _d2;
      _dydt = (_y1 - _y0)/_time_step;
      _v[0] = _y1, _dudt[0] = _dydt;
      _y.push_back(_v);
      _dy.push_back(_dudt);
      return;
   }
   if (_type==SCALAR_NL) {
      _d1 = evalF(_time-_time_step,_y0);
      _y1 = _y0 + _time_step*_d1;
      _dydt = (_y1-_y0)/_time_step;
      _y2 = _y0 = _y1;
      _d1 = _d2;
      _v[0] = _y1, _dudt[0] = _dydt;
      _y.push_back(_v);
      _dy.push_back(_dudt);
      return;
   }
   if (_type==VECTOR_NL) {
      for (size_t i=1; i<=_nb_eq; ++i)
         _vF1[i-1] = evalF(_time-_time_step,_u,i);
      _v = _u + _time_step*_vF1;
      _dudt = _vF1;
      *_w = _u = _v;
      _y.push_back(_v);
      _dy.push_back(_dudt);
      return;
   }
   _b = _time_step*_f1;
   _A0->MultAdd(-_time_step,_u,_b);
   _ls.setMatrix(_A1);
   _ls.setRHS(_b);
   _ls.setSolution(_v);
   _ls.setSolver(_s,_p);
   _ls.solve();
   _dudt = (1./_time_step)*(_v - _u);
   _v += _u;
   _y.push_back(_v);
   _dy.push_back(_dudt);
   *_w = _u = _v;
   _f1 = _f2;
}


void ODESolver::solveBackwardEuler()
{
   if (_order==2)
      throw OFELIException("In ODESolver::solveBackwardEuler(): Backward Euler scheme is "
                           "implemented for first order equations only.");
   if (_type==VECTOR_NL)
      throw OFELIException("In ODESolver::solveBackwardEuler(): Backward Euler scheme is "
                           "not implemented for a nonlinear system of ODEs.");
   if (_type==SCALAR_LINEAR) {
      _y1 = (_c1*_y0+_time_step*_d2)/(_c1+_time_step*_c0);
      _dydt = (_y1 - _y0)/_time_step;
      _y2 = _y0 = _y1;
      _v[0] = _y1, _dudt[0] = _dydt;
      _y.push_back(_v);
      _dy.push_back(_dudt);
      return;
   }
   if (_type==SCALAR_NL) {
      real_t y=_y0, f=0, df=0;
      Converged = false; theIteration = 0;
      IterationLoop {
         f = evalF(_time,y);
         df = evalDF(_time,y);
         _y1 = (_y0 + _time_step*(f-df*y))/(1. - _time_step*df);
         _y2 = _y0 = _y1;
         Converged = _iter.check(y,_y1);
      }
      _dydt = (_y1 - _y0)/_time_step;
      _v[0] = _y1, _dudt[0] = _dydt;
      _y.push_back(_v);
      _dy.push_back(_dudt);
      _y0 = _y1;
      _y1 = _y2;
      return;
   }
   _b = _f2;
   _A1->MultAdd(1./_time_step,_u,_b);
   _ls.setMatrix(_A0);
   _ls.setRHS(_b);
   _ls.setSolution(_v);
   _ls.setSolver(_s,_p);
   _A0->Axpy(1./_time_step,_A1);
   _ls.solve();
   _dudt = (1./_time_step)*(_v - _u);
   _y.push_back(_v);
   _dy.push_back(_dudt);
   *_w = _u = _v;
}


void ODESolver::solveCrankNicolson()
{
   if (_order==2)
      throw OFELIException("In ODESolver::solveCrankNicolson(): Crank-Nicolson scheme is "
                           "implemented for first order equations only.");
   if (_type==VECTOR_NL)
      throw OFELIException("In ODESolver::solveCrankNicolson(): Crank-Nicolson scheme is "
                           "not implemented for nonlinear ODEs.");
   if (_type==SCALAR_LINEAR) {
      _y1 = (_c1 - 0.5*_time_step*_c0)*_y0 + 0.5*_time_step*(_d1+_d2);
      _y1 /= _c1 + 0.5*_time_step*_c0;
      _dydt = (_y1-_y0)/_time_step;
      _y2 = _y0 = _y1;
      _d1 = _d2;
      _v[0] = _y1, _dudt[0] = _dydt;
      _y.push_back(_v);
      _dy.push_back(_dudt);
      return;
   }
   if (_type==SCALAR_NL) {
      real_t y=_y0, f0=0, f1=0, df=0;
      Converged = false; theIteration = 0;
      IterationLoop {
         f0 = evalF(_time-_time_step,_y0);
         f1 = evalF(_time,y);
         df = evalDF(_time,y);
         _y1 = (_y0 + 0.5*_time_step*(f0+f1-df*y))/(1. - 0.5*_time_step*df);
         Converged = _iter.check(y,_y1);
      }
      _v[0] = _y1, _dudt[0] = _dydt;
      _y.push_back(_v);
      _dy.push_back(_dudt);
      _y2 = _y0 = _y1;
      return;
   }
   _b = 0.5*(_f1+_f2);
   _ls.setMatrix(_A0);
   _ls.setRHS(_b);
   _ls.setSolution(_v);
   _A1->MultAdd(2./_time_step,_u,_b);
   _ls.setSolver(_s,_p);
   if (_step==1 || !_constant_matrix)
      _A0->Axpy(2./_time_step,_A1);
   _ls.solve();
   if (_constant_matrix)
      _ls.setNoFact();
   _v = 2.*_v - _u;
   _dudt = (1./_time_step)*(_v-_u);
   _y.push_back(_v);
   _dy.push_back(_dudt);
   *_w = _u = _v;
   _f1 = _f2;
}


void ODESolver::solveHeun()
{
   if (_order==2)
      throw OFELIException("In ODESolver::solveHeun(): Heun scheme is implemented "
                           "for first order equations only.");
   if (_type==SCALAR_LINEAR) {
      real_t yz = _y0 + _time_step*(_d1-_c0*_y0)/_c1;
      _y1 = _y0 + 0.5*_time_step*(_d1+_d2-_c0*(_y0+yz))/_c1;
      _dydt = (_y1-_y0)/_time_step;
      _y2 = _y0 = _y1;
      _d1 = _d2;
      _v[0] = _y1, _dudt[0] = _dydt;
      _y.push_back(_v);
      _dy.push_back(_dudt);
      return;
   }
   if (_type==SCALAR_NL) {
      _d1 = evalF(_time-_time_step,_y0);
      _d2 = evalF(_time,_y0+_time_step*_d1);
      _y1 = _y0 + 0.5*_time_step*(_d1+_d2);
      _dydt = (_y1-_y0)/_time_step;
      _y2 = _y0 = _y1;
      _d1 = _d2;
      _v[0] = _y1, _dudt[0] = _dydt;
      _y.push_back(_v);
      _dy.push_back(_dudt);
      return;
   }
   if (_type==VECTOR_NL) {
         for (size_t i=0; i<_nb_eq; ++i)
            _vF1[i] = evalF(_time-_time_step,_u,i+1);
         for (size_t i=0; i<_nb_eq; ++i)
            _vF2[i] = evalF(_time-_time_step,_u+_time_step*_vF1,i+1);
      _v = _u + 0.5*_time_step*(_vF1+_vF2);
      _dudt = (1./_time_step)*(_v-_u);
      _y.push_back(_v);
      _dy.push_back(_dudt);
      *_w = _u = _v;
      return;
   }
   _ls.setMatrix(_A1);
   _ls.setRHS(_b);
   _ls.setSolver(_s,_p);
   _b = _f1;
   _A0->MultAdd(-1.,_u,_b);
   _ls.setSolution(_v);
   _ls.solve();
   if (_constant_matrix)
      _ls.setNoFact();
   _b = _f2;
   _A0->MultAdd(-1.,_u+_time_step*_v,_b);
   _ls.setSolution(*_w);
   _ls.solve();
   _u += 0.5*_time_step*((*_w)+_v);
   _dudt = (1./_time_step)*(_v-_u);
   _y.push_back(_v);
   _dy.push_back(_dudt);
   *_w = _v = _u;
   _f1 = _f2;
}


void ODESolver::solveLeapFrog()
{
   if (_order==2)
      throw OFELIException("In ODESolver::solveLeapFrog(): Leap Frog scheme is implemented "
                           "for first order equations only.");
   if (_type==SCALAR_LINEAR) {
      if (_step==1) {
         _y1 = _y0 + _time_step*(_d0 - _c0*_y0)/_c1;
         _d1 = _d2;
         _dydt = (_y1-_y0)/_time_step;
         _v[0] = _y1, _dudt[0] = _dydt;
         _y.push_back(_v);
         _dy.push_back(_dudt);
      }
      else {
         _y2 = _y0 + 2.*_time_step*(_d1 - _c0*_y1)/_c1;
         _dydt = (_y2-_y0)/_time_step;
         _y0 = _y1; _y1 = _y2;
         _d0 = _d1, _d1 = _d2;
         _v[0] = _y2, _dudt[0] = _dydt;
         _y.push_back(_v);
         _dy.push_back(_dudt);
      }
      return;
   }
   else if (_type==SCALAR_NL) {
      if (_step==1) {
         _y1 = _y0 + _time_step*evalF(_time-_time_step,_y0);
         _dydt = (_y1-_y0)/_time_step;
         _y.push_back(_y1);
         _dy.push_back(_dydt);
         _v[0] = _y1, _dudt[0] = _dydt;
         _y.push_back(_v);
         _dy.push_back(_dudt);
      }
      else {
         _y2 = _y0 + 2*_time_step*evalF(_time-_time_step,_y1);
         _dydt = (_y2-_y0)/_time_step;
         _y0 = _y1; _y1 = _y2;
         _y.push_back(_y2);
         _dy.push_back(_dydt);
         _v[0] = _y1, _dudt[0] = _dydt;
         _y.push_back(_v);
         _dy.push_back(_dudt);
      }
      return;
   }
   else if (_type==VECTOR_NL) {
      if (_step==1) {
         for (size_t i=0; i<_nb_eq; ++i)
            _v[i] = _u[i] + _time_step*evalF(_time-_time_step,_u,i+1);
         _dudt = (1./_time_step)*(_v-_u);
         _y.push_back(_v);
         _dy.push_back(_dudt);
      }
      else {
         for (size_t i=0; i<_nb_eq; i++)
            (*_w)[i] = _u[i] + 2*_time_step*evalF(_time-_time_step,_v,i+1);
         _dudt = (1./_time_step)*(*_w-_u);
         _y.push_back(_v);
         _dy.push_back(_dudt);
         _u = _v; _v = *_w;
      }
      return;
   }
   _ls.setMatrix(_A1);
   _ls.setRHS(_b);
   _ls.setSolver(_s,_p);
   if (_step==1) {
      _ls.setSolution(_v);
      _b = _time_step*_f1;
      _A0->MultAdd(-_time_step,_u,_b);
      _ls.solve();
      _v += _u;
      _dudt = (1./_time_step)*(_v-_u);
      _y.push_back(_v);
      _dy.push_back(_dudt);
      _f1 = _f2;
   }
   else {
      _b = 2.*_time_step*_f1;
      _A0->MultAdd(-2.*_time_step,_v,_b);
      _ls.setSolution(*_w);
      _ls.solve();
      *_w += _u;
      _dudt = (1./_time_step)*(*_w-_u);
      _y.push_back(*_w);
      _dy.push_back(_dudt);
      _u = _v; _v = *_w;
      _f0 = _f1; _f1 = _f2;
   }
}


void ODESolver::solveAB2()
{
   if (_order==2)
      throw OFELIException("In ODESolver::solveAB2(): Adams-Bashforth scheme is "
                           "implemented for first order equations only.");
   if (_type==SCALAR_LINEAR) {
      if (_step==1) {
         _y1 = _y0 + _time_step*(_d0-_c0*_y0)/_c1;
         _dydt = (_y1-_y0)/_time_step;
         _v[0] = _y1, _dudt[0] = _dydt;
         _y.push_back(_v);
         _dy.push_back(_dudt);
      }
      else {
         _y2 = _y1 - 0.5*_time_step*(_c0*(3*_y1-_y0) - (3*_d1-_d0))/_c1;
         _dydt = (_y2-_y0)/_time_step;
         _y0 = _y1; _y1 = _y2;
         _d0 = _d1; _d1 = _d2;
         _v[0] = _y2, _dudt[0] = _dydt;
         _y.push_back(_v);
         _dy.push_back(_dudt);
      }
      return;
   }
   else if (_type==SCALAR_NL) {
      if (_step==1) {
         _y1 = _y0 + _time_step*evalF(_time-_time_step,_y0);
         _dydt = (_y1-_y0)/_time_step;
         _v[0] = _y1, _dudt[0] = _dydt;
         _y.push_back(_v);
         _dy.push_back(_dudt);
      }
      else {
         _d1 = evalF(_time-2*_time_step,_y0);
         _d2 = evalF(_time-_time_step,_y1);
         _y2 = _y1 + 0.5*_time_step*(3*_d2-_d1);
         _dydt = (_y2-_y1)/_time_step;
         _y0 = _y1; _y1 = _y2;
         _v[0] = _y1, _dudt[0] = _dydt;
         _y.push_back(_v);
         _dy.push_back(_dudt);
      }
      return;
   }
   else if (_type==VECTOR_NL) {
      if (_step==1) {
         for (size_t i=0; i<_nb_eq; ++i)
            _v[i] = _u[i] + _time_step*evalF(_time-_time_step,_u,i+1);
         _dudt = (1./_time_step)*(_v-_u);
         _y.push_back(_v);
         _dy.push_back(_dudt);
      }
      else {
         Vect<real_t> d1(_nb_eq), d2(_nb_eq);
         for (size_t i=0; i<_nb_eq; ++i)
            d1[i] = evalF(_time-2*_time_step,_u,i+1);
         for (size_t i=0; i<_nb_eq; ++i)
            d2[i] = evalF(_time-_time_step,_v,i+1);
         *_w = _v + 0.5*_time_step*(3.*d2-d1);
         _dudt = (1./_time_step)*(*_w-_v);
         _y.push_back(*_w);
         _dy.push_back(_dudt);
         _u = _v; _v = *_w;
      }
      return;
   }
   _ls.setMatrix(_A1);
   _ls.setRHS(_b);
   _ls.setSolver(_s,_p);
   if (_step==1) {
      _b = _time_step*_f0;
      _A0->MultAdd(-_time_step,_u,_b);
      _ls.setSolution(_v);
      _ls.solve();
      if (_constant_matrix)
         _ls.setNoFact();
      _dudt = (1./_time_step)*(_v-_u);
      _y.push_back(_v);
      _dy.push_back(_dudt);
      _v += _u;
      _f1 = _f2;
   }
   else {
      _b = 0.5*_time_step*(3.*_f1-_f0);
      _A0->MultAdd(0.5*_time_step,_u-3.*_v,_b);
      _ls.setSolution(*_w);
      _ls.solve();
      _dudt = (1./_time_step)*(*_w-_v);
      *_w += _v;
      _y.push_back(*_w);
      _dy.push_back(_dudt);
      _u = _v; _v = *_w;
      _f0 = _f1; _f1 = _f2;
   }
}


void ODESolver::solveRK4()
{
   if (_order==2)
      throw OFELIException("In ODESolver::solveRK4(): Runge-Kutta scheme is valid "
                           "for first order equations only.");
   real_t k1, k2, k3, k4;
   if (_type==SCALAR_LINEAR) {
      k1 = _d1  - _c0*_y0;
      k2 = _d01 - _c0*(_y0 + 0.5*_time_step*k1); 
      k3 = _d01 - _c0*(_y0 + 0.5*_time_step*k2);
      k4 = _d2  - _c0*(_y0 +     _time_step*k3);
      _y2 = _y1 = _y0 + OFELI_SIXTH*_time_step*(k1+2*(k2+k3)+k4)/_c1;
      _dydt = (_y2-_y0)/_time_step;
      _y0 = _y1;
      _d0 = _d1 = _d2;
      _v[0] = _y2, _dudt[0] = _dydt;
      _y.push_back(_v);
      _dy.push_back(_dudt);
   return;
   }
   else if (_type==SCALAR_NL) {
      k1 = evalF(_time-_time_step,_y0);
      k2 = evalF(_time-0.5*_time_step,_y0+0.5*_time_step*k1);
      k3 = evalF(_time-0.5*_time_step,_y0+0.5*_time_step*k2);
      k4 = evalF(_time               ,_y0+    _time_step*k3);
      _y2 = _y1 = _y0 + OFELI_SIXTH*_time_step*(k1+2*(k2+k3)+k4);
      _dydt = (_y2-_y0)/_time_step;
      _y0 = _y1;
      _v[0] = _y2, _dudt[0] = _dydt;
      _y.push_back(_v);
      _dy.push_back(_dudt);
      return;
   }
   else if (_type==VECTOR_NL) {
      Vect<real_t> z1(_nb_eq), z2(_nb_eq), z3(_nb_eq), z4(_nb_eq);
      for (size_t i=0; i<_nb_eq; ++i)
         z1[i] = evalF(_time-_time_step,_u,i+1);
      for (size_t i=0; i<_nb_eq; ++i)
         z2[i] = evalF(_time-0.5*_time_step,_u+0.5*_time_step*z1,i+1);
      for (size_t i=0; i<_nb_eq; ++i)
         z3[i] = evalF(_time-0.5*_time_step,_u+0.5*_time_step*z2,i+1);
      for (size_t i=0; i<_nb_eq; ++i)
         z4[i] = evalF(_time               ,_u+    _time_step*z3,i+1);
      (*_w) = _v = _u + OFELI_SIXTH*_time_step*(z1+2.*(z2+z3)+z4);
      _dudt = (1./_time_step)*(_v-_u);
      _y.push_back(*_w);
      _dy.push_back(_dudt);
      _u = _v;
      return;
   }
   _k1.resize(_nb_eq), _k2.resize(_nb_eq);
   _k3.resize(_nb_eq), _k4.resize(_nb_eq);
   _k1 = _f1; _k4 = _f2;
   if (_RK4_rhs)
      _k3 = _k2 = _f01;
   else
      _k3 = _k2 = 0.5*(_f1 + _f2);
   _A0->MultAdd(-1.,_u,_k1);
   _A0->MultAdd(-1.,_u+0.5*_time_step*_k1,_k2);
   _A0->MultAdd(-1.,_u+0.5*_time_step*_k2,_k3);
   _A0->MultAdd(-1.,_u+_time_step*_k3,_k4);
   _b = OFELI_SIXTH*_time_step*(_k1+2.*(_k2+_k3)+_k4);
   _ls.setMatrix(_A1);
   _ls.setRHS(_b);
   _ls.setSolver(_s,_p);
   _ls.setSolution(_v);
   _ls.solve();
   _dudt = (1./_time_step)*(_v-_u);
   _y.push_back(_v);
   _dy.push_back(_dudt);
   _v += _u;
   *_w = _u = _v;
   _f0 = _f1; _f1 = _f2;
}


void ODESolver::solveRK3_TVD()
{
   if (_order==2)
      throw OFELIException("In ODESolver::solveRK3_TVD(): Runge-Kutta scheme is valid "
                           "for first order equations only.");
   if (_type==SCALAR_LINEAR) {
//      real_t t = _time - _time_step;
      real_t y1 = _y0 + _time_step*(_d1 - _c0*_y0);
//      t += _time_step;
      real_t y2 = 0.75*_y0 + 0.25*y1 + 0.25*_time_step*(_d01 - _c0*y1);
//      t -= 0.5*_time_step;
      _y2 = _y1 = OFELI_THIRD*(_y0 + 2*y2 + 2*_time_step*(_d2 - _c0*y2));
      _y0 = _y1;
      _d0 = _d1 = _d2;
      _dydt = (_y2-_y0)/_time_step;
      _v[0] = _y2, _dudt[0] = _dydt;
      _y.push_back(_v);
      _dy.push_back(_dudt);
      return;
   }
   else if (_type==SCALAR_NL) {
      real_t y1 = _y0 + _time_step*evalF(_time-_time_step,_y0);
      real_t y2 = 0.75*_y0 + 0.25*y1 + 0.25*_time_step*evalF(_time,y1);
      _y2 = _y1 = OFELI_THIRD*(_y0 + 2*y2 + 2*_time_step*evalF(_time-0.5*_time_step,y2));
      _dydt = (_y2-_y1)/_time_step;
      _y0 = _y1;
      _v[0] = _y2, _dudt[0] = _dydt;
      _y.push_back(_v);
      _dy.push_back(_dudt);
      return;
   }
   else if (_type==VECTOR_NL) {
      Vect<real_t> z1(_nb_eq), z2(_nb_eq);
      for (size_t i=0; i<_nb_eq; ++i)
         z1[i] = _u[i] + _time_step*evalF(_time-_time_step,_u,i+1);
      for (size_t i=0; i<_nb_eq; ++i)
         z2[i] = 0.75*_u[i] + 0.25*z1[i] + 0.25*_time_step*evalF(_time,z1,i+1);
      for (size_t i=0; i<_nb_eq; ++i)
         (*_w)[i] = _v[i] = OFELI_THIRD*(_u[i] + 2*z2[i] + 2*_time_step*evalF(_time-0.5*_time_step,z2,i+1));
      _dudt = (1./_time_step)*(_v-_u);
      _y.push_back(*_w);
      _dy.push_back(_dudt);
      _u = _v;
      return;
   }
   _k1.resize(_nb_eq), _k2.resize(_nb_eq);
   _k1 = _f1;
   if (_RK4_rhs)
      _k2 = _f01;
   else
      _k3 = _k2 = 0.5*(_f1 + _f2);
   _A0->MultAdd(-1.,_u,_k1);
   _A0->MultAdd(-1.,_u+0.5*_time_step*_k1,_k2);
   _A0->MultAdd(-1.,_u+0.5*_time_step*_k2,_k3);
   _b = OFELI_SIXTH*_time_step*(_k1+2.*_k2);
   _ls.setMatrix(_A1);
   _ls.setRHS(_b);
   _ls.setSolver(_s,_p);
   _ls.setSolution(_v);
   _ls.solve();
   _dudt = 0;
   _y.push_back(_v);
   _dy.push_back(_dudt);
   _v += _u;
   *_w = _u = _v;
   _f0 = _f1; _f1 = _f2;
}


void ODESolver::solveNewmark()
{
   if (_order==1)
      throw OFELIException("In ODESolver::solveNewmark(): Newmark scheme is "
                           "valid for second order equations only.");
   if (_type==SCALAR_NL || _type==VECTOR_NL)
      throw OFELIException("In ODESolver::solveNewmark(): Newmark scheme is "
                           "not implemented for a nonlinear ODE.");
   if (_type==SCALAR_LINEAR) {
      if (_step==1)
         _ddy = (_d0 - _c1*_dy1 - _c0*_y1)/_c2;
      _y1 += _time_step*(_dy1 + _time_step*(0.5-_beta)*_ddy);
      _dy2 = _dy1 + (1-_gamma)*_time_step*_ddy;
      _y2 = _beta*_time_step*_time_step*(_d2-_c1*_dy1) + (_c2+_gamma*_time_step*_c1)*_y2;
      _y2 /= _c2 + _gamma*_time_step*_c1 + _beta*_time_step*_time_step;
      _ddy = (_y2-_y1)/(_beta*_time_step*_time_step);
      _dydt = 0;
      _y1 = _y2; _dy1 = _dy2;
      return;
   }
   if (_step==1) {
      _b = _f1;
      _A1->MultAdd(-1.,*_du,_b);
      _A0->MultAdd(-1.,_u,_b);
      LinearSolver ls0;
      DMatrix<real_t> A2(*_A2);
      ls0.setMatrix(&A2);
      ls0.setSolution(_ddu);
      ls0.setRHS(_b);
      ls0.setSolver(_s,_p);
      ls0.solve();
      _dudt = 0.;
   }
   real_t a=1./(_beta*_time_step*_time_step), b=_gamma/(_beta*_time_step);
   _v = _u + _time_step*(*_du) + (_time_step*_time_step*(0.5-_beta))*_ddu;
   *_du += ((1.-_gamma)*_time_step)*_ddu;
   _b = _f2;
   for (size_t i=1; i<=_nb_eq; i++) {
      for (size_t j=1; j<=_nb_eq; j++) {
         real_t z = a*(*_A2)(i,j) + b*(*_A1)(i,j);
         (*_A0)(i,j) += z;
         _b(i) += z*_v(j) - (*_A1)(i,j)*(*_du)(j);
      }
   }
   _ls.setMatrix(_A0);
   _ls.setRHS(_b);
   _ls.setSolver(_s,_p);
   _ls.setSolution(*_w);
   _ls.solve();
   _ddu = a*(*_w-_v);
   *_du += (_time_step*_gamma)*_ddu;
   _dudt = 0.;
   _u = *_w;
   _f1 = _f2;
}


void ODESolver::solveBDF2()
{
   if (_order==2)
      throw OFELIException("In ODESolver::solveBDF2(): BDF2 scheme is "
                           "implemented for first order equations only.");
   if (_type==SCALAR_NL || _type==VECTOR_NL)
      throw OFELIException("In ODESolver::solveBDF2(): BDF2 scheme is "
                           "not implemented for a nonlinear ODE.");
   if (_type==SCALAR_LINEAR) {
      if (_step==1) {
         _y1 = (_time_step*_d2 + _c1*_y0)/(_c1+_time_step*_c0);
         _dydt = (_y1-_y0)/_time_step;
         _y.push_back(_y1);
         _dy.push_back(_dydt);
      }
      else {
         _y2 = (2*_time_step*_d2 + _c1*(4*_y1-_y0))/(3*_c1+2*_time_step*_c0);
         _dydt = 0.5*(3*_y2-4.*_y1+_y0)/_time_step;
         _y0 = _y1; _y1 = _y2;
         _y.push_back(_y2);
         _dy.push_back(_dydt);
      }
      return;
   }
   if (_type==SCALAR_NL) {
      if (_step==1) {
         real_t y=_y0, f=0, df=0;
         Converged = false; theIteration = 0;
         IterationLoop {
            f = evalF(_time,y);
            df = evalDF(_time,y);
            _y1 = (_y0 + _time_step*(f-df*y))/(1. - _time_step*df);
            Converged = _iter.check(y,_y1);
         }
         _dydt = (_y1-_y0)/_time_step;
         _y.push_back(_y1);
         _dy.push_back(_dydt);
      }
      else {
         real_t y=_y0, f=0, df=0;
         Converged = false; theIteration = 0;
         IterationLoop {
            f = evalF(_time,y);
            df = evalDF(_time,y);
            _y2 = (4*_y1-_y0+2*_time_step*(f-df*y))/(3-2*_time_step*df);
            Converged = _iter.check(y,_y1);
         }
         _dydt = 0.5*(3*_y2-4.*_y1+_y0)/_time_step;
         _y0 = _y1; _y1 = _y2;
         _y.push_back(_y1);
         _dy.push_back(_dydt);
      }
      return;
   }
   _ls.setMatrix(_A0);
   _b = _f2;
   _ls.setRHS(_b);
   _ls.setSolver(_s,_p);
   if (_step==1) {
      _A1->MultAdd(1./_time_step,_u,_b);
      _ls.setSolution(_v);
      _A0->Axpy(1./_time_step,_A1);
      _ls.solve();
      _dudt = (1./_time_step)*(_v-_u);
      *_w = _v;
   }
   else {
      _A1->MultAdd(1./_time_step,2.*_v-0.5*_u,_b);
      _A0->Axpy(1.5/_time_step,_A1);
      _ls.setSolution(*_w);
      _ls.solve();
      _dudt = (0.5/_time_step)*(3.*(*_w)-4.*_v+_u);
      _u = _v; _v = *_w;
   }
}


Vect<real_t>& ODESolver::setF_ForwardEuler()
{
   return _vF1;
}


Vect<real_t>& ODESolver::setF_BackwardEuler()
{
   return _vF2;
}


Vect<real_t>& ODESolver::setF_CrankNicolson()
{
   _vF = 0.5*(_vF1 + _vF2);
   return _f;
}


Vect<real_t>& ODESolver::setF_Heun()
{
   if (_sstep==1)
      return _vF1;
   else {
      _vF = 0.5*(_vF1 + _vF2);
      return _f;
   }
}


Vect<real_t>& ODESolver::setF_AB2()
{
   if (_step==1)
      return _vF1;
   else
      _vF = 1.5*_vF2 - 0.5*_vF1;
   return _f;
}


Vect<real_t>& ODESolver::setF_LeapFrog()
{
   return _vF1;
}


Vect<real_t>& ODESolver::setF_RK4()
{
   if (_sstep==1)
      return _f1;
   else if (_sstep==2)
      _vF = 0.5*(_vF1+_vF2);
   else if (_sstep==3)
      _vF = 0.5*(_vF1+_vF2);
   else if (_sstep==4)
      return _vF2;
   return _vF;
}


Vect<real_t>& ODESolver::setF_RK3_TVD()
{
   if (_sstep==1)
      return _f1;
   else if (_sstep==2)
      _vF = 0.5*(_vF1+_vF2);
   else if (_sstep==3)
      _vF = 0.5*(_vF1+_vF2);
   else if (_sstep==4)
      return _vF2;
   return _vF;
}


Vect<real_t>& ODESolver::setF_Newmark()
{
   if (_step==1 && _sstep==1)
      return _vF1;
   return _vF2;
}


Vect<real_t>& ODESolver::setF_BDF2()
{
   return _vF2;
}


ostream& operator<<(ostream&         s,
                    const ODESolver& ode)
{
   s << "\nNUMERICAL SOLUTION OF A DIFFERENTIAL SYSTEM\n\n";
   s << "Number of equations: \t\t" << ode._nb_eq << endl;
   string scheme;
   if (ode._sc==FORWARD_EULER)
      scheme = "Forward Euler";
   else if (ode._sc==BACKWARD_EULER)
      scheme = "Backward Euler";
   else if (ode._sc==CRANK_NICOLSON)
      scheme = "Crank-Nicolson";
   else if (ode._sc==HEUN)
      scheme = "Heun";
   else if (ode._sc==NEWMARK)
      scheme = "Newmark";
   else if (ode._sc==LEAP_FROG)
      scheme = "LeapFrog";
   else if (ode._sc==AB2)
      scheme = "Adams-Bashforth";
   else if (ode._sc==RK4)
      scheme = "Runge-Kutta 4";
   else if (ode._sc==RK3_TVD)
      scheme = "Runge-Kutta 3 TVD";
   else if (ode._sc==BDF2)
      scheme = "BDF 2";
   else
      ;
   s << "Time integration scheme: \t" << scheme << endl;
   s << "Final time: \t\t\t" << ode._final_time << endl;
   s << "First time step value:\t\t" << ode._time_step0 << endl;
   s << "Final time step value:\t\t" << ode._time_step << endl;
   s << "Number of time steps:\t\t" << ode._step << endl;
   return s;
}

} /* namespace OFELI */
