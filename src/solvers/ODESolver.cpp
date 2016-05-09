/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzani

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

                       Implementation of class 'DESolver'

  ==============================================================================*/

#include "solvers/ODESolver.h"
using std::cout;

namespace OFELI {

ODESolver::TSPtr ODESolver::TS [] = {
   NULL,                                 // Stationary problem (not used)
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
   NULL,
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


ODESolver::ODESolver() :
           _order(0), _nb_eq(1), _step(0), _verb(1),
           _s(DIRECT_SOLVER), _p(DIAG_PREC), _a0(false), _a1(false), _a2(false),
           _constant_matrix(false), _regex(false), _explicit(false), _init(false),
           _lhs(false), _rhs(false), _time(0.), _beta(0.25), _gamma(0.5),
           _alloc_f2(false), _alloc_f01(false)
{ }


ODESolver::ODESolver(TimeScheme s,
                     real_t     time_step,
                     real_t     final_time,
                     size_t     nb_eq) :
           _order(0), _nb_eq(nb_eq), _step(0), _verb(1),
           _s(DIRECT_SOLVER), _p(DIAG_PREC), _a0(false), _a1(false), _a2(false),
           _constant_matrix(false), _regex(false), _explicit(false), _init(false),
           _lhs(false), _rhs(false), _time_step(time_step), _time(0.),
           _final_time(final_time), _beta(0.25), _gamma(0.5), _alloc_f2(false),
           _alloc_f01(false)
{
   set(s,time_step,final_time);
   if (_nb_eq>1) {
      _u.setSize(_nb_eq);
      _v.setSize(_nb_eq);
      _b.setSize(_nb_eq);
      _f0.setSize(_nb_eq);
      _vF2.setSize(_nb_eq);
      _ddu.setSize(_nb_eq);
   }
}


ODESolver::~ODESolver()
{
   if (_alloc_f2)
      delete _f2;
  if (_alloc_f01)
     delete _f01;
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
   _exc[0] = a0;
   _exc[1] = a1;
   _exc[2] = a2;
   _expF.push_back(f);
   _nb_eq = 1;
   _type = SCALAR_LINEAR;
   _regex = true;
   _lhs = true;
}


void ODESolver::setMatrices(DMatrix<real_t>& A0,
                            DMatrix<real_t>& A1)
{
   try {
      if (A0.getNbRows()!=_nb_eq || A1.getNbRows()!=_nb_eq)
         THROW_RT("setMatrices(...): Matrix size is different from system size");
   }
   CATCH_EXIT("ODESolver");
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
   try {
      if (A2.getNbRows()!=_nb_eq)
         THROW_RT("setMatrices(...): Matrix size is different from system size");
   }
   CATCH_EXIT("ODESolver");
   _A2 = &A2;
   _order = 2;
   _lhs = true;
}


void ODESolver::setF(string F)
{
   static size_t nb_eq=0;
   if (_step==1)
      nb_eq = 0;
   nb_eq++;
   try {
      if (nb_eq>_nb_eq && _step>1)
         THROW_RT("setF(string): Number of calls is larger than system size.");
   }
   CATCH_EXIT("ODESolver");
   _type = SCALAR_NL;
   if (nb_eq>1)
      _type = VECTOR_NL;
   _regex = true;
   _expF.push_back(F);
   _vF1.push_back(eval(F,_time,_y0));
   _lhs = _rhs = true;
}


void ODESolver::setRK4RHS(Vect<real_t>& f)
{
   _f01 = &f;
}


real_t ODESolver::eval(string exp,
                       real_t t,
                       real_t y)
{
   real_t v, d[2];
   d[0] = t; d[1] = y;
   int err;
   theParser.Parse(exp.c_str(),"t,y");
   try {
      v = theParser.Eval(d);
      if ((err=theParser.EvalError()))
         THROW_RT("eval(string,real_t,real_t): Illegal algebraic expression "+itos(err));
   }
   CATCH("ODESolver");
   return v;
}


real_t ODESolver::eval(string exp,
                       real_t t)
{
   real_t v;
   int err;
   theParser.Parse(exp.c_str(),"t");
   try {
      v = theParser.Eval(t);
      if ((err=theParser.EvalError()))
         THROW_RT("eval(string,real_t): Illegal algebraic expression "+itos(err));
   }
   CATCH("ODESolver");
   return v;
}


void ODESolver::set(TimeScheme s,
                    real_t     time_step,
                    real_t     final_time)
{
   _final_time = final_time;
   _time_step = _time_step0 = time_step;
   _nb_ssteps = 1;
   _explicit = false;
   if (s==FORWARD_EULER || s==HEUN || s==RK4 || s==LEAP_FROG || s==AB2)
      _explicit = true;
   if (s==HEUN)
      _nb_ssteps = 2;
   if (s==RK4)
      _nb_ssteps = 4;
   _sc = int(s);
   try {
      if (s<0 || s>10) {
         _solve = NULL;
         THROW_RT("set(...): Time integration scheme not available.");
      }
      else
         _solve = TS[int(s)];
   }
   CATCH_EXIT("ODESolver");
   _constant_matrix = false;
}


void ODESolver::setInitial(real_t u)
{
   _y0 = u;
   _nb_eq = 1;
   _init = true;
}


void ODESolver::setInitial(real_t u,
                           real_t v)
{
   _y0 = u;
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
   try {
      if (u.size()!=_nb_eq)
         THROW_RT("setInitial(Vect<real_t>): Vector size is different from system size");
   }
   CATCH_EXIT("ODESolver");
   _w = &u;
   _u = *_w;
   _order = 1;
   _init = true;
}


void ODESolver::setInitial(Vect<real_t>& u,
                           Vect<real_t>& v)
{
   setInitial(u);
   try {
      if (v.size()!=_nb_eq)
         THROW_RT("setInitial(Vect<real_t>,Vect<real_t>): Vector size is different from system size");
   }
   CATCH_EXIT("ODESolver");
   _du = &v;
   _order = 2;
   _init = true;
}


void ODESolver::setInitialRHS(Vect<real_t>& f)
{
   _f1.setSize(_nb_eq);
   _f0 = _f1 = f;
}


void ODESolver::setRHS(real_t f)
{
   _d2 = f;
   _rhs = true;
}


void ODESolver::setRHS(string f)
{
   static size_t nb_eq=0;
   if (theTime>0)
      nb_eq = 0;
   nb_eq++;
   try {
      if (nb_eq>_nb_eq && _step>1)
         THROW_RT("setRHS(string): Number of calls is larger than system size.");
   }
   CATCH_EXIT("ODESolver");
   _type = SCALAR_LINEAR;
   if (nb_eq>1)
      _type = VECTOR_LINEAR;
   _expF.push_back(f);
   _f1.push_back(eval(f,_time));
   _regex = true;
   _rhs = true;
}


void ODESolver::setRHS(Vect<real_t>& f)
{
   try {
      if (f.size()!=_nb_eq)
         THROW_RT("setRHS(Vect<real_t>): Vector size is different from system size");
   }
   CATCH_EXIT("ODESolver");
   _f2 = &f;
   if (_step==0 && _f1.size()==0)
      _f1.setSize(_nb_eq);
   _rhs = true;
}


void ODESolver::setNewmarkParameters(real_t beta,
                                     real_t gamma)
{
   _beta = beta;
   _gamma = gamma;
}


real_t ODESolver::runOneTimeStep()
{
   if (_step==0) {
      try {
         if (_init==false)
         THROW_RT("runOneTimeStep(): No initial condition has been given.");
      }
      CATCH_EXIT("ODESolver");
      try {
         if (_lhs==false && _rhs==false)
         THROW_RT("runOneTimeStep(): ODE insufficiently defined.");
      }
      CATCH_EXIT("ODESolver");
   }
   _step++;
   _time = theTime;
   if (_verb>0)
      cout << "Running time step: " << _step << ", time: " << _time << " ..." << endl;
   if (_regex) {
      if (_type==SCALAR_LINEAR) {
         _c0 = eval(_exc[0],_time);
         _c1 = eval(_exc[1],_time);
         _c2 = eval(_exc[2],_time);
         _d2 = eval(_expF[0],_time);
         _d01 = eval(_expF[0],_time-0.5*_time_step);
         _order = 2;
         if (_c2==0)
            _order = 1;
      }
      else if (_type==VECTOR_LINEAR || _type==VECTOR_NL) {
         if (_step==1) {
            _f2 = new Vect<real_t>(_nb_eq);
            _f01 = new Vect<real_t>(_nb_eq);
         }
         for (size_t i=0; i<_nb_eq; i++) {
            _f1[i] = eval(_expF[i],_time-_time_step);
            (*_f2)[i] = eval(_expF[i],_time);
            (*_f01)[i] = eval(_expF[i],_time-0.5*_time_step);
         }
      }
   }
   (this->*_solve)();
   theTimeStep = _time_step;
   return _time_step;
}


void ODESolver::run(bool opt)
{
   if (_verb>0)
      cout << "Running time integration ..." << endl;
   _constant_matrix = opt;
   TimeLoop {
      runOneTimeStep();
      theTimeStep = _time_step;
   }
}


void ODESolver::solveForwardEuler()
{
   try {
      if (_order==2)
         THROW_RT("solveForwardEuler(): Forward Euler scheme is implemented for first order equations only.");
   }
   CATCH_EXIT("ODESolver");
   if (_type==SCALAR_LINEAR) {
      _y1 = _y0 + _time_step*(_d1-_c0*_y0)/_c1;
      _y2 = _y0 = _y1;
      _d1 = _d2;
      return;
   }
   if (_type==SCALAR_NL) {
      if (_regex)
         _d1 = eval(_expF[0],_time-_time_step,_y0);
      _y1 = _y0 + _time_step*_d1;
      _y2 = _y0 = _y1;
      _d1 = _d2;
      return;
   }
   if (_type==VECTOR_NL) {
     for (size_t i=0; i<_nb_eq; i++)
         _vF1[i] = eval(_expF[i],_time-_time_step,_u[i]);
      _v = _u + _time_step*_vF1;
      *_w = _u = _v;
      return;
   }
   _b = _time_step*_f1;
   _A0->MultAdd(-_time_step,_u,_b);
   _ls.setMatrix(_A1);
   _ls.setRHS(_b);
   _ls.setSolution(_v);
   _ls.setSolver(_s,_p);
   if (_step==1 || !_constant_matrix)
      _ls.setFact();
   _ls.solve();
   if (_constant_matrix)
      _ls.setNoFact();
   _v += _u;
   *_w = _u = _v;
   _f1 = *_f2;
}


void ODESolver::solveBackwardEuler()
{
   try {
      if (_order==2)
         THROW_RT("solveBackwardEuler(): Backward Euler scheme is implemented for first order equations only.");
   }
   CATCH_EXIT("ODESolver");
   try {
      if (_type==SCALAR_NL || _type==VECTOR_NL)
         THROW_RT("solveBackwardEuler(): Backward Euler scheme is not implemented for a nonlinear ODE.");
   }
   CATCH_EXIT("ODESolver");
   if (_type==SCALAR_LINEAR) {
      _y1 = (_c1*_y0+_time_step*_d2)/(_c1+_time_step*_c0);
      _y2 = _y0 = _y1;
      return;
   }
   _b = (*_f2);
   _A1->MultAdd(1./_time_step,_u,_b);
   _ls.setMatrix(_A0);
   _ls.setRHS(_b);
   _ls.setSolution(_v);
   _ls.setSolver(_s,_p);
   if (_step==1 || !_constant_matrix) {
      _A0->Axpy(1./_time_step,_A1);
      _ls.setFact();
   }
   _ls.solve();
   if (_constant_matrix)
      _ls.setNoFact();
   *_w = _u = _v;
}


void ODESolver::solveCrankNicolson()
{
   try {
      if (_order==2)
         THROW_RT("solveCrankNicolson(): Crank-Nicolson scheme is implemented for first order equations only.");
   }
   CATCH_EXIT("ODESolver");
   try {
      if (_type==SCALAR_NL || _type==VECTOR_NL)
         THROW_RT("solveCrankNicolson(): Crank-Nicolson scheme is not implemented for nonlinear ODEs.");
   }
   CATCH_EXIT("ODESolver");
   if (_type==SCALAR_LINEAR) {
      _y1 = (_c1 - 0.5*_time_step*_c0)*_y0 + 0.5*_time_step*(_d1+_d2);
      _y1 /= _c1 + 0.5*_time_step*_c0;
      _y2 = _y0 = _y1;
      _d1 = _d2;
      return;
   }
   _b = 0.5*(_f1+(*_f2));
   _ls.setMatrix(_A0);
   _ls.setRHS(_b);
   _ls.setSolution(_v);
   _A1->MultAdd(2./_time_step,_u,_b);
   _ls.setSolver(_s,_p);
   if (_step==1 || !_constant_matrix) {
      _A0->Axpy(2./_time_step,_A1);
      _ls.setFact();
   }
   _ls.solve();
   if (_constant_matrix)
      _ls.setNoFact();
   _v = 2.*_v - _u;
   *_w = _u = _v;
   _f1 = (*_f2);
}


void ODESolver::solveHeun()
{
   try {
      if (_order==2)
         THROW_RT("solveHeun(): Heun scheme is implemented for first order equations only.");
   }
   CATCH_EXIT("ODESolver");
   if (_type==SCALAR_LINEAR) {
      real_t yz = _y0 + _time_step*(_d1-_c0*_y0)/_c1;
      _y1 = _y0 + 0.5*_time_step*(_d1+_d2-_c0*(_y0+yz))/_c1;
      _y2 = _y0 = _y1;
      _d1 = _d2;
      return;
   }
   if (_type==SCALAR_NL) {
      if (_regex) {
         _d1 = eval(_expF[0],_time-_time_step,_y0);
         _d2 = eval(_expF[0],_time,_y0+_time_step*_d1);
      }
      _y1 = _y0 + 0.5*_time_step*(_d1+_d2);
      _y2 = _y0 = _y1;
      _d1 = _d2;
      return;
   }
   if (_type==VECTOR_NL) {
      for (size_t i=0; i<_nb_eq; i++) {
         _vF1[i] = eval(_expF[i],_time-_time_step,_u[i]);
         _vF2[i] = eval(_expF[i],_time-_time_step,_u[i]+_time_step*_vF1[i]);
      }
      _v = _u + _time_step*_vF1;
      *_w = _u = _v;
      _vF1 = _vF2;
      return;
   }
   _ls.setMatrix(_A1);
   _ls.setRHS(_b);
   _ls.setSolver(_s,_p);
   _b = _f1;
   _A0->MultAdd(-1.,_u,_b);
   _ls.setSolution(_v);
   if (_step==1 || !_constant_matrix)
      _ls.setFact();
   _ls.solve();
   if (_constant_matrix)
      _ls.setNoFact();
   _b = (*_f2);
   _A0->MultAdd(-1.,_u+_time_step*_v,_b);
   _ls.setSolution(*_w);
   _ls.solve();
   _u += 0.5*_time_step*((*_w)+_v);
   *_w = _v = _u;
   _f1 = (*_f2);
}


void ODESolver::solveLeapFrog()
{
   try {
      if (_order==2)
         THROW_RT("solveLeapFrog(): Leap Frog scheme is implemented for first order equations only.");
   }
   CATCH_EXIT("DESolver");
   if (_type==SCALAR_LINEAR) {
      if (_step==1) {
         _y1 = _y0 + _time_step*(_d0 - _c0*_y0)/_c1;
         _d1 = _d2;
      }
      else {
         _y2 = _y0 + 2.*_time_step*(_d1 - _c0*_y1)/_c1;
         _y0 = _y1; _y1 = _y2;
         _d0 = _d1, _d1 = _d2;
      }
      return;
   }
   else if (_type==SCALAR_NL) {
      if (_step==1)
         _y1 = _y0 + _time_step*eval(_expF[0],_time-_time_step,_y0);
      else {
         _y2 = _y0 + 2*_time_step*eval(_expF[0],_time-_time_step,_y1);
         _y0 = _y1; _y1 = _y2;
      }
      return;
   }
   else if (_type==VECTOR_NL) {
      if (_step==1) {
         for (size_t i=0; i<_nb_eq; i++)
            _v[i] = _u[i] + _time_step*eval(_expF[i],_time-_time_step,_u[i]);
      }
      else {
         for (size_t i=0; i<_nb_eq; i++)
            (*_w)[i] = _u[i] + 2*_time_step*eval(_expF[i],_time-_time_step,_v[i]);
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
      _ls.setFact();
      _ls.solve();
      _v += _u;
      _f1 = (*_f2);
   }
   else {
      _b = 2.*_time_step*_f1;
      _A0->MultAdd(-2.*_time_step,_v,_b);
      _ls.setSolution(*_w);
      if (!_constant_matrix)
         _ls.setFact();
      else
         _ls.setNoFact();
      _ls.solve();
      *_w += _u;
      _u = _v; _v = *_w;
      _f0 = _f1; _f1 = (*_f2);
   }
}


void ODESolver::solveAB2()
{
   try {
      if (_order==2)
         THROW_RT("solveAB2(): Adams-Bashforth scheme is implemented for first order equations only.");
   }
   CATCH_EXIT("ODESolver");
   if (_type==SCALAR_LINEAR) {
      if (_step==1)
         _y1 = _y0 + _time_step*(_d0-_c0*_y0)/_c1;   
      else {
         _y2 = _y1 - 0.5*_time_step*(_c0*(3*_y1-_y0) - (3*_d1-_d0))/_c1;
         _y0 = _y1; _y1 = _y2;
         _d0 = _d1; _d1 = _d2;
      }
      return;
   }
   else if (_type==SCALAR_NL) {
      if (_step==1)
         _y1 = _y0 + _time_step*eval(_expF[0],_time-_time_step,_y0);
      else {
         _d1 = eval(_expF[0],_time-2*_time_step,_y0);
         _d2 = eval(_expF[0],_time-_time_step,_y1);
         _y2 = _y1 + 0.5*_time_step*(3*_d2-_d1);
         _y0 = _y1; _y1 = _y2;
      }
      return;
   }
   else if (_type==VECTOR_NL) {
      if (_step==1)
         for (size_t i=0; i<_nb_eq; i++)
            _v[i] = _u[i] + _time_step*eval(_expF[i],_time-_time_step,_u[i]);
      else {
         for (size_t i=0; i<_nb_eq; i++) {
            _d1 = eval(_expF[i],_time-2*_time_step,_u[i]);
            _d2 = eval(_expF[i],_time-_time_step,_v[i]);
            (*_w)[i] = _v[i] + 0.5*_time_step*(3*_d2-_d1);
         }
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
      _ls.setFact();
      _ls.solve();
      if (_constant_matrix)
         _ls.setNoFact();
      _v += _u;
      _f1 = (*_f2);
   }
   else {
      _b = 0.5*_time_step*(3.*_f1-_f0);
      _A0->MultAdd(0.5*_time_step,_u-3.*_v,_b);
      _ls.setSolution(*_w);
      _ls.solve();
      *_w += _v;
      _u = _v; _v = *_w;
      _f0 = _f1; _f1 = (*_f2);
   }
}


void ODESolver::solveRK4()
{
   try {
      if (_order==2)
         THROW_RT("solveRK4(): Runge-Kutta scheme is valid for first order equations only.");
   }
   CATCH_EXIT("ODESolver");
   real_t k1, k2, k3, k4;
   if (_type==SCALAR_LINEAR) {
      k1 = _d1  - _c0*_y0;
      k2 = _d01 - _c0*(_y0 + 0.5*_time_step*k1); 
      k3 = _d01 - _c0*(_y0 + 0.5*_time_step*k2);
      k4 = _d2  - _c0*(_y0 +     _time_step*k3);
      _y2 = _y1 = _y0 + OFELI_SIXTH*_time_step*(k1+2*(k2+k3)+k4)/_c1;
      _y0 = _y1;
      _d0 = _d1 = _d2;
      return;
   }
   else if (_type==SCALAR_NL) {
      k1 = eval(_expF[0],_time-_time_step,_y0);
      k2 = eval(_expF[0],_time-0.5*_time_step,_y0+0.5*_time_step*k1);
      k3 = eval(_expF[0],_time-0.5*_time_step,_y0+0.5*_time_step*k2);
      k4 = eval(_expF[0],_time               ,_y0+    _time_step*k3);
      _y2 = _y1 = _y0 + OFELI_SIXTH*_time_step*(k1+2*(k2+k3)+k4);
      _y0 = _y1;
      return;
   }
   else if (_type==VECTOR_NL) {
      for (size_t i=0; i<_nb_eq; i++) {
         k1 = eval(_expF[i],_time-_time_step,_u[i]);
         k2 = eval(_expF[i],_time-0.5*_time_step,_u[i]+0.5*_time_step*k1);
         k3 = eval(_expF[i],_time-0.5*_time_step,_u[i]+0.5*_time_step*k2);
         k4 = eval(_expF[i],_time               ,_u[i]+    _time_step*k3);
         (*_w)[i] = _v[i] = _u[i] + OFELI_SIXTH*_time_step*(k1+2*(k2+k3)+k4);
      }
      _u = _v;
      return;
   }
   _k1.resize(_nb_eq), _k2.resize(_nb_eq);
   _k3.resize(_nb_eq), _k4.resize(_nb_eq);
   _k1 = _f1; _k4 = (*_f2);
   if (_f01)
      _k3 = _k2 = (*_f01);
   else
      _k3 = _k2 = 0.5*(_f1 + (*_f2));
   _A0->MultAdd(-1.,_u,_k1);
   _A0->MultAdd(-1.,_u+0.5*_time_step*_k1,_k2);
   _A0->MultAdd(-1.,_u+0.5*_time_step*_k2,_k3);
   _A0->MultAdd(-1.,_u+_time_step*_k3,_k4);
   _b = OFELI_SIXTH*_time_step*(_k1+2.*(_k2+_k3)+_k4);
   _ls.setMatrix(_A1);
   _ls.setRHS(_b);
   _ls.setSolver(_s,_p);
   _ls.setSolution(_v);
   if (_step==1 || !_constant_matrix)
      _ls.setFact();
   _ls.solve();
   if (_constant_matrix)
      _ls.setNoFact();
   _v += _u;
   *_w = _u = _v;
   _f0 = _f1; _f1 = (*_f2);
}


void ODESolver::solveRK3_TVD()
{
   try {
      if (_order==2)
         THROW_RT("solveRK3_TVD(): Runge-Kutta scheme is valid for first order equations only.");
   }
   CATCH_EXIT("ODESolver");
   real_t k1, k2;
   if (_type==SCALAR_LINEAR) {
      k1 = _d1  - _c0*_y0;
      k2 = _d01 - _c0*(_y0 + 0.5*_time_step*k1); 
      _y2 = _y1 = _y0 + OFELI_THIRD*_time_step*(k1+k2)/_c1;
      _y0 = _y1;
      _d0 = _d1 = _d2;
      return;
   }
   else if (_type==SCALAR_NL) {
      k1 = eval(_expF[0],_time-_time_step,_y0);
      k2 = eval(_expF[0],_time-0.5*_time_step,_y0+0.5*_time_step*k1);
      _y2 = _y1 = _y0 + OFELI_SIXTH*_time_step*(k1+2*k2);
      _y0 = _y1;
      return;
   }
   else if (_type==VECTOR_NL) {
      for (size_t i=0; i<_nb_eq; i++) {
         k1 = eval(_expF[i],_time-_time_step,_u[i]);
         k2 = eval(_expF[i],_time-0.5*_time_step,_u[i]+0.5*_time_step*k1);
         (*_w)[i] = _v[i] = _u[i] + OFELI_SIXTH*_time_step*(k1+2*k2);
      }
      _u = _v;
      return;
   }
   _k1.resize(_nb_eq), _k2.resize(_nb_eq);
   _k1 = _f1;
   if (_f01)
      _k2 = (*_f01);
   else
      _k3 = _k2 = 0.5*(_f1 + (*_f2));
   _A0->MultAdd(-1.,_u,_k1);
   _A0->MultAdd(-1.,_u+0.5*_time_step*_k1,_k2);
   _A0->MultAdd(-1.,_u+0.5*_time_step*_k2,_k3);
   _b = OFELI_SIXTH*_time_step*(_k1+2.*_k2);
   _ls.setMatrix(_A1);
   _ls.setRHS(_b);
   _ls.setSolver(_s,_p);
   _ls.setSolution(_v);
   if (_step==1 || !_constant_matrix)
      _ls.setFact();
   _ls.solve();
   if (_constant_matrix)
      _ls.setNoFact();
   _v += _u;
   *_w = _u = _v;
   _f0 = _f1; _f1 = (*_f2);
}


void ODESolver::solveNewmark()
{
   try {
      if (_order==1)
         THROW_RT("solveNewmark(): Newmark scheme is valid for second order equations only.");
   }
   CATCH_EXIT("ODESolver");
   try {
      if (_type==SCALAR_NL)
         THROW_RT("solveNewmark(): Newmark scheme is not implemented for a nonlinear ODE.");
   }
   CATCH_EXIT("ODESolver");
   if (_type==SCALAR_LINEAR) {
      if (_step==1)
         _ddy = (_d0 - _c1*_dy1 - _c0*_y1)/_c2;
      _y1 += _time_step*(_dy1 + _time_step*(0.5-_beta)*_ddy);
      _dy2 = _dy1 + (1-_gamma)*_time_step*_ddy;
      _y2 = _beta*_time_step*_time_step*(_d2-_c1*_dy1) + (_c2+_gamma*_time_step*_c1)*_y2;
      _y2 /= _c2 + _gamma*_time_step*_c1 + _beta*_time_step*_time_step;
      _ddy = (_y2-_y1)/(_beta*_time_step*_time_step);
      _y1 = _y2; _dy1 = _dy2;
      return;
   }
   if (_step==1) {
      _b = _f1;
      _A1->MultAdd(-1.,*_du,_b);
      _A0->MultAdd(-1.,_u,_b);
      LinearSolver<real_t> ls0;
      DMatrix<real_t> A2(*_A2);
      ls0.setMatrix(&A2);
      ls0.setSolution(_ddu);
      ls0.setRHS(_b);
      ls0.setSolver(_s,_p);
      ls0.setFact();
      ls0.solve();
   }
   real_t a=1./(_beta*_time_step*_time_step), b=_gamma/(_beta*_time_step);
   _v = _u + _time_step*(*_du) + (_time_step*_time_step*(0.5-_beta))*_ddu;
   *_du += ((1.-_gamma)*_time_step)*_ddu;
   _b = (*_f2);
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
   _ls.setFact();
   _ls.solve();
   _ddu = a*(*_w-_v);
   *_du += (_time_step*_gamma)*_ddu;
   _u = *_w;
   _f1 = (*_f2);
}


void ODESolver::solveBDF2()
{
   try {
      if (_order==2)
         THROW_RT("solveBDF2(): BDF2 scheme is implemented for first order equations only.");
   }
   CATCH_EXIT("ODESolver");
   try {
      if (_type==SCALAR_NL)
         THROW_RT("solveBDF2(): BDF2 scheme is not implemented for a nonlinear ODE.");
   }
   CATCH_EXIT("ODESolver");
   if (_type==SCALAR_LINEAR) {
      if (_step==1)
         _y1 = (_time_step*_d2 + _c1*_y0)/(_c1+_time_step*_c0);
      else {
         _y2 = (2*_time_step*_d2 + _c1*(4*_y1-_y0))/(3*_c1+2*_time_step*_c0);
         _y0 = _y1; _y1 = _y2;
      }
      return;
   }
   _b = (*_f2);
   _ls.setMatrix(_A0);
   _ls.setRHS(_b);
   _ls.setSolver(_s,_p);
   if (_step==1) {
      _A1->MultAdd(1./_time_step,_u,_b);
      _ls.setSolution(_v);
      _A0->Axpy(1./_time_step,_A1);
      _ls.setFact();
      _ls.solve();
      *_w = _v;
   }
   else {
      _A1->MultAdd(1./_time_step,2.*_v-0.5*_u,_b);
      _A0->Axpy(1.5/_time_step,_A1);
      _ls.setSolution(*_w);
      if (_step==2 || !_constant_matrix)
         _ls.setFact();
      _ls.solve();
      if (_constant_matrix)
         _ls.setNoFact();
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
