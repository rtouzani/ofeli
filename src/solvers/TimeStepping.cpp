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

                       Implementation of class 'TimeStepping'

  ==============================================================================*/

#include "solvers/TimeStepping.h"
#include "solvers/LinearSolver.h"
#include "OFELIException.h"
using std::cout;

namespace OFELI {

TimeStepping::TSPtr TimeStepping::TS [] = {
   &TimeStepping::solveStationary,          // Stationary problem
   &TimeStepping::solveForwardEuler,        // Forward Euler
   &TimeStepping::solveBackwardEuler,       // Backward Euler
   &TimeStepping::solveCrankNicolson,       // Crank-Nicolson
   &TimeStepping::solveHeun,                // Heun
   &TimeStepping::solveNewmark,             // Newmark
   &TimeStepping::solveLeapFrog,            // Leap Frog
   &TimeStepping::solveAB2,                 // Adams Bashforth 2
   &TimeStepping::solveRK4,                 // Runge Kutta 4
   &TimeStepping::solveRK3_TVD,             // TVD Runge Kutta 3
   &TimeStepping::solveBDF2                 // Backward Differentiation Formula
};

TimeStepping::ASPtr TimeStepping::AS [] = {
   &TimeStepping::AssembleStationary,  
   &TimeStepping::AssembleForwardEuler,  
   &TimeStepping::AssembleBackwardEuler,
   &TimeStepping::AssembleCrankNicolson,
   &TimeStepping::AssembleHeun,
   &TimeStepping::AssembleNewmark,
   &TimeStepping::AssembleLeapFrog,
   &TimeStepping::AssembleAB2,
   &TimeStepping::AssembleRK4,
   &TimeStepping::AssembleRK3_TVD,
   &TimeStepping::AssembleBDF2
};

TimeStepping::ASSPtr TimeStepping::ASS [] = {
   &TimeStepping::SAssembleStationary,  
   &TimeStepping::SAssembleForwardEuler,  
   &TimeStepping::SAssembleBackwardEuler,
   &TimeStepping::SAssembleCrankNicolson,
   &TimeStepping::SAssembleHeun,
   &TimeStepping::SAssembleNewmark,
   &TimeStepping::SAssembleLeapFrog,
   &TimeStepping::SAssembleAB2,
   &TimeStepping::SAssembleRK4,
   &TimeStepping::SAssembleRK3_TVD,
   &TimeStepping::SAssembleBDF2
};

TimeStepping::RHSPtr TimeStepping::RHS [] = {
   &TimeStepping::setRHS_Stationary,
   &TimeStepping::setRHS_ForwardEuler,
   &TimeStepping::setRHS_BackwardEuler,
   &TimeStepping::setRHS_CrankNicolson,
   &TimeStepping::setRHS_Heun,
   &TimeStepping::setRHS_Newmark,
   &TimeStepping::setRHS_LeapFrog,
   &TimeStepping::setRHS_AB2,
   &TimeStepping::setRHS_RK4,
   &TimeStepping::setRHS_RK3_TVD,
   &TimeStepping::setRHS_BDF2
};

TimeStepping::PreSolvePtr TimeStepping::PS [] = {
   &TimeStepping::PreSolve_Stationary,
   &TimeStepping::PreSolve_ForwardEuler,
   &TimeStepping::PreSolve_BackwardEuler,
   &TimeStepping::PreSolve_CrankNicolson,
   &TimeStepping::PreSolve_Heun,
   &TimeStepping::PreSolve_Newmark,
   &TimeStepping::PreSolve_LeapFrog,
   &TimeStepping::PreSolve_AB2,
   &TimeStepping::PreSolve_RK4,
   &TimeStepping::PreSolve_RK3_TVD,
   &TimeStepping::PreSolve_BDF2
};


TimeStepping::TimeStepping() :
              _order(0), _step(0), _nb_pdes(0), _ind(0), _verb(1), _regex(false), _time_step(0.1),
              _time(0.), _final_time(1.), _beta(0.25), _gamma(0.5), _nl_toler(1.e-6),
              _max_nl_it(100)
{
   for (int i=0; i<MAX_NB_PDES; ++i) {
      _theMesh[i] = NULL;
      _theEqua[i] = NULL;
      _p[i] = IDENT_PREC;
      _A[i] = NULL;
      _constant_matrix[i] = false;
      _explicit[i] = false;
      _set_bc[i] = false;
      _s[i] = DIRECT_SOLVER;
   }
}


TimeStepping::TimeStepping(TimeScheme s,
                           real_t     time_step,
                           real_t     final_time) :
              _order(0), _step(0), _nb_pdes(0), _ind(0), _verb(1), _time_step(0.1), _time(0.),
              _final_time(1.), _beta(0.25), _gamma(0.5), _nl_toler(1.e-6), _max_nl_it(100)
{
   for (int i=0; i<MAX_NB_PDES; ++i) {
      _theMesh[i] = NULL;
      _theEqua[i] = NULL;
      _p[i] = IDENT_PREC;
      _A[i] = NULL;
      _constant_matrix[i] = false;
      _explicit[i] = false;
      _s[i] = DIRECT_SOLVER;
      _set_bc[i] = false;
   }
   set(s,time_step,final_time);
}


TimeStepping::~TimeStepping()
{
   for (int i=0; i<MAX_NB_PDES; ++i)
      if (_A[i])
         delete _A[i];
}


void TimeStepping::setPDE(AbsEqua<real_t>& eq,
                          bool             nl)
{
   _theEqua[_nb_pdes] = &eq;
   _ls = &(_theEqua[_nb_pdes]->getLinearSolver());
   _s[_nb_pdes] = _ls->getSolver();
   _theMesh[_nb_pdes] = &(_theEqua[_nb_pdes]->getMesh());
   _nb_eq[_nb_pdes] = _theMesh[_nb_pdes]->getNbEq();
   _nl[_nb_pdes] = nl;
   if (_A[_nb_pdes])
      delete _A[_nb_pdes];
   if (_explicit[_nb_pdes])
      _D[_nb_pdes].setSize(_nb_eq[_nb_pdes]);
   else if (_s[_nb_pdes]==DIRECT_SOLVER)
      _A[_nb_pdes] = new SkMatrix<real_t>(*_theMesh[_nb_pdes]);
   else
      _A[_nb_pdes] = new SpMatrix<real_t>(*_theMesh[_nb_pdes]);
   if (_sc==int(NEWMARK))
      _D[_nb_pdes].setSize(_nb_eq[_nb_pdes]);
   if (_sc==int(RK4)) {
      _k1[_nb_pdes].setSize(_theMesh[_nb_pdes]->getNbDOF());
      _k2[_nb_pdes].setSize(_theMesh[_nb_pdes]->getNbDOF());
      _k3[_nb_pdes].setSize(_theMesh[_nb_pdes]->getNbDOF());
   }
   _b[_nb_pdes].setSize(_nb_eq[_nb_pdes]);
   _vv[_nb_pdes].setSize(_theMesh[_nb_pdes]->getNbEq());
   _f[_nb_pdes].setSize(_theMesh[_nb_pdes]->getNbDOF());
   if (!_nl[_nb_pdes])
      _max_nl_it = 1;
   _nb_pdes++, _ind++;
}


void TimeStepping::setRK4RHS(Vect<real_t>& f)
{
   _f01[_nb_pdes-1] = &f;
}


void TimeStepping::set(TimeScheme s,
                       real_t     time_step,
                       real_t     final_time)
{
   _final_time = final_time;
   _time_step = _time_step0 = time_step;
   _nb_ssteps = 1;
   for (int i=0; i<MAX_NB_PDES; ++i) {
      _explicit[i] = false;
      if (s==FORWARD_EULER || s==HEUN || s==RK4 || s==RK3_TVD || s==LEAP_FROG || s==AB2)
         _explicit[i] = true;
      _constant_matrix[i] = false;
   }
   if (s==HEUN || s==NEWMARK)
      _nb_ssteps = 2;
   if (s==RK4)
      _nb_ssteps = 4;
   _sc = int(s);
   if (s<0 || s>10) {
      _solve = NULL;
      throw OFELIException("In TimeStepping::set(...): Time integration scheme not available.");
   }
   else {
      _presolve = PS[int(s)];
      _solve = TS[int(s)];
      _assemb = AS[int(s)];
      _set_rhs = RHS[int(s)];
   }
}


void TimeStepping::setLinearSolver(Iteration      s,
                                   Preconditioner p)
{
   _s[_nb_pdes-1] = s;
   _p[_nb_pdes-1] = p;
   if (_A[_nb_pdes-1])
      delete _A[_nb_pdes-1];
   if (_s[_nb_pdes-1]==DIRECT_SOLVER)
      _A[_nb_pdes-1] = new SkMatrix<real_t>(*_theMesh[_nb_pdes-1]);
   else
      _A[_nb_pdes-1] = new SpMatrix<real_t>(*_theMesh[_nb_pdes-1]);
}


void TimeStepping::setInitial(Vect<real_t>& u)
{
   _nb_eq[_nb_pdes-1] = _theMesh[_nb_pdes-1]->getNbEq();
   _nb_dof[_nb_pdes-1] = _theMesh[_nb_pdes-1]->getNbDOF();
   _theEqua[_nb_pdes-1]->setInput(INITIAL_FIELD,u);
   _v[_nb_pdes-1].setSize(u.getNx(),u.getNy(),u.getNz());
   _u[_nb_pdes-1].setSize(u.getNx(),u.getNy(),u.getNz());
   _w[_nb_pdes-1] = &u;
   _u[_nb_pdes-1] = *_w[_nb_pdes-1];
   _b[_nb_pdes-1].setSize(_nb_eq[_nb_pdes-1]);
   if (_f0[_nb_pdes-1].size()==0)
      _f0[_nb_pdes-1].setSize(u.getNx(),u.getNy(),u.getNz());
   _f1[_nb_pdes-1].setSize(u.getNx(),u.getNy(),u.getNz());
   _order = 1;
}


void TimeStepping::setInitial(Vect<real_t>& u,
                              Vect<real_t>& v)
{
   setInitial(u);
   _du[_nb_pdes-1] = &v;
   _ddu[_nb_pdes-1].setSize(u.getNx(),u.getNy());
   _order = 2;
}


void TimeStepping::setInitialRHS(Vect<real_t>& f)
{
   if (_f0[_nb_pdes-1].size()==0)
      _f0[_nb_pdes-1].setSize(f.getNx(),f.getNy(),f.getNz());
   _f0[_nb_pdes-1] = f;
}


void TimeStepping::setRHS(Vect<real_t>& f)
{
   _f2[_nb_pdes-1] = &f;
   if (_f0[_nb_pdes-1].size()==0)
      _f0[_nb_pdes-1].setSize(f.getNx(),f.getNy(),f.getNz());
   if (_f1[_nb_pdes-1].size()==0)
      _f1[_nb_pdes-1].setSize(f.getNx(),f.getNy(),f.getNz());
}


void TimeStepping::setBC(Vect<real_t>& u)
{
   _bc[_nb_pdes-1] = &u;
   _set_bc[_nb_pdes-1] = true;
}


void TimeStepping::setNewmarkParameters(real_t beta,
                                        real_t gamma)
{
   _beta = beta;
   _gamma = gamma;
}


void TimeStepping::SAssembly(const Side& sd,
                             real_t*     b,
                             real_t*     A)
{
   (this->*_sassemb)(sd,b,A);
}


void TimeStepping::Assembly(const Element& el,
                            real_t*        b,
                            real_t*        A0,
                            real_t*        A1,
                            real_t*        A2)
{
   (this->*_assemb)(el,b,A0,A1,A2);
}


real_t TimeStepping::runOneTimeStep()
{
   _step++;
   _time = theTime;
   if (_verb>0)
      cout << "Running time step: " << _step << ", time: " << _time << " ..." << endl;
   for (int i=0; i<_nb_pdes; ++i) {
      if (_bc[i])
         _theEqua[i]->setInput(BOUNDARY_CONDITION,*_bc[i]);
      if (_nl[i]) {
         for (_sstep=1; _sstep<=_nb_ssteps; _sstep++) {
            int it=1;
            real_t err=1;
            while (it<=_max_nl_it && err>_nl_toler) {
               _theEqua[i]->setInput(SOURCE,(this->*_set_rhs)());
               if (_explicit[i])
                  _D[i] = 0;
               else
                  (*_A)[i] = 0;
               _b[i] = 0;
               if (_sc==int(RK4))
                  _k1[i] = 0, _k2[i] = 0, _k3[i] = 0;
               (this->*_presolve)();
               _theEqua[i]->build(*this);
               (this->*_solve)();
               it++;
            }
         }
      }
      else {
         for (_sstep=1; _sstep<=_nb_ssteps; _sstep++) {
            _theEqua[i]->setInput(SOURCE,(this->*_set_rhs)());
            if (_explicit[i])
               _D[i] = 0;
            else
               *_A[i] = 0;
            _b[i] = 0;
            if (_sc==int(RK4))
               _k1[i] = 0, _k2[i] = 0, _k3[i] = 0;
            (this->*_presolve)();
            _theEqua[i]->build(*this);
            (this->*_solve)();
         }
      }
   }
   theTimeStep = _time_step;
   return _time_step;
}


void TimeStepping::run(bool opt)
{
   if (_verb>0)
      cout << "Running time integration ..." << endl;
   for (int i=0; i<_nb_pdes; ++i)
      _constant_matrix[i] = opt;
   TimeLoop {
      runOneTimeStep();
      theTimeStep = _time_step;
   }
}


void TimeStepping::solveStationary()
{
   for (int i=0; i<_nb_pdes; ++i) {
      _ls->setMatrix(_A[i]);
      _ls->setRHS(_b[i]);
      _ls->setSolution(_vv[i]);
      _ls->solve(_s[i],_p[i]);
      insertBC(_vv[i],_v[i]);
      *_w[i] = _u[i] = _v[i];
   }
}


void TimeStepping::solveForwardEuler()
{
   if (_order==2)
      throw OFELIException("In TimeStepping::solveForwardEuler(): Forward Euler "
                           "scheme is implemented for first order equations only.");
   for (int n=0; n<_nb_pdes; ++n) {
      for (size_t i=1; i<=_nb_eq[n]; ++i)
         _b[n](i) /= _D[n](i);
      insertBC(_b[n],_v[n]);
      *_w[n] = _u[n] = _v[n];
      _f1[n] = *_f2[n];
   }
}


void TimeStepping::solveBackwardEuler()
{
   if (_order==2)
      throw OFELIException("In TimeStepping::solveBackwardEuler(): Backward Euler "
                           "scheme is implemented for first order equations only.");
   _ls->setMatrix(_A[_ind-1]);
   _ls->setRHS(_b[_ind-1]);
   _ls->setSolution(_vv[_ind-1]);
   _ls->solve(_s[_ind-1],_p[_ind-1]);
   insertBC(_vv[_ind-1],_v[_ind-1]);
   *_w[_ind-1] = _u[_ind-1] = _v[_ind-1];
}


void TimeStepping::solveCrankNicolson()
{
   if (_order==2)
      throw OFELIException("In TimeStepping::solveCrankNicolson(): "
                           "Crank-Nicolson scheme is implemented for first order equations only.");
   _ls->setMatrix(_A[_ind-1]);
   _ls->setRHS(_b[_ind-1]);
   _ls->setSolution(_vv[_ind-1]);
   _ls->solve(_s[_ind-1],_p[_ind-1]);
   insertBC(_vv[_ind-1],_v[_ind-1]);
   *_w[_ind-1] = _u[_ind-1] = _v[_ind-1];
   _f1[_ind-1] = *_f2[_ind-1];
}


void TimeStepping::solveHeun()
{
   if (_order==2)
      throw OFELIException("In TimeStepping::solveHeun(): Heun scheme is "
                           "implemented for first order equations only.");
   for (size_t i=1; i<=_nb_eq[_ind-1]; i++)
      _b[_ind-1](i) /= _D[_ind-1](i);
   if (_sstep==1)
      insertBC(_b[_ind-1],_v[_ind-1]);
   if (_sstep==2) {
      insertBC(_b[_ind-1],*_w[_ind-1]);
      _u[_ind-1] = *_w[_ind-1];
   }
   *_w[_ind-1] = _v[_ind-1] = _u[_ind-1];
   _f1[_ind-1] = *_f2[_ind-1];
}


void TimeStepping::solveLeapFrog()
{
   if (_order==2)
      throw OFELIException("In TimeStepping::solveLeapFrog(): Leap Frog "
                           "scheme is implemented for first order equations only.");
   for (size_t i=1; i<=_nb_eq[_ind-1]; i++)
      _b[_ind-1](i) /= _D[_ind-1](i);
   if (_step==1) {
      insertBC(_b[_ind-1],_v[_ind-1]);
      _f1[_ind-1] = *_f2[_ind-1];
   }
   else {
      insertBC(_b[_ind-1],*_w[_ind-1]);
      _u[_ind-1] = _v[_ind-1]; _v[_ind-1] = *_w[_ind-1];
      _f0[_ind-1] = _f1[_ind-1];
      _f1[_ind-1] = *_f2[_ind-1];
   }
}


void TimeStepping::solveAB2()
{
   if (_order==2)
      throw OFELIException("In TimeStepping::solveAB2(): Adams-Bashforth "
                           "scheme is implemented for first order equations only.");
   for (size_t i=1; i<=_nb_eq[_ind-1]; i++)
      _b[_ind-1](i) /= _D[_ind-1](i);
   if (_step==1) {
      insertBC(_b[_ind-1],_v[_ind-1]);
      _f1[_ind-1] = *_f2[_ind-1];
   }
   else {
      insertBC(_b[_ind-1],*_w[_ind-1]);
      _u[_ind-1] = _v[_ind-1]; _v[_ind-1] = *_w[_ind-1];
      _f0[_ind-1] = _f1[_ind-1];
      _f1[_ind-1] = *_f2[_ind-1];
   }
}


void TimeStepping::solveRK4()
{
   if (_order==2)
      throw OFELIException("In TimeStepping::solveRK4(): Runge-Kutta scheme "
                           "is valid for first order equations only.");
   if (_sstep==4) {
      for (size_t i=1; i<=_nb_eq[_ind-1]; i++)
         _b[_ind-1](i) /= _D[_ind-1](i);
      insertBC(_b[_ind-1],_v[_ind-1]);
      *_w[_ind-1] = _u[_ind-1] = _v[_ind-1];
      _f0[_ind-1] = _f1[_ind-1];
      _f1[_ind-1] = *_f2[_ind-1];
      _f01[_ind-1] = NULL;
   }
   else
      return;
}


void TimeStepping::solveRK3_TVD()
{
}


void TimeStepping::solveNewmark()
{
   if (_order!=2)
      throw OFELIException("In TimeStepping::solveNewmark(): Newmark scheme "
                           "is valid for second order equations only.");
   if (_step==1 && _sstep==1) {
      for (size_t i=1; i<=_nb_eq[_ind-1]; i++)
         _b[_ind-1](i) /= _D[_ind-1](i);
      insertBC0(_b[_ind-1],_ddu[_ind-1]);
      return;
   }
   if (_sstep>1) {
      _ls->setMatrix(_A[_ind-1]);
      _ls->setRHS(_b[_ind-1]);
      _ls->setSolution(_vv[_ind-1]);
      _ls->solve(_s[_ind-1],_p[_ind-1]);
      insertBC(_vv[_ind-1],*_w[_ind-1]);
      _ddu[_ind-1] = (1./(_beta*_time_step*_time_step))*(*_w[_ind-1]-_v[_ind-1]);
      *_du[_ind-1] += (_time_step*_gamma)*_ddu[_ind-1];
      _u[_ind-1] = *_w[_ind-1];
   }
}


void TimeStepping::solveBDF2()
{
  if (_order==2)
     throw OFELIException("In TimeStepping::solveBDF2(): BDF2 scheme is "
                          "implemented for first order equations only.");
   _ls->setMatrix(_A[_ind-1]);
   _ls->setRHS(_b[_ind-1]);
   _ls->setSolution(_vv[_ind-1]);
   _ls->solve(_s[_ind-1],_p[_ind-1]);
   if (_step==1) {
      insertBC(_vv[_ind-1],_v[_ind-1]);
      *_w[_ind-1] = _v[_ind-1];
   }
   else {
      insertBC(_vv[_ind-1],*_w[_ind-1]);
      _u[_ind-1] = _v[_ind-1];
      _v[_ind-1] = *_w[_ind-1];
   }
}


void TimeStepping::AssembleStationary(const Element& el,
                                      real_t*        eb,
                                      real_t*        eA0,
                                      real_t*        eA1,
                                      real_t*        eA2)
{
   if (_set_bc[_ind-1])
      update_bc(el,Vect<real_t>(&el,*_bc[_ind-1]),eA0,eb);
   element_assembly(el,eA0,eb,_A[_ind-1],_b[_ind-1]);
}


void TimeStepping::SAssembleStationary(const Side& sd,
                                       real_t*     sb,
                                       real_t*     sA)
{
   element_assembly(sd,sA,sb,_A[_ind-1],_b[_ind-1]);
}


void TimeStepping::AssembleForwardEuler(const Element& el,
                                        real_t*        eb,
                                        real_t*        eA0,
                                        real_t*        eA1,
                                        real_t*        eA2)
{
   Vect<real_t> eu(&el,_u[_ind-1]);
   size_t k=0, n=el.getNbNodes()*el.getNbDOF();
   real_t d=1/_time_step;
   for (size_t i=0; i<n; i++) {
      for (size_t j=0; j<n; j++, k++) {
         eA1[k] *= d;
         eb[i] += (eA1[k]-eA0[k])*eu[j];
      }
   }
   if (_set_bc[_ind-1])
      update_bc_diag(el,Vect<real_t>(&el,*_bc[_ind-1]),eA1,eb);
   element_assembly(el,eA1,eb,_D[_ind-1],_b[_ind-1]);
}


void TimeStepping::SAssembleForwardEuler(const Side& sd,
                                         real_t*     sb,
                                         real_t*     sA)
{
   Vect<real_t> su(&sd,_u[_ind-1]);
   size_t k=0, n=sd.getNbNodes()*sd.getNbDOF();
   for (size_t i=0; i<n; i++) {
      for (size_t j=0; j<n; j++, k++)
         sb[i] -= sA[k]*su[j];
   }
   element_assembly(sd,sb,_b[_ind-1]);
}


void TimeStepping::AssembleBackwardEuler(const Element& el,
                                         real_t*        eb,
                                         real_t*        eA0,
                                         real_t*        eA1,
                                         real_t*        eA2)
{
   Vect<real_t> eu(&el,_u[_ind-1]);
   size_t k=0, n=el.getNbNodes()*el.getNbDOF();
   real_t d=1.0/_time_step;
   for (size_t i=0; i<n; i++) {
      for (size_t j=0; j<n; j++, k++) {
         eA0[k] += d*eA1[k];
         eb[i] += d*eA1[k]*eu[j];
      }
   }
   if (_set_bc[_ind-1])
      update_bc(el,Vect<real_t>(&el,*_bc[_ind-1]),eA0,eb);
   element_assembly(el,eA0,eb,_A[_ind-1],_b[_ind-1]);
}


void TimeStepping::SAssembleBackwardEuler(const Side& sd,
                                          real_t*     sb,
                                          real_t*     sA)
{
   Vect<real_t> su(&sd,_u[_ind-1]);
   size_t k=0, n=sd.getNbNodes()*sd.getNbDOF();
   real_t d=1.0/_time_step;
   for (size_t i=0; i<n; i++) {
      for (size_t j=0; j<n; j++, k++)
         sb[i] += d*sA[k]*su[j];
   }
   element_assembly(sd,sA,sb,_A[_ind-1],_b[_ind-1]);
}


void TimeStepping::AssembleCrankNicolson(const Element& el,
                                         real_t*        eb,
                                         real_t*        eA0,
                                         real_t*        eA1,
                                         real_t*        eA2)
{
   Vect<real_t> eu(&el,_u[_ind-1]);
   size_t k=0, n=el.getNbNodes()*el.getNbDOF();
   real_t d=1.0/_time_step;
   for (size_t i=0; i<n; i++) {
      for (size_t j=0; j<n; j++, k++) {
         real_t z=0.5*eA0[k];
         eA0[k] = d*eA1[k] + z;
         eb[i] += (d*eA1[k] - z)*eu[j];
      }
   }
   if (_set_bc[_ind-1])
      update_bc(el,Vect<real_t>(&el,*_bc[_ind-1]),eA0,eb);
   element_assembly(el,eA0,eb,_A[_ind-1],_b[_ind-1]);
}


void TimeStepping::SAssembleCrankNicolson(const Side& sd,
                                          real_t*     sb,
                                          real_t*     sA)
{
   Vect<real_t> su(&sd,_u[_ind-1]);
   size_t k=0, n=sd.getNbNodes()*sd.getNbDOF();
   for (size_t i=0; i<n; i++) {
      for (size_t j=0; j<n; j++, k++) {
         sA[k] *= 0.5;
         sb[i] -= sA[k]*su[j];
      }
   }
   element_assembly(sd,sA,sb,_A[_ind-1],_b[_ind-1]);
}


void TimeStepping::AssembleHeun(const Element& el,
                                real_t*        eb,
                                real_t*        eA0,
                                real_t*        eA1,
                                real_t*        eA2)
{
   Vect<real_t> eu(&el,_u[_ind-1]);
   size_t k=0, n=el.getNbNodes()*el.getNbDOF();
   real_t d=1.0/_time_step;
   if (_sstep==1) {
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++) {
            eA1[k] *= d;
            eb[i] += (eA1[k]-eA0[k])*eu[j];
         }
      }
   }
   if (_sstep==2) {
      Vect<real_t> ev(&el,_v[_ind-1]);
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++) {
            eA1[k] *= d;
            eb[i] += (eA1[k] - 0.5*eA0[k])*eu[j] - 0.5*eA0[k]*ev[j];
         }
      }
   }
   if (_set_bc[_ind-1])
      update_bc_diag(el,Vect<real_t>(&el,*_bc[_ind-1]),eA1,eb);
   element_assembly(el,eA1,eb,_D[_ind-1],_b[_ind-1]);
}


void TimeStepping::SAssembleHeun(const Side& sd,
                                 real_t*     sb,
                                 real_t*     sA)
{
   Vect<real_t> su(&sd,_u[_ind-1]);
   size_t k=0, n=sd.getNbNodes()*sd.getNbDOF();
   if (_sstep==1) {
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++)
            sb[i] -= sA[k]*su[j];
      }
   }
   if (_sstep==2) {
      Vect<real_t> sv(&sd,_v[_ind-1]);
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++)
            sb[i] -=  0.5*sA[k]*(su[j] + sv[j]);
      }
   }
   Element_Assembly(sd,sb,_b[_ind-1]);
}


void TimeStepping::AssembleAB2(const Element& el,
                               real_t*        eb,
                               real_t*        eA0,
                               real_t*        eA1,
                               real_t*        eA2)
{
   Vect<real_t> eu(&el,_u[_ind-1]);
   size_t k=0, n=el.getNbNodes()*el.getNbDOF();
   real_t d=0.5/_time_step;
   if (_step==1) {
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++) {
            eA1[k] *= d;
            eb[i] += (eA1[k]-eA0[k])*eu[j];
         }
      }
   }
   else {
      Vect<real_t> ev(&el,_v[_ind-1]);
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++) {
            eA1[k] *= d;
            eb[i] += (eA1[k] - 1.5*eA0[k])*ev[j] + 0.5*eA0[k]*eu[j];
         }
      }
   }
   if (_set_bc[_ind-1])
      update_bc_diag(el,Vect<real_t>(&el,*_bc[_ind-1]),eA1,eb);
   element_assembly(el,eA1,eb,_D[_ind-1],_b[_ind-1]);
}


void TimeStepping::SAssembleAB2(const Side& sd,
                                real_t*     sb,
                                real_t*     sA)
{
   Vect<real_t> su(&sd,_u[_ind-1]);
   size_t k=0, n=sd.getNbNodes()*sd.getNbDOF();
   if (_step==1) {
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++)
            sb[i] -= sA[k]*su[j];
      }
   }
   else {
      Vect<real_t> sv(&sd,_v[_ind-1]);
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++)
            sb[i] -= 1.5*sA[k]*(1.5*sv[j] - 0.5*su[j]);
      }
   }
   Element_Assembly(sd,sb,_b[_ind-1]);
}


void TimeStepping::AssembleLeapFrog(const Element& el,
                                    real_t*        eb,
                                    real_t*        eA0,
                                    real_t*        eA1,
                                    real_t*        eA2)
{
   Vect<real_t> eu(&el,_u[_ind-1]);
   size_t k=0, n=el.getNbNodes()*el.getNbDOF();
   real_t d=0.5/_time_step;
   if (_step==1) {
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++) {
            eA1[k] *= 2*d;
            eb[i] += (eA1[k]-eA0[k])*eu[j];
         }
      }
   }
   else {
      Vect<real_t> ev(&el,_v[_ind-1]);
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++) {
            eA1[k] *= d;
            eb[i] += eA1[k]*eu[j] - eA0[k]*ev[j];
         }
      }
   }
   if (_set_bc[_ind-1])
      update_bc_diag(el,Vect<real_t>(&el,*_bc[_ind-1]),eA1,eb);
   element_assembly(el,eA1,eb,_D[_ind-1],_b[_ind-1]);
}


void TimeStepping::SAssembleLeapFrog(const Side& sd,
                                     real_t*     sb,
                                     real_t*     sA)
{
   Vect<real_t> su(&sd,_u[_ind-1]);
   size_t k=0, n=sd.getNbNodes()*sd.getNbDOF();
   if (_step==1) {
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++)
            sb[i] -= sA[k]*su[j];
      }
   }
   else {
      Vect<real_t> sv(&sd,_v[_ind-1]);
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++)
            sb[i] -= sA[k]*sv[j];
      }
   }
   Element_Assembly(sd,sb,_b[_ind-1]);
}


void TimeStepping::AssembleRK4(const Element& el,
                               real_t*        eb,
                               real_t*        eA0,
                               real_t*        eA1,
                               real_t*        eA2)
{
   Vect<real_t> eu(&el,_u[_ind-1]);
   size_t m=0, n=el.getNbNodes()*el.getNbDOF();
   switch (_sstep) {

      case 1:
         {
            for (size_t i=0; i<n; i++) {
               for (size_t j=0; j<n; j++, m++)
                  eb[i] -= eA0[m]*eu[j];
            }
            element_assembly(el,eb,_k1);
         }
         break;

      case 2:
         {
            Vect<real_t> ek(&el,_k1[_ind-1]);
            for (size_t i=0; i<n; i++) {
               for (size_t j=0; j<n; j++, m++)
                  eb[i] -= eA0[m]*(eu[j] + 0.5*_time_step*ek[j]);
            }
            element_assembly(el,eb,_k2[_ind-1]);
         }
         break;

      case 3:
         {
            Vect<real_t> ek(&el,_k2[_ind-1]);
            for (size_t i=0; i<n; i++) {
               for (size_t j=0; j<n; j++, m++)
                  eb[i] -= eA0[m]*(eu[j] + 0.5*_time_step*ek[j]);
            }
            element_assembly(el,eb,_k3[_ind-1]);
         }
         break;

      case 4:
         {
            Vect<real_t> ek1(&el,_k1[_ind-1]), ek2(&el,_k2[_ind-1]), ek3(&el,_k3[_ind-1]);
            real_t d=1./_time_step;
            for (size_t i=0; i<n; i++) {
               for (size_t j=0; j<n; j++, m++)
                  eb[i] -= eA0[m]*(eu[j] + _time_step*ek3[j]);
            }
            for (size_t i=0, m=0; i<n; i++) {
               eb[i] = OFELI_SIXTH*(ek1[i]+2*ek2[i]+2*ek3[i]+eb[i]);
               for (size_t j=0; j<n; j++, m++) {
                  eA1[m] *= d;
                  eb[i] += eA1[m]*eu[j];
               }
            }
            if (_set_bc[_ind-1])
               update_bc_diag(el,Vect<real_t>(&el,*_bc[_ind-1]),eA1,eb);
            element_assembly(el,eA1,eb,_D[_ind-1],_b[_ind-1]);
         }
         break;
   }
}


void TimeStepping::AssembleRK3_TVD(const Element& el,
                                   real_t*        eb,
                                   real_t*        eA0,
                                   real_t*        eA1,
                                   real_t*        eA2)
{
}


void TimeStepping::SAssembleRK4(const Side& sd,
                                real_t*     sb,
                                real_t*     sA)
{
   Vect<real_t> su(&sd,_u[_ind-1]);
   size_t m=0, n=sd.getNbNodes()*sd.getNbDOF();
   switch (_sstep) {

      case 1:
         {
            for (size_t i=0; i<n; i++) {
               for (size_t j=0; j<n; j++, m++)
                  sb[i] -= sA[m]*su[j];
            }
            Element_Assembly(sd,sb,_k1[_ind-1]);
         }
         break;

      case 2:
         {
            Vect<real_t> sk(&sd,_k1[_ind-1]);
            for (size_t i=0; i<n; i++) {
               for (size_t j=0; j<n; j++, m++)
                  sb[i] -= sA[m]*(su[j] + 0.5*_time_step*sk[j]);
            }
            Element_Assembly(sd,sb,_k2[_ind-1]);
         }
         break;

      case 3:
         {
            Vect<real_t> sk(&sd,_k2[_ind-1]);
            for (size_t i=0; i<n; i++) {
               for (size_t j=0; j<n; j++, m++)
                  sb[i] -= sA[m]*(su[j] + 0.5*_time_step*sk[j]);
            }
            Element_Assembly(sd,sb,_k3[_ind-1]);
         }
         break;

      case 4:
         {
            Vect<real_t> sk1(&sd,_k1[_ind-1]), sk2(&sd,_k2[_ind-1]), sk3(&sd,_k3[_ind-1]);
            for (size_t i=0; i<n; i++) {
               for (size_t j=0; j<n; j++, m++)
                  sb[i] -= sA[m]*(su[j] + _time_step*sk3[j]);
            }
            Element_Assembly(sd,sb,_b[_ind-1]);
         }
         break;
   }
}


void TimeStepping::SAssembleRK3_TVD(const Side& sd,
                                    real_t*     sb,
                                    real_t*     sA)
{
}


void TimeStepping::AssembleNewmark(const Element& el,
                                   real_t*        eb,
                                   real_t*        eA0,
                                   real_t*        eA1,
                                   real_t*        eA2)
{
   size_t n=el.getNbNodes()*el.getNbDOF();
   Vect<real_t> edu(&el,*_du[_ind-1]);
   if (_step==1 && _sstep==1) {
      Vect<real_t> eu(&el,_u[_ind-1]);
      size_t k=0;
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++)
            eb[i] -= eA1[k]*edu[j] + eA0[k]*eu[j];
      }
      element_assembly(el,eA2,eb,_D[_ind-1],_b[_ind-1]);
   }
   if (_sstep>1) {
      Vect<real_t> ev(&el,_v[_ind-1]);
      real_t a=1./(_beta*_time_step*_time_step), b=_gamma/(_beta*_time_step);
      size_t k=0;
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++) {
            real_t z = a*eA2[k] + b*eA1[k];
            eA0[k] += z;
            eb[i] += z*ev[j] - eA1[k]*edu[j];
         }
      }
      if (_set_bc[_ind-1])
         update_bc(el,Vect<real_t>(&el,*_bc[_ind-1]),eA0,eb);
      element_assembly(el,eA0,eb,_A[_ind-1],_b[_ind-1]);
   }
}


void TimeStepping::SAssembleNewmark(const Side& sd,
                                    real_t*     sb,
                                    real_t*     sA)
{
   size_t n=sd.getNbNodes()*sd.getNbDOF();
   Vect<real_t> sdu(&sd,*_du[_ind-1]);
   if (_step==1 && _sstep==1) {
      Vect<real_t> su(&sd,_u[_ind-1]);
      size_t k=0;
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++)
            sb[i] -= sA[k]*su[j];
      }
      Element_Assembly(sd,sb,_b[_ind-1]);
   }
   if (_sstep>1)
      element_assembly(sd,sA,sb,_A[_ind-1],_b[_ind-1]);
}


void TimeStepping::AssembleBDF2(const Element& el,
                                real_t*        eb,
                                real_t*        eA0,
                                real_t*        eA1,
                                real_t*        eA2)
{
   Vect<real_t> eu(&el,_u[_ind-1]);
   size_t k=0, n=el.getNbNodes()*el.getNbDOF();
   real_t d=1.0/_time_step;
   if (_step==1) {
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++) {
            eA0[k] += d*eA1[k];
            eb[i] += d*eA1[k]*eu[j];
         }
      }
   }
   else {
      Vect<real_t> ev(&el,_v[_ind-1]);
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++) {
            eA0[k] += 1.5*d*eA1[k];
            eb[i] += d*eA1[k]*(2*ev[j] - 0.5*eu[j]);
         }
      }
   }
   if (_set_bc[_ind-1])
      update_bc(el,Vect<real_t>(&el,*_bc[_ind-1]),eA0,eb);
   element_assembly(el,eA0,eb,_A[_ind-1],_b[_ind-1]);
}


void TimeStepping::SAssembleBDF2(const Side& sd,
                                 real_t*     sb,
                                 real_t*     sA)
{
   element_assembly(sd,sA,sb,_A[_ind-1],_b[_ind-1]);
}


Vect<real_t>& TimeStepping::setRHS_Stationary()
{
   return *_f2[_ind-1];
}


Vect<real_t>& TimeStepping::setRHS_ForwardEuler()
{
   return _f1[_ind-1];
}


Vect<real_t>& TimeStepping::setRHS_BackwardEuler()
{
   return *_f2[_ind-1];
}


Vect<real_t>& TimeStepping::setRHS_CrankNicolson()
{
   _f[_ind-1] = 0.5*(_f1[_ind-1] + (*_f2)[_ind-1]);
   return _f[_ind-1];
}


Vect<real_t>& TimeStepping::setRHS_Heun()
{
   if (_sstep==1)
      return _f1[_ind-1];
   else {
      _f[_ind-1] = 0.5*(_f1[_ind-1] + (*_f2)[_ind-1]);
      return _f[_ind-1];
   }
}


Vect<real_t>& TimeStepping::setRHS_AB2()
{
   if (_step==1)
      return _f1[_ind-1];
   else
      _f[_ind-1] = 1.5*_f1[_ind-1] - 0.5*_f0[_ind-1];
   return _f[_ind-1];
}


Vect<real_t>& TimeStepping::setRHS_LeapFrog()
{
   return _f1[_ind-1];
}


Vect<real_t>& TimeStepping::setRHS_RK4()
{
   if (_sstep==1)
      return _f1[_ind-1];
   else if (_sstep==2)
      _f[_ind-1] = 0.5*(_f1[_ind-1]+(*_f2)[_ind-1]);
   else if (_sstep==3)
      _f[_ind-1] = 0.5*(_f1[_ind-1]+(*_f2)[_ind-1]);
   else if (_sstep==4)
      return *_f2[_ind-1];
   return _f[_ind-1];
}


Vect<real_t>& TimeStepping::setRHS_RK3_TVD()
{
   return _f[_ind-1];
}


Vect<real_t>& TimeStepping::setRHS_Newmark()
{
   if (_step==1 && _sstep==1)
      return _f1[_ind-1];
   return *_f2[_ind-1];
}


Vect<real_t>& TimeStepping::setRHS_BDF2()
{
   return *_f2[_ind-1];
}


void TimeStepping::PreSolve_Newmark()
{
   if (_step==1 && _sstep==1)
      _D[_ind-1] = 0;
   if (_sstep>1) {
      _v[_ind-1] = _u[_ind-1] + _time_step*(*_du)[_ind-1] + (_time_step*_time_step*(0.5-_beta))*_ddu[_ind-1];
      *_du[_ind-1] += ((1-_gamma)*_time_step)*_ddu[_ind-1];
   }
}


void TimeStepping::insertBC(const Vect<real_t>& b,
                            Vect<real_t>&       v)
{
   if (_bc[_ind-1])
      v.insertBC(*_theMesh[_ind-1],b,*_bc[_ind-1]);
   else
      v = b;
}


void TimeStepping::insertBC0(const Vect<real_t>& b,
                             Vect<real_t>&       v)
{
   v.insertBC(*_theMesh[_ind-1],b);
}


ostream& operator<<(ostream&            s,
                    const TimeStepping& ts)
{
   string scheme;
   if (ts._sc==STATIONARY)
      scheme = "Stationary";
   else if (ts._sc==FORWARD_EULER)
      scheme = "Forward Euler";
   else if (ts._sc==BACKWARD_EULER)
      scheme = "Backward Euler";
   else if (ts._sc==CRANK_NICOLSON)
      scheme = "Crank-Nicolson";
   else if (ts._sc==HEUN)
      scheme = "Heun";
   else if (ts._sc==NEWMARK)
      scheme = "Newmark";
   else if (ts._sc==LEAP_FROG)
      scheme = "LeapFrog";
   else if (ts._sc==AB2)
      scheme = "Adams-Bashforth";
   else if (ts._sc==RK4)
      scheme = "Runge-Kutta 4";
   else if (ts._sc==RK3_TVD)
      scheme = "Runge-Kutta 3, TVD";
   else if (ts._sc==BDF2)
      scheme = "BDF 2";
   else
      ;
   if (ts._sc==STATIONARY) {
      s << "STATIONARY ANALYSIS." << endl;
      s << "Number of equations: \t" << ts._nb_eq << endl;
   }
   else {
      s << "\nTIME STEPPING FOR A TIME DEPENDENT PROBLEM\n\n";
      for (int i=0; i<ts._nb_pdes; ++i)
         s << "PDE " << i+1 << ", Number of equations: \t\t" << ts._nb_eq[i] << endl;
      s << "Time integration scheme: \t" << scheme << endl;
      s << "Final time: \t\t\t" << ts._final_time << endl;
      s << "First time step value:\t\t" << ts._time_step0 << endl;
      s << "Final time step value:\t\t" << ts._time_step << endl;
      s << "Number of time steps:\t\t" << ts._step << endl;
   }
   return s;
}

} /* namespace OFELI */
