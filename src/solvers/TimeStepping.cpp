/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2026 Rachid Touzani

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
#include "equations/Equation_impl.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/SkMatrix_impl.h"
#include "linear_algebra/Assembly_impl.h"
#include "OFELIException.h"


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
   &TimeStepping::solveBDF2,                // Backward Differentiation Formula
   &TimeStepping::solveBuiltIn              // A specific time scheme implemented in the equation class
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
              _order(0), _step(0), _rhs_ok(0), _nb_des(0), _ind(0), _regex(false), _time_step(0.1),
              _time(0.), _final_time(1.), _beta(0.25), _gamma(0.5), _nl_toler(1.e-6),
              _max_nl_it(100)
{
}


TimeStepping::TimeStepping(TimeScheme s,
                           real_t     time_step,
                           real_t     final_time) :
              _order(0), _step(0), _rhs_ok(0), _nb_des(0), _ind(0), _time_step(0.1), _time(0.),
              _final_time(1.), _beta(0.25), _gamma(0.5), _nl_toler(1.e-6), _max_nl_it(100)
{
   set(s,time_step,final_time);
}


TimeStepping::~TimeStepping()
{
   for (int e=0; e<_nb_des; ++e)
      if (_de[e].A != nullptr)
         delete _de[e].A;
}


void TimeStepping::set(TimeScheme s,
                       real_t     time_step,
                       real_t     final_time)
{
   _final_time = final_time;
   _time_step = _time_step0 = time_step;
   _nb_ssteps = 1;
   for (int e=0; e<_nb_des; ++e) {
      _de[e].expl = false;
      if (s&FORWARD_EULER || s&HEUN || s&RK4 || s&RK3_TVD || s&LEAP_FROG || s&AB2)
         _de[e].expl = true;
      _de[e].constant_matrix = false;
   }
   if (s&HEUN || s&NEWMARK)
      _nb_ssteps = 2;
   if (s==RK4)
      _nb_ssteps = 4;
   _sc = int(s);
   if (s&FORWARD_EULER || s&BACKWARD_EULER || s&CRANK_NICOLSON || s&HEUN || s&NEWMARK ||
       s&LEAP_FROG || s&ADAMS_BASHFORTH || s&AB2 || s&RUNGE_KUTTA || s&RK4 || s&RK3_TVD ||
       s&BDF2) {
      _presolve = PS[int(s)];
      _solve = TS[int(s)];
      _assemb = AS[int(s)];
      _set_rhs = RHS[int(s)];
   }
   else {
      _solve = nullptr;
      throw OFELIException("In TimeStepping::set(...): Time integration scheme not available.");
   }
}


void TimeStepping::setPDE(Equa& eq,
                          bool  nl)
{
   DE de;
   de.constant_matrix = false;
   de.expl = false;
   de.set_bc = false;
   de.eq = &eq;
   _ls = &(de.eq->getLinearSolver());
   de.itsolver = _ls->getSolver();
   de.prec = _ls->getPreconditioner();
   de.mesh = &(eq.getMesh());
   int nb_dof=de.mesh->getNbDOF();
   de.nb_eq = de.mesh->getNbEq();
   de.nl = nl;
   de.A = nullptr;
   if (_sc==int(FORWARD_EULER) || _sc==int(HEUN) || _sc==int(RK4) || _sc==int(RK3_TVD) || _sc==int(LEAP_FROG) || _sc==int(AB2)) {
      de.expl = true;
      de.D.setSize(de.nb_eq);
   }
   if (_sc==int(NEWMARK))
      de.D.setSize(de.nb_eq);
   de.constant_matrix = false;
   if (de.itsolver==DIRECT_SOLVER)
      de.A = new SkMatrix<real_t>(*de.mesh);
   else
      de.A = new SpMatrix<real_t>(*de.mesh);
   if (_sc==int(RK4)) {
      de.k1.setSize(nb_dof);
      de.k2.setSize(nb_dof);
      de.k3.setSize(nb_dof);
   }
   de.b = new Vect<real_t>(de.nb_eq);
   de.vv.setSize(de.nb_eq);
   de.f.setSize(nb_dof);
   if (!de.nl)
      _max_nl_it = 1;
   _nb_des++, _ind++;
   _de.push_back(de);
}


void TimeStepping::setRK4RHS(Vect<real_t>& f)
{
   _de[_nb_des-1].f01 = &f;
}


void TimeStepping::setLinearSolver(Iteration      s,
                                   Preconditioner p)
{
   DE &de = _de[_nb_des-1];
   de.itsolver = s;
   de.prec = p;
   if (de.A != nullptr)
      delete de.A;
   if (de.itsolver==DIRECT_SOLVER)
      de.A = new SkMatrix<real_t>(*de.mesh);
   else
      de.A = new SpMatrix<real_t>(*de.mesh);
}


void TimeStepping::setInitial(Vect<real_t>& u)
{
   DE &de = _de[_nb_des-1];
   de.eq->setInitial(u);
   de.v.setSize(u.getNx(),u.getNy(),u.getNz());
   de.u.setSize(u.getNx(),u.getNy(),u.getNz());
   de.w = &u;
   de.u = *de.w;
   de.b = new Vect<real_t>(de.nb_eq);
   if (de.f0.size()==0)
      de.f0.setSize(u.getNx(),u.getNy(),u.getNz());
   de.f1.setSize(u.getNx(),u.getNy(),u.getNz());
   _order = 1;
}


void TimeStepping::setInitial(Vect<real_t>& u,
                              Vect<real_t>& v)
{
   setInitial(u);
   DE &de = _de[_nb_des-1];
   de.du = &v;
   de.ddu.setSize(u.getNx(),u.getNy());
   _order = 2;
}


void TimeStepping::setInitialRHS(Vect<real_t>& f)
{
   DE &de = _de[_nb_des-1];
   if (de.f0.size()==0)
      de.f0.setSize(f.getNx(),f.getNy(),f.getNz());
   de.f0 = f;
}


void TimeStepping::setRHS(Vect<real_t>& f)
{
   _rhs_ok++;
   DE &de = _de[_nb_des-1];
   de.f2 = &f;
   if (de.f0.size()==0)
      de.f0.setSize(f.getNx(),f.getNy(),f.getNz());
   if (de.f1.size()==0)
      de.f1.setSize(f.getNx(),f.getNy(),f.getNz());
}


void TimeStepping::setRHS(string exp)
{
   _ff.setMesh(*_de[_nb_des-1].mesh);
   _ff.setTime(theTime);
   _ff.set(exp);
   setRHS(_ff);
}


void TimeStepping::setBC(Vect<real_t>& u)
{
   _de[_nb_des-1].bc = &u;
   _de[_nb_des-1].set_bc = true;
}


void TimeStepping::setBC(int    code,
			 string exp)
{
   _fbc.setMesh(*_de[_nb_des-1].mesh);
   _fbc.setTime(theTime);
   _fbc.setNodeBC(code,exp);
   setBC(_fbc);
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
   Vect<real_t> f;
   _step++;
   _time = theTime;
   if (Verbosity>1)
      cout << "   Time step: " << _step << ", time: " << _time << " ..." << endl;

   for (int e=0; e<_nb_des; ++e) {
      DE &de = _de[e];
      if (_rhs_ok==0) {
         f.resize(de.mesh->getNbDOF());
         f = 0.;
         de.f2 = &f;
      }
      if (_sc==int(BUILTIN)) {
         if (de.bc != nullptr)
            de.eq->setDirichlet(*de.bc);
         de.eq->runOneTimeStep();
         break;
      }
      if (de.bc != nullptr)
         de.eq->setDirichlet(*de.bc);
      if (de.nl) {
         for (_sstep=1; _sstep<=_nb_ssteps; _sstep++) {
            int it=1;
            real_t err=1;
            while (it<=_max_nl_it && err>_nl_toler) {
               de.eq->setBodyForce((this->*_set_rhs)());
               if (de.expl)
                  de.D = 0;
               else
                  *de.A = 0;
               *de.b = 0;
               if (_sc==int(RK4))
                  de.k1 = 0, de.k2 = 0, de.k3 = 0;
               (this->*_presolve)();
               de.eq->build(*this);
               (this->*_solve)();
               it++;
            }
         }
      }
      else {
         for (_sstep=1; _sstep<=_nb_ssteps; _sstep++) {
            de.eq->setBodyForce((this->*_set_rhs)());
            if (de.expl)
               de.D = 0;
            else
               *de.A = 0;
            *de.b = 0;
            if (_sc==int(RK4))
               de.k1 = 0, de.k2 = 0, de.k3 = 0;
            (this->*_presolve)();
            de.eq->build(*this);
            (this->*_solve)();
         }
      }
   }
   theTimeStep = _time_step;
   return _time_step;
}


void TimeStepping::run(bool opt)
{
   if (Verbosity>0)
      cout << "Running time integration ..." << std::endl;
   for (int e=0; e<_nb_des; ++e)
      _de[e].constant_matrix = opt;
   TimeLoop {
      runOneTimeStep();
      theTimeStep = _time_step;
   }
}


void TimeStepping::solveStationary()
{
   for (int e=0; e<_nb_des; ++e) {
      DE &de = _de[e];
      _ls->setMatrix(de.A);
      _ls->setRHS(*de.b);
      _ls->setSolution(de.vv);
      _ls->solve(de.itsolver,de.prec);
      insertBC(de.vv,de.v);
      *de.w = de.u = de.v;
   }
}


void TimeStepping::solveForwardEuler()
{
   if (_order==2)
      throw OFELIException("In TimeStepping::solveForwardEuler(): Forward Euler "
                           "scheme is implemented for first order equations only.");
   for (int e=0; e<_nb_des; ++e) {
      DE &de = _de[e];
      for (int i=1; i<=de.nb_eq; ++i)
         (*de.b)(i) /= de.D(i);
      insertBC(*de.b,de.v);
      *de.w = de.u = de.v;
      de.f1 = *de.f2;
   }
}


void TimeStepping::solveBackwardEuler()
{
   if (_order==2)
      throw OFELIException("In TimeStepping::solveBackwardEuler(): Backward Euler "
                           "scheme is implemented for first order equations only.");
   DE &de = _de[_ind-1];
   _ls->setMatrix(de.A);
   _ls->setRHS(*de.b);
   _ls->setSolution(de.vv);
   _ls->solve(de.itsolver,de.prec);
   insertBC(_de[_ind-1].vv,de.v);
   *de.w = de.u = de.v;
}


void TimeStepping::solveCrankNicolson()
{
   if (_order==2)
      throw OFELIException("In TimeStepping::solveCrankNicolson(): "
                           "Crank-Nicolson scheme is implemented for first order equations only.");
   DE &de = _de[_ind-1];
   _ls->setMatrix(de.A);
   _ls->setRHS(*de.b);
   _ls->setSolution(de.vv);
   _ls->solve(de.itsolver,de.prec);
   insertBC(de.vv,de.v);
   *de.w = de.u = de.v;
   de.f1 = *de.f2;
}


void TimeStepping::solveHeun()
{
   if (_order==2)
      throw OFELIException("In TimeStepping::solveHeun(): Heun scheme is "
                           "implemented for first order equations only.");
   DE &de = _de[_ind-1];
   for (int i=1; i<=de.nb_eq; i++)
      (*de.b)(i) /= de.D(i);
   if (_sstep==1)
      insertBC(*de.b,de.v);
   if (_sstep==2) {
      insertBC(*de.b,*de.w);
      de.u = *de.w;
   }
   *de.w = de.v = de.u;
   de.f1 = *de.f2;
}


void TimeStepping::solveLeapFrog()
{
   if (_order==2)
      throw OFELIException("In TimeStepping::solveLeapFrog(): Leap Frog "
                           "scheme is implemented for first order equations only.");
   DE &de = _de[_ind-1];
   for (int i=1; i<=de.nb_eq; i++)
      (*de.b)(i) /= de.D(i);
   if (_step==1) {
      insertBC(*de.b,de.v);
      de.f1 = *de.f2;
   }
   else {
      insertBC(*de.b,*de.w);
      de.u = de.v;
      de.v = *de.w;
      de.f0 = de.f1;
      de.f1 = *de.f2;
   }
}


void TimeStepping::solveAB2()
{
   if (_order==2)
      throw OFELIException("In TimeStepping::solveAB2(): Adams-Bashforth "
                           "scheme is implemented for first order equations only.");
   DE &de = _de[_ind-1];
   for (int i=1; i<=de.nb_eq; i++)
      (*de.b)(i) /= de.D(i);
   if (_step==1) {
      insertBC(*de.b,de.v);
      de.f1 = *de.f2;
   }
   else {
      insertBC(*de.b,*de.w);
      de.u = de.v; de.v = *de.w;
      de.f0 = de.f1;
      de.f1 = *de.f2;
   }
}


void TimeStepping::solveRK4()
{
   if (_order==2)
      throw OFELIException("In TimeStepping::solveRK4(): Runge-Kutta scheme "
                           "is valid for first order equations only.");
   DE &de = _de[_ind-1];
   if (_sstep==4) {
      for (int i=1; i<=de.nb_eq; i++)
         (*de.b)(i) /= de.D(i);
      insertBC(*de.b,de.v);
      *de.w = de.u = de.v;
      de.f0 = de.f1;
      de.f1 = *de.f2;
      de.f01 = nullptr;
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
   DE &de = _de[_ind-1];
   if (_step==1 && _sstep==1) {
      for (int i=1; i<=de.nb_eq; ++i)
         (*de.b)(i) /= de.D(i);
      insertBC0(*de.b,de.ddu);
      return;
   }
   if (_sstep>1) {
      _ls->setMatrix(de.A);
      _ls->setRHS(*de.b);
      _ls->setSolution(de.vv);
      _ls->solve(de.itsolver,de.prec);
      insertBC(de.vv,*de.w);
      de.ddu = (1./(_beta*_time_step*_time_step))*(*de.w - de.v);
      *de.du += (_time_step*_gamma)*de.ddu;
      de.u = *de.w;
   }
}


void TimeStepping::solveBDF2()
{
   if (_order==2)
      throw OFELIException("In TimeStepping::solveBDF2(): BDF2 scheme is "
                           "implemented for first order equations only.");
   DE &de = _de[_ind-1];
   _ls->setMatrix(de.A);
   _ls->setRHS(*de.b);
   _ls->setSolution(de.vv);
   _ls->solve(de.itsolver,de.prec);
   if (_step==1) {
      insertBC(de.vv,de.v);
      *de.w = de.v;
   }
   else {
      insertBC(de.vv,*de.w);
      de.u = de.v;
      de.v = *de.w;
   }
}


void TimeStepping::solveBuiltIn()
{
}


void TimeStepping::AssembleStationary(const Element& el,
                                      real_t*        eb,
                                      real_t*        eA0,
                                      real_t*        eA1,
                                      real_t*        eA2)
{
   DE &de = _de[_ind-1];
   if (de.set_bc)
      update_bc(el,Vect<real_t>(&el,*de.bc),eA0,eb);
   element_assembly(el,eA0,eb,de.A,*de.b);
}


void TimeStepping::SAssembleStationary(const Side& sd,
                                       real_t*     sb,
                                       real_t*     sA)
{
   element_assembly(sd,sA,sb,_de[_ind-1].A,*_de[_ind-1].b);
}


void TimeStepping::AssembleForwardEuler(const Element& el,
                                        real_t*        eb,
                                        real_t*        eA0,
                                        real_t*        eA1,
                                        real_t*        eA2)
{
   DE &de = _de[_ind-1];
   Vect<real_t> eu(&el,de.u);
   size_t k=0, n=el.getNbNodes()*el.getNbDOF();
   real_t d=1/_time_step;
   for (size_t i=0; i<n; i++) {
      for (size_t j=0; j<n; j++, k++) {
         eA1[k] *= d;
         eb[i] += (eA1[k]-eA0[k])*eu[j];
      }
   }
   if (de.set_bc)
      update_bc_diag(el,Vect<real_t>(&el,*de.bc),eA1,eb);
   element_assembly(el,eA1,eb,de.D,*de.b);
}


void TimeStepping::SAssembleForwardEuler(const Side& sd,
                                         real_t*     sb,
                                         real_t*     sA)
{
   Vect<real_t> su(&sd,_de[_ind-1].u);
   size_t k=0, n=sd.getNbNodes()*sd.getNbDOF();
   for (size_t i=0; i<n; i++) {
      for (size_t j=0; j<n; j++, k++)
         sb[i] -= sA[k]*su[j];
   }
   element_assembly(sd,sb,_de[_ind-1].b);
}


void TimeStepping::AssembleBackwardEuler(const Element& el,
                                         real_t*        eb,
                                         real_t*        eA0,
                                         real_t*        eA1,
                                         real_t*        eA2)
{
   DE &de = _de[_ind-1];
   Vect<real_t> eu(&el,de.u);
   size_t k=0, n=el.getNbNodes()*el.getNbDOF();
   real_t d=1.0/_time_step;
   for (size_t i=0; i<n; i++) {
      for (size_t j=0; j<n; j++, k++) {
         eA0[k] += d*eA1[k];
         eb[i] += d*eA1[k]*eu[j];
      }
   }
   if (de.set_bc)
      update_bc(el,Vect<real_t>(&el,*de.bc),eA0,eb);
   element_assembly(el,eA0,eb,de.A,*de.b);
}


void TimeStepping::SAssembleBackwardEuler(const Side& sd,
                                          real_t*     sb,
                                          real_t*     sA)
{
   DE &de = _de[_ind-1];
   Vect<real_t> su(&sd,de.u);
   size_t k=0, n=sd.getNbNodes()*sd.getNbDOF();
   real_t d=1.0/_time_step;
   for (size_t i=0; i<n; i++) {
      for (size_t j=0; j<n; j++, k++)
         sb[i] += d*sA[k]*su[j];
   }
   element_assembly(sd,sA,sb,de.A,*de.b);
}


void TimeStepping::AssembleCrankNicolson(const Element& el,
                                         real_t*        eb,
                                         real_t*        eA0,
                                         real_t*        eA1,
                                         real_t*        eA2)
{
   DE &de = _de[_ind-1];
   Vect<real_t> eu(&el,de.u);
   size_t k=0, n=el.getNbNodes()*el.getNbDOF();
   real_t d=1.0/_time_step;
   for (size_t i=0; i<n; i++) {
      for (size_t j=0; j<n; j++, k++) {
         real_t z=0.5*eA0[k];
         eA0[k] = d*eA1[k] + z;
         eb[i] += (d*eA1[k] - z)*eu[j];
      }
   }
   if (de.set_bc)
      update_bc(el,Vect<real_t>(&el,*de.bc),eA0,eb);
   element_assembly(el,eA0,eb,de.A,*de.b);
}


void TimeStepping::SAssembleCrankNicolson(const Side& sd,
                                          real_t*     sb,
                                          real_t*     sA)
{
   DE &de = _de[_ind-1];
   Vect<real_t> su(&sd,de.u);
   size_t k=0, n=sd.getNbNodes()*sd.getNbDOF();
   for (size_t i=0; i<n; i++) {
      for (size_t j=0; j<n; j++, k++) {
         sA[k] *= 0.5;
         sb[i] -= sA[k]*su[j];
      }
   }
   element_assembly(sd,sA,sb,de.A,*de.b);
}


void TimeStepping::AssembleHeun(const Element& el,
                                real_t*        eb,
                                real_t*        eA0,
                                real_t*        eA1,
                                real_t*        eA2)
{
   DE &de = _de[_ind-1];
   Vect<real_t> eu(&el,de.u);
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
      Vect<real_t> ev(&el,de.v);
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++) {
            eA1[k] *= d;
            eb[i] += (eA1[k] - 0.5*eA0[k])*eu[j] - 0.5*eA0[k]*ev[j];
         }
      }
   }
   if (de.set_bc)
      update_bc_diag(el,Vect<real_t>(&el,*de.bc),eA1,eb);
   element_assembly(el,eA1,eb,de.D,*de.b);
}


void TimeStepping::SAssembleHeun(const Side& sd,
                                 real_t*     sb,
                                 real_t*     sA)
{
   DE &de = _de[_ind-1];
   Vect<real_t> su(&sd,de.u);
   size_t k=0, n=sd.getNbNodes()*sd.getNbDOF();
   if (_sstep==1) {
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++)
            sb[i] -= sA[k]*su[j];
      }
   }
   if (_sstep==2) {
      Vect<real_t> sv(&sd,de.v);
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++)
            sb[i] -=  0.5*sA[k]*(su[j] + sv[j]);
      }
   }
   Element_Assembly(sd,sb,*de.b);
}


void TimeStepping::AssembleAB2(const Element& el,
                               real_t*        eb,
                               real_t*        eA0,
                               real_t*        eA1,
                               real_t*        eA2)
{
   DE &de = _de[_ind-1];
   Vect<real_t> eu(&el,de.u);
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
      Vect<real_t> ev(&el,de.v);
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++) {
            eA1[k] *= d;
            eb[i] += (eA1[k] - 1.5*eA0[k])*ev[j] + 0.5*eA0[k]*eu[j];
         }
      }
   }
   if (de.set_bc)
      update_bc_diag(el,Vect<real_t>(&el,*de.bc),eA1,eb);
   element_assembly(el,eA1,eb,de.D,*de.b);
}


void TimeStepping::SAssembleAB2(const Side& sd,
                                real_t*     sb,
                                real_t*     sA)
{
   DE &de = _de[_ind-1];
   Vect<real_t> su(&sd,de.u);
   size_t k=0, n=sd.getNbNodes()*sd.getNbDOF();
   if (_step==1) {
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++)
            sb[i] -= sA[k]*su[j];
      }
   }
   else {
      Vect<real_t> sv(&sd,de.v);
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++)
            sb[i] -= 1.5*sA[k]*(1.5*sv[j] - 0.5*su[j]);
      }
   }
   Element_Assembly(sd,sb,*de.b);
}


void TimeStepping::AssembleLeapFrog(const Element& el,
                                    real_t*        eb,
                                    real_t*        eA0,
                                    real_t*        eA1,
                                    real_t*        eA2)
{
   DE &de = _de[_ind-1];
   Vect<real_t> eu(&el,de.u);
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
      Vect<real_t> ev(&el,de.v);
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++) {
            eA1[k] *= d;
            eb[i] += eA1[k]*eu[j] - eA0[k]*ev[j];
         }
      }
   }
   if (de.set_bc)
      update_bc_diag(el,Vect<real_t>(&el,*de.bc),eA1,eb);
   element_assembly(el,eA1,eb,de.D,*de.b);
}


void TimeStepping::SAssembleLeapFrog(const Side& sd,
                                     real_t*     sb,
                                     real_t*     sA)
{
   DE &de = _de[_ind-1];
   Vect<real_t> su(&sd,de.u);
   size_t k=0, n=sd.getNbNodes()*sd.getNbDOF();
   if (_step==1) {
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++)
            sb[i] -= sA[k]*su[j];
      }
   }
   else {
      Vect<real_t> sv(&sd,de.v);
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++)
            sb[i] -= sA[k]*sv[j];
      }
   }
   Element_Assembly(sd,sb,*de.b);
}


void TimeStepping::AssembleRK4(const Element& el,
                               real_t*        eb,
                               real_t*        eA0,
                               real_t*        eA1,
                               real_t*        eA2)
{
   DE &de = _de[_ind-1];
   Vect<real_t> eu(&el,de.u);
   size_t m=0, n=el.getNbNodes()*el.getNbDOF();
   switch (_sstep) {

      case 1:
         {
            for (size_t i=0; i<n; i++) {
               for (size_t j=0; j<n; j++, m++)
                  eb[i] -= eA0[m]*eu[j];
            }
            element_assembly(el,eb,de.k1);
         }
         break;

      case 2:
         {
            Vect<real_t> ek(&el,de.k1);
            for (size_t i=0; i<n; i++) {
               for (size_t j=0; j<n; j++, m++)
                  eb[i] -= eA0[m]*(eu[j] + 0.5*_time_step*ek[j]);
            }
            element_assembly(el,eb,de.k2);
         }
         break;

      case 3:
         {
            Vect<real_t> ek(&el,de.k2);
            for (size_t i=0; i<n; i++) {
               for (size_t j=0; j<n; j++, m++)
                  eb[i] -= eA0[m]*(eu[j] + 0.5*_time_step*ek[j]);
            }
            element_assembly(el,eb,de.k3);
         }
         break;

      case 4:
         {
            Vect<real_t> ek1(&el,de.k1), ek2(&el,de.k2), ek3(&el,de.k3);
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
            if (de.set_bc)
               update_bc_diag(el,Vect<real_t>(&el,*de.bc),eA1,eb);
            element_assembly(el,eA1,eb,de.D,*de.b);
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
   DE &de = _de[_ind-1];
   Vect<real_t> su(&sd,de.u);
   size_t m=0, n=sd.getNbNodes()*sd.getNbDOF();
   switch (_sstep) {

      case 1:
         {
            for (size_t i=0; i<n; i++) {
               for (size_t j=0; j<n; j++, m++)
                  sb[i] -= sA[m]*su[j];
            }
            Element_Assembly(sd,sb,de.k1);
         }
         break;

      case 2:
         {
            Vect<real_t> sk(&sd,de.k1);
            for (size_t i=0; i<n; i++) {
               for (size_t j=0; j<n; j++, m++)
                  sb[i] -= sA[m]*(su[j] + 0.5*_time_step*sk[j]);
            }
            Element_Assembly(sd,sb,de.k2);
         }
         break;

      case 3:
         {
            Vect<real_t> sk(&sd,de.k2);
            for (size_t i=0; i<n; i++) {
               for (size_t j=0; j<n; j++, m++)
                  sb[i] -= sA[m]*(su[j] + 0.5*_time_step*sk[j]);
            }
            Element_Assembly(sd,sb,de.k3);
         }
         break;

      case 4:
         {
            Vect<real_t> sk1(&sd,de.k1), sk2(&sd,de.k2), sk3(&sd,de.k3);
            for (size_t i=0; i<n; i++) {
               for (size_t j=0; j<n; j++, m++)
                  sb[i] -= sA[m]*(su[j] + _time_step*sk3[j]);
            }
            Element_Assembly(sd,sb,*de.b);
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
   DE &de = _de[_ind-1];
   size_t n=el.getNbNodes()*el.getNbDOF();
   Vect<real_t> edu(&el,*de.du);
   if (_step==1 && _sstep==1) {
      Vect<real_t> eu(&el,de.u);
      size_t k=0;
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++)
            eb[i] -= eA0[k]*eu[j];
      }
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++)
            eb[i] -= eA1[k]*edu[j];
      }
      element_assembly(el,eA2,eb,de.D,*de.b);
   }
   if (_sstep>1) {
      Vect<real_t> ev(&el,de.v);
      real_t a=1./(_beta*_time_step*_time_step), b=_gamma/(_beta*_time_step);
      size_t k=0;
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++) {
            real_t z = a*eA2[k] + b*eA1[k];
            eA0[k] += z;
            eb[i] += z*ev[j] - eA1[k]*edu[j];
         }
      }
      if (de.set_bc)
         update_bc(el,Vect<real_t>(&el,*de.bc),eA0,eb);
      element_assembly(el,eA0,eb,de.A,*de.b);
   }
}


void TimeStepping::SAssembleNewmark(const Side& sd,
                                    real_t*     sb,
                                    real_t*     sA)
{
   DE &de = _de[_ind-1];
   size_t n=sd.getNbNodes()*sd.getNbDOF();
   Vect<real_t> sdu(&sd,*de.du);
   if (_step==1 && _sstep==1) {
      Vect<real_t> su(&sd,de.u);
      size_t k=0;
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++)
            sb[i] -= sA[k]*su[j];
      }
      Element_Assembly(sd,sb,*de.b);
   }
   if (_sstep>1)
      element_assembly(sd,sA,sb,de.A,*de.b);
}


void TimeStepping::AssembleBDF2(const Element& el,
                                real_t*        eb,
                                real_t*        eA0,
                                real_t*        eA1,
                                real_t*        eA2)
{
   DE &de = _de[_ind-1];
   Vect<real_t> eu(&el,de.u);
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
      Vect<real_t> ev(&el,de.v);
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<n; j++, k++) {
            eA0[k] += 1.5*d*eA1[k];
            eb[i] += d*eA1[k]*(2*ev[j] - 0.5*eu[j]);
         }
      }
   }
   if (de.set_bc)
      update_bc(el,Vect<real_t>(&el,*de.bc),eA0,eb);
   element_assembly(el,eA0,eb,de.A,*de.b);
}


void TimeStepping::SAssembleBDF2(const Side& sd,
                                 real_t*     sb,
                                 real_t*     sA)
{
   element_assembly(sd,sA,sb,_de[_ind-1].A,*_de[_ind-1].b);
}


Vect<real_t>& TimeStepping::setRHS_Stationary()
{
   return *_de[_ind-1].f2;
}


Vect<real_t>& TimeStepping::setRHS_ForwardEuler()
{
   return _de[_ind-1].f1;
}


Vect<real_t>& TimeStepping::setRHS_BackwardEuler()
{
   return *_de[_ind-1].f2;
}


Vect<real_t>& TimeStepping::setRHS_CrankNicolson()
{
   DE &de = _de[_ind-1];
   de.f = 0.5*(de.f1 + *de.f2);
   return de.f;
}


Vect<real_t>& TimeStepping::setRHS_Heun()
{
   DE &de = _de[_ind-1];
   if (_sstep==1)
      return de.f1;
   else {
      de.f = 0.5*(de.f1 + *de.f2);
      return de.f;
   }
}


Vect<real_t>& TimeStepping::setRHS_AB2()
{
   DE &de = _de[_ind-1];
   if (_step==1)
      return de.f1;
   else
      de.f = 1.5*de.f1 - 0.5*de.f0;
   return de.f;
}


Vect<real_t>& TimeStepping::setRHS_LeapFrog()
{
   return _de[_ind-1].f1;
}


Vect<real_t>& TimeStepping::setRHS_RK4()
{
   DE &de = _de[_ind-1];
   if (_sstep==1)
      return de.f1;
   else if (_sstep==2)
      de.f = 0.5*(de.f1 + *de.f2);
   else if (_sstep==3)
      de.f = 0.5*(de.f1 + *de.f2);
   else if (_sstep==4)
      return *de.f2;
   return de.f;
}


Vect<real_t>& TimeStepping::setRHS_RK3_TVD()
{
   return _de[_ind-1].f;
}


Vect<real_t>& TimeStepping::setRHS_Newmark()
{
   if (_step==1 && _sstep==1)
      return _de[_ind-1].f1;
   return *_de[_ind-1].f2;
}


Vect<real_t>& TimeStepping::setRHS_BDF2()
{
   return *_de[_ind-1].f2;
}


void TimeStepping::PreSolve_Newmark()
{
   DE &de = _de[_ind-1];
   if (_step==1 && _sstep==1)
      de.D = 0;
   if (_sstep>1) {
      if (de.du==nullptr)
         throw OFELIException("In TimeStepping::PreSolve_Newmark(): Incorrect initial data.");
      de.v = de.u + _time_step*(*de.du) + (_time_step*_time_step*(0.5-_beta))*de.ddu;
      *de.du += ((1-_gamma)*_time_step)*de.ddu;
   }
}


void TimeStepping::insertBC(const Vect<real_t>& b,
                            Vect<real_t>&       v)
{
   DE &de = _de[_ind-1];
   if (de.bc != nullptr)
      v.insertBC(*de.mesh,b,*de.bc);
   else
      v = b;
}


void TimeStepping::insertBC0(const Vect<real_t>& b,
                             Vect<real_t>&       v)
{
   v.insertBC(*_de[_ind-1].mesh,b);
}


std::ostream& operator<<(std::ostream& s,
                         TimeStepping& ts)
{
   if (ts._sc==STATIONARY) {
      s << "STATIONARY ANALYSIS." << std::endl;
      for (int i=0; i<ts._nb_des; ++i)
         s << "Number of equations: \t" << ts._de[i].nb_eq << std::endl;
   }
   else {
      s << "\nTIME STEPPING FOR A TIME DEPENDENT PROBLEM\n\n";
      for (int i=0; i<ts._nb_des; ++i)
         s << "PDE " << i+1 << ", Number of equations: \t\t" << ts._de[i].nb_eq << std::endl;
      s << "Time integration scheme: \t" << ts._scs[ts._sc] << std::endl;
      s << "Final time: \t\t\t" << ts._final_time << std::endl;
      s << "First time step value:\t\t" << ts._time_step0 << std::endl;
      s << "Final time step value:\t\t" << ts._time_step << std::endl;
      s << "Number of time steps:\t\t" << ts._step << std::endl;
   }
   return s;
}

} /* namespace OFELI */
