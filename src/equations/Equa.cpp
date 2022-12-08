/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani

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

                  Implementation of abstract class 'Equa'

  ==============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "OFELI_Config.h"
#include "equations/Equa.h"
#include "mesh/Mesh.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "mesh/Material.h"
#include "util/Gauss.h"
#include "io/Prescription.h"
#include "linear_algebra/DMatrix_impl.h"
#include "linear_algebra/DSMatrix_impl.h"
#include "linear_algebra/SpMatrix_impl.h"
#include "linear_algebra/SkMatrix_impl.h"
#include "linear_algebra/SkSMatrix_impl.h"
#include "linear_algebra/BMatrix_impl.h"
#include "linear_algebra/TrMatrix_impl.h"
#include "io/Fct.h"


namespace OFELI {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

extern Material theMaterial;
class EigenProblemSolver;
class TimeStepping;

Equa::Equa()
     : _theMesh(nullptr), _solver(-1), _nb_fields(1), _eigen(false),
       _set_matrix(false), _set_solver(false), _analysis(STATIONARY),
       _A(nullptr), _b(nullptr), _u(nullptr), _bc(nullptr), _bf(nullptr),
       _sf(nullptr), _pf(nullptr), _v(nullptr)
{
   setTimeIntegrationParam();
   _rho_set = _Cp_set = _kappa_set = _mu_set = _sigma_set = _Mu_set = false;
   _epsilon_set = _omega_set = _beta_set = _v_set = _young_set = _poisson_set = false;
}


Equa::~Equa()
{
   if (_A!=nullptr)
      delete _A;
   if (_b!=nullptr)
      delete _b;
}


void Equa::set_rho(string exp)
{
   _rho_fct.set(exp);
   _rho_set = true;
}


void Equa::set_Cp(string exp)
{
   _Cp_fct.set(exp);
   _Cp_set = true;
}


void Equa::set_kappa(string exp)
{
   _kappa_fct.set(exp);
   _kappa_set = true;
}


void Equa::set_mu(string exp)
{
   _mu_fct.set(exp);
   _mu_set = true;
}


void Equa::set_sigma(string exp)
{
   _sigma_fct.set(exp);
   _sigma_set = true;
}


void Equa::set_Mu(string exp)
{
   _Mu_fct.set(exp);
   _Mu_set = true;
}


void Equa::set_epsilon(string exp)
{
   _epsilon_fct.set(exp);
   _epsilon_set = true;
}


void Equa::set_omega(string exp)
{
   _omega_fct.set(exp);
   _omega_set = true;
}


void Equa::set_beta(string exp)
{
   _beta_fct.set(exp);
   _beta_set = true;
}


void Equa::set_v(string exp)
{
   _v_fct.set(exp);
   _v_set = true;
}


void Equa::set_young(string exp)
{
   _young_fct.set(exp);
   _young_set = true;
}


void Equa::set_poisson(string exp)
{
   _poisson_fct.set(exp);
   _poisson_set = true;
}


bool Equa::SolverIsSet() const
{
   return _set_solver;
}


void Equa::setMesh(Mesh &m)
{
   _theMesh = &m;
   _theMesh->removeImposedDOF();
   _nb_nodes = _theMesh->getNbNodes();
   _nb_sides = _theMesh->getNbSides();
   _nb_boundary_sides = _theMesh->getNbBoundarySides();
   _nb_el = _theMesh->getNbElements();
   _nb_dof_total = _theMesh->getNbDOF();
   _nb_eq = _theMesh->getNbEq();
   setMatrixType(SPARSE|SYMMETRIC);
   setSolver(CG_SOLVER,DILU_PREC);
   _its = CG_SOLVER;
   _prec = DILU_PREC;
}


Mesh& Equa::getMesh() const
{
   return *_theMesh;
}


void Equa::build(EigenProblemSolver& e) { }


void Equa::build(TimeStepping& s) { }


int Equa::runOneTimeStep()
{
   return run(TRANSIENT_ONE_STEP);
}


int Equa::runSteadyState()
{
   return run(STATIONARY);
}


int Equa::runTransient()
{
   return run(TRANSIENT);
}


int Equa::run(Analysis   a,
              TimeScheme s)
{
   _analysis = a;
   int ret=0;
   if (_u==nullptr)
      throw OFELIException("In Equa<T>::run(Analysis,TimeScheme): No solution vector provided.");
   if (_b==nullptr)
      _b = new Vect<real_t>(_nb_eq);
   _b->clear();
   _uu.setSize(_nb_eq);
   if (a==STATIONARY) {
      build();
      ret = solveLinearSystem(*_b,_uu);
      if (_bc!=nullptr)
         _u->insertBC(*_theMesh,_uu,*_bc);
      else
         *_u = _uu;
   }
   else if (a==TRANSIENT) {
      for (_TimeInt.time=_TimeInt.init;
           _TimeInt.time<=_TimeInt.final;
           _TimeInt.time+=_TimeInt.delta, _TimeInt.step++) {
         build();
         ret = solveLinearSystem(*_b,_uu);
         if (_bc!=nullptr)
            _u->insertBC(*_theMesh,_uu,*_bc);
         else
            *_u = _uu;
      }
   }
   else if (a==TRANSIENT_ONE_STEP) {
      build();
      ret = solveLinearSystem(*_b,_uu);
      if (_bc!=nullptr)
         _u->insertBC(*_theMesh,_uu,*_bc);
      else
         *_u = _uu;
   }
   return ret;
}


void Equa::getTangent(Matrix<real_t>* Df)
{
   _Df = Df;
}


void Equa::setTransient(TimeScheme s,
                        real_t     p1,
                        real_t     p2,
                        real_t     p3)
{
   _TimeInt.scheme = s;
   _TimeInt.time_parameter1 = p1;
   _TimeInt.time_parameter2 = p2;
   _TimeInt.time_parameter3 = p3;
   if (s==FORWARD_EULER)
      _TimeInt.theta = 0;
   else if (s==BACKWARD_EULER)
      _TimeInt.theta = 1;
   else if (s==CRANK_NICOLSON)
      _TimeInt.theta = 0.5;
}


LinearSolver& Equa::getLinearSolver()
{
   return _ls;
}


#if defined (USE_PETSC)
PETScMatrix<real_t>* Equa::getMatrix() const
{
   return _A;
}
#else
Matrix<real_t>* Equa::getMatrix() const
{
   return _A;
}
#endif


void Equa::setSolver(Iteration      ls,
                     Preconditioner pc)
{
   _its = ls;
   _prec = pc;
   if (_matrix_type==TRIDIAGONAL) {
      _its = DIRECT_SOLVER;
      setMatrixType(TRIDIAGONAL);
      _ls.setSolver(_its,_prec);
      _set_solver = true;
      return;
   }
   if (_its==DIRECT_SOLVER)
      setMatrixType(SKYLINE);
   else
      setMatrixType(SPARSE);
   _ls.setSolver(_its,_prec);
   _set_solver = true;
}


void Equa::setMatrixType(int t)
{
   _matrix_type = t;
   if (_A!=nullptr)
      delete _A, _A = nullptr;
#if defined(USE_PETSC)
   _A = new PETScMatrix<real_t>(*_theMesh);
#else 
   if (_matrix_type&SPARSE)
      _A = new SpMatrix<real_t>(*_theMesh);
   else if (_matrix_type&(DENSE&(DENSE|SYMMETRIC)))
      _A = new DSMatrix<real_t>(*_theMesh);
   else if (_matrix_type&DENSE)
      _A = new DMatrix<real_t>(*_theMesh);
   else if (_matrix_type&(SKYLINE&(SKYLINE|SYMMETRIC)))
      _A = new SkSMatrix<real_t>(*_theMesh);
   else if (_matrix_type&SKYLINE)
      _A = new SkMatrix<real_t>(*_theMesh);
   else if (_matrix_type&TRIDIAGONAL)
      _A = new TrMatrix<real_t>(_theMesh->getNbEq());
   else if (_matrix_type&BAND)
      _A = new BMatrix<real_t>;
#endif
   _ls.setMatrix(_A);
   _set_matrix = true;
}


void Equa::setWithConvection(int f)
{
   f = 0;
}


void Equa::setTerms(PDE_Terms t)
{
   _terms = t;
}


#if defined(USE_PETSC)
int Equa::solveLinearSystem(PETScMatrix<real_t>* A,
                            PETScVect<real_t>&   b,
                            PETScVect<real_t>&   x)
#else
int Equa::solveLinearSystem(Matrix<real_t>* A,
                            Vect<real_t>&   b,
                            Vect<real_t>&   x)
#endif
{
   if (_matrix_type&(SPARSE&(SPARSE|SYMMETRIC)) && _its==DIRECT_SOLVER)
      throw OFELIException("In Equa::solveLinearSystem(A,b,x): "
                           "Linear solver and matrix storage format are incompatible.");
   _ls.setMatrix(A);
   _ls.setRHS(b);
   _ls.setSolution(x);
   int ret = _ls.solve();
   return ret;
}

void Equa::LinearSystemInfo()
{
   if (_matrix_type&(SKYLINE&(SKYLINE|SYMMETRIC)))
      cout << "Matrix is stored in skyline format." << endl;
   if (_its==DIRECT_SOLVER)
      cout << "Linear system solver: Direct solver" << endl;
   else if (_its==CG_SOLVER)
      cout << "Linear system solver: Conjugate Gradient method." << endl;
   else if (_its==CGS_SOLVER)
      cout << "Linear system solver: Squared Conjugate Gradient method." << endl;
   else if (_its==BICG_SOLVER)
      cout << "Linear system solver: Bi-Conjugate Gradient method." << endl;
   else if (_its==BICG_STAB_SOLVER)
      cout << "Linear system solver: Bi-Conjugate Stabilized Gradient method." << endl;
   else if (_its==GMRES_SOLVER)
      cout << "Linear system solver: GMRES method." << endl;
   if (_its==DIRECT_SOLVER)
      return;
   if (_prec==IDENT_PREC)
      cout << "Preconditioner: Identity." << endl;
   else if (_prec==DIAG_PREC)
      cout << "Preconditioner: Diagonal." << endl;
   else if (_prec==DILU_PREC)
      cout << "Preconditioner: Diagonal Incomplete LU factorization." << endl;
   else if (_prec==ILU_PREC)
      cout << "Preconditioner: Incomplete LU factorization." << endl;
   else if (_prec==SSOR_PREC)
      cout << "Preconditioner: Symmetric Successive Over Relaxation." << endl;
}


#if defined(USE_PETSC)
int Equa::solveLinearSystem(PETScVect<real_t>& b,
                            PETScVect<real_t>& x)
#else
int Equa::solveLinearSystem(Vect<real_t>& b,
                                   Vect<real_t>& x)
#endif
{
   return solveLinearSystem(_A,b,x);
}


void Equa::setAnalysis(Analysis a)
{
   _analysis = a;
}


void Equa::setTimeIntegrationParam()
{
   _TimeInt.step = 0;
   _TimeInt.init = 0.;
   _TimeInt.delta = theTimeStep;
   _TimeInt.time = theTime;
   _TimeInt.final = theFinalTime;
}
    

void Equa::setTimeIndex(size_t step)
{
   _TimeInt.step = step;
}


void Equa::setInitTime(real_t t)
{
   _TimeInt.init = t;
}


void Equa::setTimeStep(real_t t)
{
   _TimeInt.delta = t;
}


real_t Equa::getTimeStep() const
{
   return _TimeInt.delta;
}


void Equa::setTime(real_t t)
{
   _TimeInt.time = t;
}


void Equa::setFinalTime(real_t t)
{
   _TimeInt.final = t;
}


size_t Equa::getNbFields() const
{
   return _nb_fields;
}


void Equa::setTimeIntegration(TimeScheme s)
{
   _TimeInt.scheme = s;
}


TimeScheme Equa::getTimeIntegration() const
{
   return _TimeInt.scheme;
}


string Equa::getEquationName() const
{
   return _equation_name;
}


string Equa::getFiniteElementType() const
{
   return _finite_element;
}


#if defined(USE_PETSC)
void Equa::setInput(EqDataType     opt,
                    PETScVect<real_t>& u)
#else
void Equa::setInput(EqDataType opt,
                    Vect<real_t>&  u)
#endif
{
   if (opt==INITIAL_FIELD || opt==SOLUTION)
      _u = &u;
   else if (opt==BOUNDARY_CONDITION)
      _bc = &u;
   else if (opt==SOURCE || opt==BODY_FORCE)
      _bf = &u;
   else if (opt==FLUX || opt==TRACTION || opt==BOUNDARY_FORCE)
      _sf = &u;
   else if (opt==POINT_FORCE)
      _pf = &u;
   else
      ;
}


void Equa::set(Prescription& p)
{
   _prescription = &p;
}


void Equa::setTolerance(real_t toler)
{
   _toler = toler;
}


#if defined(USE_PETSC)
void Equa::setMatrix(PETScMatrix<real_t> &A)
{
   _A = &A;
   _matrix_type = SKYLINE|SYMMETRIC;
   _A->setMesh(*_theMesh);
}
#else
void Equa::setMatrix(SkSMatrix<real_t> &A)
{
   _A = &A;
   _matrix_type = SKYLINE|SYMMETRIC;
   _A->setMesh(*_theMesh);
}


void Equa::setMatrix(SkMatrix<real_t> &A)
{
   _A = &A;
   _matrix_type = SKYLINE;
   _A->setMesh(*_theMesh);
}


void Equa::setMatrix(SpMatrix<real_t> &A)
{
   _A = &A;
   _matrix_type = SPARSE;
   _A->setMesh(*_theMesh);
}
#endif


bool Equa::isConstantMatrix() const
{
   return _constant_matrix;
}


bool Equa::isConstantMesh() const
{
   return _constant_mesh;
}


void Equa::setConstantMatrix()
{
   _constant_matrix = true;
}


void Equa::setConstantMesh()
{
   _constant_mesh = true;
}


void Equa::setTerms(int opt)
{
   _terms = opt;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */
