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

                  Implementation of abstract class 'AbsEqua'

  ==============================================================================*/

#ifndef __ABS_EQUA_IMPL_H
#define __ABS_EQUA_IMPL_H

#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "OFELI_Config.h"
#include "equations/AbsEqua.h"
#include "mesh/Mesh.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "mesh/Material.h"
#include "solvers/LinearSolver.h"
#include "util/Gauss.h"
#include "io/Prescription.h"
#include "solvers/TimeStepping.h"
#include "solvers/EigenProblemSolver.h"
#include "linear_algebra/DMatrix_impl.h"
#include "linear_algebra/DSMatrix_impl.h"
#include "linear_algebra/SpMatrix_impl.h"
#include "linear_algebra/SkMatrix_impl.h"
#include "linear_algebra/SkSMatrix_impl.h"
#include "linear_algebra/BMatrix_impl.h"
#include "linear_algebra/TrMatrix_impl.h"

namespace OFELI {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

extern Material theMaterial;


template<class T_>
AbsEqua<T_>::AbsEqua()
            : _theMesh(nullptr), _solver(-1), _nb_fields(1), _eigen(false),
              _set_matrix(false), _set_solver(false), _analysis(STATIONARY),
              _A(nullptr), _b(nullptr), _u(nullptr),
              _bc(nullptr), _bf(nullptr), _sf(nullptr), _pf(nullptr), _v(nullptr)
{
   setTimeIntegrationParam();
}


template<class T_>
AbsEqua<T_>::~AbsEqua()
{
   if (_A!=nullptr)
      delete _A;
   if (_b!=nullptr)
      delete _b;
}


template<class T_>
bool AbsEqua<T_>::SolverIsSet() const
{
   return _set_solver;
}


template<class T_>
void AbsEqua<T_>::setMesh(Mesh &m)
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
}


template<class T_>
Mesh& AbsEqua<T_>::getMesh() const
{
   return *_theMesh;
}


template<class T_>
void AbsEqua<T_>::build(EigenProblemSolver& e) { }


template<class T_>
void AbsEqua<T_>::build(TimeStepping& s) { }


template<class T_>
int AbsEqua<T_>::runOneTimeStep()
{
   return run(TRANSIENT_ONE_STEP);
}


template<class T_>
int AbsEqua<T_>::runSteadyState()
{
   return run(STATIONARY);
}


template<class T_>
int AbsEqua<T_>::runTransient()
{
   return run(TRANSIENT);
}


template<class T_>
int AbsEqua<T_>::run(Analysis   a,
                     TimeScheme s)
{
   _analysis = a;
   int ret=0;
   if (_u==nullptr)
      throw OFELIException("In AbsEqua<T>::run(Analysis,TimeScheme): No solution vector provided.");
   if (_b==nullptr)
      _b = new Vect<T_>(_nb_eq);
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


template<class T_>
void AbsEqua<T_>::getTangent(Matrix<T_>* Df)
{
   _Df = Df;
}


template<class T_>
void AbsEqua<T_>::setTransient(TimeScheme s,
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


template<class T_>
LinearSolver<T_>& AbsEqua<T_>::getLinearSolver()
{
   return _ls;
}


#if defined (USE_PETSC)
template<class T_>
PETScMatrix<T_>* AbsEqua<T_>::getMatrix() const
{
   return _A;
}
#else
template<class T_>
Matrix<T_>* AbsEqua<T_>::getMatrix() const
{
   return _A;
}
#endif


template<class T_>
void AbsEqua<T_>::setSolver(Iteration      ls,
                            Preconditioner pc)
{
   if (ls==DIRECT_SOLVER) {
   }
   if (ls==DIRECT_SOLVER && _matrix_type==SPARSE)
      throw OFELIException("In AbsEqua::setSolver(Iteration,Preconditioner): "
                           "Choices of solver and storage modes are incompatible.");
   if (!_set_matrix)
      setMatrixType(SPARSE);
   _ls.setSolver(ls,pc);
   _set_solver = true;
}


template<class T_>
void AbsEqua<T_>::setLinearSolver(Iteration      ls,
                                  Preconditioner pc)
{
   setSolver(ls,pc);
}


template<class T_>
void AbsEqua<T_>::setMatrixType(int t)
{
   _matrix_type = t;
   if (_A!=nullptr)
      delete _A, _A = nullptr;
#if defined(USE_PETSC)
   _A = new PETScMatrix<T_>(*_theMesh);
#else 
   if (_matrix_type&SPARSE)
      _A = new SpMatrix<T_>(*_theMesh);
   else if (_matrix_type&(DENSE&(DENSE|SYMMETRIC)))
      _A = new DSMatrix<T_>(*_theMesh);
   else if (_matrix_type&DENSE)
      _A = new DMatrix<T_>(*_theMesh);
   else if (_matrix_type&(SKYLINE&(SKYLINE|SYMMETRIC)))
      _A = new SkSMatrix<T_>(*_theMesh);
   else if (_matrix_type&SKYLINE)
      _A = new SkMatrix<T_>(*_theMesh);
   else if (_matrix_type&TRIDIAGONAL)
      _A = new TrMatrix<T_>(_theMesh->getNbEq());
   else if (_matrix_type&BAND)
      _A = new BMatrix<T_>;
#endif
   _ls.setMatrix(_A);
   _set_matrix = true;
}


template<class T_>
void AbsEqua<T_>::setWithConvection(int f)
{
   f = 0;
}


template<class T_>
void AbsEqua<T_>::setTerms(PDE_Terms t)
{
   _terms = t;
}


template<class T_>
#if defined(USE_PETSC)
int AbsEqua<T_>::solveLinearSystem(PETScMatrix<T_>* A,
                                   PETScVect<T_>&   b,
                                   PETScVect<T_>&   x)
#else
int AbsEqua<T_>::solveLinearSystem(Matrix<T_>* A,
                                   Vect<T_>&   b,
                                   Vect<T_>&   x)
#endif
{
   _ls.setMatrix(A);
   _ls.setRHS(b);
   _ls.setSolution(x);
   return _ls.solve();
}


template<class T_>
#if defined(USE_PETSC)
int AbsEqua<T_>::solveLinearSystem(PETScVect<T_>& b,
                                   PETScVect<T_>& x)
#else
int AbsEqua<T_>::solveLinearSystem(Vect<T_>& b,
                                   Vect<T_>& x)
#endif
{
   return solveLinearSystem(_A,b,x);
}


template<class T_>
void AbsEqua<T_>::setAnalysis(Analysis a)
{
   _analysis = a;
}


template<class T_>
void AbsEqua<T_>::setTimeIntegrationParam()
{
   _TimeInt.step = 0;
   _TimeInt.init = 0.;
   _TimeInt.delta = theTimeStep;
   _TimeInt.time = theTime;
   _TimeInt.final = theFinalTime;
}
    

template<class T_>
void AbsEqua<T_>::setTimeIndex(size_t step)
{
   _TimeInt.step = step;
}


template<class T_>
void AbsEqua<T_>::setInitTime(real_t t)
{
   _TimeInt.init = t;
}


template<class T_>
void AbsEqua<T_>::setTimeStep(real_t t)
{
   _TimeInt.delta = t;
}


template<class T_>
real_t AbsEqua<T_>::getTimeStep() const
{
   return _TimeInt.delta;
}


template<class T_>
void AbsEqua<T_>::setTime(real_t t)
{
   _TimeInt.time = t;
}


template<class T_>
void AbsEqua<T_>::setFinalTime(real_t t)
{
   _TimeInt.final = t;
}


template<class T_>
size_t AbsEqua<T_>::getNbFields() const
{
   return _nb_fields;
}


template<class T_>
void AbsEqua<T_>::setTimeIntegration(TimeScheme s)
{
   _TimeInt.scheme = s;
}


template<class T_>
TimeScheme AbsEqua<T_>::getTimeIntegration() const
{
   return _TimeInt.scheme;
}


template<class T_>
string AbsEqua<T_>::getEquationName() const
{
   return _equation_name;
}


template<class T_>
string AbsEqua<T_>::getFiniteElementType() const
{
   return _finite_element;
}


template<class T_>
#if defined(USE_PETSC)
void AbsEqua<T_>::setInput(EqDataType     opt,
                           PETScVect<T_>& u)
#else
void AbsEqua<T_>::setInput(EqDataType opt,
                           Vect<T_>&  u)
#endif
{
   if (opt==INITIAL_FIELD || opt==SOLUTION)
      _u = &u;
   else if (opt==BOUNDARY_CONDITION)
      _bc = &u;
   else if (opt==SOURCE || opt==BODY_FORCE)
      _bf = &u;
   else if (opt==FLUX || opt==TRACTION || opt==BOUNDARY_FORCE) {
      _sf = &u;
   }
   else if (opt==POINT_FORCE)
      _pf = &u;
   else
      ;
}


template<class T_>
void AbsEqua<T_>::set(Prescription& p)
{
   _prescription = &p;
}


template<class T_>
void AbsEqua<T_>::setTolerance(real_t toler)
{
   _toler = toler;
}


template<class T_>
#if defined(USE_PETSC)
void AbsEqua<T_>::setMatrix(PETScMatrix<T_> &A)
{
   _A = &A;
   _matrix_type = SKYLINE|SYMMETRIC;
   _A->setMesh(*_theMesh);
}
#else
void AbsEqua<T_>::setMatrix(SkSMatrix<T_> &A)
{
   _A = &A;
   _matrix_type = SKYLINE|SYMMETRIC;
   _A->setMesh(*_theMesh);
}


template<class T_>
void AbsEqua<T_>::setMatrix(SkMatrix<T_> &A)
{
   _A = &A;
   _matrix_type = SKYLINE;
   _A->setMesh(*_theMesh);
}


template<class T_>
void AbsEqua<T_>::setMatrix(SpMatrix<T_> &A)
{
   _A = &A;
   _matrix_type = SPARSE;
   _A->setMesh(*_theMesh);
}
#endif


template<class T_>
bool AbsEqua<T_>::isConstantMatrix() const
{
   return _constant_matrix;
}


template<class T_>
bool AbsEqua<T_>::isConstantMesh() const
{
   return _constant_mesh;
}


template<class T_>
void AbsEqua<T_>::setConstantMatrix()
{
   _constant_matrix = true;
}


template<class T_>
void AbsEqua<T_>::setConstantMesh()
{
   _constant_mesh = true;
}


template<class T_>
void AbsEqua<T_>::setTerms(int opt)
{
   _terms = opt;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */

#endif
