/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2017 Rachid Touzani

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

#include "equations/AbsEqua.h"
#include "solvers/EigenProblemSolver.h"

namespace OFELI {

EigenProblemSolver::EigenProblemSolver()
                   : _theEqua(NULL), _theMesh(NULL), _K(NULL), _M(NULL),
                     _K_alloc(false), _M_alloc(false), _lM_alloc(false), _diag(false),
                     _max_it(100), _dim(1), _nb_eq(0), _nb_eigv(1), _opt(1),
                     _epsv(1.e-8), _epsj(1.e-10)
{ }


EigenProblemSolver::EigenProblemSolver(DSMatrix<real_t>& A,
                                       int               n)
                   : _theEqua(NULL), _theMesh(NULL), _K(&A), _M(NULL),
                     _K_alloc(false), _M_alloc(false), _lM_alloc(true), _diag(true),
                     _max_it(100), _dim(1), _opt(1), _epsv(1.e-8), _epsj(1.e-10)
{
   _nb_eq = _nb_eigv = _K->getNbRows();
   if (n>0)
      _nb_eq = _nb_eigv = n;
   _lM = new Vect<real_t>(_nb_eq);
   *_lM = 1;
   run(_nb_eigv);
}


EigenProblemSolver::EigenProblemSolver(DSMatrix<real_t>& A,
                                       Vect<real_t>&     ev,
                                       int               n)
                   : _theEqua(NULL), _theMesh(NULL), _K(&A), _M(NULL),
                     _K_alloc(false), _M_alloc(false), _lM_alloc(true), _diag(true),
                     _max_it(100), _dim(1), _opt(1), _epsv(1.e-8), _epsj(1.e-10)
{
   _nb_eq = _nb_eigv = _K->getNbRows();
   if (n>0)
      _nb_eq = _nb_eigv = n;
   _lM = new Vect<real_t>(_nb_eq);
   *_lM = 1;
   run(n);
   ev.setSize(_nb_eigv);
   ev = _eigv;
}


EigenProblemSolver::EigenProblemSolver(AbsEqua<real_t>& eq,
                                       bool             lumped)
                   : _K(NULL), _M(NULL),
                     _K_alloc(true), _M_alloc(false), _lM_alloc(true), _diag(true),
                     _max_it(100), _dim(1), _nb_eq(0), _nb_eigv(1), _opt(1),
                     _epsv(1.e-8), _epsj(1.e-10)
{
   setPDE(eq,lumped);
}


EigenProblemSolver::EigenProblemSolver(SkSMatrix<real_t>& K,
                                       SkSMatrix<real_t>& M,
                                       int                n)
                   : _K_alloc(false), _M_alloc(false), _lM_alloc(false), _diag(false),
                     _max_it(100), _dim(1), _opt(1), _epsv(1.e-8), _epsj(1.e-10)
{
   setMatrix(K,M);
   if (n==0)
      _nb_eigv = n;
}


EigenProblemSolver::EigenProblemSolver(SkSMatrix<real_t>& K,
                                       Vect<real_t>&      M,
                                       int                n)
                   : _K_alloc(false), _M_alloc(false), _lM_alloc(false), _diag(true),
                     _max_it(100), _dim(1), _opt(1), _epsv(1.e-8), _epsj(1.e-10)
{
   setMatrix(K,M);
   if (n==0)
      _nb_eigv = n;
}


EigenProblemSolver::~EigenProblemSolver()
{
   if (_K_alloc)
      delete _K;
   if (_M_alloc)
      delete _M;
   if (_lM_alloc)
      delete _lM;
}


void EigenProblemSolver::setMatrix(SkSMatrix<real_t>& K,
                                   SkSMatrix<real_t>& M)
{
   _K = &K;
   _M = &M;
   _nb_eq = _K->getNbRows();
   _diag = false;
}


void EigenProblemSolver::setMatrix(SkSMatrix<real_t>& K,
                                   Vect<real_t>&      M)
{
   _K = &K;
   _lM = &M;
   _nb_eq = _K->getNbRows();
   _diag = true;
}


void EigenProblemSolver::setMatrix(DSMatrix<real_t>& K)
{
   _K = &K;
   _lM = new Vect<real_t>(_nb_eq);
   _lM_alloc = true;
   *_lM = 1;
   _diag = true;
}


void EigenProblemSolver::setPDE(AbsEqua<real_t>& eq,
                                bool             lumped)
{
   _theEqua = &eq;
   _theMesh = &(_theEqua->getMesh());
   _theMesh->removeImposedDOF();
   _theMesh->NumberEquations();
   _nb_eq = _theMesh->getNbEq();
   _K = new SkSMatrix<real_t>(*_theMesh);
   if (lumped)
      _lM = new Vect<real_t>(_nb_eq);
   else {
      _diag = _lM_alloc = false, _M_alloc = true;
      _M = new SkSMatrix<real_t>(_nb_eq);
   }
}


int EigenProblemSolver::run(int nb)
{
   if (_theEqua)
      _theEqua->build(*this);
   _nb_eigv = nb;
   if (nb==0)
      _nb_eigv = _nb_eq;
   return runSubSpace(_nb_eigv);
}


void EigenProblemSolver::Assembly(const Element& el,
                                  real_t*        eK,
                                  real_t*        eM)
{
   element_assembly(el,eK,_K);
   if (_diag)
      element_assembly(el,eM,_lM);
   else
      element_assembly(el,eM,_M);
}


void EigenProblemSolver::SAssembly(const Side& sd,
                                   real_t*     sK)
{
   element_assembly(sd,sK,_K);
}


int EigenProblemSolver::runSubSpace(size_t nb_eigv,
                                    size_t ss_dim)
{
   int ret=0;
   _nb_eigv = nb_eigv;
   if (_nb_eigv>_nb_eq)
      _nb_eigv = _nb_eq;
   _dim = ss_dim;
   if (_dim==0)
      _dim = _nb_eigv + 1;
   if (_dim > _nb_eq)
      _dim = _nb_eq;
   Vect<real_t> old_eigv(_dim), wv(_nb_eq), ww(_dim);
   DMatrix<real_t> wm(_dim);
   _pK.setSize(_dim);
   _pM.setSize(_dim);
   _ev.setSize(_dim,_nb_eq);
   _eigv.setSize(_dim);
   _rc.setSize(_dim);

// Is automatic generation required for initial subspace
   if (_opt)
      ret = init();
   else {
      for (size_t i=1; i<=size_t(_dim); i++) {
         Mxv(_ev.getRow(i),wv);
         _ev.setRow(i,wv);
      }
   }

// Initialization
   int err=0;

// Factorize K
   _K->Factor();

// Starting iterations
   size_t it = 0;
   while (it<_max_it) {
      it++;

//    Find projection of *_K and store in _pK
      for (size_t i=1; i<=_dim; i++) {
         _K->solve(_ev.getRow(i),wv);
         for (size_t j=i; j<=_dim; j++)
            _pK(i,j) = (_ev.getRow(j),wv);
         _ev.setRow(i,wv);
      }

//    Find projection of *_M and store in _pM
      for (size_t i=1; i<=_dim; i++) {
         Mxv(_ev.getRow(i),wv);
         for (size_t j=i; j<=size_t(_dim); j++)
            _pM(i,j) = (_ev.getRow(j),wv);
         if (err==0)
            _ev.setRow(i,wv);
      }

//    Solve eigensystem of subspace operators
      runJacobi(wm);

//    Arrange eigenvalues in ascending order
      int ii, jj, nnend, istp=int(_dim)-1;
      do {
         nnend=istp; istp=0; ii=1;
         for (size_t i=1; i<=size_t(nnend); i++) {
            jj = ii + int(_dim) - int(i) - 1;
            if (_eigv(i+1)<_eigv(i)) {
               istp = int(i);
               Swap(_eigv(i),_eigv(i+1));
               Swap(_pM(jj),_pM(ii));
               for (size_t j=1; j<=_dim; j++)
                  Swap(wm(i,j),wm(i+1,j));
               ii = jj;
            }
         }
      } while (istp>1);

//    Calculate M x (approximate eigenvectors) (err=0)
//    or final eigenvector approximations (err>0)
      for (size_t i=1; i<=_nb_eq; i++) {
         _ev.getColumn(i,ww);
         for (size_t j=1; j<=_dim; j++)
            _ev(j,i) = (ww,wm.getRow(j));
      }
      if (err>0) {
         _max_it = it;
         return 0;
      }

//    Check for convergence of eigenvalues
      for (size_t i=1; i<=_dim; i++)
         _rc(i) = fabs(_eigv(i)-old_eigv(i))/_eigv(i);
      if (_rc.getNormMax()<_epsv)
         err = 1;
      else {
         if (it<_max_it)
            old_eigv = _eigv;
         else
            err = 2;
      }
   }
   return ret;
}


void EigenProblemSolver::Mxv(const Vect<real_t>& b,
                             Vect<real_t>&       c)
{
   if (_diag) {
      for (size_t i=1; i<=b.size(); i++)
         c(i) = (*_lM)(i)*b(i);
      return;
   }
   else
      _M->Mult(b,c);
}


int EigenProblemSolver::runJacobi(DMatrix<real_t>& wm)
{
   const size_t max_it = 12;
   real_t ajj, aii, ab, det, ca, cg;
   real_t aj, bj, ak, bk, xj, xk, bb, eps;
   size_t i, j, k;
   Vect<real_t> wv(_dim);

// Initialization
   for (i=1; i<=_dim; i++) {
      if (_pK(i,i)<=0. || _pM(i,i)<=0.)
         return 1;
      _eigv(i) = wv(i) = _pK(i,i)/_pM(i,i);
   }
   wm = 1;
   if (_dim==1)
      return 0;

// Initialize iterations
   int nswp = 0;

// Beginning of loop
   do {
      nswp++;

//    Check if present off-diagonal element is large enough to require zeroing
      eps = pow(0.01,2*nswp);
      for (i=1; i<=_dim-1; i++) {
         for (j=i+1; j<=_dim; j++) {
            real_t epsa = _pK(i,j)*_pK(i,j)/(_pK(i,i)*_pK(j,j));
            real_t epsb = _pM(i,j)*_pM(i,j)/(_pM(i,i)*_pM(j,j));
            if (epsa>=eps || epsb>=eps) {

//             If zeroing is required, calculate the rotation matrix elements ca and cg
               ajj = _pK(j,j)*_pM(i,j) - _pM(j,j)*_pK(i,j);
               aii = _pK(i,i)*_pM(i,j) - _pM(i,i)*_pK(i,j);
               ab  = _pK(i,i)*_pM(j,j) - _pK(j,j)*_pM(i,i);
               det = 0.25*ab*ab + ajj*aii;
               if (det < 0)
                  return 1;
               real_t sqdt = sqrt(det);
               real_t det1=0.5*ab+sqdt, det2=0.5*ab-sqdt;
               real_t den = det1;
               if (fabs(det2) > fabs(det1))
                  den = det2;
               if (den == 0) {
                  ca = 0.;
                  cg = -_pK(i,j)/_pK(j,j);
               }
               else {
                  ca =  ajj/den;
                  cg = -aii/den;
               }

//             Find rotation of column i and column j above the diagonal
               if (_dim!=2) {
                  if (i>1) {
                     for (k=1; k<=i-1; k++) {
                        aj = _pK(k,i); bj = _pM(k,i);
                        ak = _pK(k,j); bk = _pM(k,j);
                        _pK(k,i) = aj + cg*ak;
                        _pM(k,i) = bj + cg*bk;
                        _pK(k,j) = ak + ca*aj;
                        _pM(k,j) = bk + ca*bj;
                     }
                  }

//                Find rotation of row i and row j starting from column j+1
                  if (j<_dim) {
                     for (k=j+1; k<=_dim; k++) {
                        aj = _pK(i,k); bj = _pM(i,k);
                        ak = _pK(j,k); bk = _pM(j,k);
                        _pK(i,k) = aj + cg*ak;
                        _pM(i,k) = bj + cg*bk;
                        _pK(j,k) = ak + ca*aj;
                        _pM(j,k) = bk + ca*bj;
                     }
                  }

//                Find rotation of row i between diagonal and
//                column j, and column j between row i and diagonal
                  if (i+1 <= j-1) {
                     for (k=i+1; k<=j-1; k++) {
                        aj = _pK(i,k); bj = _pM(i,k);
                        ak = _pK(k,j); bk = _pM(k,j);
                        _pK(i,k) = aj + cg*ak;
                        _pM(i,k) = bj + cg*bk;
                        _pK(k,j) = ak + ca*aj;
                        _pM(k,j) = bk + ca*bj;
                     }
                  }
	       }

//             Find rotation of elements (i,i) and (j,j) and zero element (i,j)
               ak = _pK(j,j); bk = _pM(j,j);
               _pK(j,j) = ak + ca*(2*_pK(i,j) + ca*_pK(i,i));
               _pM(j,j) = bk + ca*(2*_pM(i,j) + ca*_pM(i,i));
               _pK(i,i) = _pK(i,i) + cg*(2*_pK(i,j) + cg*ak);
               _pM(i,i) = _pM(i,i) + cg*(2*_pM(i,j) + cg*bk);
               _pK(i,j) = 0; _pM(i,j) = 0;

//             Update the eigenvector matrix after each rotation
               for (k=1; k<=_dim; k++) {
                  xj = wm(i,k);
                  xk = wm(j,k);
                  wm(i,k) = xj + cg*xk;
                  wm(j,k) = xk + ca*xj;
               }
            }
         }
      }

//    Update eigenvalues after each sweep
      for (i=1; i<=_dim; i++) {
         if (_pK(i,i)<=0. || _pM(i,i)<=0.)
            return 1;
         _eigv(i) = _pK(i,i)/_pM(i,i);
      }

//    Check for convergence
      for (i=0; i<_dim; i++)
         if (fabs(_eigv[i]-wv[i]) > _epsj*wv[i])
            goto L240;

//    Check all off diagonal elements to see if another sweep is required
      eps = _epsj*_epsj;
      for (i=1; i<=_dim-1; i++) {
         for (j=i+1; j<=_dim; j++) {
            real_t epsa = _pK(i,j)*_pK(i,j)/(_pK(i,i)*_pK(j,j)),
                   epsb = _pM(i,j)*_pM(i,j)/(_pM(i,i)*_pM(j,j));
            if (epsa>=eps || epsb>=eps)
               goto L240;
         }
      }

//    Fill out bottom triangle of resultant matrix and scale eigenvector
L210: for (i=1; i<=_dim; i++) {
         bb = sqrt(_pM(i,i));
         for (j=1; j<=_dim; j++)
            wm(i,j) /= bb;
      }
      return 0;

//    Update work array and start new sweep if allowed
L240: wv = _eigv;
   }
   while (nswp < int(max_it));

// Jacobi tolerance not satisfied
   goto L210;
}


int EigenProblemSolver::init()
{
   size_t ind=0, ncon=0;
   Vect<real_t> w(_nb_eq);
   int nums=int(_nb_eq/_dim);
   if (_diag) {
      for (size_t i=1; i<=_nb_eq; i++) {
         _ev(1,i) = (*_lM)(i);
         if ((*_lM)(i) > 0.)
            ncon++;
         w(i) = (*_lM)(i)/(*_K)(i,i);
      }
   }
   else {
      for (size_t i=1; i<=_nb_eq; i++) {
         _ev(1,i) = (*_M)(i,i);
         if ((*_M)(i,i) > 0.)
            ncon++;
         w(i) = (*_M)(i,i)/(*_K)(i,i);
      }
   }

// Check degrees of freedom
   if (_dim > ncon)
      return 1;
   for (size_t i=1; i<=_nb_eq; i++)
      for (size_t j=2; j<=_dim; j++)
         _ev(j,i) = 0.;

// Find lowest dividers and place them in ascending order
   int ldif=int(_nb_eq)-nums;
   for (size_t j=2; j<=_dim; j++) {
      real_t s = 0.;
      for (size_t i=1; i<=size_t(ldif); i++) {
         if (w(i) >= s) {
            s = w(i);
            ind = i;
         }
      }
      for (size_t i=size_t(ldif); i<=_nb_eq; i++) {
         if (w(i) > s) {
            s = w(i);
            ind = i;
         }
      }
      w(ind) = 0;
      ldif -= nums;
      _ev(j,ind) = 1;
   }
   return 0;
}


int EigenProblemSolver::checkSturm(int& nb_found,
                                   int& nb_lost)
{
   size_t ncon=0, err=0, i=0;
   real_t eps=0.01, eigu, eigl, eigenb;

   for (i=0; i<_dim; i++)
      if (_rc[i]<_epsv)
         ncon++;
   if (ncon==0)
      return 1;

// Find upper bound on cluster around reigv(ncon)
   if (ncon!=_dim) {
      eigu = _eigv(ncon)*(1.+eps);
      for (i=ncon+1; i<=_dim; i++) {
         eigl = _eigv(i)*(1.-eps);
         if (eigu<eigl)
            break;
      }
//    Register number of eigenvalues less than upper bound
      ncon = i;
   }
   eigenb = _eigv(ncon)*(1+eps);

// Shift K by eigenb and factorize
   if (_diag)
      for (i=1; i<=_nb_eq; i++)
         (*_K)(i,i) -= eigenb*(*_lM)(i);
   else
      for (i=1; i<=_nb_eq; i++)
         (*_K)(i,i) -= eigenb*(*_M)(i,i);
   _K->Factor();

// Count negative elements on diagonal
   int negd = 0;
   for (i=1; i<=_nb_eq; i++)
      if ((*_K)(i,i)<0)
         negd++;

// Give nb. of eigenvalues found and nb of missing.
   nb_found = negd;
   nb_lost = negd - ncon;
   return err;
}


real_t EigenProblemSolver::getEigenValue(int n) const
{
   return _eigv(n);
}


void EigenProblemSolver::getEigenVector(int           n,
                                        Vect<real_t>& v) const
{
   if (n>int(_nb_eigv))
      return;
   if (_theEqua) {
      v.resize(_theMesh->getNbDOF());
      Vect<real_t> w(_nb_eq);
      w = _ev.getRow(n);
      v.insertBC(*_theMesh,w);
   }
   else {
      v.resize(_nb_eq);
      v = _ev.getRow(n);
   }
}


ostream& operator<<(ostream&                  s,
                    const EigenProblemSolver& es)
{
   s << "\nEIGEN PROBLEM SOLVER\n\n";
   s << "Number of equations: \t\t" << es._nb_eq << endl;
   return s;
}

} /* namespace OFELI */
