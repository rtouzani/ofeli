/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani

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

                      Implementation of class 'EigenProblemSolver'

  ==============================================================================*/

#include "solvers/EigenProblemSolver.h"
#include "linear_algebra/DMatrix_impl.h"
#include "linear_algebra/DSMatrix_impl.h"
#include "linear_algebra/SkSMatrix_impl.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/Assembly_impl.h"
#include "equations/Equa.h"

using std::abs;

namespace OFELI {

EigenProblemSolver::EigenProblemSolver()
                   : _theEqua(nullptr), _theMesh(nullptr), _K(nullptr), _M(nullptr),
                     _K_alloc(false), _M_alloc(false), _lM_alloc(false), _diag(false),
                     _max_it(100), _dim(1), _nb_eq(0), _nb_eigv(1), _opt(1),
                     _epsv(1.e-8), _epsj(1.e-10), _sym(false), _eigv(false)
{ }


EigenProblemSolver::EigenProblemSolver(DMatrix<real_t>& A,
                                       bool             eigv)
                   : _theEqua(nullptr), _theMesh(nullptr), _K(&A), _M(nullptr),
                     _K_alloc(false), _M_alloc(false), _lM_alloc(false), _diag(true),
                     _max_it(100), _dim(1), _opt(1), _epsv(1.e-8), _epsj(1.e-10),
                     _sym(false), _eigv(eigv)
{
   _nb_eq = _nb_eigv = _K->getNbRows();
   _method = QR;
   _evr.setSize(_nb_eq);
   _evi.setSize(_nb_eq);
   run(_nb_eigv);
}


EigenProblemSolver::EigenProblemSolver(DSMatrix<real_t>& A,
                                       int               n)
                   : _theEqua(nullptr), _theMesh(nullptr), _K(&A), _M(nullptr),
                     _K_alloc(false), _M_alloc(false), _lM_alloc(true), _diag(true),
                     _max_it(100), _dim(1), _opt(1), _epsv(1.e-8), _epsj(1.e-10),
                     _sym(true), _eigv(false)
{
   _nb_eq = _nb_eigv = _K->getNbRows();
   if (n>0)
      _nb_eigv = n;
   _lM = new Vect<real_t>(_nb_eq);
   *_lM = 1;
   _method = SUBSPACE;
   run(_nb_eigv);
}


EigenProblemSolver::EigenProblemSolver(DSMatrix<real_t>& A,
                                       Vect<real_t>&     ev,
                                       int               n)
                   : _theEqua(nullptr), _theMesh(nullptr), _K(&A), _M(nullptr),
                     _K_alloc(false), _M_alloc(false), _lM_alloc(true), _diag(true),
                     _max_it(100), _dim(1), _opt(1), _epsv(1.e-8), _epsj(1.e-10),
                     _sym(true), _eigv(false)
{
   _nb_eq = _nb_eigv = _K->getNbRows();
   _nb_eigv = n;
   if (n==0)
      _nb_eigv = _nb_eq;
   _lM = new Vect<real_t>(_nb_eq);
   *_lM = 1;
   _method = SUBSPACE;
   run(n);
   ev.setSize(_nb_eigv);
   ev = _evr;
}


EigenProblemSolver::EigenProblemSolver(DMatrix<real_t>& A,
                                       Vect<real_t>&    evr,
                                       Vect<real_t>&    evi,
                                       bool             eigv)
                   : _theEqua(nullptr), _theMesh(nullptr), _K(&A), _M(nullptr),
                     _K_alloc(false), _M_alloc(false), _lM_alloc(false), _diag(true),
                     _max_it(50), _dim(1), _opt(1), _epsv(1.e-8), _epsj(1.e-10),
                     _sym(false), _eigv(eigv)
{
   _nb_eq = _nb_eigv = _K->getNbRows();
   _method = QR;
   run(_nb_eigv);
   evr.setSize(_nb_eigv);
   evi.setSize(_nb_eigv);
   evr = _evr;
   evi = _evi;
}


EigenProblemSolver::EigenProblemSolver(Equa& eq,
                                       bool  lumped)
                   : _K(nullptr), _M(nullptr),
                     _K_alloc(true), _M_alloc(false), _lM_alloc(true), _diag(true),
                     _max_it(100), _dim(1), _nb_eq(0), _nb_eigv(1), _opt(1),
                     _epsv(1.e-8), _epsj(1.e-10), _sym(false), _eigv(false)
{
   _method = SUBSPACE;
   setPDE(eq,lumped);
}


EigenProblemSolver::EigenProblemSolver(SkSMatrix<real_t>& K,
                                       SkSMatrix<real_t>& M,
                                       int                n)
                   : _K_alloc(false), _M_alloc(false), _lM_alloc(false), _diag(false),
                     _max_it(100), _dim(1), _opt(1), _epsv(1.e-8), _epsj(1.e-10),
                     _sym(true), _eigv(false)
{
   setMatrix(K,M);
   _nb_eigv = n;
   if (n==0)
      _nb_eigv = _nb_eq;
}


EigenProblemSolver::EigenProblemSolver(SkSMatrix<real_t>& K,
                                       Vect<real_t>&      M,
                                       int                n)
                   : _K_alloc(false), _M_alloc(false), _lM_alloc(false), _diag(true),
                     _max_it(100), _dim(1), _opt(1), _epsv(1.e-8), _epsj(1.e-10),
                     _sym(true), _eigv(false)
{
   setMatrix(K,M);
   _nb_eigv = n;
   if (n==0)
      _nb_eigv = _nb_eq;
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


void EigenProblemSolver::setNbEigv(int n)
{
    _nb_eigv = n;
}


void EigenProblemSolver::setEigenVectors()
{
   _eigv = true;
}


void EigenProblemSolver::set(EigenMethod m,
                             bool        sym)
{
   _method = m;
   _sym = sym;
   if (m==SUBSPACE && sym==false)
      throw OFELIException("In EigenProblemSolver::set(..): This method is not valid for unsymmetric matrices");
}


void EigenProblemSolver::setMatrix(Matrix<real_t>* A)
{
   _K = A;
   _nb_eq = _K->getNbRows();
   _diag = false;
}


void EigenProblemSolver::setMatrix(Matrix<real_t>* K,
                                   Matrix<real_t>* M)
{
   _K = K;
   _M = M;
   _nb_eq = _K->getNbRows();
   _diag = false;
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
   _method = SUBSPACE;
}


void EigenProblemSolver::setMatrix(DSMatrix<real_t>& K)
{
   _K = &K;
   _lM = new Vect<real_t>(_nb_eq);
   _lM_alloc = true;
   *_lM = 1;
   _diag = true;
   _method = SUBSPACE;
}


void EigenProblemSolver::setMatrix(DMatrix<real_t>& K)
{
   _K = &K;
   _lM = new Vect<real_t>(_nb_eq);
   _lM_alloc = true;
   *_lM = 1;
   _diag = true;
   _method = QR;
}


void EigenProblemSolver::setPDE(Equa& eq,
                                bool  lumped)
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
      _M = new SkSMatrix<real_t>(*_theMesh);
   }
}


int EigenProblemSolver::run(int  nb,
                            bool opt)
{
   if (_theEqua)
      _theEqua->build(*this);
   _nb_eigv = nb;
   if (nb==0)
      _nb_eigv = _nb_eq;
   if (_method==SUBSPACE)
      return runSubSpace(_nb_eigv);
   else if (_method==QR) {
      _evr.setSize(_nb_eq);
      _evi.setSize(_nb_eq);
      _evec.setSize(_nb_eq,_nb_eq);
      eigen(0,1);
      norm_2();
   }
   return 0;
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
   _nb_eigv = std::min(nb_eigv,_nb_eq);
   _dim = std::max(ss_dim,_nb_eigv+1);
   _dim = std::min(_dim,_nb_eq);
   Vect<real_t> old_eigv(_dim), wv(_nb_eq), ww(_dim);
   DMatrix<real_t> wm(_dim);
   _pK.setSize(_dim);
   _pM.setSize(_dim);
   _ev.setSize(_dim,_nb_eq);
   _evr.setSize(_dim);
   _evi.setSize(_dim);
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
         _K->solve(_ev.getRow(i),wv,false);
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
         for (int i=1; i<=nnend; i++) {
            jj = ii + int(_dim) - i - 1;
            if (_evr(i+1)<_evr(i)) {
               istp = i;
               Swap(_evr(i),_evr(i+1));
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
         _rc(i) = fabs(_evr(i)-old_eigv(i))/_evr(i);
      if (_rc.getNormMax()<_epsv)
         err = 1;
      else {
         if (it<_max_it)
            old_eigv = _evr;
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
      _evr(i) = wv(i) = _pK(i,i)/_pM(i,i);
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
               real_t det1 = 0.5*ab+sqdt, det2 = 0.5*ab-sqdt;
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
         _evr(i) = _pK(i,i)/_pM(i,i);
      }

//    Check for convergence
      for (i=0; i<_dim; i++)
         if (fabs(_evr[i]-wv[i]) > _epsj*wv[i])
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
L240: wv = _evr;
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
      eigu = _evr(ncon)*(1.+eps);
      for (i=ncon+1; i<=_dim; i++) {
         eigl = _evr(i)*(1.-eps);
         if (eigu<eigl)
            break;
      }
//    Register number of eigenvalues less than upper bound
      ncon = i;
   }
   eigenb = _evr(ncon)*(1+eps);

// Shift K by eigenb and factorize
   if (_diag)
      for (i=1; i<=_nb_eq; i++)
         _K->add(i,i,-eigenb*(*_lM)(i));
   else
      for (i=1; i<=_nb_eq; i++)
         _K->add(i,i,-eigenb*(*_M)(i,i));
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


real_t EigenProblemSolver::getEigenValue(int n,
                                         int i) const
{
   if (n<=0 || n>_nb_eigv)
      throw OFELIException("In EigenProblemSolver::getEigenValue(n,i): Illegal value of n");
   if (i!=1 && i!=2)
      throw OFELIException("In EigenProblemSolver::getEigenValue(n,i): Illegal value of i");
   if (i==1)
      return _evr(n);
   if (i==2)
      return _evi(n);
   return 0.0;
}


void EigenProblemSolver::getEigenVector(int           n,
                                        Vect<real_t>& v) const
{
   if (_method!=SUBSPACE)
      throw OFELIException("In EigenProblemSolver::getEigenVector(n,v): This function "
                           "is available for the subspace method only");
   if (n>int(_nb_eigv))
      return;
   if (_theEqua) {
      v.setSize(_theMesh->getNbDOF());
      Vect<real_t> w(_nb_eq);
      w = _ev.getRow(n);
      v.insertBC(*_theMesh,w);
   }
   else {
      v.setSize(_nb_eq);
      v = _ev.getRow(n);
   }
}


void EigenProblemSolver::getEigenVector(int           n,
                                        Vect<real_t>& vr,
                                        Vect<real_t>& vi) const
{
   if (_method!=QR)
      throw OFELIException("In EigenProblemSolver::getEigenVector(n,vr,vi): This function "
                           "is available for the QR method only");
   if (n>int(_nb_eigv))
      return;
   if (_theEqua) {
      vr.setSize(_theMesh->getNbDOF());
      Vect<real_t> w(_nb_eq);
      w = _ev.getRow(n);
      vr.insertBC(*_theMesh,w);
   }
   else {
      vr.setSize(_nb_eq);
      vi.setSize(_nb_eq);
      for (int i=1; i<=_nb_eq; ++i)
         vr(i) = _evec(i,n);
      if (_evi(n)>0.0) {
         for (int i=1; i<=_nb_eq; ++i)
            vi(i) = _evec(i,n+1);
      }
      else if (_evi(n)<0.0) {
         for (int i=1; i<=_nb_eq; ++i)
            vi(i) = -_evec(i,n+1);
      }
   }
}


void EigenProblemSolver::balance(Vect<real_t>& scale)
{
   const real_t b2 = 2.0;
   bool iter=false;
   int k=_nb_eq-1, m=0;
   do {
      iter = false;
      for (int j=k; j>=0; j--) {
         real_t r=0.0;
         for (int i=0; i<=k; i++) {
            if (i!=j)
               r += abs((*_K)(j+1,i+1));
         }
         if (r==0.0) {
            scale[k] = j;
            if (j!=k) {
               for (int i=0; i<=k; i++)
                  Swap((*_K)(i+1,j+1),(*_K)(i+1,k+1));
               for (int i=m; i<_nb_eq; i++)
                  Swap((*_K)(j+1,i+1),(*_K)(k+1,i+1));
            }
            k--;
            iter = true;
         }
      }
   }
   while (iter);

   do {
      iter = false;
      for (int j=m; j<=k; j++) {
         real_t c=0.0;
         for (int i=m; i<=k; i++) {
            if (i!=j)
               c += abs((*_K)(i+1,j+1));
         }
         if (c==0.0) {
            scale[m] = j;
            if (j!=m) {
               for (int i=1; i<=k+1; i++)
                  Swap((*_K)(i,j+1),(*_K)(i,m+1));
               for (int i=m+1; i<=_nb_eq; i++)
                  Swap((*_K)(j+1,i),(*_K)(m+1,i));
            }
            m++;
            iter = true;
         }
      }
   }
   while (iter);

   _low = m, _high = k;
   for (int i=m; i<=k; i++)
      scale[i] = 1.0;

   do {
      iter = false;
      for (int i=m; i<=k; i++) {
         real_t c=0.0, r=0.0;
         for (int j=m; j<=k; j++) {
            if (j!=i) {
               c += abs((*_K)(j+1,i+1));
               r += abs((*_K)(i+1,j+1));
            }
         }
         real_t g = r/b2, f = 1.0, s = c + r;

         while (c < g)
            f *= b2, c *= b2;
         g = r*b2;
         while (c >= g)
            f /= b2, c /= b2;
         if ((c+r)/f < 0.95*s) {
            g = 1.0/f;
            scale[i] *= f;
            iter = true;
            for (int j=m; j<_nb_eq; j++)
               (*_K)(i+1,j+1) *= g;
            for (int j=0; j<=k; j++)
               (*_K)(j+1,i+1) *= f;
         }
      }
   }
   while (iter);
}


void EigenProblemSolver::balance_back(const Vect<real_t>& scale)
{
   for (int i=_low; i<=_high; i++) {
      real_t s = scale[i];
      for (int j=1; j<=_nb_eq; j++)
         _evec(i+1,j) *= s;
   }
   for (int i=_low-1; i>=0; i--) {
      int k = scale[i];
      if (k!=i)
         for (int j=1; j<=_nb_eq; j++)
            Swap(_evec(i+1,j),_evec(k+1,j));
   }
   for (int i=_high+1; i<_nb_eq; i++) {
      int k = scale[i];
      if (k!=i)
         for (int j=1; j<=_nb_eq; j++)
            Swap(_evec(i+1,j),_evec(k+1,j));
   }
}


void EigenProblemSolver::elmhes(Vect<int>& perm)
{
   for (int m=_low+1; m<_high; m++) {
      int i = m;
      real_t x = 0.0;
      for (int j=m; j<=_high; j++) {
         if (abs((*_K)(j+1,m)) > abs(x)) {
            x = (*_K)(j+1,m);
            i = j;
         }
      }

      perm[m] = i;
      if (i != m) {
         for (int j=m; j<=_nb_eq; j++)
            Swap((*_K)(i+1,j),(*_K)(m+1,j));
         for (int j=0; j<=_high; j++)
            Swap((*_K)(j+1,i+1),(*_K)(j+1,m+1));
      }

      if (x != 0.0) {
         for (int i=m+1; i<=_high; i++) {
            real_t y = (*_K)(i+1,m);
            if (y != 0.0) {
               y = (*_K)(i+1,m) = y / x;
               for (int j=m; j<_nb_eq; j++)
                  (*_K)(i+1,j+1) -= y * (*_K)(m+1,j+1);
               for (int j=0; j<=_high; j++)
                  (*_K)(j+1,m+1) += y * (*_K)(j+1,i+1);
            }
         }
      }
   }
}


void EigenProblemSolver::elmtrans(const Vect<int>& perm)
{
   for (int i=1; i<=_nb_eq; i++) {
      for (int k=1; k<=_nb_eq; k++)
         _evec(i,k) = 0.0;
      _evec(i,i) = 1.0;
   }
   for (int i=_high-1; i>_low; i--) {
      int j = perm[i];
      for (int k=i+1; k<=_high; k++)
         _evec(k+1,i+1) = (*_K)(k+1,i);
      if (i!=j) {
         for (int k=i; k<=_high; k++) {
            _evec(i+1,k+1) = _evec(j+1,k+1);
            _evec(j+1,k+1) = 0.0;
         }
         _evec(j+1,i+1) = 1.0;
      }
   }
}


void EigenProblemSolver::orthes(Vect<real_t>& d)
{
   const real_t eps = 128.0*DBL_EPSILON;
   real_t s=0.0, x=0.0;
   for (int m=_low+1; m<_high; m++) {
      real_t y=0.0;
      for (int i=_high; i>=m; i--) {
         x = (*_K)(i+1,m);
         d[i] = x;
         y += x*x;
      }
      if (y <= eps)
         s = 0.0;
      else {
         s = (x>=0.0) ? -sqrt(y) : sqrt(y);
         y -= x*s;
         d[m] =  x - s;
         for (int j=m; j<_nb_eq; j++) {
            real_t x = 0.0;
            for (int i=_high; i>=m; i--)
               x += d[i] * (*_K)(i+1,j+1);
            x /= y;
            for (int i=m; i<=_high; i++)
               (*_K)(i+1,j+1) -= x*d[i];
         }

         for (int i=0; i<=_high; i++) {
            real_t x = 0.0;
            for (int j=_high; j>=m; j--)
               x += d[j]*(*_K)(i+1,j+1);
            x /= y;
            for (int j=m; j<=_high; j++)
               (*_K)(i+1,j+1) -= x*d[j];
         }
      }
      (*_K)(m+1,m) = s;
   }
}


void EigenProblemSolver::orttrans(Vect<real_t>& d)
{
   for (int i=1; i<=_nb_eq; ++i) {
      for (int j=1; j<=_nb_eq; ++j)
         _evec(i,j) = 0.0;
      _evec(i,i) = 1.0;
   }
   for (int m=_high-1; m>_low; m--) {
      real_t y = (*_K)(m+1,m);
      if (y != 0.0) {
         y *= d[m];
         for (int i=m+1; i<=_high; i++)
            d[i] = (*_K)(i+1,m);
         for (int j=m; j<=_high; j++) {
            real_t x = 0.0;
            for (int i=m; i<=_high; i++)
               x += d[i]*_evec(i+1,j+1);
            x /= y;
            for (int i=m; i<=_high; i++)
               _evec(i+1,j+1) += x*d[i];
         }
      }
   }
}


int EigenProblemSolver::hqrvec()
{
   real_t norm=0.0;
   for (int i=0; i<_nb_eq; i++) {
      for (int j=i; j<_nb_eq; j++)
         norm += abs((*_K)(i+1,j+1));
   }
   if (norm==0.0)
      return 1;

   real_t z=0.0, r=0.0, w=0.0, s=0.0, t=0.0;
   for (int en=_nb_eq-1; en>=0; en--) {
      real_t p = _evr[en], q = _evi[en];
      int na = en - 1;
      if (q==0.0) {
         int m = en;
         (*_K)(en+1,en+1) = 1.0;
         for (int i=na; i>=0; i--) {
            w = (*_K)(i+1,i+1)-p, r = (*_K)(i+1,en+1);
            for (int j=m; j<=na; j++)
               r += (*_K)(i+1,j+1) * (*_K)(j+1,en+1);
            if (_evi[i] < 0.0)
               z = w, s = r;
            else {
               m = i;
               if (_evi[i]==0.0)
                  (*_K)(i+1,en+1) = -r/((w!=0.0) ? (w) : (DBL_EPSILON*norm));
               else {
                  real_t x = (*_K)(i+1,i+2), y = (*_K)(i+2,i+1);
                  q = SQR(_evi[i]-p) + SQR(_evi[i]);
                  (*_K)(i+1,en+1) = t = (x*s - z*r)/q;
                  (*_K)(i+2,en+1) = ((abs(x)>abs(z)) ? (-r-w*t)/x : (-s-y*t)/z);
               }
            }
         }
      }

      else if (q < 0.0) {
         int m = na;
         if (abs((*_K)(en+1,na+1)) > abs((*_K)(na+1,en+1))) {
            (*_K)(na+1,na+1) = - ((*_K)(en+1,en+1) - p) / (*_K)(en+1,na+1);
            (*_K)(na+1,en+1) = - q/(*_K)(en+1,na+1);
         } else {
            complex_t c = complex_t(-(*_K)(na+1,en+1),0.0)/complex_t((*_K)(na+1,na+1)-p,q);
            (*_K)(na+1,na+1) = c.real(), (*_K)(na+1,en+1) = c.imag();
         }
         (*_K)(en+1,na+1) = 1.0, (*_K)(en+1,en+1) = 0.0;
         for (int i=na-1; i>=0; i--) {
            real_t w = (*_K)(i+1,i+1) - p;
            real_t ra = (*_K)(i+1,en+1), sa = 0.0;
            for (int j=m; j<=na; j++) {
               ra += (*_K)(i+1,j+1) * (*_K)(j+1,na+1);
               sa += (*_K)(i+1,j+1) * (*_K)(j+1,en+1);
            }
            if (_evi[i]<0.0)
               z = w, r = ra, s = sa;
            else {
               m = i;
               if (_evi[i]==0.0) {
                  complex_t c = complex_t(-ra,-sa)/complex_t(w,q);
                  (*_K)(i+1,na+1) = c.real(); (*_K)(i+1,en+1) = c.imag();
               }
               else {
                  real_t x = (*_K)(i+1,i+2), y = (*_K)(i+2,i+1);
                  real_t vr = SQR(_evr[i]-p) + SQR(_evi[i]) - SQR(q);
                  real_t vi = 2.0*q*(_evr[i] - p);
                  if (vr==0.0 && vi==0.0)
                     vr = DBL_EPSILON*norm*(abs(w) + abs(q) + abs(x) + abs(y) + abs(z));
                  complex_t c=complex_t(x*r-z*ra+q*sa,x*s-z*sa-q*ra)/complex_t(vr,vi);
                  (*_K)(i+1,na+1) = c.real(), (*_K)(i+1,en+1) = c.imag();
                  if (abs(x) > abs(z)+abs(q)) {
                     (*_K)(i+2,na+1) = (-ra - w*(*_K)(i+1,na+1) + q*(*_K)(i+1,en+1))/x;
                     (*_K)(i+2,en+1) = (-sa - w*(*_K)(i+1,en+1) - q*(*_K)(i+1,na+1))/x;
                  }
                  else {
                     complex_t c=complex_t(-r-y*(*_K)(i+1,na+1),-s-y*(*_K)(i+1,en+1))/complex_t(z,q);
                     (*_K)(i+2,na+1) = c.real(), (*_K)(i+2,en+1) = c.imag();
                  }
               }
            }
         }
      }
   }

   for (int i=0; i<_nb_eq; i++)
      if (i<_low || i>_high)
         for (int k=i+1; k<_nb_eq; k++)
            _evec(i+1,k+1) = (*_K)(i+1,k+1);

   for (int j=_nb_eq-1; j>=_low; j--) {
      int m = (j<=_high) ? j : _high;
      if (_evi[j]<0.0) {
         for (int l=j-1, i=_low; i<=_high; i++) {
            real_t y=0.0, z=0.0;
            for (int k=_low; k<=m; k++) {
               y += _evec(i+1,k+1) * (*_K)(k+1,l+1);
               z += _evec(i+1,k+1) * (*_K)(k+1,j+1);
            }
            _evec(i+1,l+1) = y, _evec(i+1,j+1) = z;
         }
      }
      else {
         if (_evi[j]==0.0) {
            for (int i=_low; i<=_high; i++) {
               real_t z=0.0;
               for (int k=_low; k<=m; k++)
                  z += _evec(i+1,k+1) * (*_K)(k+1,j+1);
               _evec(i+1,j+1) = z;
            }
         }
      }
   }
   return 0;
}


int EigenProblemSolver::hqr2(Vect<int>& cnt)
{
   int  na, iter, l;
   real_t p=0.0, q=0.0, r=0.0, s, w, x, y, z;

   for (int i=0; i<_nb_eq; i++) {
      if (i<_low || i>_high) {
         _evr[i] = (*_K)(i+1,i+1);
         _evi[i] = 0.0;
         cnt[i] = 0;
      }
   }

   int en=_high;
   real_t t=0.0;
   while (en>=_low) {
      iter = 0;
      na = en - 1;

      for ( ; ; ) {
         for (l=en; l>_low; l--)
            if (abs((*_K)(l+1,l))<=DBL_EPSILON*(abs((*_K)(l,l)) + abs((*_K)(l+1,l+1))))
               break;
         x = (*_K)(en+1,en+1);
         if (l==en) {
            _evr[en] = (*_K)(en+1,en+1) = x + t;
            _evi[en] = 0.0;
            cnt[en--] = iter;
            break;
         }
         y = (*_K)(na+1,na+1);
         w = (*_K)(en+1,na+1)*(*_K)(na+1,en+1);

         if (l==na) {
            p = (y-x)*0.5;
            q = p*p + w;
            z = sqrt(abs(q));
            x = (*_K)(en+1,en+1) = x + t;
            (*_K)(na+1,na+1) = y + t;
            cnt[en] = -iter;
            cnt[na] = iter;
            if (q >= 0.0) {
               z = (p<0.0) ? (p-z) : (p+z);
               _evr[na] = x + z;
               _evr[en] = s = x - w/z;
               _evi[na] = _evi[en] = 0.0;
               x = (*_K)(en+1,na+1);
               r = sqrt(x*x + z*z);

               if (_eigv) {
                  p = x/r, q = z/r;
                  for (int j=na; j<_nb_eq; j++) {
                     z = (*_K)(na+1,j+1);
                     (*_K)(na+1,j+1) = q*z + p*(*_K)(en+1,j+1);
                     (*_K)(en+1,j+1) = q*(*_K)(en+1,j+1) - p*z;
                  }

                  for (int i=0; i<=en; i++) {
                     z = (*_K)(i+1,na+1);
                     (*_K)(i+1,na+1) = q*z + p*(*_K)(i+1,en+1);
                     (*_K)(i+1,en+1) = q*(*_K)(i+1,en+1) - p*z;
                  }

                  for (int i=_low; i<=_high; i++) {
                     z = _evec(i+1,na+1);
                     _evec(i+1,na+1) = q*z + p*_evec(i+1,en+1);
                     _evec(i+1,en+1) = q*_evec(i+1,en+1) - p*z;
                  }
               }
            }
            else {
               _evr[na] = _evr[en] = x + p;
               _evi[na] =  z;
               _evi[en] = -z;
            }
            en -= 2;
            break;
         }

         if (iter >= _max_it) {
            cnt[en] = _max_it + 1;
            return en;
         }

         if (iter!=0 && iter%10==0) {
            t += x;
            for (int i=_low; i<=en; i++)
               (*_K)(i+1,i+1) -= x;
            s = abs((*_K)(en+1,na+1)) + abs((*_K)(na+1,en-1));
            x = y = 0.75*s;
            w = -0.4375*s*s;
         }
         iter++;

         int m=0;
         for (m=en-2; m>=l; m--) {
            z = (*_K)(m+1,m+1);
            r = x - z;
            s = y - z;
            p = (r*s - w)/(*_K)(m+2,m+1) + (*_K)(m+1,m+2);
            q = (*_K)(m+2,m+2) - z - r - s;
            r = (*_K)(m+3,m+2);
            s = abs(p) + abs(q) + abs(r);
            p /= s, q /= s, r /= s;
            if (m==l)
               break;
            if (abs((*_K)(m+1,m))*(abs(q)+abs(r)) <=
               DBL_EPSILON*abs(p)*(abs((*_K)(m,m))+abs(z)+abs((*_K)(m+2,m+2))) )
               break;
         }

         for (int i=m+2; i<=en; i++)
            (*_K)(i+1,i-1) = 0.0;
         for (int i=m+3; i<=en; i++)
            (*_K)(i+1,i-2) = 0.0;

         for (int k=m; k<=na; k++) {
            if (k != m) {
               p = (*_K)(k+1,k), q = (*_K)(k+2,k);
               r = (k != na) ? (*_K)(k+3,k) : 0.0;
               x = abs(p) + abs(q) + abs(r);
               if (x==0.0)
                  continue;
               p /= x, q /= x, r /= x;
            }
            s = sqrt(p*p + q*q + r*r);
            (p < 0.0) ? s = -s : 0;

            if (k!=m)
               (*_K)(k+1,k) = -s*x;
            else if (l!=m)
               (*_K)(k+1,k) = -(*_K)(k+1,k);
            p += s;
            x = p/s, y = q/s, z = r/s;
            q /= p, r /= p;

            for (int j=k; j<_nb_eq; j++) {
               p = (*_K)(k+1,j+1) + q*(*_K)(k+2,j+1);
               if (k!=na) {
                  p += r*(*_K)(k+3,j+1);
                  (*_K)(k+3,j+1) -= p*z;
               }
               (*_K)(k+2,j+1) -= p*y, (*_K)(k+1,j+1) -= p*x;
            }

            int j = (k+3<en) ? (k+3) : en;
            for (int i=0; i<=j; i++) {
               p = x*(*_K)(i+1,k+1) + y*(*_K)(i+1,k+2);
               if (k != na) {
                  p += z*(*_K)(i+1,k+3);
                  (*_K)(i+1,k+3) -= p*r;
               }
               (*_K)(i+1,k+2) -= p*q;
               (*_K)(i+1,k+1) -= p;
            }

            if (_eigv) {
               for (int i=_low; i<=_high; i++) {
                  p = x*_evec(i+1,k+1) + y*_evec(i+1,k+2);
                  if (k != na) {
                     p += z*_evec(i+1,k+3);
                     _evec(i+1,k+3) -= p*r;
                  }
                  _evec(i+1,k+2) -= p*q;
                  _evec(i+1,k+1) -= p;
               }
            }
         }
      }
   }
   if (_eigv && hqrvec())
      return 99;
   return 0;
}


void EigenProblemSolver::norm_2()
{
   for (int j=1; j<=_nb_eq; ++j) {
      real_t xn = 0.;
      if (_evi(j)==0.0) {
         for (int i=1; i<=_nb_eq; ++i)
            xn += SQR(_evec(i,j));
         xn = sqrt(xn);
         for (int i=1; i<=_nb_eq; ++i)
            _evec(i,j) /= xn;
      }
      else {
         for (int i=1; i<=_nb_eq; ++i)
            xn += SQR(_evec(i,j)) + SQR(_evec(i,j+1));
         xn = sqrt(xn);
         for (int i=1; i<=_nb_eq; ++i) {
            _evec(i,j  ) /= xn;
            _evec(i,j+1) /= xn;
         }
         j++;
      }
   }
}


int EigenProblemSolver::eigen(int ortho,
                              int ev_norm)
{
   if (_nb_eq<1)
      return 1;

   if (_nb_eq==1) {
      if (_eigv)
         _evec(1,1) = 1.0;
      _evr(1) = (*_K)(1,1);
      _evi(1) = 0.0;
      return 0;
   }

   _low = _high = 0;
   Vect<real_t> scale(_nb_eq), d(_nb_eq);
   balance(scale);

   Vect<int> cnt(_nb_eq);
   cnt = 0;
   int ret = 0;
   if (ortho)
      orthes(d);
   else
      elmhes(cnt);
   if (ret)
      return ret;
   if (_eigv) {
      if (ortho)
         orttrans(d);
      else
         elmtrans(cnt);
   }
   ret = hqr2(cnt);
   if (ret)
      return ret;
   if (_eigv) {
      balance_back(scale);
      if (ret)
         return ret;
  }
  return ret;
}


ostream& operator<<(ostream&                  s,
                    const EigenProblemSolver& es)
{
   s << "\nEIGEN PROBLEM SOLVER\n\n";
   s << "Number of equations:         " << setw(6) << es._nb_eq << endl;
   s << "Number of found eigenvalues: " << setw(6) << es._nb_eigv << endl;
   s << "Eigenvalue extraction method: ";
   if (es._method==SUBSPACE)
      s << "Subspace iteration (Symmetric matrix)" << endl;
   else if (es._method==QR)
      s << "QR method (General matrix)" << endl;
   if (Verbosity>2) {
      Vect<real_t> vr(es._nb_eq), vi(es._nb_eq);
      s.setf(ios::scientific);
      for (int i=1; i<=es._nb_eq; i++) {
         s << "\nEigenvalue #" << i << ": " 
           << out_complex(es._evr(i),es._evi(i)) << endl;
         if (Verbosity>5 && es._eigv) {
            s << "Eigen vector:" << endl;
            es.getEigenVector(i,vr,vi);
            for (int j=1; j<=es._nb_eq; j++)
               s << setw(8) << j << " " << setprecision(8) 
                 << setw(8) << out_complex(vr(j),vi(j)) << endl;
         }
      }
   }
   return s;
}

} /* namespace OFELI */
