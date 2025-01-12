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

                       Implementation of class LaplaceDG2DP1

  ==============================================================================*/


#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"
#include "solvers/LinearSolver.h"
#include "post/Reconstruction.h"
#include "equations/laplace/LaplaceDG2DP1.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/LocalMatrix_impl.h"

namespace OFELI {

LaplaceDG2DP1::LaplaceDG2DP1(Mesh&         ms,
                             Vect<real_t>& bf,
                             Vect<real_t>& bc,
                             Vect<real_t>& sf,
                             Vect<real_t>& u)
              : DG(ms,1)
{
   setMatrixType(SPARSE);
   setSolver(CG_SOLVER,DILU_PREC);
   _u = &u;
   _bf = &bf;
   _bc = &bc;
   _sf = &sf;
   _K = 1;
   _sigma = 100.;
   _eps = -1.;
}


LaplaceDG2DP1::~LaplaceDG2DP1()
{
}


void LaplaceDG2DP1::set(real_t sigma,
                        real_t eps)
{
   _sigma = sigma;
   _eps = eps;
}


void LaplaceDG2DP1::set(const LocalMatrix<real_t,2,2>& K)
{
   _K = K;
}


void LaplaceDG2DP1::set(Element& el)
{
   _theElement = &el;
   Triang3 t(_theElement);
   _el_geo.area = t.getArea();
   _dSh = t.DSh();
   _ne = _theElement->n();
}

  /*
void LaplaceDG2DP1::setSide(size_t k)
{
   _theSide = _theElement->getPtrSide(k);
   Point<real_t> N = _theElement->getUnitNormal(k);
   _el_geo.length = Line2(_theSide).getLength();
   _is(1) = (*_theSide)(1)->n(), _is(2) = (*_theSide)(2)->n();
   Triang3 T(_theElement);
   _el_geo.area = T.getArea();
   LocalVect<size_t,2> ls;
   _ls1(1) = _g2l[_ne-1][_is(1)-1];
   _ls1(2) = _g2l[_ne-1][_is(2)-1];
   _ls1(3) = 6 - _ls1(1) - _ls1(2);
   LocalVect<Point<real_t>,3> s;
   s(1).x = 0, s(1).y = 0;
   s(2).x = 1, s(2).y = 0;
   s(3).x = 0, s(3).y = 1;
   for (size_t i=1; i<=3; i++) {
      _z(i) = T.Sh(i,s(_ls1(1))) + T.Sh(i,s(_ls1(2)));
      _F1(i) = N.x*(_K(1,1)*T.DSh(i).x + _K(1,2)*T.DSh(i).y) +
               N.y*(_K(2,1)*T.DSh(i).x + _K(2,2)*T.DSh(i).y);
   }
   if (_theSide->isOnBoundary()==false) {
      Element *el = _theSide->getOtherNeighborElement(_theElement);
      _nf = el->n();
      _ls2(1) = _g2l[_nf-1][_is(1)-1];
      _ls2(2) = _g2l[_nf-1][_is(2)-1];
      _ls2(3) = 6 - _ls2(1) - _ls2(2);
      Triang3 t(el);
      for (size_t i=1; i<=3; i++)
         _F2(i) = N.x*(_K(1,1)*t.DSh(i).x + _K(1,2)*t.DSh(i).y) +
                  N.y*(_K(2,1)*t.DSh(i).x + _K(2,2)*t.DSh(i).y);
   }
}


void LaplaceDG2DP1::build()
{
   *_A = 0;
   *_b = 0;
   MESH_EL {
      set(TheElement);
      for (size_t i=1; i<=3; i++) {
         for (size_t j=1; j<=3; j++) {
            (*_A)(II(i),II(j)) += _el_geo.area*(_K(1,1)*_dSh(i).x*_dSh(j).x +
                                                _K(1,2)*_dSh(i).x*_dSh(j).y +
                                                _K(2,1)*_dSh(i).y*_dSh(j).x +
                                                _K(2,2)*_dSh(i).y*_dSh(j).y);
         }
         (*_b)(II(i)) += OFELI_THIRD*_el_geo.area*(*_bf)((*_theElement)(i)->n());
      }

//    Loop on element sides
      for (size_t k=1; k<=3; k++) {
         setSide(k);
         real_t c = 0.5*_el_geo.length;

//       Internal sides
         if (_theSide->isOnBoundary()==false) {
            for (size_t i=1; i<=3; i++) {
               size_t ii = _ls1(i);
               (*_A)(II(ii),II(ii)) += 0.5*_sigma*_z(ii);
               (*_A)(II(ii),IJ(_ls2(i))) -= 0.5*_sigma*_z(ii);
               for (size_t j=1; j<=3; j++) {
                  size_t jj = _ls1(j), kk = _ls2(j);
                  (*_A)(II(ii),II(jj)) -= 0.5*c*(_F1(jj)*_z(ii) - _eps*_F1(ii)*_z(jj));
                  (*_A)(II(ii),IJ(kk)) -= 0.5*c*(_F2(kk)*_z(ii) + _eps*_F1(ii)*_z(jj));
               }
            }
         }

//       Dirichlet boundary condition sides
         if ((*_theSide)(1)->getCode(1)>0 && (*_theSide)(2)->getCode(1)>0) {
            real_t g = c*((*_bc)(_is(1))+(*_bc)(_is(2)));
            for (size_t i=1; i<=2; i++) {
               size_t ii = _ls1(i);
               (*_b)(II(ii)) += _eps*_F1(ii)*g + 0.5*_z(ii)*_sigma*(*_bc)(_is(i));
               (*_A)(II(ii),II(ii)) += 0.5*_sigma*_z(ii);
               for (size_t j=1; j<=2; j++) {
                  size_t jj = _ls1(j);
                  (*_A)(II(ii),II(jj)) -= c*(_F1(jj)*_z(ii)-_eps*_F1(ii)*_z(jj));
               }
            }
         }
 
//       Neumann boundary condition sides
         if (_theSide->getCode(1)>0) {
            for (size_t i=1; i<=3; i++)
               (*_b)(II(_ls1(i))) += c*(*_sf)(_theSide->n());
         }
      }
   }
}


int LaplaceDG2DP1::run()
{
   build();
   return solve();
}


int LaplaceDG2DP1::solve()
{
   LinearSolver<double> ls(1000,1.e-6);
   int nb_it = 0;
   if (_eps==-1)
      nb_it = ls.solve(*_A,*_b,*_u,CG_SOLVER,DILU_PREC);
   else
      nb_it = ls.solve(*_A,*_b,*_u,GMRES_SOLVER,DILU_PREC);
   return nb_it;
}


void LaplaceDG2DP1::Smooth(Vect<real_t>& u)
{
   Reconstruction r(*_theMesh);
   r.DP1toP1(*_u,u);
}
  */
} /* namespace OFELI */
