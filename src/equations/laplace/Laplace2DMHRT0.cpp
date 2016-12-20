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

                      Class Laplace2DMHRT0: 2-D Laplace equation
        using 3-Node Mixed Hybrid Lowest degree Raviart-Thomas Finite element

  ==============================================================================*/

#include "equations/laplace/Laplace2DMHRT0.h"

namespace OFELI {

Laplace2DMHRT0::Laplace2DMHRT0(Mesh&             ms,
                               SpMatrix<real_t>& A,
                               Vect<real_t>&     b)
{
   _theMesh = &ms;
   _A = &A;
   _u = &b;
   _K.set(IDENTITY);
   _IK = _K;
   _bf_given = _bc_given = _sf_given = false;
}


void Laplace2DMHRT0::Set(const Element* el)
{
   _nb_dof = 1;
   Init(el);
   Triang3 tr(_theElement);
   _area = tr.getArea();
   _sd1 = _theElement->getPtrSide(1);
   _sd2 = _theElement->getPtrSide(2);
   _sd3 = _theElement->getPtrSide(3);
   getNormals(tr);
   eMat = 0;
   eRHS = 0;
}


void Laplace2DMHRT0::setDiffusivity(const LocalMatrix<real_t,2,2>& K)
{
   _K = K;
   _IK = K;
   _K.Invert(_IK);
   _K = K;
}


void Laplace2DMHRT0::build()
{
   MESH_EL {
      Set(theElement);
      LM_LHS();
      LM_RHS();
      side_assembly(*theElement,eMat,_A);
      side_assembly(*theElement,eRHS,*_u);
   }
}


int Laplace2DMHRT0::solve(Vect<real_t>& u)
{
   Vect<real_t> x(_A->size());
   LinearSolver<real_t> ls(*_A,*_u,x);
   ls.setTolerance(1.e-7);
   int nb_it = ls.solve(CG_SOLVER,ILU_PREC);
   u.insertBC(*_theMesh,x,*_bc);
   return nb_it;
}


void Laplace2DMHRT0::LM_LHS()
{
   for (size_t i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eMat(i,j) = getIP(_K,_n(i),_n(j))/_area;
}


void Laplace2DMHRT0::LM_RHS()
{
// Sources
   real_t ff=0.;
   if (_bf_given)
      ff = 0.5*(*_bf)(_theElement->n());

// Contribution of Dirichlet Boundary Condition
   Point<real_t> g;
   if (_bc_given) {
      if (_sd1->getCode(1)>0)
         g += (*_bc)(_sd1->n())*getP(_K,_n(1))/_area;
      if (_sd2->getCode(1)>0)
         g += (*_bc)(_sd2->n())*getP(_K,_n(2))/_area;
      if (_sd3->getCode(1)>0)
         g += (*_bc)(_sd3->n())*getP(_K,_n(3))/_area;
   }

// Contribution of Neumann Boundary Condition
   LocalVect<real_t,3> h;
   if (_sf_given) {
      if (_sd1->getCode(1)<0)
         h(1) = (*_sf)(_sd1->n())*Line2(_sd1).getLength();
      if (_sd2->getCode(1)<0)
         h(2) = (*_sf)(_sd2->n())*Line2(_sd2).getLength();
      if (_sd3->getCode(1)<0)
         h(3) = (*_sf)(_sd3->n())*Line2(_sd3).getLength();
   }

// Right-Hand Side
   Point<real_t> d = getP(_IK,_c);
   eRHS(1) = ff*(_n(1)*_ce(1) - getIP(_K,d,_n(1))) - g*_n(1) + h(1);
   eRHS(2) = ff*(_n(2)*_ce(2) - getIP(_K,d,_n(2))) - g*_n(2) + h(2);
   eRHS(3) = ff*(_n(3)*_ce(3) - getIP(_K,d,_n(3))) - g*_n(3) + h(3);
}


void Laplace2DMHRT0::Post(const Vect<real_t>&         lambda,
                          const Vect<real_t>&         f,
                                Vect<real_t>&         v,
                                Vect<Point<real_t> >& p,
                                Vect<real_t>&         u)
{
   _vv.resize(_theMesh->getNbNodes(),0);
   MESH_EL {
      size_t n = TheElement.n();
      Set(theElement);
      real_t l1=lambda(_sd1->n()), l2=lambda(_sd2->n()), l3=lambda(_sd3->n());
      SolAtNodes(theElement,lambda,v);
      Point<real_t> s = l1*_n(1) + l2*_n(2) + l3*_n(3);
      Point<real_t> d = getP(_IK,_c);
      p(n) = getP(_K,0.5*f(n)*d + s/_area);
      real_t t = _n(1)*_ce(1)*l1 + _n(2)*_ce(2)*l2 + _n(3)*_ce(3)*l3;
      real_t k = OFELI_THIRD*_area*(getIP(_IK,_ce(1),_ce(1))+getIP(_IK,_ce(2),_ce(2))+getIP(_IK,_ce(3),_ce(3)));
      u.set(n,0.5*(0.5*f(n)*k + t)/_area - 0.5*d*p(n));
      p.add(n,-0.5*f(n)*_c);
   }
}


void Laplace2DMHRT0::SolAtNodes(      Element*      el,
                                const Vect<real_t>& lambda,
                                      Vect<real_t>& v)
{
   size_t n1=(*_sd1)(1)->n(), n2=(*_sd1)(2)->n();
   size_t n3 = NotOnSide(el,n1,n2);
   if (_vv(n1)==0)
      v(n1) += lambda(_sd1->n());
   if (_vv(n2)==0)
      v(n2) += lambda(_sd1->n());
   if (_vv(n3)==0)
      v(n3) -= lambda(_sd1->n());
   n1 = _sd2->getNodeLabel(1);
   n2 = _sd2->getNodeLabel(2);
   n3 = NotOnSide(el,n1,n2);
   if (_vv(n1)==0)
      v(n1) += lambda(_sd2->n());
   if (_vv(n2)==0)
      v(n2) += lambda(_sd2->n());
   if (_vv(n3)==0)
      v(n3) -= lambda(_sd2->n());
   n1 = _sd3->getNodeLabel(1);
   n2 = _sd3->getNodeLabel(2);
   n3 = NotOnSide(el,n1,n2);
   if (_vv(n1)==0)
      v(n1) += lambda(_sd3->n());
   if (_vv(n2)==0)
      v(n2) += lambda(_sd3->n());
   if (_vv(n3)==0)
      v(n3) -= lambda(_sd3->n());
   _vv(n1) = _vv(n2) = _vv(n3) = true;
}


size_t Laplace2DMHRT0::NotOnSide(Element* el,
                                 size_t   n1,
                                 size_t   n2)
{
   if ((el->getNodeLabel(1)==n1&&el->getNodeLabel(2)==n2) || (el->getNodeLabel(2)==n1&&el->getNodeLabel(1)==n2))
      return el->getNodeLabel(3);
   else if ((el->getNodeLabel(1)==n1&&el->getNodeLabel(3)==n2) || (el->getNodeLabel(3)==n1&&el->getNodeLabel(1)==n2))
      return el->getNodeLabel(2);
   else if ((el->getNodeLabel(2)==n1&&el->getNodeLabel(3)==n2) || (el->getNodeLabel(3)==n1&&el->getNodeLabel(2)==n2))
      return el->getNodeLabel(1);
   else
      return 0;
}


void Laplace2DMHRT0::getNormals(const Triang3& tr)
{
   _c = tr.getCenter();
   Point<real_t> x11, x12, x21, x22, x31, x32;
   x11 = (*_sd1)(1)->getCoord();
   x12 = (*_sd1)(2)->getCoord();
   x21 = (*_sd2)(1)->getCoord();
   x22 = (*_sd2)(2)->getCoord();
   x31 = (*_sd3)(1)->getCoord();
   x32 = (*_sd3)(2)->getCoord();
   _ce(1) = 0.5*(x11 + x12);
   _ce(2) = 0.5*(x21 + x22);
   _ce(3) = 0.5*(x31 + x32);

   Point<real_t> tangent;
   tangent = x11 - x12;
   _n(1).x = -tangent.y;
   _n(1).y =  tangent.x;
   if (_n(1)*(x11-_c) < 0)
      _n(1) = -_n(1);
   tangent  = x21 - x22;
   _n(2).x = -tangent.y;
   _n(2).y =  tangent.x;
   if (_n(2)*(x21-_c) < 0)
      _n(2) = -_n(2);
   tangent  = x31 - x32;
   _n(3).x = -tangent.y;
   _n(3).y =  tangent.x;
   if (_n(3)*(x31-_c) < 0)
      _n(3) = -_n(3);
}


real_t Laplace2DMHRT0::getIP(const LocalMatrix<real_t,2,2>& A,
                             const Point<real_t>&           a,
                             const Point<real_t>&           b)
{
   return (A(1,1)*b.x + A(1,2)*b.y)*a.x + (A(2,1)*b.x + A(2,2)*b.y)*a.y;
}


Point<real_t> Laplace2DMHRT0::getP(const LocalMatrix<real_t,2,2>& A,
                                   const Point<real_t>&           a)
{
   return Point<real_t>(A(1,1)*a.x+A(1,2)*a.y,A(2,1)*a.x+A(2,2)*a.y);
}

} /* namespace OFELI */
