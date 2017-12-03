/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2018 Rachid Touzani

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

                      Class Laplace2DFVT : Laplace Equation
              using classical triangular finite volume scheme in 2-D

  ==============================================================================*/


#include "equations/laplace/Laplace2DFVT.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"

using std::cout;

namespace OFELI {

Laplace2DFVT::Laplace2DFVT(Mesh&         ms,
                           Vect<real_t>& b,
                           Vect<real_t>& u)
{
   _theMesh = &ms;
   _A = new SpMatrix<real_t>(_theMesh->getNbEq());
   _b = &b;
   _u = &u;
   _in_place = true;
}


Laplace2DFVT::Laplace2DFVT(Mesh&             ms,
                           SpMatrix<real_t>& A,
                           Vect<real_t>&     b)
{
   _theMesh = &ms;
   _A = &A;
   _b = &b;
   _in_place = false;
}

Laplace2DFVT::~Laplace2DFVT()
{
   if (_in_place)
      delete _A;
}


int Laplace2DFVT::checkDelaunay(int verb)
{
   int ret=0, rb=0;
   if (verb>1)
      cout << "Checking for Delaunay mesh ..." << endl;
   MESH_EL {
      Triang3 tr(theElement);
      _xc1 = tr.getCircumcenter();
      if (tr.isStrictlyIn(_xc1)==false) {
         if (theElement->isOnBoundary()) {
            if (verb>1) {
               cout << "Element on boundary " << theElement->n();
               cout << ": Circumcenter does not lie in the interior of the triangle !" << endl;
            }
            rb++;
            ret++;
         }
         else {
            if (verb>1) {
               cout << "Interior Element    " << theElement->n();
               cout << ": Circumcenter does not lie in the interior of the triangle !" << endl;
            }
            ret++;
         }
      }
   }
   if (ret && verb) {
      cout << "Number of non regular triangles on boundary: " << rb  << endl;
      cout << "Total number of non regular triangles:       " << ret << endl;
   }
   return ret;
}


void Laplace2DFVT::build(const Vect<real_t>& f)
{
   MESH_SD {
      Element *e1=theSide->getNeighborElement(1), *e2=theSide->getNeighborElement(2);
      Line2 ln(theSide);
      _ll = ln.getLength();
      if (theSide->isOnBoundary()) {
         _xc2 = ln.getCenter();
         if (e1) {
            Triang3 tr(e1);
            _m1 = tr.getArea();
            _xc1 = tr.getCircumcenter();
         }
         if (e2) {
            Triang3 tr(e2);
            _m1 = tr.getArea();
            _xc1 = tr.getCircumcenter();
         }
      }
      else {
         Triang3 tr1(e1), tr2(e2);
         _m1 = tr1.getArea();
         _xc1 = tr1.getCircumcenter();
         _m2 = tr2.getArea();
         _xc2 = tr2.getCircumcenter();
      }
      _dd = Distance(_xc1,_xc2);
      LHS(e1,e2);
   }
   RHS(f);
}


int Laplace2DFVT::run(const Vect<real_t>& f)
{
   build(f);
   LinearSolver<real_t> ls(*_A,*_b,*_u);
   ls.setTolerance(1.e-7);
   int nb_it = ls.solve(CG_SOLVER,ILU_PREC);
   return nb_it;
}


void Laplace2DFVT::LHS(const Element* e1,
                       const Element* e2)
{
   size_t i, j;
   if (e1 && e2) {
      i=e1->n(), j=e2->n();
      _A->add(i,i, _ll/_dd);
      _A->add(j,j, _ll/_dd);
      _A->add(i,j,-_ll/_dd);
      _A->add(j,i,-_ll/_dd);
   }
   else if (e1)
      _A->add(e1->n(),e1->n(),_ll/_dd);
   else
      _A->add(e2->n(),e2->n(),_ll/_dd);
}


void Laplace2DFVT::RHS(const Vect<real_t>& f)
{
   MESH_EL
      _b->set(theElementLabel,Triang3(theElement).getArea()*f(theElementLabel));
}

} /* namespace OFELI */
