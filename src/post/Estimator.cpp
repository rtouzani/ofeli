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

                        Implementation of class 'Estimator'

  ==============================================================================*/


#include "post/Estimator.h"
#include "shape_functions/Triang3.h"

namespace OFELI {

Estimator::Estimator(Mesh&         m,
                     Vect<real_t>& e)
{
   _mesh = &m;
   _est = &e;
   _est_type = ESTIM_ZZ;
   _est->setMesh(*_mesh,1,ELEMENT_DOF);
   _nb_el = _mesh->getNbElements();
   _nb_nd = _mesh->getNbNodes();
}


void Estimator::setType(EstimatorType t)
{
   _est_type = t;
   if (_est_type==ESTIM_ZZ)
      _est->setSize(_nb_el);
   else if (_est_type==ESTIM_ND_JUMP) {
      _mesh->getAllSides();
      _nb_sd = _mesh->getNbSides();
      _N.setSize(_nb_sd);
      _est->setSize(_nb_sd);
      mesh_sides(*_mesh) {
         Point<real_t> T(The_side(1)->getCoord()-The_side(2)->getCoord());
         _N(side_label) = Point<real_t>(T.y,-T.x);
      }
   }
}


void Estimator::setError(const Vect<real_t>& u)
{
   if (_est_type==ESTIM_ZZ) {
      Vect<real_t> M(_mesh->getNbNodes());
      Vect<Point<real_t> > b(*_mesh);
      the_element = (*_mesh)(1);
      try {
         if (The_element.getShape()==TRIANGLE && The_element.getNbNodes()==3)
            elementT3_ZZ(u, M, b);
         else
            THROW_RT("This element is not implemented for Estimator calculation.");
      }
      CATCH("Estimator");
      _average = 0;
      mesh_elements(*_mesh)
         _average += (*_est)(element_label);
      _average /= _mesh->getNbElements();
   }
   else if (_est_type==ESTIM_ND_JUMP) {
      elementT3_ND_JUMP(u);
      _average = 0;
      mesh_sides(*_mesh)
         _average += (*_est)(side_label);
      _average /= _nb_sd;
      mesh_sides(*_mesh)
         (*_est)(side_label) /= _average;
   }
}


void Estimator::elementT3_ND_JUMP(const Vect<real_t>& u)
{
   mesh_sides(*_mesh) {
      (*_est)(side_label) = 0;
      Element *el1 = The_side.getNeighborElement(1),
              *el2 = The_side.getNeighborElement(2);
      real_t u1 = u((*el1)(1)->n()), u2 = u((*el1)(2)->n()), u3 = u((*el1)(3)->n());
      Triang3 tr(el1);
      real_t dudn1 = (_N(side_label),u1*tr.DSh(1)+u2*tr.DSh(2)+u3*tr.DSh(3));
      if (el2) {
         Triang3 tr(el2);
         u1 = u((*el2)(1)->n()), u2 = u((*el2)(2)->n()), u3 = u((*el2)(3)->n());
         real_t dudn2 = (_N(side_label),u1*tr.DSh(1)+u2*tr.DSh(2)+u3*tr.DSh(3));
         (*_est)(side_label) = fabs(dudn1 - dudn2);
      }
   }
}


void Estimator::elementT3_ZZ(const Vect<real_t>&   u,
                             Vect<real_t>&         M,
                             Vect<Point<real_t> >& b)
{
   Vect<Point<real_t> > Du(*_mesh,0,ELEMENT_DOF);
   size_t nb=u.getNbDOF();
   try {
      if (nb==0)
         THROW_RT("This procedure is not allowed with non constant nb of DOF per node.");
   }
   CATCH("Estimator");
   mesh_elements(*_mesh) {
      Triang3 tr(the_element);
      real_t c = tr.getArea()*OFELI_THIRD;
      for (size_t k=1; k<=nb; k++) {
         Du(element_label,k) =  tr.DSh(1)*u(The_element(1)->n(),k)
                              + tr.DSh(2)*u(The_element(2)->n(),k)
                              + tr.DSh(3)*u(The_element(3)->n(),k);
      }
      for (size_t i=1; i<=3; i++) {
         size_t n = The_element(i)->n();
         M(n) += c;
         for (size_t k=1; k<=nb; k++)
            b(n,k) += c*Du(element_label,k);
      }
   }

   mesh_nodes(*_mesh) {
      for (size_t k=1; k<=nb; k++)
         b(node_label,k) /= M(node_label);
   }
  
   Point<real_t> g[3];
   mesh_elements(*_mesh) {
      Triang3 tr(the_element);
      size_t n = element_label;
      size_t n1 = The_element(1)->n(),
             n2 = The_element(2)->n(),
             n3 = The_element(3)->n();
      for (size_t k=1; k<=nb; k++) {
         g[0] = 0.5*(b(n2,k) + b(n3,k));
         g[1] = 0.5*(b(n3,k) + b(n1,k));
         g[2] = 0.5*(b(n1,k) + b(n2,k));
         real_t c = tr.getArea()*OFELI_THIRD;
         for (size_t i=0; i<3; i++) {
            (*_est)(n) += c*((Du(n,k).x-g[i].x)*(Du(n,k).x-g[i].x) +
                             (Du(n,k).y-g[i].y)*(Du(n,k).y-g[i].y));
         }
      }
   }
   mesh_elements(*_mesh)
      (*_est)(element_label) = sqrt((*_est)(element_label)/_mesh->getNbElements());
}


ostream& operator<<(ostream&         s,
                    const Estimator& r)
{
   s << "LOCAL ERROR INFORMATION" << endl << endl;
   s.setf(ios::scientific);
   s << "Errors in elements" << endl << endl;
   mesh_elements(r.getMesh()) {
      s << setw(6) << element_label << "   ";
      s << setprecision(8) << setw(18) << (*r._est)(element_label) << endl;
   }
   s << endl << "Average Error: " << r.getAverage() << endl;
   s << "Relative Errors in Elements" << endl << endl;
   mesh_nodes(r.getMesh()) {
      s << setw(6) << element_label << "   ";
      s << setprecision(8) << setw(18) << ((*r._est)(element_label)-r.getAverage())/r.getAverage() << endl;
   }
   return s;
}

} /* namespace OFELI */
