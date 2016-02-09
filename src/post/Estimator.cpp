/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzani

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

Estimator::Estimator(Mesh& m)
{
   _mesh = &m;
   Err.setMesh(*_mesh);
}


void Estimator::setError(const Vect<real_t>& u)
{
   Vect<real_t> M(_mesh->getNbNodes());
   Vect<Point<real_t> > b(*_mesh);
   the_element = (*_mesh)(1);
   try {
      if (The_element.getShape()==TRIANGLE && The_element.getNbNodes()==3)
         elementT3(u, M, b);
      else
         THROW_RT("This element is not implemented for Estimator calculation.");
   }
   CATCH("Estimator");
   _average = 0;
   mesh_elements(*_mesh)
      _average += Err(element_label);
   _average /= _mesh->getNbElements();
}


void Estimator::elementT3(const Vect<real_t>&         u,
                                Vect<real_t>&         M,
                                Vect<Point<real_t> >& b)
{
   Vect<Point<real_t> > Du(*_mesh,0,ELEMENT_DOF);
   size_t k, n, nb=u.getNbDOF();
   try {
      if (nb==0)
         THROW_RT("This procedure is not allowed with non constant nb of DOF per node.");
   }
   CATCH("Estimator");
   mesh_elements(*_mesh) {
      Triang3 tr(the_element);
      real_t c = tr.getArea()*OFELI_THIRD;
      for (k=1; k<=nb; k++) {
         Du.set(element_label,k,   tr.DSh(1)*u(The_element(1)->n(),k)
                                 + tr.DSh(2)*u(The_element(2)->n(),k)
                                 + tr.DSh(3)*u(The_element(3)->n(),k));
      }
      for (size_t i=1; i<=3; i++) {
         n = The_element(i)->n();
         M(n) += c;
         for (k=1; k<=nb; k++)
            b.add(n,k,c*Du(element_label,k));
      }
   }

   mesh_nodes(*_mesh) {
      for (size_t k=1; k<=nb; k++)
          b.set(node_label,k,b(node_label,k)/M(node_label));
   }

   Point<real_t> g[3];
   mesh_nodes(*_mesh) {
      Triang3 tr(the_element);
      n = element_label;
      for (size_t k=1; k<=nb; k++) {
         g[0] = 0.5*(b(2,k) + b(3,k));
         g[1] = 0.5*(b(3,k) + b(1,k));
         g[2] = 0.5*(b(1,k) + b(2,k));
         real_t c = tr.getArea()*OFELI_THIRD;
         for (size_t i=0; i<3; i++) {
            Err(n) += c*((Du(n).x-g[i].x)*(Du(n).x-g[i].x) +
                         (Du(n).y-g[i].y)*(Du(n).y-g[i].y));
         }
      }
   }
   mesh_elements(*_mesh)
      Err.set(element_label,sqrt(Err(element_label)/nb));
}


ostream& operator<<(      ostream&   s,
                    const Estimator& r)
{
   s << "LOCAL ERROR INFORMATION" << endl << endl;
   s.setf(ios::scientific);
   s << "Errors in elements" << endl << endl;
   mesh_elements(r.getMesh()) {
      s << setw(6) << element_label << "   ";
      s << setprecision(8) << setw(18) << r.Err(element_label) << endl;
   }
   s << endl << "Average Error : " << r.getAverage() << endl;
   s << "Relative Errors in Elements" << endl << endl;
   mesh_nodes(r.getMesh()) {
      s << setw(6) << element_label << "   ";
      s << setprecision(8) << setw(18) << (r.Err(element_label)-r.getAverage())/r.getAverage() << endl;
   }
   return s;
}

} /* namespace OFELI */
