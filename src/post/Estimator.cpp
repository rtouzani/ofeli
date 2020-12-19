/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

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
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

Estimator::Estimator(Mesh& m)
{
   _mesh = &m;
   _est_type = ESTIM_ZZ;
   _nb_el = _mesh->getNbElements();
   _nb_nd = _mesh->getNbNodes();
   _nb_dof = (*_mesh)[1]->getNbDOF();
}


void Estimator::setType(EstimatorType t)
{
   _est_type = t;
   if (_est_type==ESTIM_ZZ) {
      _nd_I.setSize(_nb_nd);
      _el_I.setSize(_nb_el);
   }
   else if (_est_type==ESTIM_ND_JUMP) {
      _mesh->getAllSides();
      _nb_sd = _mesh->getNbSides();
      _N.setSize(_nb_sd);
      _sd_I.setSize(_nb_sd);
      side_loop(_mesh) {
         Point<real_t> T(The_side(1)->getCoord()-The_side(2)->getCoord());
         _N(side_label) = Point<real_t>(T.y,-T.x);
      }
   }
}


void Estimator::setSolution(const Vect<real_t>& u)
{
   if (_est_type==ESTIM_ZZ) {
      Vect<real_t> M(_nb_nd);
      Vect<Point<real_t> > b(*_mesh);
      the_element = (*_mesh)(1);
      if (The_element.getShape()==TRIANGLE && The_element.getNbNodes()==3)
         elementT3_ZZ(u,b);
      else
         throw OFELIException("Estimator::This element is not implemented for Estimator calculation.");
      _average = 0;
      element_loop(_mesh)
         _average += _el_I(element_label);
      _average /= _nb_el;
   }
   else if (_est_type==ESTIM_ND_JUMP) {
      elementT3_ND_JUMP(u);
      _average = 0;
      side_loop(_mesh)
         _average += _sd_I(side_label);
      _average /= _nb_sd;
      side_loop(_mesh)
         _sd_I(side_label) /= _average;
   }
}


void Estimator::getNodeWiseIndex(Vect<real_t>& I)
{
   I.setSize(_nb_nd,_nb_dof);
   I = _nd_I;
   real_t emax = I.getNormMax();
   for (size_t i=0; i<_nb_nd; i++)
      I(i) /= emax;
}


void Estimator::getElementWiseIndex(Vect<real_t>& I)
{
   I.setSize(_nb_el,_nb_dof);
   I = _el_I;
}


void Estimator::getSideWiseIndex(Vect<real_t>& I)
{
   I.setSize(_nb_sd,_nb_dof);
   I = _sd_I;
}


void Estimator::elementT3_ND_JUMP(const Vect<real_t>& u)
{
   _sd_I.setSize(_nb_sd);
   side_loop(_mesh) {
      _sd_I(side_label) = 0;
      Element *el1 = The_side.getNeighborElement(1),
              *el2 = The_side.getNeighborElement(2);
      real_t u1 = u((*el1)(1)->n()), u2 = u((*el1)(2)->n()), u3 = u((*el1)(3)->n());
      Triang3 tr(el1);
      std::vector<Point<real_t> > dsh = tr.DSh();
      real_t dudn1 = (_N(side_label),u1*dsh[0]+u2*dsh[1]+u3*dsh[1]);
      if (el2) {
         Triang3 tr(el2);
         std::vector<Point<real_t> > dsh = tr.DSh();
         u1 = u((*el2)(1)->n()), u2 = u((*el2)(2)->n()), u3 = u((*el2)(3)->n());
         real_t dudn2 = (_N(side_label),u1*dsh[0]+u2*dsh[1]+u3*dsh[2]);
         _sd_I(side_label) = fabs(dudn1 - dudn2);
      }
   }
}


void Estimator::elementT3_ZZ(const Vect<real_t>&   u,
                             Vect<Point<real_t> >& b)
{
   _el_I.setSize(_nb_el);
   _nd_I.setSize(_nb_nd);
   Vect<real_t> E(_nb_el), M(_nb_nd);
   M = 0;
   Vect<Point<real_t> > Du(_nb_el,_nb_dof);
   if (_nb_dof==0)
      throw OFELIException("Estimator::This procedure is not allowed with"
                           " non constant nb of DOF per node.");
   element_loop(_mesh) {
      Triang3 tr(the_element);
      std::vector<Point<real_t> > dsh = tr.DSh();
      real_t c = tr.getArea()*OFELI_THIRD;
      for (size_t k=1; k<=_nb_dof; k++) {
        Du(element_label,k) =  dsh[0]*u(The_element(1)->n(),k)
                             + dsh[1]*u(The_element(2)->n(),k)
                             + dsh[2]*u(The_element(3)->n(),k);
      }
      for (size_t i=1; i<=3; i++) {
         M(The_element(i)->n()) += c;
         for (size_t k=1; k<=_nb_dof; k++)
            b(The_element(i)->n(),k) += c*Du(element_label,k);
      }
   }

   node_loop(_mesh) {
      for (size_t k=1; k<=_nb_dof; k++)
         b(node_label,k) /= M(node_label);
   }
  
   Point<real_t> g[3];
   element_loop(_mesh) {
      Triang3 tr(the_element);
      size_t n = element_label;
      size_t n1 = The_element(1)->n(),
             n2 = The_element(2)->n(),
             n3 = The_element(3)->n();
      for (size_t k=1; k<=_nb_dof; k++) {
         g[0] = 0.5*(b(n2,k) + b(n3,k));
         g[1] = 0.5*(b(n3,k) + b(n1,k));
         g[2] = 0.5*(b(n1,k) + b(n2,k));
         real_t c = tr.getArea()*OFELI_THIRD;
         for (size_t i=0; i<3; i++) {
            _el_I(n) += c*((Du(n,k).x-g[i].x)*(Du(n,k).x-g[i].x) +
                           (Du(n,k).y-g[i].y)*(Du(n,k).y-g[i].y));
         }
      }
   }

   /*   _nd_I = 0;
   mesh_elements(*_mesh) {
      _el_I(element_label) = sqrt(_el_I(element_label)/_nb_el);
      Triang3 tr(the_element);
      real_t c = tr.getArea()*OFELI_THIRD;
      _nd_I(The_element(1)->n()) += c*_el_I(element_label);
      _nd_I(The_element(2)->n()) += c*_el_I(element_label);
      _nd_I(The_element(3)->n()) += c*_el_I(element_label);
   }
   mesh_nodes(*_mesh)
      _nd_I(node_label) /= M(node_label);
   */

}


ostream& operator<<(ostream&         s,
                    const Estimator& r)
{
   s << "LOCAL ERROR INFORMATION" << endl << endl;
   s.setf(ios::scientific);
   s << "Errors in elements" << endl << endl;
   element_loop(&(r.getMesh())) {
      s << setw(6) << element_label << "   ";
      s << setprecision(8) << setw(18) << r._el_I(element_label) << endl;
   }
   s << endl << "Average Error: " << r.getAverage() << endl;
   s << "Relative Errors in Elements" << endl << endl;
   node_loop(&(r.getMesh())) {
      s << setw(6) << element_label << "   ";
      s << setprecision(8) << setw(18) << (r._el_I(element_label)-r.getAverage())/r.getAverage() << endl;
   }
   return s;
}

} /* namespace OFELI */
