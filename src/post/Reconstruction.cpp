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

                    Implementation of class 'Reconstruction'

  ==============================================================================*/


#include "post/Reconstruction.h"
#include "shape_functions/Line2.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Quad4.h"
#include "shape_functions/Tetra4.h"
#include "linear_algebra/Vect_impl.h"
#include "util/util.h"

namespace OFELI {

Reconstruction::Reconstruction(const Mesh& ms)
{
   setMesh(ms);
}


void Reconstruction::setMesh(const Mesh& ms)
{
  _theMesh = &ms;
}

void Reconstruction::P0toP1(const Vect<real_t>& u,
                            Vect<real_t>&       v)
{
   size_t nb_dof = u.getNbDOF(), dim = _theMesh->getDim();
   v.setSize(_theMesh->getNbNodes(),nb_dof);
   v = 0;
   _M.setSize(_theMesh->getNbNodes());
   _M = 0;
   MESH_EL {
      size_t n = element_label;
      if (dim==1 && The_element.getShape()==LINE) {
         real_t a = Line2(the_element).getLength();
         size_t n1=The_element(1)->n(), n2=The_element(2)->n();
         _M(n1) += a; _M(n2) += a;
         for (size_t k=1; k<=nb_dof; k++) {
            v(n1,k) += a*u(n,k);
            v(n2,k) += a*u(n,k);
         }
      }
      else if (dim==2 && The_element.getShape()==TRIANGLE) {
         real_t a = Triang3(the_element).getArea();
         size_t n1=The_element(1)->n(), n2=The_element(2)->n(), n3=The_element(3)->n();
         _M(n1) += a; _M(n2) += a; _M(n3) += a;
         for (size_t k=1; k<=nb_dof; k++) {
            v(n1,k) += a*u(n,k);
            v(n2,k) += a*u(n,k);
            v(n3,k) += a*u(n,k);
         }
      }
      else if (dim==2 && The_element.getShape()==QUADRILATERAL) {
         Quad4 q(the_element);
         q.setLocal(Point<real_t>(0.,0.));
         real_t a = fabs(q.getDet());
         size_t n1=The_element(1)->n(), n2=The_element(2)->n(),
                n3=The_element(3)->n(), n4=The_element(4)->n();
         _M(n1) += a; _M(n2) += a;
         _M(n3) += a; _M(n4) += a;
         for (size_t k=1; k<=nb_dof; k++) {
            v(n1,k) += a*u(n,k);
            v(n2,k) += a*u(n,k);
            v(n3,k) += a*u(n,k);
            v(n4,k) += a*u(n,k);
         }
      }
      else if (dim==3 && The_element.getShape()==TETRAHEDRON) {
         real_t a = fabs(Tetra4(the_element).getVolume());
         size_t n1=The_element(1)->n(), n2=The_element(2)->n(),
                n3=The_element(3)->n(), n4=The_element(4)->n();
         _M(n1) += a; _M(n2) += a;
         _M(n3) += a; _M(n4) += a;
         for (size_t k=1; k<=nb_dof; k++) {
            v(n1,k) += a*u(n,k);
            v(n2,k) += a*u(n,k);
            v(n3,k) += a*u(n,k);
            v(n4,k) += a*u(n,k);
         }
      }
      else
         throw OFELIException("Reconstruction::P0toP1(...): Not valid for element: " +
                              itos(element_label));
   }

   MESH_ND {
      for (size_t k=1; k<=nb_dof; k++)
         v(node_label,k) /= _M(node_label);
   }
}


void Reconstruction::DP1toP1(const Vect<real_t>& u,
                             Vect<real_t>&       v)
{
   _M.setSize(_theMesh->getNbNodes());
   _M = 0;
   v = 0;
   MESH_EL {
      size_t n = element_label;
      if (_theMesh->getDim()==2 && TheElement.getShape()==TRIANGLE) {
         real_t a = Triang3(the_element).getArea();
         _M(The_element(1)->n()) += a;
         _M(The_element(2)->n()) += a;
         _M(The_element(3)->n()) += a;
         v(The_element(1)->n()) += 0.25*a*(2*u(n,1) + u(n,2) + u(n,3));
         v(The_element(2)->n()) += 0.25*a*(2*u(n,2) + u(n,3) + u(n,1));
         v(The_element(3)->n()) += 0.25*a*(2*u(n,3) + u(n,1) + u(n,2));
      }
      else
         throw OFELIException("Reconstruction::DP1toP1(...): Not valid for element: " + itos(n));
   }
   for (size_t i=1; i<=_theMesh->getNbNodes(); i++)
      v(i) /= _M(i);
}

} /* namespace OFELI */
