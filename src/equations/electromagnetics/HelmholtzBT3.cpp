/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2024 Rachid Touzani

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

                         Implementation of class HelmholtzBT3
                       for Helmholtz Equation in a Bounded Domain
                         using 3-node triangular finite element

  ==============================================================================*/


#include "equations/electromagnetics/HelmholtzBT3.h"
#include "linear_algebra/Vect_impl.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"

namespace OFELI {

HelmholtzBT3::HelmholtzBT3()
             : Equation<3,6,2,4>()
{
}


HelmholtzBT3::HelmholtzBT3(Mesh& ms)
             : Equation<3,6,2,4>(ms)
{
   if (Equa::_nb_dof!=2)
      throw OFELIException("In HelmholtzBT3::HelmholtzBT3(..): Nodes must have "
                            "2 degrees of freedom because of complex-valued formulation.");
   _omega = 1.;
   _equation_name = "Helmholtz equation in a bounded domain";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _theMesh->removeImposedDOF();
   setMatrixType(SKYLINE);
   setSolver(DIRECT_SOLVER);
}


HelmholtzBT3::HelmholtzBT3(Mesh&         ms,
                           Vect<real_t>& u)
             : Equation<3,6,2,4>(ms,u)
{
   if (Equa::_nb_dof!=2)
      throw OFELIException("In HelmholtzBT3::HelmholtzBT3(..): Nodes must have "
                            "2 degrees of freedom because of complex-valued formulation.");
   _omega = 1.;
   _equation_name = "Helmholtz equation in a bounded domain";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   _theMesh->removeImposedDOF();
   setMatrixType(SKYLINE);
   setSolver(DIRECT_SOLVER);
}


HelmholtzBT3::~HelmholtzBT3()
{
}


void HelmholtzBT3::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   setMaterial();
   Triang3 tr(_theElement);
   _el_geo.area = tr.getArea();
   _el_geo.center = tr.getCenter();
   _dSh = tr.DSh();
   ElementNodeCoordinates();
   ElementNodeVector(*_u,_eu);
   if (_omega_set)
      _omega = _omega_fct(_el_geo.center,0.);
   eMat = 0.;
   eRHS = 0.;
}


void HelmholtzBT3::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Line2 ln(sd);
   SideNodeCoordinates();
   _el_geo.length = ln.getLength();
   sMat = 0.;
   sRHS = 0.;
}


void HelmholtzBT3::LHS()
{
   real_t c = _omega*_el_geo.area*OFELI_TWELVETH;
   for (size_t i=1; i<=3; i++) {
      Point<real_t> a = _el_geo.area*_dSh[i-1];
      for (size_t j=1; j<=3; j++) {
         eMat(2*i-1,2*j-1) = (a,_dSh[j-1]) - c;
         eMat(2*i  ,2*j  ) = (a,_dSh[j-1]) - c;
      }
      eMat(2*i-1,2*i-1) -= c;
      eMat(2*i  ,2*i  ) -= c;
   }
   eA0 = eMat;
}


void HelmholtzBT3::BodyRHS(Vect<real_t>& f)
{
   real_t c = OFELI_THIRD*_el_geo.area;
   for (size_t i=1; i<=3; i++) {
      eRHS(2*i-1) += c*f(2*(*_theElement)(i)->n()-1);
      eRHS(2*i  ) += c*f(2*(*_theElement)(i)->n()  );
   }
}


void HelmholtzBT3::BoundaryRHS(Vect<real_t>& f)
{
   if (_theSide->getCode(1)>0) {
      sRHS(1) += 0.5*_el_geo.length*f(2*(*_theSide)(1)->n()-1);
      sRHS(2) += 0.5*_el_geo.length*f(2*(*_theSide)(1)->n()  );
      sRHS(3) += 0.5*_el_geo.length*f(2*(*_theSide)(2)->n()-1);
      sRHS(4) += 0.5*_el_geo.length*f(2*(*_theSide)(2)->n()  );
   }
}


void HelmholtzBT3::build()
{
   MESH_EL {
      set(the_element);
      LHS();
      Equa::_A->Assembly(The_element,eMat.get());
      if (Equa::_bf!=nullptr)
         BodyRHS(*Equa::_bf);
      if (Equa::_bc!=nullptr)
         this->updateBC(*_theElement,*Equa::_bc);
      Equa::_b->Assembly(The_element,eRHS.get());
   }
   MESH_SD {
      set(the_side);
      if (Equa::_sf!=nullptr)
         BoundaryRHS(*_sf);
      Equa::_b->Assembly(The_side,sRHS.get());
   }
}

} /* namespace OFELI */
