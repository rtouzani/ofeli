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

                        Implementation of class NS2DT3BT3
      for 2-D Navier-Stokes equations using P1-Bubble/P1 (Mini) finite element

  ==============================================================================*/


#include "equations/fluid/NS2DT3BT3.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"

namespace OFELI {

NS2DT3BT3::NS2DT3BT3(Element* el)
{
   _equation_name = "Navier-Stokes";
   _finite_element = "2-D, P1+Bubble/P1";
   set(el);
   Misc();
}


NS2DT3BT3::NS2DT3BT3(      Element*      el,
                     const Vect<real_t>& u,
                     const real_t&       time)
{
   _equation_name = "Navier-Stokes";
   _finite_element = "2-D, P1+Bubble/P1";
   set(el);
   ePrev[0] = u(_theElement->getNodeLabel(1));
   ePrev[1] = u(_theElement->getNodeLabel(2));
   ePrev[2] = u(_theElement->getNodeLabel(3));
   Misc();
}


NS2DT3BT3::NS2DT3BT3(Side* sd)
{
   _equation_name = "Navier-Stokes";
   _finite_element = "2-D, P1+Bubble/P1";
   set(sd);
}


NS2DT3BT3::NS2DT3BT3(Side*               sd,
                     const Vect<real_t>& u,
                     const real_t&       time)
{
   _equation_name = "Navier-Stokes";
   _finite_element = "2-D, P1+Bubble/P1";
   set(sd);
   SideVector(u);
}


NS2DT3BT3::NS2DT3BT3(Mesh&         mesh,
                     Vect<real_t>& u,
                     real_t        Re)
          : Equation<real_t,3,9,2,6>(mesh)
{
   _equation_name = "Navier-Stokes";
   _finite_element = "2-D, P1+Bubble/P1";
   _step = 1;
   _u = &u;
   _Re = Re;
   setInput(SOLUTION,u);
   _bc_given = _bf_given = _sf_given = false;
}


NS2DT3BT3::~NS2DT3BT3() { }


void NS2DT3BT3::set(const Element* el)
{
   _nb_dof = 3;
   Init(el);
   setMaterial();
   Triang3 tr(_theElement);
   _area = tr.getArea();
   ElementNodeCoordinates();
   _dSh(1) = tr.DSh(1);
   _dSh(2) = tr.DSh(2);
   _dSh(3) = tr.DSh(3);
   eMat = 0;
   eRHS = 0;
}


void NS2DT3BT3::set(const Side* sd)
{
   _nb_dof = 3;
   Init(sd);
   Line2 ln(_theSide);
   _length = ln.getLength();
   SideNodeCoordinates();
   _dSh(1) = ln.DSh(1);
   _dSh(2) = ln.DSh(2);
   sMat = 0;
   sRHS = 0;
}


int NS2DT3BT3::run()
{
   MESH_EL {
      set(theElement);
      Viscous();
      PressureGradient();
      updateBC(*_theElement,*_bc);
      ElementAssembly(_A);
      ElementAssembly(*_b);
   }
   //   int ret = SolveLinearSystem(*_A,*_b,_uu);
   _u->insertBC(*_theMesh,_uu,*_bc);
   return 0;
}


void NS2DT3BT3::Misc()
{
   Point<real_t> x32=_x[2]-_x[1], x13=_x[0]-_x[2], x21=_x[1]-_x[0];
   real_t a = 0.25*_visc/_area;

   _aa(1,1) = a*x32.NNorm();
   _aa(2,1) = -a*(x13.NNorm() + x13*x21);
   _aa(1,2) = _aa(1,2);
   _aa(3,1) = -a*(x21.NNorm() + x13*x21);
   _aa(1,3) = _aa(3,1);
   _aa(2,2) = a*x13.NNorm();
   _aa(3,2) = a*x13*x21;
   _aa(2,3) = _aa(3,2);
   _aa(3,3) = a*x21.NNorm();
   _aa(4,4) = a*(x13.NNorm() + x21.NNorm() + x13*x21)/90.;

   _bb(1,1) = -OFELI_SIXTH*x32.y; 
   _bb(2,1) = _bb(1,1); 
   _bb(3,1) = _bb(1,1); 
   _bb(1,4) = 0.1*OFELI_TWELVETH*x32.y;
   _bb(1,2) = -OFELI_SIXTH*x13.y; 
   _bb(2,2) = _bb(1,2); 
   _bb(3,2) = -OFELI_SIXTH*x13.y; 
   _bb(2,4) = 0.1*OFELI_TWELVETH*x13.y;
   _bb(1,3) = -OFELI_SIXTH*x21.y; 
   _bb(2,3) = -OFELI_SIXTH*x21.y; 
   _bb(3,3) = -OFELI_SIXTH*x21.y; 
   _bb(3,4) = 0.1*OFELI_TWELVETH*x21.y;

   _cc(1,1) = OFELI_SIXTH*x32.x;
   _cc(2,1) = _cc(1,1);
   _cc(3,1) = _cc(1,1);
   _cc(1,2) = OFELI_SIXTH*x13.x;
   _cc(2,2) = _cc(1,2);
   _cc(3,2) = _cc(1,2);
   _cc(1,3) = OFELI_SIXTH*x21.x;
   _cc(2,3) = _cc(1,3);
   _cc(3,3) = _cc(1,3);
   _cc(1,4) = -0.1*OFELI_TWELVETH*x32.x;
   _cc(2,4) = -0.1*OFELI_TWELVETH*x13.x;
   _cc(3,4) = -0.1*OFELI_TWELVETH*x21.x;
}


void NS2DT3BT3::LMass(real_t coef)
{
   real_t c = coef*_dens*OFELI_THIRD*_area;
   for (size_t i=1; i<=3; i++) {
      eMat(3*i-2,3*i-2) += c;
      eMat(3*i-1,3*i-1) += c;
      eMat(3*i  ,3*i  ) += c;
      eRHS(3*i-2) += c*ePrev(3*i-2);
      eRHS(3*i-1) += c*ePrev(3*i-1);
      eRHS(3*i  ) += c*ePrev(3*i  );
   }
}


void NS2DT3BT3::Viscous(real_t coef)
{
   for (size_t i=0; i<3; i++) {
      eMat(3*i+1,1) = _aa(1,i+1)*coef;
      eMat(3*i+2,2) = _aa(1,i+1)*coef;
      eMat(3*i+1,4) = _aa(2,i+1)*coef;
      eMat(3*i+2,5) = _aa(2,i+1)*coef;
      eMat(3*i+1,7) = _aa(3,i+1)*coef;
      eMat(3*i+2,8) = _aa(3,i+1)*coef;
   }
}


void NS2DT3BT3::PressureGradient(real_t coef)
{
   for (size_t i=0; i<3; i++) {
      eMat(3*i+1,3) = -_bb(1,i+1);
      eMat(3*i+1,6) = -_bb(2,i+1);
      eMat(3*i+1,9) = -_bb(3,i+1);
      eMat(3*i+2,3) = -_cc(1,i+1);
      eMat(3*i+2,6) = -_cc(2,i+1);
      eMat(3*i+2,9) = -_cc(3,i+1);
      eMat(3*i+3,1) =  _bb(i+1,1);
      eMat(3*i+3,2) =  _cc(i+1,1);
      eMat(3*i+3,3) = (_bb(i+1,4)*_bb(1,4) + _cc(i+1,4)*_cc(1,4))/_aa(4,4);
      eMat(3*i+3,4) =  _bb(i+1,2);
      eMat(3*i+3,5) =  _cc(i+1,2);
      eMat(3*i+3,6) = (_bb(i+1,4)*_bb(2,4) + _cc(i+1,4)*_cc(2,4))/_aa(4,4);
      eMat(3*i+3,7) =  _bb(i+1,3);
      eMat(3*i+3,8) =  _cc(i+1,3);
      eMat(3*i+3,9) = (_bb(i+1,4)*_bb(3,4) + _cc(i+1,4)*_cc(3,4))/_aa(4,4);
   }
   for (size_t i=0; i<3; i++) {
      eMat(3*i+1,3) = -_bb(1,i+1);
      eMat(3*i+1,6) = -_bb(2,i+1);
      eMat(3*i+1,9) = -_bb(3,i+1);
      eMat(3*i+2,3) = -_cc(1,i+1);
      eMat(3*i+2,6) = -_cc(2,i+1);
      eMat(3*i+2,9) = -_cc(3,i+1);
      eMat(3*i+3,1) =  _bb(i+1,1);
      eMat(3*i+3,2) =  _cc(i+1,1);
      eMat(3*i+3,3) = (_bb(i+1,4)*_bb(1,4) + _cc(i+1,4)*_cc(1,4))/_aa(4,4);
      eMat(3*i+3,4) =  _bb(i+1,2);
      eMat(3*i+3,5) =  _cc(i+1,2);
      eMat(3*i+3,6) = (_bb(i+1,4)*_bb(2,4) + _cc(i+1,4)*_cc(2,4))/_aa(4,4);
      eMat(3*i+3,7) =  _bb(i+1,3);
      eMat(3*i+3,8) =  _cc(i+1,3);
      eMat(3*i+3,9) = (_bb(i+1,4)*_bb(3,4) + _cc(i+1,4)*_cc(3,4))/_aa(4,4);
   }
}


void NS2DT3BT3::RHS_Convection(real_t coef)
{
   Point<real_t> du = ePrev(1)*_dSh(1) + ePrev(3)*_dSh(2) + ePrev(5)*_dSh(3);
   Point<real_t> dv = ePrev(2)*_dSh(1) + ePrev(4)*_dSh(2) + ePrev(6)*_dSh(3);
   real_t c = coef*OFELI_THIRD*_area;
   for (size_t i=1; i<=3; i++) {
      eRHS(3*i-2) -= c*(ePrev(2*i-1)*du.x + ePrev(2*i)*du.y);
      eRHS(3*i-1) -= c*(ePrev(2*i-1)*dv.x + ePrev(2*i)*dv.y);
   }
}


void NS2DT3BT3::BodyRHS(UserData<real_t>& ud)
{
   for (size_t i=1; i<=3; i++) {
      real_t fx = ud.BodyForce(_x[i-1],_time,1);
      real_t fy = ud.BodyForce(_x[i-1],_time,2);
      eRHS(3*i-2) += OFELI_THIRD*_area*fx;
      eRHS(3*i-1) += OFELI_THIRD*_area*fy;
   }
}


void NS2DT3BT3::BoundaryRHS(UserData<real_t>& ud)
{
   for (size_t i=1; i<=2; i++) {
      real_t fx = ud.SurfaceForce(_x[i-1],_theSide->getCode(1),_time,1);
      real_t fy = ud.SurfaceForce(_x[i-1],_theSide->getCode(2),_time,2);
      sRHS(3*i-2) += 0.5*_length*fx;
      sRHS(3*i-1) += 0.5*_length*fy;
   }
}

} /* namespace OFELI */
