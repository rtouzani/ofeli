/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2019 Rachid Touzani

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

                         Implementation of class Bar2DL2

  ==============================================================================*/


#include "equations/solid/Bar2DL2.h"
#include "shape_functions/Line2.h"

namespace OFELI {

Bar2DL2::Bar2DL2(Element* el,
                 real_t   section)
{
   set(el);
   _section = section;
}


void Bar2DL2::set(const Element* el)
{
   _nb_dof = 2;
   Init(el);
   setMaterial();
   Line2 ln(_theElement);
   _length = ln.getLength();
   _center = ln.getCenter();
   ElementNodeCoordinates();
   real_t alpha = atan((_x[1].y-_x[0].y)/(_x[1].x-_x[0].x));
   real_t c = cos(alpha), s = sin(alpha);
   _cc = c*c;
   _ss = s*s;
   _sc = c*s;
   eMat = 0;
   eRHS = 0;
}


void Bar2DL2::LMassToLHS(real_t coef)
{
   real_t c = 0.5*_rho*coef*_length;
   eMat(1,1) += c; eMat(2,2) += c;
   eMat(3,3) += c; eMat(4,4) += c;
}


void Bar2DL2::MassToLHS(real_t coef)
{
   real_t c = coef*_rho*_length*OFELI_SIXTH;
   real_t cc = c*_cc, ss = c*_ss, sc = c*_sc;
   eMat(1,1) += 2*cc; eMat(1,2) += 2*sc; eMat(1,3) +=   cc; eMat(1,4) +=   sc;
   eMat(2,1) += 2*sc; eMat(2,2) += 2*ss; eMat(2,3) +=   sc; eMat(2,4) +=   ss;
   eMat(3,1) +=   cc; eMat(3,2) +=   sc; eMat(3,3) += 2*cc; eMat(3,4) += 2*sc;
   eMat(4,1) +=   sc; eMat(4,2) +=   ss; eMat(4,3) += 2*sc; eMat(4,4) += 2*ss;
}


void Bar2DL2::LMassToRHS(real_t coef)
{
   real_t c = 0.5*_length*_rho*coef;
   eRHS(1) += c*_cc*ePrev(1);
   eRHS(2) += c*_ss*ePrev(2);
   eRHS(3) += c*_cc*ePrev(3);
   eRHS(4) += c*_ss*ePrev(4);
}


void Bar2DL2::MassToRHS(real_t coef)
{
   real_t c = coef*_rho*_length*OFELI_SIXTH;
   real_t cc = c*_cc, ss = c*_ss, sc = c*_sc;
   eRHS(1) += 2*cc*ePrev(1) + 2*sc*ePrev(2) +   cc*ePrev(3) +   sc*ePrev(4);
   eRHS(2) += 2*sc*ePrev(1) + 2*ss*ePrev(2) +   sc*ePrev(3) +   ss*ePrev(4);
   eRHS(3) +=   cc*ePrev(1) +   sc*ePrev(2) + 2*cc*ePrev(3) + 2*sc*ePrev(4);
   eRHS(4) +=   sc*ePrev(1) +   ss*ePrev(2) + 2*sc*ePrev(3) + 2*ss*ePrev(4);
}


void Bar2DL2::Stiffness(real_t coef)
{
   real_t c = _E*_section*coef/_length;
   real_t cc = c*_cc, ss = c*_ss, sc = c*_sc;
   eMat(1,1) += cc; eMat(1,2) += sc; eMat(1,3) -= cc; eMat(1,4) -= sc;
   eMat(2,1) += sc; eMat(2,2) += ss; eMat(2,3) -= sc; eMat(2,4) -= ss;
   eMat(3,1) -= cc; eMat(3,2) -= sc; eMat(3,3) += cc; eMat(3,4) += sc;
   eMat(4,1) -= sc; eMat(4,2) -= ss; eMat(4,3) += sc; eMat(4,4) += ss;
}


void Bar2DL2::BodyRHS(UserData<real_t>& ud)
{
   eRHS(1) += 0.5*_length*ud.BodyForce(_x[0],_time,1);
   eRHS(2) += 0.5*_length*ud.BodyForce(_x[0],_time,2);
   eRHS(3) += 0.5*_length*ud.BodyForce(_x[1],_time,1);
   eRHS(4) += 0.5*_length*ud.BodyForce(_x[1],_time,2);
}


int Bar2DL2::run()
{
   int ret=0;
   _b = new Vect<real_t>(_theMesh->getNbEq());
   if (_analysis==STEADY_STATE) {
      build();
      _A.Factor();
      _A.Solve(*_b,_uu);
      _u->insertBC(*_theMesh,_uu,*_bc);
   }
   else {
      for (_time=0.; _time<=_final_time; _time+=_time_step)
         _time_step = runOneTimeStep();
   }
   delete _b;
   return ret;
}


int Bar2DL2::runOneTimeStep()
{
   build();
   _A.Factor();
   _A.Solve(*_b,_uu);
   _u->insertBC(*_theMesh,_uu,*_bc);
   return 0;
}


void Bar2DL2::build()
{
   _A = 0;
   MESH_EL {
      set(theElement);
      ElementVector(*_u);
      if (_time_scheme==STEADY_STATE) {
         if (_terms&DEVIATORIC)
            Deviator();
         if (_terms&DILATATION)
            Dilatation();
      }
      else if (_time_scheme==TRANSIENT) {
         if (_terms&LUMPED_MASS)
            setLumpedMass();
         if (_terms&MASS)
            setMass();
         if (_terms&DEVIATORIC)
            setDeviator();
         if (_terms&DILATATION)
            setDilatation();
      }
      ElementAssembly(_A);
      ElementAssembly(*_b);
   }
}


real_t Bar2DL2::Stress() const
{
   return (_E*_section*(ePrev(2)-ePrev(1))/_length);
}


void Bar2DL2::getStresses(const Vect<real_t>& u,
                                Vect<real_t>& s)
{
   real_t d11, d12, d21, d22;
   mesh_elements(*_theMesh) {
      set(the_element);
      d11 = u(The_element(1)->n(),1);
      d12 = u(The_element(1)->n(),2);
      d21 = u(The_element(2)->n(),1);
      d22 = u(The_element(2)->n(),2);
      s(element_label,1) = _E*_section*(d21-d11)/_length;
      s(element_label,2) = _E*_section*(d22-d12)/_length;
   }
}


void Bar2DL2::buildEigen(SkSMatrix<real_t>& K,
                         SkSMatrix<real_t>& M)
{
   mesh_elements(*_theMesh) {
      set(the_element);
      real_t c = _E*_section/_length;
      real_t cc = c*_cc, ss = c*_ss, sc = c*_sc;
      eMat(1,1) += cc; eMat(1,2) += sc; eMat(1,3) -= cc; eMat(1,4) -= sc;
      eMat(2,1) += sc; eMat(2,2) += ss; eMat(2,3) -= sc; eMat(2,4) -= ss;
      eMat(3,1) -= cc; eMat(3,2) -= sc; eMat(3,3) += cc; eMat(3,4) += sc;
      eMat(4,1) -= sc; eMat(4,2) -= ss; eMat(4,3) += sc; eMat(4,4) += ss;
      ElementAssembly(K);
   }
   mesh_elements(*_theMesh) {
      set(the_element);
      real_t c = _rho*_length*OFELI_SIXTH;
      real_t cc = c*_cc, ss = c*_ss, sc = c*_sc;
      eMat(1,1) += 2*cc; eMat(1,2) += 2*sc; eMat(1,3) +=   cc; eMat(1,4) +=   sc;
      eMat(2,1) += 2*sc; eMat(2,2) += 2*ss; eMat(2,3) +=   sc; eMat(2,4) +=   ss;
      eMat(3,1) +=   cc; eMat(3,2) +=   sc; eMat(3,3) += 2*cc; eMat(3,4) += 2*sc;
      eMat(4,1) +=   sc; eMat(4,2) +=   ss; eMat(4,3) += 2*sc; eMat(4,4) += 2*ss;
      ElementAssembly(M);
   }
}


void Bar2DL2::buildEigen(SkSMatrix<real_t>& K,
                         Vect<real_t>&      M)
{
   MESH_EL {
      set(theElement);
      real_t c = _E*_section/_length;
      real_t cc = c*_cc, ss = c*_ss, sc = c*_sc;
      eMat(1,1) += cc; eMat(1,2) += sc; eMat(1,3) -= cc; eMat(1,4) -= sc;
      eMat(2,1) += sc; eMat(2,2) += ss; eMat(2,3) -= sc; eMat(2,4) -= ss;
      eMat(3,1) -= cc; eMat(3,2) -= sc; eMat(3,3) += cc; eMat(3,4) += sc;
      eMat(4,1) -= sc; eMat(4,2) -= ss; eMat(4,3) += sc; eMat(4,4) += ss;
      c = 0.5*_rho*_length;
      eRHS(1) += c; eRHS(2) += c;
      eRHS(3) += c; eRHS(4) += c;
      ElementAssembly(K);
      ElementAssembly(M);
   }
}

} /* namespace OFELI */
