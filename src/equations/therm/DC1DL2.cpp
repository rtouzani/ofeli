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

	           Class DC1DL2 : Diffusion-Convection Element
                  using 2-Node Line Finite element in two dimensions

  ==============================================================================*/


#include "equations/therm/DC1DL2.h"
#include "shape_functions/Line2.h"

namespace OFELI {

DC1DL2::DC1DL2()
{
   _equation_name = "Diffusion/Convection";
   _finite_element = "1-D, 2-Node Lines (P1)";
}


DC1DL2::DC1DL2(const Element* el)
{
   set(el);
   _equation_name = "Diffusion/Convection";
   _finite_element = "1-D, 2-Node Lines (P1)";
}


DC1DL2::DC1DL2(const Element*      el,
               const Vect<real_t>& u,
                     real_t        time)
{
   set(el);
   _time = time;
   ElementVector(u);
   setMaterial();
   _equation_name = "Diffusion/Convection";
   _finite_element = "1-D, 2-Node Lines (P1)";
}


DC1DL2::DC1DL2(const Element*      el,
               const Vect<real_t>& u,
                     real_t        time,
                     real_t        deltat,
                     int           scheme)
{
   set(el);
   _time = time;
   ElementVector(u);
   setMaterial();
   _time_step = deltat;
   setTimeIntegration(scheme);
   _equation_name = "Diffusion/Convection";
   _finite_element = "1-D, 2-Node Lines (P1)";
}


DC1DL2::~DC1DL2() { }


void DC1DL2::set(const Element* el)
{
   _nb_dof = 1;
   Init(el);
   setMaterial();
   Line2 ln(_theElement);
   _center = ln.getCenter();
   _dSh[0] = ln.DSh(1);
   _dSh(2) = ln.DSh(2);
   ElementNodeCoordinates();
   eMat = 0;
   eRHS = 0;
}


void DC1DL2::setInput(EqDataType    opt,
                      Vect<real_t>& u)
{
   AbsEqua<real_t>::setInput(opt,u);
   if (opt==VELOCITY_FIELD)
      _vel = &u;
}


void DC1DL2::LCapacityToLHS(real_t coef)
{
   real_t c = 0.5*_length*_rhocp*coef;
   eMat(1,1) += c;
   eMat(2,2) += c;
}


void DC1DL2::LCapacityToRHS(real_t coef)
{
   real_t c = 0.5*_length*_rhocp*coef;
   eRHS(1) += c*ePrev(1);
   eRHS(2) += c*ePrev(2);
}


void DC1DL2::CapacityToLHS(real_t coef)
{
   real_t c = OFELI_SIXTH*_length*_rhocp*coef;
   eMat(1,1) += 2*c;
   eMat(2,2) += 2*c;
   eMat(1,2) +=   c; 
   eMat(2,1) +=   c;
}


void DC1DL2::CapacityToRHS(real_t coef)
{
   real_t c = OFELI_SIXTH*_length*_rhocp*coef;
   eRHS(1) += c*(2*ePrev(1) + ePrev(2));
   eRHS(2) += c*(2*ePrev(2) + ePrev(1));
}


void DC1DL2::Diffusion(real_t coef)
{
   eMat(1,1) += coef*_diff*_length*_dSh[0].x*_dSh[0].x;
   eMat(1,2) += coef*_diff*_length*_dSh[0].x*_dSh[1].x;
   eMat(2,1) += coef*_diff*_length*_dSh[1].x*_dSh[0].x;
   eMat(2,2) += coef*_diff*_length*_dSh[1].x*_dSh[1].x;
}


void DC1DL2::DiffusionToRHS(real_t coef)
{
   real_t u = _dSh[0].x*ePrev(1) + _dSh[1].x*ePrev(2);
   eRHS(1) -= coef*_diff*_length*_dSh[0].x*u;
   eRHS(2) -= coef*_diff*_length*_dSh[1].x*u;
}


void DC1DL2::Convection(const real_t& v,
                              real_t  coef)
{
   eMat(1,1) += 0.5*_length*coef*v*_dSh[0].x;
   eMat(1,2) += 0.5*_length*coef*v*_dSh[1].x;
   eMat(2,1) += 0.5*_length*coef*v*_dSh[0].x;
   eMat(2,2) += 0.5*_length*coef*v*_dSh[1].x;
}


void DC1DL2::Convection(real_t coef)
{
   real_t c = 0.25*_length*coef*((*_vel)((*_theElement)(1)->n()) + (*_vel)((*_theElement)(2)->n()));
   eMat(1,1) += c*_dSh[0].x;
   eMat(1,2) += c*_dSh[1].x;
   eMat(2,1) += c*_dSh[0].x;
   eMat(2,2) += c*_dSh[1].x;
}


void DC1DL2::Convection(const Vect<real_t>& v,
                              real_t        coef)
{
   real_t c = 0.25*_length*coef*(v((*_theElement)(1)->n()) + v((*_theElement)(2)->n()));
   eMat(1,1) += c*_dSh[0].x;
   eMat(1,2) += c*_dSh[1].x;
   eMat(2,1) += c*_dSh[0].x;
   eMat(2,2) += c*_dSh[1].x;
}


void DC1DL2::ConvectionToRHS(const real_t& v,
                                   real_t  coef)
{
   real_t u = _dSh[0].x*ePrev(1) + _dSh[1].x*ePrev(2);
   eRHS(1) -= 0.5*_length*coef*v*u;
   eRHS(2) -= 0.5*_length*coef*v*u;
}


void DC1DL2::ConvectionToRHS(real_t coef)
{
   real_t u = _dSh[0].x*ePrev(1) + _dSh[1].x*ePrev(2);
   eRHS(1) -= 0.25*_length*coef*((*_vel)((*_theElement)(1)->n())+(*_vel)((*_theElement)(2)->n()))*u;
   eRHS(2) -= 0.25*_length*coef*((*_vel)((*_theElement)(1)->n())+(*_vel)((*_theElement)(2)->n()))*u;
}


void DC1DL2::BodyRHS(UserData<real_t>& ud,
                     real_t            coef)
{
   eRHS(1) +=  0.5*coef*_length*ud.BodyForce(_x[0],_time);
   eRHS(2) +=  0.5*coef*_length*ud.BodyForce(_x[1],_time);
}


void DC1DL2::BodyRHS(const Vect<real_t>& bf,
                           int           opt)
{
   if (opt==LOCAL_ARRAY) {
      eRHS(1) += 0.5*bf(1)*_length;
      eRHS(2) += 0.5*bf(2)*_length;
   }
   else {
     eRHS(1) += 0.5*_length*bf((*_theElement)(1)->n());
     eRHS(2) += 0.5*_length*bf((*_theElement)(2)->n());
   }
}


real_t DC1DL2::Flux() const
{
   return _diff*(ePrev(1)*_dSh(1).x + ePrev(2)*_dSh(2).x);
}


void DC1DL2::build()
{
   _A = 0;
   if (_time_scheme==FORWARD_EULER)
      _theta = 0;
   else if (_time_scheme==BACKWARD_EULER)
      _theta = 1;
   else if (_time_scheme==CRANK_NICOLSON)
      _theta = 0.5;
   MESH_EL {
      set(theElement);
      ElementVector(_uu);
      if (_terms&CAPACITY)
         setCapacity();
      if (_terms&LUMPED_CAPACITY)
         setLumpedCapacity();
      if (_terms&DIFFUSION)
         setDiffusion();
      if (_terms&CONVECTION)
         setConvection();
      updateBC(*_theElement,*_bc);
      ElementAssembly(_A);
      ElementAssembly(*_b);
   }
}

} /* namespace OFELI */
