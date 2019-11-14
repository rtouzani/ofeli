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

                        Implementation of class NSP2DQ41
             for 2-D Navier-Stokes equations using penalty formulation and
                               Q1/P0 finite element

  ==============================================================================*/


#include "equations/fluid/NSP2DQ41.h"
#include "post/Reconstruction.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/SkSMatrix_impl.h"

namespace OFELI {


NSP2DQ41::NSP2DQ41(Mesh& ms)
         : Equation<real_t,4,8,2,4>(ms), _p(nullptr), _quad(nullptr), _ln(nullptr),
           _penal(1.e07)
{
  //   setMatrixType(SPARSE);
  //   setSolver(CG_SOLVER,DILU_PREC);
   setMatrixType(SKYLINE|SYMMETRIC);
   setSolver(DIRECT_SOLVER);
   _equation_name = "Navier-Stokes";
   _finite_element = "2-D, 4-Node Quadrilaterals (Q1/P0), Penalty";
   _dens = _visc = 0.;
   _terms = VISCOSITY|DILATATION|BODY_FORCE|TRACTION;
}


NSP2DQ41::NSP2DQ41(Mesh&         ms,
                   Vect<real_t>& u)
         : Equation<real_t,4,8,2,4>(ms,u), _p(nullptr), _quad(nullptr), _ln(nullptr),
           _penal(1.e07)
{
   //   setMatrixType(SPARSE);
   //   setSolver(CG_SOLVER,DILU_PREC);
   setMatrixType(SKYLINE|SYMMETRIC);
   setSolver(DIRECT_SOLVER);
   _equation_name = "Navier-Stokes";
   _finite_element = "2-D, 4-Node Quadrilaterals (Q1/P0), Penalty";
   _dens = _visc = 0.;
   _terms = VISCOSITY|DILATATION|BODY_FORCE|TRACTION;
   _terms = MASS|VISCOSITY|DILATATION;
}


NSP2DQ41::~NSP2DQ41()
{
   if (_quad != nullptr)
      delete _quad;
   if (_ln != nullptr)
      delete _ln;
}


void NSP2DQ41::set(const Element *el)
{
   _theElement = el, _theSide = nullptr;
   setMaterial();
   if (_quad != nullptr)
      delete _quad, _quad = nullptr;
   if (_ln != nullptr)
      delete _ln, _ln = nullptr;
   _dens = _visc = 1.;
   _quad = new Quad4(_theElement);
   _quad->atGauss(2,_sh,_dSh,_wg);
   _xl[0] = _yl[0] = _yl[1] = _xl[3] = -1.;
   _xl[1] = _xl[2] = _yl[2] = _yl[3] =  1.;
   ElementNodeCoordinates();
   eMat = 0;
   eRHS = 0;
}


void NSP2DQ41::set(const Side *sd)
{
   _theElement = nullptr, _theSide = sd;
   if (_quad != nullptr)
      delete _quad, _quad = nullptr;
   if (_ln != nullptr)
      delete _ln, _ln = nullptr;
   _ln = new Line2(_theSide);
   _xl[0] = _yl[0] = _yl[1] = _xl[3] = -1.;
   _xl[1] = _xl[2] = _yl[2] = _yl[3] =  1.;
   SideNodeCoordinates();
   _dSh = _ln->DSh();
   sMat = 0;
   sRHS = 0;
}


void NSP2DQ41::setInput(EqDataType    opt,
                        Vect<real_t>& u)
{
   AbsEqua<real_t>::setInput(opt,u);
   if (opt==PRESSURE_FIELD)
      _p = &u;
}


void NSP2DQ41::LMass(real_t coef)
{
   for (size_t i=0; i<4; i++) {
      _quad->setLocal(Point<real_t>(_xl[i],_yl[i]));
      real_t c = _dens*coef*_quad->getDet();
      eMat(2*i+1,2*i+1) += c;
      eMat(2*i+2,2*i+2) += c;
      eRHS(2*i+1) += c*_eu(2*i+1);
      eRHS(2*i+2) += c*_eu(2*i+2);
   }
}


void NSP2DQ41::Mass(real_t coef)
{
   coef = 1;
   std::cerr << "Sorry, consistent mass matrix is not implemented !\n";
}


void NSP2DQ41::Viscous(real_t coef)
{
   for (size_t k=0; k<4; ++k) {
      real_t c = coef*_visc*_wg[k];
      for (size_t i=0; i<4; i++) {
         Point<real_t> a = c*_dSh[4*i+k];
         for (size_t j=0; j<4; j++) {
            Point<real_t> b = _dSh[4*j+k];
            eMat(2*i+1,2*j+1) += 2*a.x*b.x + a.y*b.y;
            eMat(2*i+1,2*j+2) += a.y*b.x;
            eMat(2*i+2,2*j+1) += a.x*b.y;
            eMat(2*i+2,2*j+2) += 2*a.y*b.y + a.x*b.x;
         }
      }
   }
}


void NSP2DQ41::RHS_Viscous(real_t coef)
{
   coef = 1;
   std::cerr << "Sorry, not yet implemented !\n";
}


void NSP2DQ41::Penal(real_t coef)
{
   _quad->atGauss(1,_sh,_dSh,_wg);
   for (size_t i=0; i<4; i++) {
      Point<real_t> a = coef*_wg[0]*_dSh[i];
      for (size_t j=0; j<4; j++) {
         eMat(2*i+1,2*j+1) += a.x*_dSh[j].x;
         eMat(2*i+1,2*j+2) += a.x*_dSh[j].y;
         eMat(2*i+2,2*j+1) += a.y*_dSh[j].x;
         eMat(2*i+2,2*j+2) += a.y*_dSh[j].y;
      }
   }
}


void NSP2DQ41::RHS_Convection(real_t coef)
{
   real_t ug[2];
   Point<real_t> dug[2];
   for (size_t k=0; k<4; ++k) {
      for (size_t j=0; j<2; j++) {
         ug[j] = dug[j].x = dug[j].y = 0.;
         for (size_t i=1; i<=4; i++) {
            size_t ii = 2*i+j-1;
            ug[j]    += _sh[4*(i-1)+k]*_eu(ii);
            dug[j].x += _dSh[4*(i-1)+k].x*_eu(ii);
            dug[j].y += _dSh[4*(i-1)+k].y*_eu(ii);
         }
      }
      real_t t1 = ug[0]*dug[0].x + ug[1]*dug[0].y;
      real_t t2 = ug[0]*dug[1].x + ug[1]*dug[1].y;
      for (size_t i=1; i<=4; ++i) {
         real_t c = coef*_dens*_ln->getLength()*_sh[4*(i-1)+k];
         eRHS(2*i-1) -= c*t1;
         eRHS(2*i  ) -= c*t2;
      }
   }
}


void NSP2DQ41::LHS1_Convection(real_t coef)
{
   real_t ug[2];
   for (size_t k=0; k<4; ++k) {
      for (size_t j=0; j<2; j++) {
         ug[j] = 0.;
         for (size_t i=1; i<=4; ++i)
            ug[j] += _sh[4*(i-1)+k]*_eu(2*i+j-1);
      }
      for (size_t ii=1; ii<=4; ++ii) {
         for (size_t jj=1; jj<=4; ++jj) {
            real_t c = coef*_dens*_wg[k]*(ug[0]*_dSh[4*(jj-1)+k].x + ug[1]*_dSh[4*(jj-1)+k].y);
            eMat(2*ii-1,2*jj-1) += c*_sh[4*(ii-1)+k];
            eMat(2*ii  ,2*jj  ) += c*_sh[4*(ii-1)+k];
         }
      }
   }
}


void NSP2DQ41::LHS2_Convection(real_t coef)
{
   Point<real_t>  dug[2];
   for (size_t k=0; k<4; ++k) {
      for (size_t j=0; j<2; j++) {
         dug[j].x = dug[j].y = 0.;
         for (size_t i=1; i<=4; i++) {
            dug[j].x += _dSh[4*(i-1)+k].x*_eu(2*i+j-1)*coef;
            dug[j].y += _dSh[4*(i-1)+k].y*_eu(2*i+j-1)*coef;
         }
      }
      for (size_t ii=1; ii<=4; ii++) {
         for (size_t jj=1; jj<=4; jj++) {
            real_t c = _dens*_quad->getDet()*_sh[4*(ii-1)+k]*_sh[4*(jj-1)+k];
            eMat(2*ii-1,2*jj-1) += c*dug[0].x;
            eMat(2*ii-1,2*jj  ) += c*dug[0].y;
            eMat(2*ii  ,2*jj-1) += c*dug[1].x;
            eMat(2*ii  ,2*jj  ) += c*dug[1].y;
         }
      }
   }
}


void NSP2DQ41::BodyRHS(Vect<real_t>& f)
{
   for (size_t k=0; k<4; ++k) {
      real_t fx=0., fy=0.;
      for (size_t j=0; j<4; ++j) {
         fx += _sh[4*j+k]*f((*_theElement)(j+1)->n(),1);
         fy += _sh[4*j+k]*f((*_theElement)(j+1)->n(),2);
      }
      for (size_t i=0; i<4; ++i) {
         eRHS(2*i+1) += _wg[k]*fx*_sh[4*i+k];
         eRHS(2*i+2) += _wg[k]*fy*_sh[4*i+k];
      }
   }
}


void NSP2DQ41::BoundaryRHS(Vect<real_t>& f)
{
}


void NSP2DQ41::Periodic(real_t coef)
{
   for (size_t i=1; i<=2; i++) {
      real_t c = 0.5*_ln->getLength()*coef;
      if (_theSide->getCode(1) == PERIODIC_A)
         sMat(2*i-1,2*i-1) += c;
      else if (_theSide->getCode(1) == PERIODIC_B)
         sMat(2*i-1,2*i-1) -= c;
      if (_theSide->getCode(2) == PERIODIC_A)
         sMat(2*i  ,2*i  ) += c;
      else if (_theSide->getCode(2) == PERIODIC_B)
         sMat(2*i  ,2*i  ) -= c;
   }
}


void NSP2DQ41::getPressure()
{
   if (_p==nullptr)
      return;
   Vect<real_t> ep(_nb_el);
   mesh_elements(*_theMesh) {
      set(the_element);
      real_t pr = 0.;
      for (size_t k=0; k<4; ++k) {
         for (size_t i=0; i<4; i++)
            pr += _quad->getDet() * (_dSh[4*i+k].x*(*_u)(The_element(i+1)->n(),1)
                                  +  _dSh[4*i+k].y*(*_u)(The_element(i+1)->n(),2));
      }
      ep(The_element.n()) = -pr*_penal;
   }
   _p->setSize(_nb_nodes);
   _p->setName("Pressure");
   _p->setTime(theTime);
   Reconstruction rr(*_theMesh);
   rr.P0toP1(ep,*_p);
}


void NSP2DQ41::build()
{
   if (_u==nullptr)
      throw OFELIException("In NSP2DQ41::build: No solution vector given.");
   _A->clear();
   _b->clear();
   mesh_elements(*_theMesh) {
      set(the_element);
      ElementNodeVector(*_u,_eu);
      if (_terms&MASS)
         LMass(1./theTimeStep);
      if (_terms&DILATATION)
         Penal(1.e04);
      if (_terms&VISCOSITY)
         Viscous();
      if (_terms&CONVECTION)
         RHS_Convection();
      if (_terms&BODY_FORCE && AbsEqua<real_t>::_bf!=nullptr)
         BodyRHS(*AbsEqua<real_t>::_bf);
      Equation<real_t,4,8,2,4>::updateBC(The_element,AbsEqua<real_t>::_bc);
      AbsEqua<real_t>::_A->Assembly(The_element,eMat.get());
      AbsEqua<real_t>::_b->Assembly(The_element,eRHS.get());
   }
   mesh_sides(*_theMesh) {
      set(the_side);
      ElementSideVector(*_u,_su);
      if (_terms&TRACTION && AbsEqua<real_t>::_sf!=nullptr)
         BoundaryRHS(*AbsEqua<real_t>::_sf);
      AbsEqua<real_t>::_b->Assembly(The_side,sRHS.get());
   }
}


int NSP2DQ41::runOneTimeStep()
{
   build();
   int ret = solveLinearSystem(AbsEqua<real_t>::_A,*_b,_uu);
   if (_bc == nullptr)
      _u->insertBC(*_theMesh,_uu);
   else
      _u->insertBC(*_theMesh,_uu,*_bc);
   getPressure();
   _u->setTime(theTime);
   return ret;
}

} /* namespace OFELI */
