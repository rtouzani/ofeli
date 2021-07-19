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

                       Implementation of class TINS2DT3S

  ==============================================================================*/

#include "equations/Equation_impl.h"
#include "equations/fluid/TINS2DT3S.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/LocalVect_impl.h"
#include "linear_algebra/Assembly_impl.h"
#include <algorithm>

using std::cout;

namespace OFELI {

TINS2DT3S::TINS2DT3S()
{
   _Re = 0.;
}


TINS2DT3S::TINS2DT3S(Mesh& ms)
          : Equation<3,6,2,4>(ms), _constant_matrix(false), _p(nullptr)
{
   _TimeInt.step = 0;
   init();
   _Re = 0.;
}


TINS2DT3S::TINS2DT3S(Mesh&         ms,
                     Vect<real_t>& u)
          : Equation<3,6,2,4>(ms,u), _constant_matrix(false), _p(nullptr)
{
   _TimeInt.step = 0;
   init();
   _Re = 0.;
}


TINS2DT3S::~TINS2DT3S() { }


void TINS2DT3S::init()
{
   if (Verbosity>2)
      cout << "Initializing Navier-Stokes equations settings ..." << endl;
   _MM.setSize(_nb_nodes);
   _PM.setExtendedGraph();
   _PM.setMesh(*_theMesh,1);
   for (size_t i=0; i<=_PM.size(); i++)
      _row_ptr.push_back(_PM.getRowPtr(i+1));
   for (size_t i=0; i<_PM.getLength(); i++)
      _col_ind.push_back(_PM.getColInd(i+1));
   _bp.setSize(_nb_nodes);
   PressureMatrix();
#if !defined(USE_EIGEN)
   _PP.setType(DILU_PREC);
   //   _PP.setMatrix(_PM);
#endif
   _q.setSize(_nb_nodes);
   _c.setSize(_nb_nodes,2);
   _cfl = 0;
}


void TINS2DT3S::setInput(EqDataType    opt,
                         Vect<real_t>& u)
{
   Equa::setInput(opt,u);
   if (opt==PRESSURE_FIELD) {
      _p = &u;
      getPressure();
      updateVelocity();
   }
}


void TINS2DT3S::set(Element* el)
{
   _theElement = el, _theSide = nullptr;
   _mu = _rho = 1.;
   if (_Re==0.) {
      setMaterial();
      _Re = 1.0;
   }
   _ne = element_label;
   Triang3 tr(el);
   _el_geo.det = tr.getDet();
   _el_geo.center = tr.getCenter();
   _cr = tr.getCircumRadius();
   _dSh = tr.DSh();
   _cd = OFELI_SIXTH*_el_geo.det;
   for (size_t i=0; i<3; ++i)
      _en[i] = (*el)(i+1)->n();
   if (_rho_set)
      _rho = _rho_fct(_el_geo.center,_TimeInt.time);
   if (_mu_set)
      _mu = _mu_fct(_el_geo.center,_TimeInt.time);
   if (_beta_set)
      _beta = _beta_fct(_el_geo.center,_TimeInt.time);
   eMat = 0; eRHS = 0;
}


void TINS2DT3S::set(Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Line2 ln(sd);
   SideNodeCoordinates();
   SideNodeVector(*_u,_su);
   _el_geo.length = ln.getLength();
   sMat = 0; sRHS = 0;
}


void TINS2DT3S::build()
{
   _TimeInt.step++;
   getMomentum();
   getPressure();
   updateVelocity();
}


int TINS2DT3S::runOneTimeStep()
{
   if (_TimeInt.step==0) {
      if (_u==nullptr)
         throw OFELIException("In TINS2DT3S::runOneTimeStep(): No solution vector (initial condition) given.");
      if (_p==nullptr)
         throw OFELIException("In TINS2DT3S::runOneTimeStep(): No pressure vector given.");
      PressureMatrix();
   }
   build();
   return 0;
}


void TINS2DT3S::ElementVelocityMatrix()
{
   real_t a = _rho*_cd/_TimeInt.delta;
   real_t c = 0.25*_el_geo.det*_mu/_Re;
   for (size_t i=0; i<3; ++i) {
      for (size_t j=0; j<3; ++j) {
         eMat(2*i+1,2*j+1) += c*(2*_dSh[i].x*_dSh[j].x + _dSh[i].y*_dSh[j].y);
         eMat(2*i+1,2*j+2) += c*_dSh[i].y*_dSh[j].x;
         eMat(2*i+2,2*j+1) += c*_dSh[i].x*_dSh[j].y;
         eMat(2*i+2,2*j+2) += c*(2*_dSh[i].y*_dSh[j].y + _dSh[i].x*_dSh[j].x);
      }
      eMat(2*i+1,2*i+1) += a;
      eMat(2*i+2,2*i+2) += a;
   }
}


void TINS2DT3S::SideVelocityMatrix()
{
   const real_t penal=1.e20;
   for (size_t i=1; i<=2; i++) {
      real_t c = 0.5*_el_geo.length*penal;
      if (The_side.getCode(1)==PERIODIC_A)
         sMat(2*i-1,2*i-1) += c;
      else if (The_side.getCode(1) == PERIODIC_B)
         eMat(2*i-1,2*i-1) -= c;
      if (The_side.getCode(2)==PERIODIC_A)
         sMat(2*i  ,2*i  ) += c;
      else if (The_side.getCode(2)==PERIODIC_B)
         sMat(2*i  ,2*i  ) -= c;
   }
}


void TINS2DT3S::PressureMatrix()
{
   if (Verbosity>2)
      cout << "Calculating pressure matrix ..." << endl;
   SpMatrix<real_t> Dx(1,*_theMesh,0), Dy(1,*_theMesh,0);
   MESH_EL {
      set(the_element);
      real_t a=0.5*_TimeInt.delta*_el_geo.det;
      for (size_t i=0; i<3; i++) {
         _MM(_en[i]) += _cd/_TimeInt.delta;
         for (size_t j=0; j<3; j++) {
            _PM.add(_en[i],_en[j],a*(_dSh[i],_dSh[j]));
            Dx.add(_en[i],_en[j],_dSh[j].x*_el_geo.det);
            Dy.add(_en[i],_en[j],_dSh[j].y*_el_geo.det);
         }
      }
   }
   for (size_t i=1; i<=_PM.size(); i++) {
      for (size_t jj=_row_ptr[i-1]; jj<_row_ptr[i]; jj++) {
         size_t j=_col_ind[jj];
         for (size_t kk=_row_ptr[i-1]; kk<_row_ptr[i]; kk++) {
            size_t k=_col_ind[kk];
            real_t d=1./(36*_MM(k));
            for (size_t ll=_row_ptr[j-1]; ll<_row_ptr[j]; ll++) {
               if (k==_col_ind[ll]) {
                  _PM.add(i,j,d*(Dx(i,k)*Dx(j,k)+Dy(i,k)*Dy(j,k)));
                  break;
               }
            }
         }
      }
   }
}


void TINS2DT3S::getMomentum()
{
   if (Verbosity>2)
      cout << "Solving momentum equations ..." << endl;
   LocalVect<real_t,6> ce;
   _A->clear();
   _b->clear();
   size_t j=0;
   MESH_ND {
      if (The_node.getCode(1)==0)
         (*_b)[j++] = 0.5*_c(node_label,1);
      if (The_node.getCode(2)==0)
         (*_b)[j++] = 0.5*_c(node_label,2);
   }

// Loop over elements
   MESH_EL {
      set(the_element);
      ElementNodeVector(*_u,_eu);

//    cfl in element
      real_t d1=(_eu(1)+_eu(3)+_eu(5))/3., d2=(_eu(2)+_eu(4)+_eu(6))/3.;
      _cfl = std::max(_cfl,_TimeInt.delta*sqrt(d1*d1+d2*d2)/_cr);

//    Element matrix
      ElementVelocityMatrix();

//    Mass (R.H.S.)
      real_t a = _rho*_cd/_TimeInt.delta;
      for (size_t i=1; i<=3; ++i) {
         eRHS(2*i-1) = a*_eu(2*i-1);
         eRHS(2*i  ) = a*_eu(2*i  );
      }

//    Viscous Term (R.H.S.)
      a = 0.25*_mu*_el_geo.det/_Re;
      for (size_t i=0; i<3; i++) {
         for (size_t j=0; j<3; j++) {
            eRHS(2*i+1) -= a*(_dSh[i].y*_dSh[j].x*_eu[2*j+1] + (2*_dSh[i].x*_dSh[j].x+_dSh[i].y*_dSh[j].y)*_eu[2*j  ]);
            eRHS(2*i+2) -= a*(_dSh[i].x*_dSh[j].y*_eu[2*j  ] + (2*_dSh[i].y*_dSh[j].y+_dSh[i].x*_dSh[j].x)*_eu[2*j+1]);
         }
      }

//    Convection
      a = 0.25*_cd*_rho;
      Point<real_t> du = _dSh[0]*_eu[0] + _dSh[1]*_eu[2] + _dSh[2]*_eu[4],
                    dv = _dSh[0]*_eu[1] + _dSh[1]*_eu[3] + _dSh[2]*_eu[5];
      for (size_t i=0; i<3; i++) {
         ce[2*i  ] = a*((d1 + _eu[2*i])*du.x + (d2 + _eu[2*i+1])*du.y);
         ce[2*i+1] = a*((d1 + _eu[2*i])*dv.x + (d2 + _eu[2*i+1])*dv.y);
      }
      Axpy(-1.5,ce,eRHS);
      Assembly(The_element,ce,_c);

//    Pressure Gradient
      Point<real_t> dp = _cd*((*_p)(_en[0])*_dSh[0]+(*_p)(_en[1])*_dSh[1]+(*_p)(_en[2])*_dSh[2]);
      for (size_t i=1; i<=3; ++i) {
         eRHS(2*i-1) -= dp.x;
         eRHS(2*i  ) -= dp.y;
      }

//    Body Force
      if (_bf!=nullptr) {
         for (size_t i=0; i<3; ++i) {
            eRHS(2*i+1) += _cd*(*_bf)(2*_en[i]-1);
            eRHS(2*i+2) += _cd*(*_bf)(2*_en[i]  );
         }
      }

//    Boundary conditions
      Equa_Fluid<3,6,2,4>::updateBC(The_element,_bc);

//    Assembly of matrix and R.H.S.
      Equa::_A->Assembly(The_element,eMat.get());
      Equa::_b->Assembly(The_element,eRHS.get());
   }

// Solve the linear system
#ifdef USE_EIGEN
   LinearSolver<real_t> ls(1000,toler,1);
   int nb_it = ls.solve(*_A,*_b,_uu,CG_SOLVER,ILU_PREC);
#else
   int nb_it = solveLinearSystem(_A,*_b,_uu);
#endif
   if (_bc!=nullptr)
      _u->insertBC(_uu,*_bc);
   if (Verbosity>3)
      cout << "Nb. of CG iterations for Momentum: " << nb_it << endl;
}


int TINS2DT3S::getPressure()
{
   _bp.clear();
   if (Verbosity>2)
      cout << "Solving pressure equation ..." << endl;
   MESH_EL {
      set(the_element);
      ElementNodeVector(*_u,_eu);
      real_t d = _cd*(_dSh[0].x*_eu[0] + _dSh[0].y*_eu[1] +
                      _dSh[1].x*_eu[2] + _dSh[1].y*_eu[3] +
                      _dSh[2].x*_eu[4] + _dSh[2].y*_eu[5]);
      for (size_t i=0; i<3; ++i)
         _bp(_en[i]) -= d;
   }
   real_t toler=1.e-7;
#ifdef USE_EIGEN
   LinearSolver<real_t> ls(1000,toler,1);
   int nb_it = ls.solve(_PM,_bp,_q,CG_SOLVER,ILU_PREC);
#else
   _PP.setMatrix(_PM);
   int nb_it = CG(_PM,_PP,_bp,_q,1000,toler);
#endif
   if (Verbosity>3)
      cout << "Nb. of CG iterations for pressure: " << nb_it << endl;
   *_p += 2.*_q;
   return nb_it;
}


void TINS2DT3S::updateVelocity()
{
   if (Verbosity>2)
      cout << "Updating velocity ..." << endl;
   MESH_EL {
      set(the_element);
      ElementNodeVector(*_u,_eu);
      LocalVect<real_t,3> qe(the_element,_q,1);
      Point<real_t> dp = _cd*(qe[0]*_dSh[0]+qe[1]*_dSh[1]+qe[2]*_dSh[2]);
      for (size_t i=0; i<3; ++i) {
         size_t n = _en[i];
         if ((*_theMesh)[n]->getCode(1)==0)
            (*_u)(n,1) -= dp.x/_MM(n);
         if ((*_theMesh)[n]->getCode(2)==0)
            (*_u)(n,2) -= dp.y/_MM(n);
      }
   }
}

} /* namespace OFELI */
