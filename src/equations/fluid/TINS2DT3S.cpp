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

                       Implementation of class TINS2DT3S

  ==============================================================================*/

#include "equations/fluid/TINS2DT3S.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"
#include <algorithm>

using std::cout;

namespace OFELI {

TINS2DT3S::TINS2DT3S() : _Re(0.)
{
}


TINS2DT3S::TINS2DT3S(Mesh&         mesh,
                     Vect<real_t>& u,
                     Vect<real_t>& p,
                     real_t&       ts,
                     real_t        Re)
          :  _constant_matrix(true), _Re(Re), _p(&p)
{
   _step = 0;
   _u = &u;
   init(mesh,ts);
   _bc_given = _bf_given = _sf_given = false;
   _dens = 0.;
}


TINS2DT3S::~TINS2DT3S() { }


void TINS2DT3S::init(Mesh&  mesh,
                     real_t ts)
{
   if (_verbose>2)
      cout << "Initializing Navier-Stokes equations settings ..." << endl;
   _theMesh = &mesh;
   _time_step = ts;
   _theMesh->removeImposedDOF();
   _VM.setMesh(*_theMesh);
   _MM.setSize(_theMesh->getNbNodes());
   _PM.setExtendedGraph();
   _PM.setMesh(*_theMesh,1);
   for (size_t i=0; i<=_PM.size(); i++)
      _row_ptr.push_back(_PM.getRowPtr(i+1));
   for (size_t i=0; i<_PM.getLength(); i++)
      _col_ind.push_back(_PM.getColInd(i+1));
   PressureMatrix();
#if !defined(USE_EIGEN)
   _PP.setType(DILU_PREC);
   _PP.setMatrix(_PM);
#endif
   _uu.setSize(_theMesh->getNbEq());
   _b.setSize(_theMesh->getNbEq());
   _q.setSize(_theMesh->getNbNodes());
   _c.setSize(_theMesh->getNbNodes(),2);
   _cfl = 0;
}


void TINS2DT3S::setInput(EqDataType    opt,
                         Vect<real_t>& u)
{
   AbsEqua<real_t>::setInput(opt,u);
   if (opt==INITIAL_FIELD) {
      _u = &u;
      getPressure();
      updateVelocity();
   }
}


void TINS2DT3S::set(Element* el)
{
   Init(el);
   _visc = _dens = 1.;
   if (_Re==0.) {
      setMaterial();
      _Re = 1;
   }
   _ne = element_label;
   Triang3 tr(el);
   _center = tr.getCenter();
   _det = tr.getDet();
   _cr = tr.getCircumRadius();
   for (size_t i=0; i<3; ++i) {
      _dSh[i] = tr.DSh(i+1);
      _en[i] = (*el)(i+1)->n();
   }
   eMat = 0; eRHS = 0;
}


void TINS2DT3S::set(Side* sd)
{
   Init(sd);
   Line2 ln(sd);
   SideNodeCoordinates();
   _center = ln.getCenter();
   _length = ln.getLength();
   sMat = 0; sRHS = 0;
}


void TINS2DT3S::build()
{
   _step++;
   getMomentum();
   getPressure();
   updateVelocity();
}


int TINS2DT3S::runOneTimeStep()
{
   build();
   return 0;
}


void TINS2DT3S::ElementVelocityMatrix()
{
   real_t a=0.25*_det*_visc/_Re;
   for (size_t i=1; i<=3; i++) {
      for (size_t j=1; j<=3; j++) {
         eMat(2*i-1,2*j-1) = a*(2*_dSh(i).x*_dSh(j).x + _dSh(i).y*_dSh(j).y);
         eMat(2*i-1,2*j  ) = a*_dSh(i).y*_dSh(j).x;
         eMat(2*i  ,2*j-1) = a*_dSh(i).x*_dSh(j).y;
         eMat(2*i  ,2*j  ) = a*(2*_dSh(i).y*_dSh(j).y + _dSh(i).x*_dSh(j).x);
      }
   }
   a = _dens*_det/(6.*_time_step);
   for (size_t i=1; i<=3; i++) {
      eMat(2*i-1,2*i-1) += a;
      eMat(2*i  ,2*i  ) += a;
   }
}


void TINS2DT3S::SideVelocityMatrix()
{
   const real_t penal=1.e20;
   for (size_t i=1; i<=2; i++) {
      real_t c = 0.5*_length*penal;
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
   if (_verbose>2)
      cout << "Calculating pressure matrix ..." << endl;
   SpMatrix<real_t> Dx(1,*_theMesh,0), Dy(1,*_theMesh,0);
   mesh_elements(*_theMesh) {
      set(the_element);
      real_t z=_det/(6*_time_step), a=0.5*_time_step*_det;
      for (size_t i=0; i<3; i++) {
         _MM(_en[i]) += z;
         for (size_t j=0; j<3; j++) {
            _PM.add(_en[i],_en[j],a*(_dSh[i]*_dSh[j]));
            Dx.add(_en[i],_en[j],_dSh[j].x*_det);
            Dy.add(_en[i],_en[j],_dSh[j].y*_det);
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
   if (_verbose>2)
      cout << "Solving momentum equations ..." << endl;
   LocalVect<real_t,6> ce;
   if (_step==1 || _constant_matrix==false)
      _VM = 0;
   _b = 0;
   size_t j=0;
   _c = 0;
   mesh_nodes(*_theMesh) {
      if (The_node.getCode(1)==0)
         _b[j++] = 0.5*_c(node_label,1);
      if (The_node.getCode(2)==0)
         _b[j++] = 0.5*_c(node_label,2);
   }

// Loop over elements
   mesh_elements(*_theMesh) {
      set(the_element);
      LocalVect<real_t,6> ue(the_element,*_u);
      LocalVect<real_t,3> pe(the_element,*_p,1);

//    cfl in element
      real_t d1=(ue[0]+ue[2]+ue[4])/3., d2=(ue[1]+ue[3]+ue[5])/3.;
      _cfl = std::max(_cfl,_time_step*sqrt(d1*d1+d2*d2)/_cr);

//    Element matrix
      if (_step==1 || _constant_matrix==false)
         ElementVelocityMatrix();

//    Mass (R.H.S.)
      real_t cc=_dens*_det/(6.*_time_step);
      for (size_t i=0; i<3; ++i) {
         eRHS(2*i+1) = cc*ue[2*i  ];
         eRHS(2*i+2) = cc*ue[2*i+1];
      }

//    Viscous Term (R.H.S.)
      cc = 0.25*_visc*_det/_Re;
      for (size_t i=0; i<3; i++) {
         for (size_t j=0; j<3; j++) {
            eRHS(2*i+1) -= cc*(_dSh[i].y*_dSh[j].x*ue[2*j+1] + (2*_dSh[i].x*_dSh[j].x+_dSh[i].y*_dSh[j].y)*ue[2*j  ]);
            eRHS(2*i+2) -= cc*(_dSh[i].x*_dSh[j].y*ue[2*j  ] + (2*_dSh[i].y*_dSh[j].y+_dSh[i].x*_dSh[j].x)*ue[2*j+1]);
         }
      }

//    Convection
      cc = _det*_dens/24.;
      Point<real_t> du = _dSh[0]*ue[0] + _dSh[1]*ue[2] + _dSh[2]*ue[4],
                    dv = _dSh[0]*ue[1] + _dSh[1]*ue[3] + _dSh[2]*ue[5];
      for (size_t i=0; i<3; i++) {
         ce[2*i  ] = cc*((d1 + ue[2*i])*du.x + (d2 + ue[2*i+1])*du.y);
         ce[2*i+1] = cc*((d1 + ue[2*i])*dv.x + (d2 + ue[2*i+1])*dv.y);
      }
      Axpy(-1.5,ce,eRHS);
      Assembly(The_element,ce,_c);

//    Pressure Gradient
      Point<real_t> dp = _det/6.*(pe[0]*_dSh[0]+pe[1]*_dSh[1]+pe[2]*_dSh[2]);
      for (size_t i=0; i<3; ++i) {
         eRHS(2*i+1) -= dp.x;
         eRHS(2*i+2) -= dp.y;
      }

//    Body Force
      if (_bf_given) {
         for (size_t i=0; i<3; ++i) {
            eRHS(2*i+1) += (*_bf)(2*_en[i]-1)*_det/6.;
            eRHS(2*i+2) += (*_bf)(2*_en[i]  )*_det/6.;
         }
      }
      
//    Boundary conditions
      Equa_Fluid<real_t,3,6,2,4>::updateBC(The_element,*_bc);

//    Assembly of matrix and R.H.S.
      if (_step==1 || _constant_matrix==false)
         Assembly(The_element,eMat,_VM);
      Assembly(The_element,eRHS,_b);
   }

// Solve the linear system
#ifdef USE_EIGEN
   LinearSolver<real_t> ls(1000,toler,1);
   int nb_it = ls.solve(_VM,_b,_uu,CG_SOLVER,ILU_PREC);
#else
   LinearSolver<real_t> ls(_VM,_b,_uu);
   ls.setTolerance(_toler);
   int nb_it = ls.solve(CG_SOLVER,DILU_PREC);
#endif
   if (_bc_given)
      _u->insertBC(_uu,*_bc);
   if (_verbose>1)
      cout << "Nb. of CG iterations for Momentum: " << nb_it << endl;
}


int TINS2DT3S::getPressure()
{
   Vect<real_t> b(_theMesh->getNbNodes());
   if (_verbose>2)
      cout << "Solving pressure equation ..." << endl;
   mesh_elements(*_theMesh) {
      set(the_element);
      LocalVect<real_t,6> ue(the_element,*_u);
      double d = _det/6.*(_dSh[0].x*ue[0] + _dSh[0].y*ue[1] +
                          _dSh[1].x*ue[2] + _dSh[1].y*ue[3] +
                          _dSh[2].x*ue[4] + _dSh[2].y*ue[5]);
      for (size_t i=0; i<3; ++i)
         b(_en[i]) -= d;
   }
   real_t toler=1.e-7;
#ifdef USE_EIGEN
   LinearSolver<double> ls(1000,toler,1);
   int nb_it = ls.solve(_PM,b,_q,CG_SOLVER,ILU_PREC);
#else
   int nb_it = CG(_PM,_PP,b,_q,1000,toler,0);
#endif
   if (_verbose>1)
      cout << "Nb. of CG iterations for pressure: " << nb_it << endl;
   *_p += 2.*_q;
   return nb_it;
}


void TINS2DT3S::updateVelocity()
{
   if (_verbose>2)
      cout << "Updating velocity ..." << endl;
   mesh_elements(*_theMesh) {
      set(the_element);
      LocalVect<real_t,3> qe(the_element,_q,1);
      Point<real_t> dp = _det/6.*(qe[0]*_dSh[0]+qe[1]*_dSh[1]+qe[2]*_dSh[2]);
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
