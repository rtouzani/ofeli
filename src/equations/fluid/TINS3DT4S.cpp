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

                       Implementation of class TINS3DT4S

  ==============================================================================*/

#include "equations/fluid/TINS3DT4S.h"
#include "shape_functions/Tetra4.h"
#include "shape_functions/Triang3.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/LocalVect_impl.h"
#include "linear_algebra/SpMatrix_impl.h"
#include <algorithm>

using std::cout;

namespace OFELI {

TINS3DT4S::TINS3DT4S()
{
   _Re = 0.;
}


TINS3DT4S::TINS3DT4S(Mesh& ms)
          : Equation<real_t,4,12,3,9>(ms), _constant_matrix(true)
{
   init();
   _dens = 0.;
   _Re = 0.;
}


TINS3DT4S::TINS3DT4S(Mesh&         ms,
                     Vect<real_t>& u)
          : Equation<real_t,4,12,3,9>(ms,u), _constant_matrix(true)
{
   init();
   _dens = 0.;
   _Re = 0.;
}


TINS3DT4S::~TINS3DT4S() { }


void TINS3DT4S::init()
{
   if (Verbosity>2)
      cout << "Initializing Navier-Stokes equations settings ..." << endl;
   _theMesh->removeImposedDOF();
   _VM.setMesh(*_theMesh);
   _MM.setSize(_nb_nodes);
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
   _q.setSize(_nb_nodes);
   _c.setSize(_nb_nodes,2);
   _cfl = 0;
}


void TINS3DT4S::setInput(EqDataType    opt,
                         Vect<real_t>& u)
{
   AbsEqua<real_t>::setInput(opt,u);
   if (opt==PRESSURE_FIELD) {
      _p = &u;
      getPressure();
      updateVelocity();
   }
}


void TINS3DT4S::set(Element* el)
{
   _theElement = el, _theSide = nullptr;
   _visc = _dens = 1.;
   if (_Re==0.) {
      setMaterial();
      _Re = 1;
   }
   _ne = element_label;
   Tetra4 te(el);
   _el_geo.det = te.getDet();
   _cr = te.getMaxEdgeLength();
   _c24 = _el_geo.det/24.;
   _vol = _el_geo.det/6.;
   _dSh = te.DSh();
   for (size_t i=0; i<4; ++i)
      _en[i] = (*el)(i+1)->n();
   eMat = 0; eRHS = 0;
}


void TINS3DT4S::build()
{
   _TimeInt.step++;
   getMomentum();
   getPressure();
   updateVelocity();
}


void TINS3DT4S::ElementVelocityMatrix()
{
   real_t c=_vol*_visc/_Re;
   for (size_t j=1; j<=4; j++) {
      Point<real_t> db=c*_dSh[j-1];
      for (size_t i=1; i<=4; i++) {
         Point<real_t> a=_dSh[i-1];
         eMat(3*i-2,3*j-2) += 2*a.x*db.x + a.z*db.z + a.y*db.y;
         eMat(3*i-2,3*j-1) += a.y*db.x;
         eMat(3*i-2,3*j  ) += a.z*db.x;
         eMat(3*i-1,3*j-2) += a.x*db.y;
         eMat(3*i-1,3*j-1) += 2*a.y*db.y + a.z*db.z + a.x*db.x;
         eMat(3*i-1,3*j  ) += a.z*db.y;
         eMat(3*i  ,3*j-2) += a.x*db.z;
         eMat(3*i  ,3*j-1) += a.y*db.z;
         eMat(3*i  ,3*j  ) += 2*a.z*db.z + a.y*db.y + a.x*db.x;
      }
   }
   c = _c24*_dens;
   for (size_t i=1; i<=4; i++) {
      eMat(3*i-2,3*i-2) += c;
      eMat(3*i-1,3*i-1) += c;
      eMat(3*i  ,3*i  ) += c;
   }
}


void TINS3DT4S::PressureMatrix()
{
   if (Verbosity>2)
      cout << "Calculating pressure matrix ..." << endl;
   SpMatrix<real_t> Dx(1,*_theMesh,0), Dy(1,*_theMesh,0);
   MESH_EL {
      set(the_element);
      real_t z=_c24/_TimeInt.delta, a=_vol*_TimeInt.delta;
      for (size_t i=0; i<4; i++) {
         _MM(_en[i]) += z;
         for (size_t j=0; j<4; j++) {
            _PM.add(_en[i],_en[j],a*(_dSh[i]*_dSh[j]));
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


void TINS3DT4S::getMomentum()
{
   if (Verbosity>2)
      cout << "Solving momentum equations ..." << endl;
   LocalVect<real_t,12> ce;
   if (_TimeInt.step==1 || _constant_matrix==false)
      _VM = 0;
   _b = 0;
   size_t j=0;
   _c = 0;
   MESH_ND {
      if (The_node.getCode(1)==0)
         _b[j++] = 0.5*_c(node_label,1);
      if (The_node.getCode(2)==0)
         _b[j++] = 0.5*_c(node_label,2);
   }

// Loop over elements
   MESH_EL {
      set(the_element);
      LocalVect<real_t,12> ue(the_element,*_u);
      LocalVect<real_t,4> pe(the_element,*_p,1);

//    cfl in element
      real_t d1=0.25*(ue[0]+ue[3]+ue[6]+ue[ 9]),
             d2=0.25*(ue[1]+ue[4]+ue[7]+ue[10]),
             d3=0.25*(ue[2]+ue[5]+ue[8]+ue[11]);
      _cfl = std::max(_cfl,_TimeInt.delta*sqrt(d1*d1+d2*d2+d3*d3)/_cr);

//    Element matrix
      if (_TimeInt.step==1 || _constant_matrix==false)
         ElementVelocityMatrix();

//    Mass (R.H.S.)
      real_t cc=_c24*_dens/_TimeInt.delta;
      for (size_t i=0; i<3; ++i) {
         eRHS(2*i+1) = cc*ue[2*i  ];
         eRHS(2*i+2) = cc*ue[2*i+1];
      }

//    Viscous Term (R.H.S.)
      cc = 0.5*_vol*_visc/_Re;
      for (size_t i=0; i<4; i++) {
         for (size_t j=0; j<4; j++) {
            eRHS(3*i+1) -= cc*(_dSh[i].y*_dSh[j].x*ue[2*j+1] + (2*_dSh[i].x*_dSh[j].x+_dSh[i].y*_dSh[j].y)*ue[2*j  ]);
            eRHS(3*i+2) -= cc*(_dSh[i].x*_dSh[j].y*ue[2*j  ] + (2*_dSh[i].y*_dSh[j].y+_dSh[i].x*_dSh[j].x)*ue[2*j+1]);
            eRHS(3*i+3) -= cc*(_dSh[i].x*_dSh[j].y*ue[2*j  ] + (2*_dSh[i].y*_dSh[j].y+_dSh[i].x*_dSh[j].x)*ue[2*j+1]);
         }
      }

//    Convection
      Point<real_t> du = _dSh[0]*ue[0] + _dSh[1]*ue[3] + _dSh[2]*ue[6] + _dSh[3]*ue[ 9],
                    dv = _dSh[0]*ue[1] + _dSh[1]*ue[4] + _dSh[2]*ue[7] + _dSh[2]*ue[10],
                    dw = _dSh[0]*ue[2] + _dSh[1]*ue[5] + _dSh[2]*ue[8] + _dSh[2]*ue[11];
      for (size_t i=0; i<4; i++) {
         ce[3*i  ] = _c24*_dens*((d1 + ue[3*i])*du.x + (d2 + ue[3*i+1])*du.y + (d3 + ue[3*i+2])*du.z);
         ce[3*i+1] = _c24*_dens*((d1 + ue[3*i])*dv.x + (d2 + ue[3*i+1])*dv.y + (d3 + ue[3*i+2])*dv.z);
         ce[3*i+2] = _c24*_dens*((d1 + ue[3*i])*dw.x + (d2 + ue[3*i+1])*dw.y + (d3 + ue[3*i+2])*dw.z);
      }
      Axpy(-1.5,ce,eRHS);
      Assembly(The_element,ce,_c);

//    Pressure Gradient
      Point<real_t> dp = _c24*(pe[0]*_dSh[0]+pe[1]*_dSh[1]+pe[2]*_dSh[2]+pe[3]*_dSh[3]);
      for (size_t i=0; i<4; ++i) {
         eRHS(3*i+1) -= dp.x;
         eRHS(3*i+2) -= dp.y;
         eRHS(3*i+3) -= dp.z;
      }

//    Body Force
      if (_bf!=nullptr) {
         for (size_t i=0; i<4; ++i) {
            eRHS(3*i+1) += (*_bf)(3*_en[i]+1)*_c24;
            eRHS(3*i+2) += (*_bf)(3*_en[i]+2)*_c24;
            eRHS(3*i+3) += (*_bf)(3*_en[i]+3)*_c24;
         }
      }

//    Boundary conditions
      Equa_Fluid<real_t,4,12,3,9>::updateBC(The_element,*_bc);

//    Assembly of matrix and R.H.S.
      if (_TimeInt.step==1 || _constant_matrix==false)
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
   if (_bc!=nullptr)
      _u->insertBC(_uu,*_bc);
   if (Verbosity>1)
      cout << "Nb. of CG iterations for Momentum: " << nb_it << endl;
}


int TINS3DT4S::getPressure()
{
   Vect<real_t> b(_theMesh->getNbNodes());
   if (Verbosity>2)
      cout << "Solving pressure equation ..." << endl;
   MESH_EL {
      set(the_element);
      LocalVect<real_t,6> ue(the_element,*_u);
      double d = _c24*(_dSh[0].x*ue[0] + _dSh[0].y*ue[1] +
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
   int nb_it = CG(_PM,_PP,b,_q,1000,toler);
#endif
   if (Verbosity>1)
      cout << "Nb. of CG iterations for pressure: " << nb_it << endl;
   *_p += 2.*_q;
   return nb_it;
}


void TINS3DT4S::updateVelocity()
{
   if (Verbosity>2)
      cout << "Updating velocity ..." << endl;
   MESH_EL {
      set(the_element);
      LocalVect<real_t,3> qe(the_element,_q,1);
      Point<real_t> dp = _c24*(qe[0]*_dSh[0]+qe[1]*_dSh[1]+qe[2]*_dSh[2]+qe[3]*_dSh[3]);
      for (size_t i=0; i<4; ++i) {
         size_t n = _en[i];
         if ((*_theMesh)[n]->getCode(1)==0)
            (*_u)(n,1) -= dp.x/_MM(n);
         if ((*_theMesh)[n]->getCode(2)==0)
            (*_u)(n,2) -= dp.y/_MM(n);
         if ((*_theMesh)[n]->getCode(3)==0)
            (*_u)(n,3) -= dp.z/_MM(n);
      }
   }
}

} /* namespace OFELI */
