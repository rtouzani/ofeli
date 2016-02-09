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

                       Implementation of class TINS2DT3B

  ==============================================================================*/


#include "equations/fluid/TINS2DT3B.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"
#include <algorithm>

using std::cout;

namespace OFELI {

TINS2DT3B::TINS2DT3B()
{
   _Re = 0;
   _b = NULL;
}


TINS2DT3B::TINS2DT3B(Mesh&         mesh,
                     Vect<real_t>& u,
                     Vect<real_t>& p,
                     real_t&       ts,
                     real_t        Re)
{
   _step = 1;
   _b = NULL;
   _u = &u;
   _p = &p;
   _Re = Re;
   initEquation(mesh,ts);
   _bc_given = _bf_given = _sf_given = false;
}


TINS2DT3B::~TINS2DT3B()
{
   if (_b)
      delete _b;
}


void TINS2DT3B::initEquation(Mesh&  mesh,
                             real_t ts)
{
   if (_verbose>2)
      cout << "Initializing Navier-Stokes equations settings ..." << endl;
   _theMesh = &mesh;
   _time_step = ts;
   _theMesh->NumberEquations();
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
   _VM.setMesh(*_theMesh);
   VelocityMatrix();
#if !defined(USE_EIGEN)
   _PV.setType(DILU_PREC);
   _PV.setMatrix(_VM);
#endif
   _uu.setSize(_theMesh->getNbEq());
   _b = new Vect<real_t>(_theMesh->getNbEq());
   _ub.setMesh(*_theMesh,2,ELEMENT_DOF);
   _q.setMesh(*_theMesh,1,NODE_DOF);
   _c.setMesh(*_theMesh,2,NODE_DOF);
   _cfl = 0;
}


void TINS2DT3B::setInput(EqDataType    opt,
                         Vect<real_t>& u)
{
   AbsEqua<real_t>::setInput(opt,u);
   if (opt==INITIAL_FIELD) {
      _u = &u;
      mesh_elements(*_theMesh) {
         size_t n1=The_element(1)->n(),
                n2=The_element(2)->n(),
                n3=The_element(3)->n();
         _ub(element_label,1) = OFELI_THIRD*((*_u)(n1,1)+(*_u)(n2,1)+(*_u)(n3,1));
         _ub(element_label,2) = OFELI_THIRD*((*_u)(n1,2)+(*_u)(n2,2)+(*_u)(n3,2));
      }
      getPressure();
      updateVelocity();
   }
}


void TINS2DT3B::set(Element* el)
{
   Init(el);
   _visc = _dens = 1.;
   if (_Re==0.) {
      setMaterial();
      _Re = 1;
   }
   Triang3 tr(el);
   _center = tr.getCenter();
   _det = tr.getDet();
   ElementNodeCoordinates();
   _dSh(1) = tr.DSh(1);
   _dSh(2) = tr.DSh(2);
   _dSh(3) = tr.DSh(3);
   eMat = 0; eRHS = 0;
}


void TINS2DT3B::set(Side* sd)
{
   Init(sd);
   Line2 ln(sd);
   SideNodeCoordinates();
   _center = ln.getCenter();
   _length = ln.getLength();
   sMat = 0; sRHS = 0;
}


void TINS2DT3B::build()
{
   _step++;
   getMomentum();
   getPressure();
   updateVelocity();
}


int TINS2DT3B::runOneTimeStep()
{
   build();
   return 0;
}


void TINS2DT3B::VelocityMatrix()
{
   const real_t penal=1.e20;
   if (_verbose>2)
      cout << "Calculating velocity matrices ..." << endl;
   Point<real_t> z;
   MESH_EL {
      set(theElement);
      real_t c=OFELI_SIXTH*_dens*_det/_time_step, aa=0.25*_det*_visc/_Re;
      for (size_t i=0; i<3; i++) {
         z = aa*_dSh[i];
         for (size_t j=i; j<3; j++) {
            eMat(2*i+1,2*j+1) = 2*z.x*_dSh[j].x + z.y*_dSh[j].y;
            eMat(2*i+1,2*j+2) = z.y*_dSh[j].x;
            eMat(2*i+2,2*j+1) = z.x*_dSh[j].y;
            eMat(2*i+2,2*j+2) = 2*z.y*_dSh[j].y + z.x*_dSh[j].x;
         }
      }
      eMat(1,1) += c; eMat(2,2) += c; eMat(3,3) += c;
      eMat(4,4) += c; eMat(5,5) += c; eMat(6,6) += c;
      eMat.Symmetrize();
      ElementAssembly(TheElement,_VM);
   }

   MESH_SD {
      set(theSide);
      for (size_t i=1; i<=2; i++) {
         real_t c = 0.5*_length*penal;
         if (TheSide.getCode(1)==PERIODIC_A)
            sMat(2*i-1,2*i-1) += c;
         else if (TheSide.getCode(1) == PERIODIC_B)
            eMat(2*i-1,2*i-1) -= c;
         if (TheSide.getCode(2)==PERIODIC_A)
            sMat(2*i  ,2*i  ) += c;
         else if (TheSide.getCode(2)==PERIODIC_B)
            sMat(2*i  ,2*i  ) -= c;
      }
      SideAssembly(TheSide,_VM);
   }
}


void TINS2DT3B::PressureMatrix()
{
   if (_verbose>2)
      cout << "Calculating pressure matrix ..." << endl;
   size_t size=_PM.size();
   SpMatrix<real_t> Dx(1,*_theMesh), Dy(1,*_theMesh);
   mesh_elements(*_theMesh) {
      set(the_element);
      real_t z=OFELI_SIXTH*_det/_time_step, a=0.35*_time_step*_det;
      for (size_t i=0; i<3; i++) {
         size_t ii=The_element(i+1)->n();
         _MM(ii) += z;
         for (size_t j=0; j<3; j++) {
            size_t jj=The_element(j+1)->n();
            _PM.add(ii,jj,a*(_dSh[i]*_dSh[j]));
            Dx.add(ii,jj,_dSh[j].x*_det);
            Dy.add(ii,jj,_dSh[j].y*_det);
         }
      }
   }
   for (size_t i=1; i<=size; i++) {
      for (size_t jj=_row_ptr[i-1]; jj<_row_ptr[i]; jj++) {
         size_t j=_col_ind[jj];
         for (size_t kk=_row_ptr[i-1]; kk<_row_ptr[i]; kk++) {
            size_t k=_col_ind[kk];
            real_t d=OFELI_SIXTH*OFELI_SIXTH/_MM(k);
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


void TINS2DT3B::getMomentum()
{
   if (_verbose>2)
      cout << "Solving momentum equations ..." << endl;
   size_t i;
   real_t ax, ay, bm, bv;
   LocalVect<real_t,6> b, ce, us;

   size_t j=0;
   *_b = 0;
   mesh_nodes(*_theMesh) {
      if (the_node->getCode(1)==0)
         (*_b)[j++] = 0.5*_c(node_label,1);
      if (the_node->getCode(2)==0)
         (*_b)[j++] = 0.5*_c(node_label,2);
   }
   _c = 0;

// Loop over elements
   mesh_elements(*_theMesh) {
      size_t n=element_label;
      set(the_element);
      LocalVect<real_t,6> ue(the_element,*_u);
      LocalVect<real_t,3> pe(the_element,*_p,1);
      us[0] = 0.5*(ue[2]+ue[4]); us[1] = 0.5*(ue[3]+ue[5]);
      us[2] = 0.5*(ue[4]+ue[0]); us[3] = 0.5*(ue[5]+ue[1]);
      us[4] = 0.5*(ue[0]+ue[2]); us[5] = 0.5*(ue[1]+ue[3]);

//    cfl in element
      real_t d1 = ue[0] + ue[2] + ue[4];
      real_t d2 = ue[1] + ue[3] + ue[5];
      _cfl = std::max(_cfl,_time_step*sqrt(2*(d1*d1+d2*d2)/(9*_det)));

//    Mass (R.H.S.)
      real_t cc=OFELI_SIXTH*_dens*_det/_time_step;
      b = cc*ue;
      cc *= 0.5;
      b[0] = cc*(ue[0]+us[0]); b[1] = cc*(ue[1]+us[1]);
      b[2] = cc*(ue[2]+us[2]); b[3] = cc*(ue[3]+us[3]);
      b[4] = cc*(ue[4]+us[4]); b[5] = cc*(ue[5]+us[5]);

//    Viscous Term (R.H.S.)
      cc = 0.25*_visc*_det/_Re;
      for (i=0; i<3; i++) {
         ax = cc*_dSh[i].x;
         ay = cc*_dSh[i].y;
         for (j=0; j<3; j++) {
            b[2*i  ] -= ay*_dSh[j].x*ue[2*j+1] + (2*ax*_dSh[j].x+ay*_dSh[j].y)*ue[2*j  ];
            b[2*i+1] -= ax*_dSh[j].y*ue[2*j  ] + (2*ay*_dSh[j].y+ax*_dSh[j].x)*ue[2*j+1];
         }
      }

//    Convection
      cc = 0.5*OFELI_TWELVETH*_det*_dens;
      Point<real_t> du = _dSh[0]*ue[0] + _dSh[1]*ue[2] + _dSh[2]*ue[4];
      Point<real_t> dv = _dSh[0]*ue[1] + _dSh[1]*ue[3] + _dSh[2]*ue[5];
      for (i=0; i<3; i++) {
         ce[2*i  ] = cc*((d1 + ue[2*i])*du.x + (d2 + ue[2*i+1])*du.y);
         ce[2*i+1] = cc*((d1 + ue[2*i])*dv.x + (d2 + ue[2*i+1])*dv.y);
      }
      Axpy(-1.5,ce,b);
      element_assembly(The_element,ce,_c);

//    Pressure Gradient
      cc = OFELI_SIXTH*_det*(pe(1)+pe(2)+pe(3));
      for (i=0; i<3; i++) {
         b[2*i  ] += cc*_dSh[i].x;
         b[2*i+1] += cc*_dSh[i].y;
      }

//    Body Force
      if (_bf_given) {
         for (i=1; i<=3; i++) {
            b(2*i-1) += (*_bf)(2*the_element->getNodeLabel(i)-1)*_det*OFELI_SIXTH;
            b(2*i  ) += (*_bf)(2*the_element->getNodeLabel(i)  )*_det*OFELI_SIXTH;
         }
      }

//    Bubble velocity
      Point<real_t> dbp=-0.225*_det*(_dSh[0]*pe(1)+_dSh[1]*pe(2)+_dSh[2]*pe(3));
      bm = 0.144642857*_det*_dens/_time_step;
      bv = 1.0125*_visc/_Re*_det*(_dSh[0]*_dSh[0]+_dSh[1]*_dSh[1]+_dSh[2]*_dSh[2]);
      _ub(n,1) = (dbp.x + _ub(2*n-1)*(bm-bv))/(bm+bv);
      _ub(n,2) = (dbp.y + _ub(2*n  )*(bm-bv))/(bm+bv);

//    Assembly of R.H.S.
      element_assembly(The_element,b,*_b);
   }

// Tractions
   if (_sf_given) {
      LocalVect<real_t,4> ss;
      mesh_sides(*_theMesh) {
         set(the_side);
         ss[0] = 0.5*_length*(*_sf)(The_side(1)->n(),1);
         ss[1] = 0.5*_length*(*_sf)(The_side(1)->n(),2);
         ss[2] = 0.5*_length*(*_sf)(The_side(2)->n(),1);
         ss[3] = 0.5*_length*(*_sf)(The_side(2)->n(),2);
         side_assembly(The_side,ss,*_b);
      }
   }

// Solve the linear system
   real_t toler=1.e-7;
#ifdef USE_EIGEN
   LinearSolver<double> ls(1000,toler,1);
   int nb_it = ls.solve(_VM,*_b,_uu,CG_SOLVER,ILU_PREC);
#else
   LinearSolver<real_t> ls(_VM,*_b,_uu);
   ls.setTolerance(toler);
   int nb_it = ls.solve(CG_SOLVER,DILU_PREC);
#endif
   if (_bc_given)
      _u->insertBC(_uu,*_bc);
   if (_verbose>1)
      cout << "Nb. of CG iterations for Momentum: " << nb_it << endl;
}


int TINS2DT3B::getPressure()
{
   Vect<real_t> b(*_theMesh,1,NODE_DOF);
   if (_verbose>2)
      cout << "Solving pressure equation ..." << endl;
   mesh_elements(*_theMesh) {
      set(the_element);
      size_t n=element_label;
      LocalVect<real_t,6> ue(the_element,*_u);
      real_t d=OFELI_SIXTH*_det*(_dSh[0].x*ue[0] + _dSh[0].y*ue[1] +
                                 _dSh[1].x*ue[2] + _dSh[1].y*ue[3] +
                                 _dSh[2].x*ue[4] + _dSh[2].y*ue[5]);
      b(The_element(1)->n()) += 0.225*_det*(_dSh[0].x*_ub(n,1)+_dSh[0].y*_ub(n,2)) - d;
      b(The_element(2)->n()) += 0.225*_det*(_dSh[1].x*_ub(n,1)+_dSh[1].y*_ub(n,2)) - d;
      b(The_element(3)->n()) += 0.225*_det*(_dSh[2].x*_ub(n,1)+_dSh[2].y*_ub(n,2)) - d;
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
   *_p += 2.0*_q;
   return nb_it;
}


void TINS2DT3B::updateVelocity()
{
   if (_verbose>2)
      cout << "Updating velocity ..." << endl;
   real_t z = 14.0/9.0*_time_step;
   mesh_elements(*_theMesh) {
      set(the_element);
      LocalVect<real_t,3> qe(the_element,_q,1);
      _ub(element_label,1) -= z*(_dSh[0].x*qe(1)+_dSh[1].x*qe(2)+_dSh[2].x*qe(3));
      _ub(element_label,2) -= z*(_dSh[0].y*qe(1)+_dSh[1].y*qe(2)+_dSh[2].y*qe(3));
      real_t d = OFELI_SIXTH*_det*(qe(1) + qe(2) + qe(3));
      for (size_t i=0; i<3; i++) {
         size_t n = The_element(i+1)->n();
         if ((*_theMesh)[n]->getCode(1)==0)
            (*_u)(n,1) += d*_dSh[i].x/_MM(n);
         if ((*_theMesh)[n]->getCode(2)==0)
            (*_u)(n,2) += d*_dSh[i].y/_MM(n);
      }
   }
}


void TINS2DT3B::setBuyoancy()
{
}

} /* namespace OFELI */
