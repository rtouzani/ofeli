/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2025 Rachid Touzani

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

                            Implementation of Class ICPG2DT
                
       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/

#include "equations/cl/ICPG2DT.h"
#include "mesh/Mesh.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "linear_algebra/LocalVect_impl.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/Point.h"
#include "shape_functions/Triang3.h"
#include <algorithm>

using std::min;
using std::max;
using std::to_string;

namespace OFELI {

ICPG2DT::e_fptr ICPG2DT::fsolver [] = {
                                       &ICPG2DT::RiemannSolver_ROE,
                                       &ICPG2DT::RiemannSolver_VFROE,
                                       &ICPG2DT::RiemannSolver_LF,
                                       &ICPG2DT::RiemannSolver_RUSANOV,
                                       &ICPG2DT::RiemannSolver_HLL,
                                       &ICPG2DT::RiemannSolver_HLLC
                                      };


ICPG2DT::ICPG2DT(Mesh& ms)
        : Muscl2DT(ms)
{
   init();
   _init_alloc = true;
   _r = new Vect<real_t>(*_theMesh,ELEMENT_DOF,1);
   _v = new Vect<real_t>(*_theMesh,ELEMENT_DOF,2);
   _p = new Vect<real_t>(*_theMesh,ELEMENT_DOF,1);
   _min_density = 1.e-9;
}


ICPG2DT::ICPG2DT(Mesh&         ms,
                 Vect<real_t>& r,
                 Vect<real_t>& v,
                 Vect<real_t>& p)
        : Muscl2DT(ms)
{
   _r = &r;
   _v = &v;
   _p = &p;
   _init_alloc = false;
   _min_density = 1.e-9;
   init();
}


void ICPG2DT::init()
{
   _Lr.setSize(_nb_sides);
   _Lv.setSize(_nb_sides,2);
   _Lp.setSize(_nb_sides);
   _Rr.setSize(_nb_sides);
   _Rv.setSize(_nb_sides,2);
   _Rp.setSize(_nb_sides);
   _Fr.setSize(_nb_sides);
   _Frv.setSize(_nb_sides,2);
   _FE.setSize(_nb_sides);
   setGamma(1.4);
   _ReferenceLength = _MinimumEdgeLength;
   setCFL(0.2);
   setTimeStep(0.);
   mySolver = fsolver[0];
}


ICPG2DT::~ICPG2DT()
{
   if (_init_alloc) {
      delete _r;
      delete _v;
      delete _p;
   }
}


void ICPG2DT::setSolver(SolverType s)
{
   if (s > 6)
      mySolver = nullptr;
   else
      mySolver = fsolver[int(s)];
}


void ICPG2DT::setInitialConditionShockTube(const LocalVect<real_t,4>& BcL,
                                           const LocalVect<real_t,4>& BcR,
                                           real_t                     x0)
{
// Two differents zones (for shock tube purpose)
   MESH_EL {
      size_t n = element_label;
//    Defining two zones
      if (Triang3(the_element).getCenter().x<x0) {
         _r->set(n,BcL(1));
         _v->set(n,1,BcL(2));
         _v->set(n,2,BcL(3));
         _p->set(n,BcL(4));
      }
      else {
         _r->set(n,BcR(1));
         _v->set(n,1,BcR(2));
         _v->set(n,2,BcR(3));
         _p->set(n,BcR(4));
      }
   }
}


void ICPG2DT::setInitialCondition(const LocalVect<real_t,4>& u)
{
   MESH_EL {
      size_t n = element_label;
      _r->set(n,u(1));
      _v->set(n,1,u(2));
      _v->set(n,2,u(3));
      _p->set(n,u(4));
   }
}


void ICPG2DT::setReconstruction()
{
   if (Muscl2DT::setReconstruction(*_r,_Lr,_Rr,1))
      throw OFELIException("ICPG2DT::setReconstruction(): Reconstruction of rho failed");
   if (Muscl2DT::setReconstruction(*_v,_Lv,_Rv,1))
      throw OFELIException("ICPG2DT::setReconstruction(): Reconstruction of u failed");
   if (Muscl2DT::setReconstruction(*_v,_Lv,_Rv,2))
      throw OFELIException("ICPG2DT::setReconstruction(): Reconstruction of v failed");
   if (Muscl2DT::setReconstruction(*_p,_Lp,_Rp,1))
      throw OFELIException("ICPG2DT::setReconstruction(): Reconstruction of p failed");
}


real_t ICPG2DT::RiemannSolver_ROE(int s)
{
   real_t lambda_entropy_left, lambda_entropy_right;          // entropy correction
   real_t rho_star, u_star, p_star, c_star;                   // entropy correction state
   real_t alpha1, alpha2, alpha4;                             // amplitude of the waves
   real_t sqrt_rho_L=sqrt(_Lr(s)), sqrt_rho_R=sqrt(_Rr(s));
   real_t one_over_sqrt = 1./(sqrt_rho_L+sqrt_rho_R);         // normalization factor of the Roe average
   real_t one_over_rho_L = 1./_Lr(s);                         // coefficient for Roe average
   real_t one_over_rho_R = 1./_Rr(s);
   sqrt_rho_L *= one_over_sqrt;                               // normalized coefficient
   sqrt_rho_R *= one_over_sqrt;                               // same !

// Averages are hence normalized !!
   if (_Lr(s)<_min_density)
      throw OFELIException("ICPG2DT::RiemannSolver_ROE(int): Left minimum density rho = " + to_string(_Lr(s)));
   if (_Rr(s)<_min_density)
      throw OFELIException("ICPG2DT::RiemannSolver_ROE(int): Right minimum density rho = " + to_string(_Rr(s)));

// get data
   _tLr = _Lr(s); _tRr = _Rr(s);
   _tLv[0] = _Lv(s,1); _tLv[1] = _Lv(s,2); 
   _tRv[0] = _Rv(s,1); _tRv[1] = _Rv(s,2);
   _Lrv[0] = _tLr*_tLv[0]; _Lrv[1] = _tLr*_tLv[1];
   _Rrv[0] = _tRr*_tRv[0]; _Rrv[1] = _tRr*_tRv[1];

// Compute pressure, internal energy and sound speed
   _tLp = _Lp(s);
   _Le = _tLp/(_Gamma-1.) + 0.5*_tLr*(_tLv[0]*_tLv[0]+_tLv[1]*_tLv[1]);
   if (_tLp<0.)
      cerr << "WARNING: Gp < 0" << endl;
   _LH = (_Le+_tLp)*one_over_rho_L;
   _tRp = _Rp(s);
   _Re = _tRp/(_Gamma-1.) + 0.5*_tRr*(_tRv[0]*_tRv[0]+_tRv[1]*_tRv[1]);
   if (_tRp<0.)
      cerr << "WARNING: Dp < 0" << endl;
   _RH = (_Re+_tRp)*one_over_rho_R;

// compute Roe average
   Point<real_t> Roe_u(sqrt_rho_L*_tLv[0] + sqrt_rho_R*_tRv[0],
                       sqrt_rho_L*_tLv[1] + sqrt_rho_R*_tRv[1]);
   real_t Roe_H = sqrt_rho_L*_LH + sqrt_rho_R*_RH;
   real_t Roe_c = (_Gamma-1.)*(Roe_H-0.5*Roe_u.NNorm());   //  c*c !!!!
   if (Roe_c<0)
      throw OFELIException("ICPG2DT::RiemannSolver_ROE(int): Roe sound speed negative in point "
                           + to_string(s) + ", on boundary (0/1) "
                           + to_string(_theMesh->getPtrSide(s)->isOnBoundary())
                           + ": "+ to_string(Triang3(_theMesh->getPtrSide(s)).getCenter().x) + "  " 
                           + to_string(Triang3(_theMesh->getPtrSide(s)).getCenter().y));
   Roe_c = sqrt(Roe_c);

// compute eigenvalues
   real_t lambda1 = Roe_u.x - Roe_c, lambda2 = Roe_u.x, lambda4 = Roe_u.x + Roe_c;

// calculate flux
   if (lambda1>0.) {
      _Fr(s)    = _Lrv[0];
      _Frv(s,1) = _Lrv[0]*_tLv[0] + _tLp;
      _Frv(s,2) = _Lrv[1]*_tLv[0];
      _FE(s)    = _tLv[0]*(_Le+_tLp);
      return (Roe_u.x+Roe_c);
   }

   alpha2  = (_tRr-_tLr)*(Roe_H-Roe_u.x*Roe_u.x) + (_Rrv[0]-_Lrv[0])*Roe_u.x -
             (_Re-_Le) + ((_Rrv[0]-_Lrv[0])-(_tRr-_tLr)*Roe_u.y)*Roe_u.y;
   alpha2 *= (_Gamma-1)/(Roe_c*Roe_c);
   alpha1  = (_tRr-_tLr)*(Roe_c+Roe_u.y)-(_Rrv[0]-_Lrv[0]) - alpha2*Roe_c;
   alpha1 *= 0.5/Roe_c;

   if (lambda2>=0.) { 
//    Entropy fix if necessary (Toro 1997, p. 338)
      lambda_entropy_left = _tLv[0] - _Lc;
      rho_star = _tLr + alpha1;
      u_star = (_Lrv[0]+alpha1*(Roe_u.y-Roe_c))/rho_star;
      p_star = (_Gamma-1)*(_Le+alpha1*(Roe_H-Roe_u.x*Roe_c)-0.5*rho_star*u_star*u_star);
      if (p_star<0)
         throw OFELIException("ICPG2DT::RiemannSolver_ROE(int): ERROR: negative pressure"
                              " in left entropy fix method");
      c_star = sqrt(_Gamma*p_star/rho_star);
      lambda_entropy_right = u_star - c_star;
      if (lambda_entropy_left<0. && 0.<lambda_entropy_right) //we have to fix the entropy
         lambda1 = lambda_entropy_left*(lambda_entropy_right-lambda1)/(lambda_entropy_right-lambda_entropy_left);
      _Fr(s)    = _Lrv[0] + alpha1*lambda1;
      _Frv(s,1) = _Lrv[0]*_tLv[0] + _tLp + alpha1*lambda1*(Roe_u.x-Roe_c);
      _Frv(s,2) = _Lrv[0]*_tLv[1] + alpha1*lambda1*Roe_u[1];
      _FE(s)    = _tLv[0]*(_Le+_tLp) + alpha1*lambda1*(Roe_H-Roe_u.x*Roe_c);
      return (Roe_u.x+Roe_c);
   }

   alpha4 = (_tRr-_tLr) - alpha2 - alpha1;
   if (lambda4>0.) {
//    Fix entropy if necessary (Toro 1997, p. 338)
      lambda_entropy_right = _tRv[0] + _Rc;
      rho_star = _tRr - alpha4;
      u_star = (_Rrv[0]-alpha4*(Roe_u.y+Roe_c))/rho_star;
      p_star = (_Gamma-1)*(_Re-alpha4*(Roe_H+Roe_u.x*Roe_c)-0.5*rho_star*u_star*u_star);
      if (p_star<0) {
         cerr << "ERROR: negative pressure in right entropy fix method" << endl;
         exit(2);
      }
      c_star = sqrt(_Gamma*p_star/rho_star);
      lambda_entropy_left = u_star + c_star;
      if (lambda_entropy_left<0. && 0.<lambda_entropy_right) // we have to fix the entropy
         lambda4 = lambda_entropy_right*(lambda4-lambda_entropy_left)/(lambda_entropy_right-lambda_entropy_left);
      _Fr.set(s,_Rrv[0] - alpha4*lambda4);
      _Frv.set(s,1,_Rrv[0]*_tRv[0] + _tRp - alpha4*lambda4*(Roe_u.x+Roe_c));
      _Frv.set(s,2,_Rrv[0]*_tRv[1] - alpha4*lambda4*Roe_u.y);
      _FE.set(s,_tRv[0]*(_Re+_tRp) - alpha4*lambda4*(Roe_H+Roe_u.x*Roe_c));
      return (-Roe_u.x+Roe_c);
   }

// lambda4>0
   _Fr.set(s,_Rrv[0]);
   _Frv.set(s,1,_Rrv[0]*_tRv[0] + _tRp);
   _Frv.set(s,2,_Rrv[0]*_tRv[1]);
   _FE.set(s,_tRv[0]*(_Re+_tRp));
   return (-Roe_u.x+Roe_c);
}


real_t ICPG2DT::RiemannSolver_VFROE(int s)
{
   real_t lambda_entropy_left, lambda_entropy_right;           // entropy correction
   real_t rho_star, u_star, p_star, c_star;                    // entropy correction state
   real_t alpha1, alpha2, alpha4;                              // amplitude of the waves
   real_t sqrt_rho_L=sqrt(_Lr(s)), sqrt_rho_R=sqrt(_Rr(s));
   real_t one_over_sqrt=1./(sqrt_rho_L+sqrt_rho_R);            // normalization factor of the Roe average
   real_t one_over_rho_L = 1./_Lr(s);                          // coefficient for Roe average
   real_t one_over_rho_R = 1./_Rr(s);                          // same !
   sqrt_rho_L *= one_over_sqrt;                                // normalized coefficient
   sqrt_rho_R *= one_over_sqrt;                                // same !

// Averages are normalized
   if (_Lr(s)<_min_density) {
      cerr << "ERROR: left minimum density rho = " << _Lr(s) << endl;
      exit(2);
   }
   if (_Rr(s)<_min_density) {
      cerr << "ERROR: right minimum density rho = " << _Rr(s) << endl;
      exit(2);
   }
   _tLr = _Lr(s); _tRr = _Rr(s);

// compute velocities
   _tLv[0] = _Lv(s,1); _tLv[1] = _Lv(s,2);
   _tRv[0] = _Rv(s,1); _tRv[1] = _Rv(s,2);
   _Lrv[0] = _tLr*_tLv[0]; _Lrv[1] = _tLr*_tLv[1];
   _Rrv[0] = _tRr*_tRv[0]; _Rrv[1] = _tRr*_tRv[1];

// compute pressure and internal energy and sound speed
   _tLp = _Lp(s);
   _Le = _tLp/(_Gamma-1.) + 0.5*_tLr*(_tLv[0]*_tLv[0]+_tLv[1]*_tLv[1]);
   if (_tLp<0.)
      cerr << "WARNING: Gp < 0" << endl;
   _LH = (_Le+_tLp)*one_over_rho_L;

   _tRp = _Rp(s);
   _Re = _tRp/(_Gamma-1.) + 0.5*_tRr*(_tRv[0]*_tRv[0]+_tRv[1]*_tRv[1]);
   if (_tRp<0.)
      cerr << "WARNING: DP < 0" << endl;
   _RH = (_Re+_tRp)*one_over_rho_R; 

// compute Roe average
   Point<real_t> Roe_u(0.5*(_tLv[0]+_tRv[0]),0.5*(_tLv[1]+_tRv[1]));
   real_t Roe_H = 0.5*(_LH+_RH);
   real_t Roe_c = (_Gamma-1.)*(Roe_H-0.5*Roe_u.NNorm());   //  c*c !!!!

   if (Roe_c<0) {
      cerr << "ERROR: Roe sound speed negative " << Roe_c
           << "       in point " << s << ", On boundary (0/1) ? " 
           << _theMesh->getPtrSide(s)->isOnBoundary()
           << ", x = " << Triang3(_theMesh->getPtrSide(s)).getCenter().x 
           << ", y = " << Triang3(_theMesh->getPtrSide(s)).getCenter().y << endl;
      exit(2);
   }
   Roe_c = sqrt(Roe_c);

// compute eigenvalues
   real_t lambda1 = Roe_u.x - Roe_c, lambda2 = Roe_u.x, lambda4 = Roe_u.x + Roe_c;

// compute flux
   if (0.<lambda1) {
      _Fr.set(s,_Lrv[0]);
      _Frv.set(s,1,_Lrv[0]*_tLv[0] + _tLp);
      _Frv.set(s,2,_Lrv[1]*_tLv[0]);
      _FE.set(s,_tLv[0]*(_Le+_tLp));
      return (Roe_u.x+Roe_c);
   }

   alpha2 = (_tRr-_tLr)*(Roe_H-Roe_u.x*Roe_u.x) + (_Rrv[1]-_Lrv[1])*Roe_u.x - (_Re-_Le) + 
			((_Rrv[1]-_Lrv[1])-(_tRr-_tLr)*Roe_u.y)*Roe_u.y;
   alpha2 *= (_Gamma-1)/(Roe_c*Roe_c);
   alpha1 = (_tRr-_tLr)*(Roe_c+Roe_u.x)-(_Rrv[0]-_Lrv[0])-alpha2*Roe_c;
   alpha1 *= 0.5/Roe_c;

   if (lambda2>=0.) { 
//    entropy fix if necessary (Toro 1997, p. 338)
      lambda_entropy_left = _tLv[0] - _Lc;
      rho_star = _tLr + alpha1;
      u_star = (_Lrv[0]+alpha1*(Roe_u.x-Roe_c))/rho_star;
      p_star = (_Gamma-1)*(_Le+alpha1*(Roe_H-Roe_u.x*Roe_c)-0.5*rho_star*u_star*u_star);
      if (p_star<0) {
         cerr << "ERROR: negative pressure in left entropy fix method" << endl;
         exit(2);
      }
      c_star = sqrt(_Gamma*p_star/rho_star);
      lambda_entropy_right = u_star - c_star;
      if (lambda_entropy_left<0. && 0.<lambda_entropy_right)
         lambda1 = lambda_entropy_left*(lambda_entropy_right-lambda1)/(lambda_entropy_right-lambda_entropy_left);
      _Fr.set(s,_Lrv[0] + alpha1*lambda1);
      _Frv.set(s,1,_Lrv[0]*_tLv[0] + _tLp + alpha1*lambda1*(Roe_u.x-Roe_c));
      _Frv.set(s,2,_Lrv[0]*_tLv[1] + alpha1*lambda1*Roe_u.y);
      _FE.set(s,_tLv[0]*(_Le+_tLp) + alpha1*lambda1*(Roe_H-Roe_u.x*Roe_c));
      return (Roe_u.x+Roe_c);
   }

   alpha4 = (_tRr-_tLr) - alpha2 - alpha1;
   if (0.<lambda4) {  // Fix entropy if necessary (Toro 1997, p. 338)
      lambda_entropy_right = _tRv[0] + _Rc;
      rho_star = _tRr - alpha4;
      u_star = (_Rrv[0]-alpha4*(Roe_u.x+Roe_c))/rho_star;
      p_star = (_Gamma-1)*(_Re-alpha4*(Roe_H+Roe_u.x*Roe_c)-0.5*rho_star*u_star*u_star);
      if (p_star<0)
         throw OFELIException("ICPG2DT::RiemannSolver_VFROE(int): Negative pressure in right "
                              "entropy fix method");
      c_star = sqrt(_Gamma*p_star/rho_star);
      lambda_entropy_left = u_star + c_star;
      if (lambda_entropy_left<0. && 0.<lambda_entropy_right) // we have to fix the entropy
         lambda4 = lambda_entropy_right*(lambda4-lambda_entropy_left)/(lambda_entropy_right-lambda_entropy_left);
      _Fr.set(s,_Rrv[0] - alpha4*lambda4);
      _Frv.set(s,1,_Rrv[0]*_tRv[0] + _tRp - alpha4*lambda4*(Roe_u.x+Roe_c));
      _Frv.set(s,2, _Rrv[0]*_tRv[1] - alpha4*lambda4*Roe_u.y);
      _FE.set(s,_tRv[0]*(_Re+_tRp) - alpha4*lambda4*(Roe_H+Roe_u.x*Roe_c));
      return (-Roe_u.x+Roe_c);
   }
   _Fr.set(s,_Rrv[0]);
   _Frv.set(s,1,_Rrv[0]*_tRv[0] + _tRp);
   _Frv.set(s,2,_Rrv[0]*_tRv[1]);
   _FE.set(s,_tRv[0]*(_Re+_tRp));
   return (-Roe_u.x+Roe_c);
}


real_t ICPG2DT::RiemannSolver_LF(int s)
{
   real_t one_over_rho_G = 1./_Lr(s), one_over_rho_D = 1./_Rr(s);
   if (_Lr(s)<_min_density) {
      if (_verbose)
         cerr << "WARNING: left minimum density rho = " << _Lr(s) << endl;
      _Lr(s) = _min_density;
   }
   if (_Rr(s)<_min_density) {
      if (_verbose)
         cerr << "WARNING: right minimum density rho = " << _Rr(s) << endl;
      _Rr(s) = _min_density;
   }

// get data
   _tLr    = _Lr(s); _tRr    = _Rr(s);
   _tLv[0] = _Lv(s,1); _tLv[1] = _Lv(s,2);
   _tRv[0] = _Rv(s,1); _tRv[1] = _Rv(s,2);
   _Lrv[0] = _tLr*_tLv[0];
   _Lrv[1] = _tLr*_tLv[1];
   _Rrv[0] = _tRr*_tRv[0];
   _Rrv[1] = _tRr*_tRv[1];

// compute pressure, internal energy and sound speed
   _tLp = _Lp(s);
   _Le = _tLp/(_Gamma-1.) + 0.5*_tLr*(_tLv[0]*_tLv[0]+_tLv[1]*_tLv[1]);
   if (_tLp<0.)
      cerr << "WARNING: Gp<0" << endl;
   _LH = (_Le+_tLp)*one_over_rho_G; 

   _tRp = _Rp(s);
   _Re = _tRp/(_Gamma-1.) + 0.5*_tRr*(_tRv[0]*_tRv[0]+_tRv[1]*_tRv[1]);
   if (_tRp<0.)
      cerr << "WARNING: Dp<0" << endl;
   _RH = (_Re+_tRp)*one_over_rho_D; 

// use VFRoe method to compute eigenvalues
// compute VFRoe average
   Point<real_t> Roe_u(0.5*(_tLv[0]+_tRv[0]),0.5*(_tLv[1]+_tRv[1]));
   real_t Roe_H = 0.5*(_LH+_RH);
   real_t Roe_c = (_Gamma-1.)*(Roe_H-0.5*Roe_u.NNorm()); // Here c*c

   if (Roe_c<0) {
      cerr << "ERROR: Roe sound speed negative " << Roe_c
           << "at point " << s << ", On boundary ? (0/1)" << _theMesh->getPtrSide(s)->isOnBoundary() << endl;
      exit(2);
   }
   Roe_c = sqrt(Roe_c);

// compute LF flux
   _Fr.set(s,0.5*(_Lrv[0]+_Rrv[0]));
   _Frv.set(s,1,0.5*(_Lrv[0]*_tLv[0]+_tLp+_Rrv[0]*_tRv[0]+_tRp));
   _Frv.set(s,2,0.5*(_Lrv[0]*_tLv[1]+_Rrv[0]*_tRv[1]));
   _FE.set(s,0.5*(_tLv[0]*(_Le+_tLp)+_tRv[0]*(_Re+_tRp)));
   if (_TimeStep) {
      _Fr.add(s,-0.5*_ReferenceLength/_TimeStep*(_tRr-_tLr));
      _Frv.add(s,1,-0.5*_ReferenceLength/_TimeStep*(_Rrv[0]-_Lrv[0]));
      _Frv.add(s,2,-0.5*_ReferenceLength/_TimeStep*(_Rrv[1]-_Lrv[1]));
      _FE.add(s,-0.5*_ReferenceLength/_TimeStep*(_Re-_Le));
   }
   return (_ReferenceLength/_TimeStep);
}


real_t ICPG2DT::RiemannSolver_RUSANOV(int s)
{
   real_t Splus, one_over_rho_G = 1./_Lr(s), one_over_rho_D = 1./_Rr(s);
   if (_Lr(s)<_min_density) {
      if (_verbose)
         cerr << "WARNING: Left minimum density rho = " << _Lr(s) << endl;
      _Lr.set(s,_min_density);
   }
   if (_Rr(s)<_min_density) {
      if (_verbose)
         cerr << "WARNING: Right minimum density rho = " << _Rr(s) << endl;
      _Rr.set(s,_min_density);
   }

// get data
   _tLr = _Lr(s); _tRr = _Rr(s);
   _tLv[0] = _Lv(s,1); _tLv[1] = _Lv(s,2);
   _tRv[0] = _Rv(s,1); _tRv[1] = _Rv(s,2);
   _Lrv[0] = _tLr*_tLv[0]; _Lrv[1] = _tLr*_tLv[1];
   _Rrv[0] = _tRr*_tRv[0]; _Rrv[1] = _tRr*_tRv[1];

// Compute pressure, internal energy and sound speed
   _tLp = _Lp(s);
   _Le  = _tLp/(_Gamma-1.) + 0.5*_tLr*(_tLv[0]*_tLv[0]+_tLv[1]*_tLv[1]);
   if (_tLp<0.)
      cerr << "WARNING: Gp < 0" << endl;
   _LH  = (_Le+_tLp)*one_over_rho_G;
   _tRp = _Rp(s);
   _Re  = _tRp/(_Gamma-1.) + 0.5*_tRr*(_tRv[0]*_tRv[0]+_tRv[1]*_tRv[1]);
   if (_tRp<0.)
      cerr << "WARNING: Dp < 0" << endl;
   _RH = (_Re+_tRp)*one_over_rho_D;

// compute S+
   _Lc = (_Gamma-1.)*(_LH-0.5*(_tLv[0]*_tLv[0]+_tLv[1]*_tLv[1]));
   _Rc = (_Gamma-1.)*(_RH-0.5*(_tRv[0]*_tRv[0]+_tRv[1]*_tRv[1]));

   if (_Lc<0. || _Rc<0.) {
      cerr << "ERROR: sound speed negative " << _Lc << "  " << _Rc << "       at point "
           << s << ", on boundary (0/1) ? "
           << _theMesh->getPtrSide(s)->isOnBoundary() << endl;
      exit(2);
   }
   _Lc = sqrt(_Lc); _Rc = sqrt(_Rc);
   Splus = max(fabs(_tLv[0]-_Lc),fabs(_tLv[0]+_Lc));
   Splus = max(Splus,fabs(_tRv[0]-_Rc));
   Splus = max(Splus,fabs(_tRv[0]+_Rc));

// Compute Rusanov flux
   _Fr.set(s,0.5*(_Lrv[0]+_Rrv[0])-0.5*Splus*(_tRr-_tLr));
   _Frv.set(s,1,0.5*(_Lrv[0]*_tLv[0]+_tLp+_Rrv[0]*_tRv[0]+_tRp) - 0.5*Splus*(_Rrv[0]-_Lrv[0]));
   _Frv.set(s,2,0.5*(_Lrv[0]*_tLv[1] + _Rrv[0]*_tRv[1]) - 0.5*Splus*(_Rrv[1]-_Lrv[1]));
   _FE.set(s,0.5*(_tLv[0]*(_Le+_tLp) + _tRv[0]*(_Re+_tRp)) - 0.5*Splus*(_Re-_Le));
   return Splus;
}


real_t ICPG2DT::RiemannSolver_HLL(int s)
{
   real_t one_over_rho_L = 1./_Lr(s), one_over_rho_R = 1./_Rr(s);
   real_t SL, SR, Splus, C1, C2, C3;
   if (_Lr(s)<_min_density) {
      if (_verbose)
         cerr << "WARNING: Left minimum density rho = " << _Lr(s) << endl;
      _Lr.set(s,_min_density);
   }
   if (_Rr(s)<_min_density) {
      if (_verbose)
         cerr << "WARNING: Right minimum density rho = " << _Rr(s) << endl;
      _Rr.set(s,_min_density);
   }
   
// get data
   _tLr = _Lr(s); _tRr = _Rr(s);
   _tLv[0] = _Lv(s,1); _tLv[1] = _Lv(s,2);
   _tRv[0] = _Rv(s,1); _tRv[1] = _Rv(s,2);
   _Lrv[0] = _tLr*_tLv[0]; _Lrv[1] = _tLr*_tLv[1];
   _Rrv[0] = _tRr*_tRv[0]; _Rrv[1] = _tRr*_tRv[1];

// compute pressure, internal energy and sound speed
   _tLp = _Lp(s);
   _Le = _tLp/(_Gamma-1.) + 0.5*_tLr*(_tLv[0]*_tLv[0]+_tLv[1]*_tLv[1]);
   if (_tLp<0.)
      cerr << "WARNING: Gp < 0" << endl;
   _LH = (_Le+_tLp)*one_over_rho_L;
   _tRp = _Rp(s);
   _Re = _tRp/(_Gamma-1.) + 0.5*_tRr*(_tRv[0]*_tRv[0]+_tRv[1]*_tRv[1]);
   if (_tRp<0.)
      cerr << "WARNING: Dp < 0" << endl;
   _RH = (_Re+_tRp)*one_over_rho_R;

// compute S+
   _Lc = (_Gamma-1.)*(_LH-(_tLv[0]*_tLv[0]+_tLv[1]*_tLv[1])*0.5); // Here c*c
   _Rc = (_Gamma-1.)*(_RH-(_tRv[0]*_tRv[0]+_tRv[1]*_tRv[1])*0.5); // Here c*c
   if (_Lc<0. || _Rc<0.) {
      cerr << "ERROR: sound speed negative " << _Lc << "  " << _Rc
           << "on point " << s << ", on boundary (0/1) ? "
           << _theMesh->getPtrSide(s)->isOnBoundary() << endl;
      exit(2);
   }
   _Lc = sqrt(_Lc); _Rc = sqrt(_Rc);

   Splus = max(fabs(_tLv[0] - _Lc),fabs(_tLv[0] + _Lc));
   Splus = max(Splus,fabs(_tRv[0] - _Rc));
   Splus = max(Splus,fabs(_tRv[0] + _Rc));
   SL = -Splus;
   SR = Splus;
   C1 =  SR/(SR-SL);
   C2 = -SL/(SR-SL);
   C3 =  SL*SR/(SR-SL);

// compute HLL flux
   if (SL>=0.) { // Flux = FG
      _Fr.set(s,_Lrv[0]);
      _Frv.set(s,1,_Lrv[0]*_tLv[0] + _tLp);
      _Frv.set(s,2,_Lrv[0]*_tLv[1]);
      _FE.set(s,_tLv[0]*(_Le+_tLp));
      return Splus;
   }
   if (SR<=0.) { // Flux = FD
      _Fr.set(s,_Rrv[0]);
      _Frv.set(s,1,_Rrv[0]*_tRv[0] + _tRp);
      _Frv.set(s,2,_Rrv[0]*_tRv[1]);
      _FE.set(s,_tRv[0]*(_Re+_tRp));
   }

// (SL<=0.) && (SR >=0.)) => Flux = FHLL
   _Fr.set(s,C1*_Lrv[0] + C2*_Rrv[0] + C3*(_tRr-_tLr));
   _Frv.set(s,1,C1*(_Lrv[0]*_tLv[0]+_tLp) + C2*(_Rrv[0]*_tRv[0]+_tRp) + C3*(_Rrv[0]-_Lrv[0]));
   _Frv.set(s,2,C1*(_Lrv[0]*_tLv[1]) + C2*(_Rrv[0]*_tRv[1]) + C3*(_Rrv[1]-_Lrv[1]));
   _FE.set(s,C1*(_tLv[0]*(_Le+_tLp)) + C2*(_tRv[0]*(_Re+_tRp)) + C3*(_Re-_Le));
   return Splus;
}


real_t ICPG2DT::RiemannSolver_HLLC(int s)
{
   real_t one_over_rho_G = 1./_Lr(s), one_over_rho_D = 1./_Rr(s);
   real_t SL, SR, Splus, Sstar, AL, AR;
   real_t rhoGstar, rhoDstar, uGstar, uDstar, vGstar, vDstar, EGstar, EDstar;

   if (_Lr(s)<_min_density) {
      cerr << "WARNING: Left minimum density rho = " << _Lr(s) << endl;
      _Lr.set(s,_min_density);
   }
   if (_Rr(s)<_min_density) {
      cerr << "WARNING: right minimum density rho = " << _Rr(s) << endl;
      _Rr.set(s,_min_density);
   }

// get data
   _tLr = _Lr(s); _tRr = _Rr(s);
   _tLv[0] = _Lv(s,1); _tLv[1] = _Lv(s,2);
   _tRv[0] = _Rv(s,1); _tRv[1] = _Rv(s,2);
   _Lrv[0] = _tLr*_tLv[0]; _Lrv[1] = _tLr*_tLv[1];
   _Rrv[0] = _tRr*_tRv[0]; _Rrv[1] = _tRr*_tRv[1];

// compute pressure, internal energy and sound speed
   _tLp = _Lp(s);
   _Le  = _tLp/(_Gamma-1.) + 0.5*_tLr*(_tLv[0]*_tLv[0]+_tLv[1]*_tLv[1]);
   _LH  = (_Le+_tLp)*one_over_rho_G; 
   _tRp = _Rp(s);
   _Re  = _tRp/(_Gamma-1.) + 0.5*_tRr*(_tRv[0]*_tRv[0]+_tRv[1]*_tRv[1]);
   _RH  = (_Re+_tRp)*one_over_rho_D;

// compute S+
   _Lc = (_Gamma-1.)*(_LH-(_tLv[0]*_tLv[0]+_tLv[1]*_tLv[1])*0.5);
   _Rc = (_Gamma-1.)*(_RH-(_tRv[0]*_tRv[0]+_tRv[1]*_tRv[1])*0.5);

   if (_Lc<0. || _Rc<0.) {
      cerr << "WARNING: sound speed negative " << _Lc << "  " << _Rc
           << "at point " << s << ", on boundary (0/1) ? "
           << _theMesh->getPtrSide(s)->isOnBoundary() << endl;
      exit(2);
   }
   _Lc = sqrt(_Lc);
   _Rc = sqrt(_Rc);
   Splus = max(fabs(_tLv[0]-_Lc),fabs(_tLv[0] + _Lc));
   Splus = max(Splus,fabs(_tRv[0] - _Rc));
   Splus = max(Splus,fabs(_tRv[0] + _Rc));
   SL = -Splus;
   SR = Splus;
   Sstar = (_tRp - _tLp + _tLv[0]*(SL-_tLv[0]) - _Rrv[0]*(SR-_tRv[0]))/(_tLr*(SL-_tLv[0]) - _tRr*(SR-_tRv[0])) ;
   AL = _tLr*(SL-_tLv[0])/(SL-Sstar);
   AR = _tRr*(SR-_tRv[0])/(SR-Sstar);

// compute star region
   rhoGstar = AL; rhoDstar = AR;
   uGstar = AL * Sstar;
   uDstar = AR*Sstar;
   vGstar = AL * _tLv[1];
   vDstar = AR * _tRv[1];
   EGstar = AL*(_Le*one_over_rho_G + (Sstar-_tLv[0])*(Sstar + _tLp/(_tLr*(SL-_tLv[0]))));
   EDstar = AR*(_Re*one_over_rho_D + (Sstar-_tRv[0])*(Sstar + _tRp/(_tRr*(SR-_tRv[0]))));

// compute HLLC flux
   if (SL>=0.) {     // Flux = FG
      _Fr.set(s,_Lrv[0]);
      _Frv.set(s,1,_Lrv[0]*_tLv[0] + _tLp);
      _Frv.set(s,2,_Lrv[0]*_tLv[1]);
      _FE.set(s,_tLv[0]*(_Le+_tLp));
      return Splus;
   }
   if (SR<=0.) {     // Flux = FD
      _Fr.set(s,_Rrv[0]);
      _Frv.set(s,1,_Rrv[0]*_tRv[0] + _tRp);
      _Frv.set(s,2,_Rrv[0]*_tRv[1]);
      _FE.set(s,_tRv[0]*(_Re+_tRp));
      return Splus;
   }
   if (Sstar >=0.) { // Flux = FHLLCG
      _Fr.set(s,_Lrv[0] + SL*(rhoGstar-_Lr(s)));
      _Frv.set(s,1,_Lrv[0]*_tLv[0] + _tLp + SL*(uGstar-_Lrv[0]));
      _Frv.set(s,2,_Lrv[0]*_tLv[1] + SL*(vGstar-_Lrv[1]));
      _FE.set(s,_tLv[0]*(_Le+_tLp) + SL*(EGstar-_Le));
   }
   else {             // Flux = FHLLCD
      _Fr.set(s,_Rrv[0]+ SR*(rhoDstar-_Rr(s)));
      _Frv.set(s,1,_Rrv[0]*_tRv[0] + _tRp + SR*(uDstar-_Rrv[0]));
      _Frv.set(s,2,_Rrv[0]*_tRv[1] + SR*(vDstar-_Rrv[1]));
      _FE.set(s,_tRv[0]*(_Re+_tRp) + SR*(EDstar-_Re));
   }
   return Splus;
}


real_t ICPG2DT::getFlux()
{
// The solver assumes that boundary condition
// i.e. _Lv at boundary, has already been defined  !!!
   real_t lambda_max=1.e-10, lambda, fu, fv;
   MESH_SD {
      size_t s = side_label;
      Point<real_t> n = _n(s);
      _tLv[0] = _Lv(s,1); _tLv[1] = _Lv(s,2);
      _tRv[0] = _Rv(s,1); _tRv[1] = _Rv(s,2);

//    change orientation
      _Lv.set(s,1, n.x*_tLv[0] + n.y*_tLv[1]);
      _Lv.set(s,2,-n.y*_tLv[0] + n.x*_tLv[1]);
      _Rv.set(s,1, n.x*_tRv[0] + n.y*_tRv[1]);
      _Rv.set(s,2,-n.y*_tRv[0] + n.x*_tRv[1]);

//    compute the flux with the new orientation
      assert(mySolver);       // just in case
 //   lambda_max = max((this->*mySolver)(int(s)),lambda_max);
      lambda = (this->*mySolver)(int(s));
      if (lambda>lambda_max)
         lambda_max = lambda;

//    get back to initial orientation
      fu = _Frv(s,1); fv = _Frv(s,2);
      _Frv.set(s,1,n.x*fu - n.y*fv);
      _Frv.set(s,2,n.y*fu + n.x*fv);

//    RESTORE VALUES: really usefull????
      _Lv.set(s,1,_tLv[0]);
      _Lv.set(s,2,_tLv[1]);
      _Rv.set(s,1,_tRv[0]);
      _Rv.set(s,2,_tRv[1]);
   }
   return lambda_max;
}


void ICPG2DT::setBC(const Side&  sd,
                    real_t       a)
{
   size_t s = sd.n();
   _Rr(s) = _Lr(s);
   real_t tmp = (a-1)*(_Lv(s,1)*_n(s).x + _Lv(s,2)*_n(s).y);
   _Rv.set(s,1,_Lv(s,1) + tmp*_n(s).x);
   _Rv.set(s,2,_Lv(s,2) + tmp*_n(s).y);
   _Rp.set(s,_Lp(s));
}


void ICPG2DT::setBC(int    code,
                    real_t a)
{
   MESH_SD 
      if (The_side.isOnBoundary() && The_side.getCode(1)==code)
         setBC(The_side,a);
}


void ICPG2DT::setBC(real_t a)
{
   MESH_SD
      if (The_side.isOnBoundary())
         setBC(The_side,a);
}


void ICPG2DT::setBC(const Side&                sd,
                    const LocalVect<real_t,4>& u)
{
   size_t s = sd.n();
   _Rr.set(s,u[0]);
   _Rv.set(s,1,u[1]);
   _Rv.set(s,2,u[2]);
   _Rp.set(s,u[3]);
}


void ICPG2DT::setBC(int                        code,
                    const LocalVect<real_t,4>& u)
{
   MESH_SD
      if (The_side.isOnBoundary() && The_side.getCode()==code)
         setBC(The_side,u);
}


void ICPG2DT::setBC(const LocalVect<real_t,4>& u)
{
   MESH_SD
      if (The_side.isOnBoundary())
         setBC(The_side,u);
}


void ICPG2DT::setInOutFlowBC(const Side&                sd,
                             const LocalVect<real_t,4>& u)
{
   size_t s = sd.n();
   real_t Gc  = sqrt(_Gamma*_Lp(s)/_Lr(s));
   real_t GUn = _Lv(s,1)*_n(s).x + _Lv(s,2)*_n(s).y;
   int choice = (GUn>=0.) + 2*(Gc>=1.);

   switch (choice) {

      case 0:                          // GU.n < 0 and Gc < 1: subsonic inflow
         _Rr.set(s,_Lr(s));            // getting from interior
         _Rv.set(s,1,u[1]);   
         _Rv.set(s,2,u[2]);
         _Rp.set(s,u[3]);
         break;

      case 1:                          // GU.n >=0 and Gc < 1: subsonic outflow
         _Rr.set(s,_Lr(s));
         _Rv.set(s,1,_Lv(s,1));
         _Rv.set(s,2,_Lv(s,2));
         _Rp.set(s,u[3]);              // imposing pressure
         break;

      case 2:                          // GU.n < 0 and Gc > 1: supersonic inflow
         _Rr.set(s,u[0]);              // imposing density
         _Rv.set(s,1,u[1]);              // because exterior normal
         _Rv.set(s,2,u[2]);
         _Rp.set(s,u[3]);
         break;

      case 3:                          //  GU.n > 0 and Gc > 1: supersonic outflow
         _Rr.set(s,_Lr(s));
         _Rv.set(s,1,_Lv(s,1));
         _Rv.set(s,2,_Lv(s,2));
         _Rp.set(s,_Lp(s));            // getting from interior
         break;

      default :
         cerr << "ERROR in setInOutflowBC" << endl;
         exit(1);
   }
}


void ICPG2DT::setInOutFlowBC(int                        code,
                             const LocalVect<real_t,4>& u)
{
   MESH_SD
     if (The_side.isOnBoundary() && The_side.getCode()==code)
         setInOutFlowBC(The_side,u);
}


void ICPG2DT::setInOutFlowBC(const LocalVect<real_t,4>& u)
{
   MESH_SD
      if (The_side.isOnBoundary())
         setInOutFlowBC(The_side,u);
}

// Variable flux contains fluxes as computed from the Riemman problem with 
// getNeighborElement(1) as left element and getNeighborElement(2) as right element.
// If getNeighborElement(2) returns NULL, we are on the boundary and we prescribe a symmetry 
// condition


void ICPG2DT::fromPrimalToConservative()
{
   real_t G = 1./(_Gamma-1.);
   MESH_EL {
      size_t n = element_label;
      real_t rho = (*_r)(n);
      _p->set(n,(*_p)(n)*G + 0.5*rho*((*_v)(n,1)*(*_v)(n,1)+(*_v)(n,2)*(*_v)(n,2)));
      _v->set(n,1,(*_v)(n,1)*rho);
      _v->set(n,2,(*_v)(n,2)*rho);
   }
}


void ICPG2DT::fromConservativeToPrimal()
{
   MESH_EL {
      size_t n = element_label;
      _v->set(n,1,(*_v)(n,1)/(*_r)(n));
      _v->set(n,2,(*_v)(n,2)/(*_r)(n));
      _p->set(n,(_Gamma-1.)*((*_p)(n) - 0.5*(*_r)(n)*((*_v)(n,1)*(*_v)(n,1)+(*_v)(n,2)*(*_v)(n,2))));
   }
}


void ICPG2DT::forward()
{
   fromPrimalToConservative();

// WARNING: u is now r*u and p is now e
   MESH_SD {
      size_t s = side_label;
      Element *el = The_side.getNeighborElement(1), 
              *er = The_side.getNeighborElement(2);
      size_t nl = el->n();
      real_t rate = _Lrate(s)*_TimeStep;
      _r->add(nl,-_Fr(s)*rate);
      if ((*_r)(nl)<0.)
         cerr << "WARNING: negative density\nFr = " << _Fr(s)*rate << endl;
      _v->add(nl,1,-_Frv(s,1)*rate);
      _v->add(nl,2,-_Frv(s,2)*rate);
      _p->add(nl,-_FE(s)*rate);
      if (The_side.isOnBoundary()==0) {
         size_t nr = er->n();
         rate = _Rrate(s)*_TimeStep;
         _r->add(nr,_Fr(s)*rate);
         _v->add(nr,1,_Frv(s,1)*rate);
         _v->add(nr,2,_Frv(s,2)*rate);
         _p->add(nr,_FE(s)*rate);
      }
   }

// We now come back to primal variables
   fromConservativeToPrimal();
}


void ICPG2DT::Forward(const Vect<real_t>& Flux,
                      Vect<real_t>&       Field)
{
   MESH_SD {
      size_t ns = side_label;
      Field.add(The_side.getNeighborElement(1)->n(),-_Lrate(ns)*_TimeStep*Flux(ns));
      if (The_side.isOnBoundary()==0)
         Field.add(The_side.getNeighborElement(2)->n(),_Rrate(ns)*_TimeStep*Flux(ns));
   }
}


void ICPG2DT::getMomentum(Vect<real_t>& m) const
{
   m.setMesh(*_theMesh,ELEMENT_DOF,2);
   MESH_EL {
      size_t n = element_label;
      m.set(n,1,(*_r)(n)*(*_v)(n,1));
      m.set(n,2,(*_r)(n)*(*_v)(n,2));
   }
}


void ICPG2DT::getInternalEnergy(Vect<real_t>& e) const
{
   e.setMesh(*_theMesh,ELEMENT_DOF,1);
   MESH_EL
      e.set(element_label,(*_p)(element_label)*(_Gamma-1.));
}


void ICPG2DT::getTotalEnergy(Vect<real_t>& e) const
{
   e.setMesh(*_theMesh,ELEMENT_DOF,1);
   Vect<real_t> tmp_e;
   getInternalEnergy(tmp_e);
   MESH_EL {
      size_t n = element_label;
      e.set(n,tmp_e(n) + 0.5*(*_r)(n)*(*_v)(n,1)*(*_v)(n,1) + (*_v)(n,2)*(*_v)(n,2));
   }
}


void ICPG2DT::getSoundSpeed(Vect<real_t>& s) const
{
   s.setMesh(*_theMesh,ELEMENT_DOF,1);
   MESH_EL {
      size_t n = element_label;
      s.set(n,sqrt(_Gamma*(*_p)(n)/(*_r)(n)));
   }
}


void ICPG2DT::getMach(Vect<real_t>& m) const
{
   m.setMesh(*_theMesh,ELEMENT_DOF,2);
   Vect<real_t> tmp;
   getSoundSpeed(tmp);
   MESH_EL {
      size_t n = element_label;
      m.set(n,1,(*_v)(n,1)/tmp(n));
      m.set(n,2,(*_v)(n,2)/tmp(n));
   }
}


real_t ICPG2DT::runOneTimeStep()
{
   real_t v = getFlux();
   _TimeStep = _ReferenceLength/v*_CFL;
   forward();
   return _TimeStep;
}

} /*  namespace OFELI */
