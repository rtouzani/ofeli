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

                            Implementation of Class ICPG3DT
                
       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/

#include "equations/cl/ICPG3DT.h"
#include "mesh/Mesh.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "linear_algebra/LocalMatrix_impl.h"
#include "linear_algebra/LocalVect_impl.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/Point.h"
#include "shape_functions/Tetra4.h"
#include "shape_functions/Triang3.h"
#include <algorithm>

using std::min;
using std::max;

using std::cout;

namespace OFELI {

ICPG3DT::e_fptr ICPG3DT::fsolver [] = {	
            &ICPG3DT::RiemannSolver_ROE,
            &ICPG3DT::RiemannSolver_VFROE,
            &ICPG3DT::RiemannSolver_LF,
            &ICPG3DT::RiemannSolver_RUSANOV,
            &ICPG3DT::RiemannSolver_HLL,
            &ICPG3DT::RiemannSolver_HLLC
           };


ICPG3DT::ICPG3DT(Mesh& ms) : Muscl3DT(ms)
{
   Init();
   _init_alloc = true;
   _r = new Vect<real_t>(*_theMesh,ELEMENT_DOF,1);
   _v = new Vect<real_t>(*_theMesh,ELEMENT_DOF,3);
   _p = new Vect<real_t>(*_theMesh,ELEMENT_DOF,1);
   _min_density = 1.e-9;
}


ICPG3DT::ICPG3DT(Mesh&         ms,
                 Vect<real_t>& r,
                 Vect<real_t>& v,
                 Vect<real_t>& p)
        : Muscl3DT(ms)
{
   _r = &r;
   _v = &v;
   _p = &p;
   _init_alloc = false;
   _min_density = 1.e-9;
   Init();
}


void ICPG3DT::Init()
{
   _Gr.setMesh(*_theMesh,SIDE_DOF,1);
   _Gv.setMesh(*_theMesh,SIDE_DOF,3);
   _Gp.setMesh(*_theMesh,SIDE_DOF,1);

   _Dr.setMesh(*_theMesh,SIDE_DOF,1);
   _Dv.setMesh(*_theMesh,SIDE_DOF,3);
   _Dp.setMesh(*_theMesh,SIDE_DOF,1);

   _Fr.setMesh(*_theMesh,SIDE_DOF,1);
   _Frv.setMesh(*_theMesh,SIDE_DOF,3);
   _FE.setMesh(*_theMesh,SIDE_DOF,1);

   setGamma(1.4);
   _ReferenceLength = _MinimumEdgeLength;
   setCFL(0.2);
   setTimeStep(0.);
   mySolver = fsolver[0];
}


ICPG3DT::~ICPG3DT()
{
   if (_init_alloc) {
      delete _r;
      delete _v;
      delete _p;
   }
}


void ICPG3DT::setInitialConditionShockTube(const LocalVect<real_t,5>& BcG,
                                           const LocalVect<real_t,5>& BcD,
                                           real_t                     x0)
{
// Two differents zones (for shock tube purpose)
   MESH_EL {
      size_t n = element_label;
//    defining two zones
      if (Tetra4(the_element).getCenter().x<x0){
         (*_r)(n)   = BcG[0];
         (*_v)(n,1) = BcG[1];
         (*_v)(n,2) = BcG[2];
         (*_v)(n,3) = BcG[3];
         (*_p)(n)   = BcG[4];
      }
      else {
         (*_r)(n)   = BcD[0];
         (*_v)(n,1) = BcD[1];
         (*_v)(n,2) = BcD[2];
         (*_v)(n,3) = BcD[3];
         (*_p)(n)   = BcD[4];
      }
   }
}


void ICPG3DT::setInitialCondition(const LocalVect<real_t,5>& u)
{
   MESH_EL {
      size_t n = element_label;
      (*_r)(n)   = u[0];
      (*_v)(n,1) = u[1];
      (*_v)(n,2) = u[2];
      (*_v)(n,3) = u[3];
      (*_p)(n)   = u[4];
   }
}


void ICPG3DT::setReconstruction()
{
    if (Muscl3DT::setReconstruction(*_r,_Gr,_Dr,1)) {
       cout << "ERROR : reconstruction of rho failed" << endl;
       exit(3);
    }
    if (Muscl3DT::setReconstruction(*_v,_Gv,_Dv,1)) {
       cout << "ERROR : reconstruction of u failed" << endl;
       exit(3);
    }
    if (Muscl3DT::setReconstruction(*_v,_Gv,_Dv,2)) {
       cout << "ERROR : reconstruction of v failed" << endl;
       exit(3);
    }
    if (Muscl3DT::setReconstruction(*_v,_Gv,_Dv,3)) {
       cout << "ERROR : reconstruction of w failed" << endl;
       exit(3);
    }
    if (Muscl3DT::setReconstruction(*_p,_Gp,_Dp,1)) {
       cout << "ERROR : reconstruction of p failed" << endl;
       exit(3);
    }
}


/*====== solve Riemann problem for a Ox direction side (Roe method) ===========*/
real_t ICPG3DT::RiemannSolver_ROE(int s)
{
// setSolvers(s);
   real_t lambda_entropy_left, lambda_entropy_right;   // entropy correction
   real_t rho_star, u_star, p_star, c_star;            // entropy correction state
   real_t alpha1, alpha2, alpha5;                      // amplitude of the waves
   real_t sqrt_rho_G=sqrt(_Gr(s)), sqrt_rho_D=sqrt(_Dr(s));
   real_t one_over_sqrt = 1./(sqrt_rho_G+sqrt_rho_D);
   real_t one_over_rho_G=1./_Gr(s), one_over_rho_D=1./_Dr(s);
   sqrt_rho_G *= one_over_sqrt;
   sqrt_rho_D *= one_over_sqrt;

   if (_Gr(s)<_min_density) {
      if (_verbose)
         cerr << "WARNING: left minimum density rho = " << _Gr(s) << endl;
      exit(2);
   }
   if (_Dr(s)<_min_density) {
      if (_verbose)
         cerr << "WARNING: right  minimum density rho = " << _Dr(s) << endl;
      exit(2);
   }
   _tGr = _Gr(s);
   _tDr = _Dr(s);

// compute velocity
   _tGv[0] = _Gv(s,1);
   _tGv[1] = _Gv(s,2);
   _tGv[2] = _Gv(s,3);
   _tDv[0] = _Dv(s,1);
   _tDv[1] = _Dv(s,2);
   _tDv[2] = _Dv(s,3);
   _Grv[0] = _tGr*_tGv[0];
   _Grv[1] = _tGr*_tGv[1];
   _Grv[2] = _tGr*_tGv[2];
   _Drv[0] = _tDr*_tDv[0];
   _Drv[1] = _tDr*_tDv[1];
   _Drv[2] = _tDr*_tDv[2];

// compute pressure, internal energy and sound speed
   _tGp = _Gp(s);
   _GE = _tGp/(_Gamma-1.) + 0.5*_tGr*(_tGv[0]*_tGv[0] + _tGv[1]*_tGv[1] + _tGv[2]*_tGv[2]);
   if (_tGp<0.)
      cout << "WARNING: Gp < 0" << endl;
   _GH = (_GE+_tGp)*one_over_rho_G; 
   _tDp = _Dp(s);
   _DE = _tDp/(_Gamma-1.) + 0.5*_tDr*(_tDv[0]*_tDv[0] + _tDv[1]*_tDv[1] + _tDv[2]*_tDv[2]);
   if (_tDp<0.)
      cout << "WARNING: Dp < 0" << endl;
   _DH = (_DE+_tDp)*one_over_rho_D; 

// compute Roe average
   Point<real_t> Roe_u(sqrt_rho_G*_tGv[0] + sqrt_rho_D*_tDv[0],
                       sqrt_rho_G*_tGv[1] + sqrt_rho_D*_tDv[1],
                       sqrt_rho_G*_tGv[2] + sqrt_rho_D*_tDv[2]);
   real_t Roe_H = sqrt_rho_G*_GH + sqrt_rho_D*_DH;
   real_t Roe_c = (_Gamma-1.)*(Roe_H-0.5*Roe_u.NNorm());

   if (Roe_c<0) {
      if (_verbose) {
         cerr << "WARNING: Roe sound speed negative: " << Roe_c << " at point " << s
              << ", on boundary ? 0/1 " << _theMesh->getPtrSide(s)->isOnBoundary() << endl
              << ", x = " << Triang3(_theMesh->getPtrSide(s)).getCenter().x
              << ", y = " << Triang3(_theMesh->getPtrSide(s)).getCenter().y
              << ", z = " << Triang3(_theMesh->getPtrSide(s)).getCenter().z << endl;
      }
      exit(2);
   }
   Roe_c = sqrt(Roe_c);

// Compute eigenvalues
   real_t lambda1 = Roe_u.x - Roe_c, lambda2 = Roe_u.x, lambda5 = Roe_u.x + Roe_c;

// Compute flux
   if (0.<lambda1) {
      _Fr(s) = _Grv[0];
      _Frv(s,1) = _Grv[0]*_tGv[0] + _tGp;
      _Frv(s,2) = _Grv[1]*_tGv[0];
      _Frv(s,3) = _Grv[2]*_tGv[0];
      _FE(s) = _tGv[0]*(_GE+_tGp);
      return (Roe_u.x+Roe_c);
   }

   alpha2 = (_tDr-_tGr)*(Roe_H-Roe_u.x*Roe_u.x) + (_Drv[1]-_Grv[1])*Roe_u.x - (_DE-_GE)
            + ((_Drv[1]-_Grv[1]) - (_tDr-_tGr)*Roe_u.y)*Roe_u.y
            + ((_Drv[2]-_Grv[2]) - (_tDr-_tGr)*Roe_u.z)*Roe_u.z;
   alpha2 *= (_Gamma-1)/(Roe_c*Roe_c);
   alpha1  = (_tDr-_tGr)*(Roe_c+Roe_u.x)-(_Drv[0]-_Grv[0]) - alpha2*Roe_c;
   alpha1 *= 0.5/Roe_c;

   if (0.<=lambda2) {
//    entropy fix if necessary (Toro  1997 p. 338)
      lambda_entropy_left = _tGv[0] - _Gc;
      rho_star = _tGr + alpha1;
      u_star = (_Grv[0]+alpha1*(Roe_u.x-Roe_c))/rho_star;
      p_star = (_Gamma-1)*(_GE+alpha1*(Roe_H-Roe_u.x*Roe_c)-0.5*rho_star*u_star*u_star);
      if (p_star<0) {
         if (_verbose)
            cerr << "Warning: negative pressure in left entropy fix method." << endl;
         exit(2);
      }
      c_star = sqrt(_Gamma*p_star/rho_star);
      lambda_entropy_right = u_star - c_star;
      if (lambda_entropy_left<0. && 0.<lambda_entropy_right)
         lambda1 = lambda_entropy_left*(lambda_entropy_right-lambda1)/(lambda_entropy_right-lambda_entropy_left);
      _Fr(s) = _Grv[0] + alpha1*lambda1;
      _Frv(s,1) = _Grv[0]*_tGv[0] + _tGp + alpha1*lambda1*(Roe_u.x-Roe_c);
      _Frv(s,2) = _Grv[0]*_tGv[1] + alpha1*lambda1*Roe_u.y;
      _Frv(s,3) = _Grv[0]*_tGv[2] + alpha1*lambda1*Roe_u.z;
      _FE(s) = _tGv[0]*(_GE+_tGp) + alpha1*lambda1*(Roe_H-Roe_u.x*Roe_c);
      return (Roe_u.x+Roe_c);
   }
   alpha5 = (_tDr-_tGr) - alpha2 - alpha1;

   if (0.<lambda5) {
      lambda_entropy_right = _tDv[0] + _Dc;
      rho_star = _tDr - alpha5;
      u_star = (_Drv[0] - alpha5*(Roe_u.x+Roe_c))/rho_star;
      p_star = (_Gamma-1)*(_DE-alpha5*(Roe_H+Roe_u.x*Roe_c)-0.5*rho_star*u_star*u_star);
      if (p_star<0) {
         if (_verbose)
            cerr << "Warning: negative pressure in right entropy fixed method." << endl;
         exit(2);
      }
      c_star = sqrt(_Gamma*p_star/rho_star);
      lambda_entropy_left = u_star + c_star;
      if (lambda_entropy_left<0. && 0.<lambda_entropy_right)
         lambda5 = lambda_entropy_right*(lambda5-lambda_entropy_left)/(lambda_entropy_right-lambda_entropy_left);
      _Fr(s) = _Drv[0] - alpha5*lambda5;
      _Frv(s,1) = _Drv[0]*_tDv[0] + _tDp - alpha5*lambda5*(Roe_u.x+Roe_c);
      _Frv(s,2) = _Drv[0]*_tDv[1] - alpha5*lambda5*Roe_u.y;
      _Frv(s,3) = _Drv[0]*_tDv[2] - alpha5*lambda5*Roe_u.z;
      _FE(s) = _tDv[0]*(_DE+_tDp) - alpha5*lambda5*(Roe_H+Roe_u.x*Roe_c);
      return (-Roe_u.x+Roe_c);
    }

   _Fr(s)    = _Drv[0];
   _Frv(s,1) = _Drv[0]*_tDv[0] + _tDp;
   _Frv(s,2) = _Drv[0]*_tDv[1];
   _Frv(s,3) = _Drv[0]*_tDv[2];
   _FE(s)    = _tDv[0]*(_DE+_tDp);
   return (-Roe_u.x+Roe_c);
}


real_t ICPG3DT::RiemannSolver_VFROE(int s)
{
   real_t lambda_entropy_left, lambda_entropy_right;   // entropy correction
   real_t rho_star, u_star, p_star, c_star;            // entropy correction state
   real_t alpha1, alpha2, alpha5;                      // amplitude of the waves
   real_t sqrt_rho_G=sqrt(_Gr(s)), sqrt_rho_D=sqrt(_Dr(s));
   real_t one_over_sqrt = 1./(sqrt_rho_G+sqrt_rho_D);	// normalization factor of the Roe average
   real_t one_over_rho_G=1./_Gr(s), one_over_rho_D=1./_Dr(s);
   sqrt_rho_G *= one_over_sqrt;                       // normalized coefficients
   sqrt_rho_D *= one_over_sqrt;

   if (_Gr(s)<_min_density) {
      if (_verbose)
         cerr << "WARNING: Left minimum density rho = " << _Gr(s) << endl;
      exit(2);
   }
   if (_Dr(s)<_min_density) {
      if (_verbose)
         cerr << "WARNING: Right minimum density rho = " << _Dr(s) << endl;
      exit(2);
   }
   _tGr = _Gr(s); _tDr = _Dr(s);

// Compute velocity
   _tGv[0] = _Gv(s,1); _tGv[1] = _Gv(s,2); _tGv[2] = _Gv(s,3);
   _tDv[0] = _Dv(s,1); _tDv[1] = _Dv(s,2); _tDv[2] = _Dv(s,3);

   _Grv[0] = _tGr*_tGv[0]; _Grv[1] = _tGr*_tGv[1]; _Grv[2] = _tGr*_tGv[2];
   _Drv[0] = _tDr*_tDv[0]; _Drv[1] = _tDr*_tDv[1]; _Drv[2] = _tDr*_tDv[2];

// Compute pressure, internal energy and sound speed
   _tGp = _Gp(s);
   _GE = _tGp/(_Gamma-1.) + 0.5*_tGr*(_tGv[0]*_tGv[0]+_tGv[1]*_tGv[1]+_tGv[2]*_tGv[2]);
   if (_tGp<0.)
      cout << "WARNING, Gp < 0" << endl;
   _GH = (_GE+_tGp)*one_over_rho_G;

   _tDp = _Dp(s);
   _DE = _tDp/(_Gamma-1.) + 0.5*_tDr*(_tDv[0]*_tDv[0]+_tDv[1]*_tDv[1]+_tDv[2]*_tDv[2]);
   if (_tDp<0.)
      cout << "WARNING: Dp < 0" << endl;
   _DH = (_DE+_tDp)*one_over_rho_D; 

// Compute Roe average
   Point<real_t> Roe_u(0.5*(_tGv[0]+_tDv[0]),0.5*(_tGv[1]+_tDv[1]),0.5*(_tGv[2]+_tDv[2]));
   real_t Roe_H = 0.5*(_GH+_DH);
   real_t Roe_c = (_Gamma-1.)*(Roe_H-0.5*Roe_u.NNorm());

   if (Roe_c<0) {
      if (_verbose)
         cerr << "WARNING: Roe sound speed negative: " << Roe_c << ", at point " << s
		      << " on boundary ? (0/1) " << _theMesh->getPtrSide(s)->isOnBoundary()
			  << ", x = " << Triang3(_theMesh->getPtrSide(s)).getCenter().x
			  << ", y = " << Triang3(_theMesh->getPtrSide(s)).getCenter().y
			  << ", z = " << Triang3(_theMesh->getPtrSide(s)).getCenter().z << endl;
      exit(2);
   }
   Roe_c = sqrt(Roe_c);

// Compute eigenvalues
   real_t lambda1 = Roe_u.x - Roe_c, lambda2 = Roe_u.x, lambda5 = Roe_u.x + Roe_c;

// Compute flux
   if (0.<lambda1) {
      _Fr(s) = _Grv[0];
      _Frv(s,1) = _Grv[0]*_tGv[0] + _tGp;
      _Frv(s,2) = _Grv[1]*_tGv[0];
      _Frv(s,3) = _Grv[2]*_tGv[0];
      _FE(s) = _tGv[0]*(_GE+_tGp);
      return (Roe_u.x+Roe_c);
    }

   alpha2 = (_tDr-_tGr)*(Roe_H-Roe_u.x*Roe_u.x)+(_Drv[1]-_Grv[1])*Roe_u.x - (_DE-_GE)
            + ((_Drv[1]-_Grv[1]) - (_tDr-_tGr)*Roe_u.y)*Roe_u.y
            + ((_Drv[2]-_Grv[2]) - (_tDr-_tGr)*Roe_u.z)*Roe_u.z;
   alpha2 *= (_Gamma-1)/(Roe_c*Roe_c);
   alpha1  = (_tDr-_tGr)*(Roe_c+Roe_u.x)-(_Drv[0]-_Grv[0]) - alpha2*Roe_c;
   alpha1 *= 0.5/Roe_c;

   if (lambda2 >= 0) {
//    entropy fix if necessary (Toro 1997 p. 338)
      lambda_entropy_left = _tGv[0] - _Gc;
      rho_star = _tGr + alpha1;
      u_star = (_Grv[0] + alpha1*(Roe_u.x-Roe_c) )/rho_star;
      p_star = (_Gamma-1)*(_GE+alpha1*(Roe_H-Roe_u.x*Roe_c)-0.5*rho_star*u_star*u_star);
      if (p_star<0) {
         if (_verbose)
            cerr << "Warning: Negative pressure in left entropy fix method." << endl;
         exit(2);
      }
      c_star = sqrt(_Gamma*p_star/rho_star);
      lambda_entropy_right = u_star - c_star;
      if (lambda_entropy_left<0. && 0.<lambda_entropy_right)
         lambda1 = lambda_entropy_left*(lambda_entropy_right-lambda1)/(lambda_entropy_right-lambda_entropy_left);
      _Fr(s)    = _Grv[0] + alpha1*lambda1;
      _Frv(s,1) = _Grv[0]*_tGv[0] + _tGp + alpha1*lambda1*(Roe_u.x-Roe_c);
      _Frv(s,2) = _Grv[0]*_tGv[1] + alpha1*lambda1*Roe_u.y;
      _Frv(s,3) = _Grv[0]*_tGv[2] + alpha1*lambda1*Roe_u.z;
      _FE(s)    = _tGv[0]*(_GE+_tGp) + alpha1*lambda1*(Roe_H-Roe_u.x*Roe_c);
      return (Roe_u.x+Roe_c);
   }
   alpha5 = (_tDr-_tGr) - alpha2 - alpha1;

   if (lambda5>0) {
//    entropy fix if necessary (Toro  1997 p. 338)
      lambda_entropy_right = _tDv[0] + _Dc;
      rho_star = _tDr - alpha5;
      u_star = (_Drv[0]-alpha5*(Roe_u.x+Roe_c))/rho_star;
      p_star = (_Gamma-1)*(_DE-alpha5*(Roe_H+Roe_u.x*Roe_c)-0.5*rho_star*u_star*u_star);
      if (p_star<0) {
         if (_verbose)
            cerr << "Warning: negative pressure in right entropy fixed method." << endl;
         exit(2);
      }
      c_star = sqrt(_Gamma*p_star/rho_star);
      lambda_entropy_left = u_star + c_star;
      if (lambda_entropy_left<0. && 0.<lambda_entropy_right) // we have to fix the entropy
         lambda5 = lambda_entropy_right*(lambda5-lambda_entropy_left)/(lambda_entropy_right-lambda_entropy_left);
	   _Fr(s)    = _Drv[0] - alpha5*lambda5;
	   _Frv(s,1) = _Drv[0]*_tDv[0] + _tDp - alpha5*lambda5*(Roe_u.x+Roe_c);
	   _Frv(s,2) = _Drv[0]*_tDv[1] - alpha5*lambda5*Roe_u.y;
	   _Frv(s,3) = _Drv[0]*_tDv[2] - alpha5*lambda5*Roe_u.z;
	   _FE(s)    = _tDv[0]*(_DE+_tDp) - alpha5*lambda5*(Roe_H+Roe_u.x*Roe_c);
        return (-Roe_u.x+Roe_c);
    }

   _Fr(s) = _Drv[0];
   _Frv(s,1) = _Drv[0]*_tDv[0] + _tDp;
   _Frv(s,2) = _Drv[0]*_tDv[1];
   _Frv(s,3) = _Drv[0]*_tDv[2];
   _FE(s) = _tDv[0]*(_DE+_tDp);
   return (-Roe_u.x+Roe_c);
}


real_t ICPG3DT::RiemannSolver_LF(int s)
{
   real_t one_over_rho_G = 1./_Gr(s), one_over_rho_D = 1./_Dr(s);
   if (_Gr(s)<_min_density) {
      if (_verbose)
         cerr << "WARNING: Left minimum density rho = " << _Gr(s) << endl;
      exit(2);
   }
   if (_Dr(s)<_min_density) {
      if (_verbose)
         cerr << "WARNING: Right minimum density rho = " << _Dr(s) << endl;
      exit(2);
   }
   _tGr = _Gr(s);
   _tDr = _Dr(s);

// Compute velocities
   _tGv[0] = _Gv(s,1); _tGv[1] = _Gv(s,2); _tGv[2] = _Gv(s,3);
   _tDv[0] = _Dv(s,1); _tDv[1] = _Dv(s,2); _tDv[2] = _Dv(s,3);
   _Grv[0] = _tGr*_tGv[0]; _Grv[1] = _tGr*_tGv[1]; _Grv[2] = _tGr*_tGv[2];
   _Drv[0] = _tDr*_tDv[0]; _Drv[1] = _tDr*_tDv[1]; _Drv[2] = _tDr*_tDv[2];

// Compute pressure, internal energy and sound speed
   _tGp = _Gp(s);
   _GE = _tGp/(_Gamma-1.) + 0.5*_tGr*(_tGv[0]*_tGv[0]+_tGv[1]*_tGv[1]+_tGv[2]*_tGv[2]);
   if (_tGp<0.)
      cerr << "WARNING: Gp < 0" << endl;
   _GH = (_GE+_tGp)*one_over_rho_G; 

   _tDp = _Dp(s);
   _DE = _tDp/(_Gamma-1.) + 0.5*_tDr*(_tDv[0]*_tDv[0]+_tDv[1]*_tDv[1]+_tDv[2]*_tDv[2]);
   if (_tDp<0.)
      cerr << "WARNING: DP<0" << endl;
   _DH = (_DE+_tDp)*one_over_rho_D;

// Use VFRoe method to compute eigenvalues
// Compute VFRoe average
   Point<real_t> Roe_u(0.5*(_tGv[0]+_tDv[0]),0.5*(_tGv[1]+_tDv[1]),0.5*(_tGv[2]+_tDv[2]));
   real_t Roe_H = 0.5*(_GH+_DH);
   real_t Roe_c = (_Gamma-1.)*(Roe_H-0.5*Roe_u.NNorm());

   if (Roe_c<0) {
      if (_verbose)
         cerr << "WARNING: Roe sound speed negative: " << Roe_c << ", at point " << s
              << " on boundary ? (0/1) " << _theMesh->getPtrSide(s)->isOnBoundary()
              << ", x = " << Triang3(_theMesh->getPtrSide(s)).getCenter().x
              << ", y = " << Triang3(_theMesh->getPtrSide(s)).getCenter().y
              << ", z = " << Triang3(_theMesh->getPtrSide(s)).getCenter().z << endl;
      exit(2);
   }
   Roe_c = sqrt(Roe_c);

// compute LF flux
   _Fr(s)    = 0.5*(_Grv[0] + _Drv[0]);
   _Frv(s,1) = 0.5*(_Grv[0]*_tGv[0] + _tGp + _Drv[0]*_tDv[0] + _tDp);
   _Frv(s,2) = 0.5*(_Grv[0]*_tGv[1] + _Drv[0]*_tDv[1]);
   _Frv(s,3) = 0.5*(_Grv[0]*_tGv[2] + _Drv[0]*_tDv[2]);
   _FE(s)    = 0.5*(_tGv[0]*(_GE+_tGp) + _tDv[0]*(_DE+_tDp));

   if (_TimeStep) {
      _Fr(s)    -= 0.5*_ReferenceLength/_TimeStep*(_tDr-_tGr);
      _Frv(s,1) -= 0.5*_ReferenceLength/_TimeStep*(_Drv[0]-_Grv[0]);
      _Frv(s,2) -= 0.5*_ReferenceLength/_TimeStep*(_Drv[1]-_Grv[1]);
      _Frv(s,3) -= 0.5*_ReferenceLength/_TimeStep*(_Drv[2]-_Grv[2]);
      _FE(s)    -= 0.5*_ReferenceLength/_TimeStep*(_DE-_GE);
   }
   return max(fabs(Roe_u.x-Roe_c),fabs(Roe_u.x+Roe_c));
}


real_t ICPG3DT::RiemannSolver_RUSANOV(int s)
{
   real_t one_over_rho_G=1./_Gr(s), one_over_rho_D=1./_Dr(s), Splus;
   if (_Gr(s)<_min_density) {
      if (_verbose)
         cerr << "WARNING: Left minimum density rho = " << _Gr(s) << endl;
      exit(2);
   }
   if (_Dr(s)<_min_density) {
      if (_verbose)
         cerr << "WARNING: Right minimum density rho = " << _Dr(s) << endl;
      exit(2);
   }
   _tGr = _Gr(s);
   _tDr = _Dr(s);

// Compute velocitiy
   _tGv[0] = _Gv(s,1); _tGv[1] = _Gv(s,2); _tGv[2] = _Gv(s,3);
   _tDv[0] = _Dv(s,1); _tDv[1] = _Dv(s,2); _tDv[2] = _Dv(s,3);

   _Grv[0] = _tGr*_tGv[0]; _Grv[1] = _tGr*_tGv[1]; _Grv[2] = _tGr*_tGv[2];
   _Drv[0] = _tDr*_tDv[0]; _Drv[1] = _tDr*_tDv[1]; _Drv[2] = _tDr*_tDv[2];

// Compute pressure, internal energy and sound speed
   _tGp = _Gp(s);
   _GE = _tGp/(_Gamma-1.) + 0.5*_tGr*(_tGv[0]*_tGv[0]+_tGv[1]*_tGv[1]+_tGv[2]*_tGv[2]);
   if (_tGp<0.)
      cerr << "WARNING: Gp < 0" << endl;
   _GH = (_GE+_tGp)*one_over_rho_G; 

   _tDp = _Dp(s);
   _DE = _tDp/(_Gamma-1.) + 0.5*_tDr*(_tDv[0]*_tDv[0]+_tDv[1]*_tDv[1]+_tDv[2]*_tDv[2]);
   if (_tDp<0.)
      cerr << "WARNING: Dp < 0" << endl;
   _DH = (_DE+_tDp)*one_over_rho_D;

// Compute S+
   _Gc = (_Gamma-1.)*(_GH-(_tGv[0]*_tGv[0]+_tGv[1]*_tGv[1]+_tGv[2]*_tGv[2])*0.5); // Here c*c
   _Dc = (_Gamma-1.)*(_DH-(_tDv[0]*_tDv[0]+_tDv[1]*_tDv[1]+_tDv[2]*_tDv[2])*0.5); // Here c*c

   if (_Gc<0. || _Dc<0.) {
      if (_verbose)
         cerr << "WARNING: Sound speed negative: " << _Dc << " " << _Gc << ", at point " << s
              << " on boundary ? (0/1) " << _theMesh->getPtrSide(s)->isOnBoundary()
              << ", x = " << Triang3(_theMesh->getPtrSide(s)).getCenter().x
              << ", y = " << Triang3(_theMesh->getPtrSide(s)).getCenter().y
              << ", z = " << Triang3(_theMesh->getPtrSide(s)).getCenter().z << endl;
      exit(2);
   }
   _Gc = sqrt(_Gc); _Dc = sqrt(_Dc);
   Splus = max(fabs(_tGv[0]-_Gc),fabs(_tGv[0]+_Gc));
   Splus = max(Splus,fabs(_tDv[0]-_Dc));
   Splus = max(Splus,fabs(_tDv[0]+_Dc));

//  Compute Rusanov flux
   _Fr(s)    = 0.5*(_Grv[0] + _Drv[0]);
   _Frv(s,1) = 0.5*(_Grv[0]*_tGv[0] + _tGp + _Drv[0]*_tDv[0] + _tDp);
   _Frv(s,2) = 0.5*(_Grv[0]*_tGv[1] + _Drv[0]*_tDv[1]);
   _Frv(s,3) = 0.5*(_Grv[0]*_tGv[2] + _Drv[0]*_tDv[2]);
   _FE(s)    = 0.5*(_tGv[0]*(_GE+_tGp) + _tDv[0]*(_DE+_tDp));

    _Fr(s)    -= 0.5*Splus*(_tDr-_tGr);
    _Frv(s,1) -= 0.5*Splus*(_Drv[0]-_Grv[0]);
    _Frv(s,2) -= 0.5*Splus*(_Drv[1]-_Grv[1]);
    _Frv(s,3) -= 0.5*Splus*(_Drv[2]-_Grv[2]);
    _FE(s)    -= 0.5*Splus*(_DE-_GE);
    return Splus;
}


real_t ICPG3DT::RiemannSolver_HLL(int s)
{
   real_t one_over_rho_G=1./_Gr(s), one_over_rho_D=1./_Dr(s);
   real_t SL, SR, Splus, C1, C2, C3;

   if (_Gr(s)<_min_density) {
      if (_verbose)
         cerr << "WARNING: Left minimum density rho = " << _Gr(s) << endl;
      _Gr(s) = _min_density;
   }
   if (_Dr(s)<_min_density) {
      if (_verbose)
         cerr << "WARNING: Right minimum density rho = " << _Dr(s) << endl;
      _Dr(s) = _min_density;
   }
   _tGr = _Gr(s); _tDr = _Dr(s);

// Compute velocity
   _tGv[0] = _Gv(s,1); _tGv[1] = _Gv(s,2); _tGv[2] = _Gv(s,3);
   _tDv[0] = _Dv(s,1); _tDv[1] = _Dv(s,2); _tDv[2] = _Dv(s,3);

   _Grv[0] = _tGr*_tGv[0]; _Grv[1] = _tGr*_tGv[1]; _Grv[2] = _tGr*_tGv[2];
   _Drv[0] = _tDr*_tDv[0]; _Drv[1] = _tDr*_tDv[1]; _Drv[2] = _tDr*_tDv[2];

// Compute pressure, internal energy and sound speed
   _tGp = _Gp(s);
   _GE = _tGp/(_Gamma-1.) + 0.5*_tGr*(_tGv[0]*_tGv[0]+_tGv[1]*_tGv[1]+_tGv[2]*_tGv[2]);
   if (_tGp<0.)
      cerr << "WARNING: Gp < 0" << endl;
   _GH = (_GE+_tGp)*one_over_rho_G; 

   _tDp = _Dp(s);
   _DE = _tDp/(_Gamma-1.) + 0.5*_tDr*(_tDv[0]*_tDv[0]+_tDv[1]*_tDv[1]+_tDv[2]*_tDv[2]);
   if (_tDp<0.)
      cerr << "WARNING: Dp < 0" << endl;
   _DH = (_DE+_tDp)*one_over_rho_D;

// compute S+
   _Gc = (_Gamma-1.)*(_GH-(_tGv[0]*_tGv[0]+_tGv[1]*_tGv[1]+_tGv[2]*_tGv[2])*0.5); // Here c*c
   _Dc = (_Gamma-1.)*(_DH-(_tDv[0]*_tDv[0]+_tDv[1]*_tDv[1]+_tDv[2]*_tDv[2])*0.5); // Here c*c

   if (_Gc<0. || _Dc<0.) {
      if (_verbose)
         cerr << "WARNING: Sound speed negative: " << _Dc << " " << _Gc << ", at point " << s
              << " on boundary ? (0/1) " << _theMesh->getPtrSide(s)->isOnBoundary()
              << ", x = " << Triang3(_theMesh->getPtrSide(s)).getCenter().x
              << ", y = " << Triang3(_theMesh->getPtrSide(s)).getCenter().y
              << ", z = " << Triang3(_theMesh->getPtrSide(s)).getCenter().z << endl;
      exit(2);
   }
   _Gc = sqrt(_Gc);
   _Dc = sqrt(_Dc);

   Splus = max(fabs(_tGv[0] - _Gc),fabs(_tGv[0] + _Gc));
   Splus = max(Splus,fabs(_tDv[0] - _Dc));
   Splus = max(Splus,fabs(_tDv[0] + _Dc));
   SL = -Splus;
   SR = Splus;
   C1 =  SR/(SR-SL);
   C2 = -SL/(SR-SL);
   C3 = SL*SR/(SR-SL);
// See other choices of SR, SL

// compute HLL flux
   if (SL>=0.) { // Flux = FG
      _Fr(s)    = _Grv[0];
      _Frv(s,1) = _Grv[0]*_tGv[0] + _tGp;
      _Frv(s,2) = _Grv[0]*_tGv[1];
      _Frv(s,3) = _Grv[0]*_tGv[2] ;
      _FE(s)    = _tGv[0]*(_GE+_tGp);
   }
   else if ( (SL<=0.) && (SR >=0.)) { // Flux = FHLL
      _Fr(s)    = C1*_Grv[0]                + C2*_Drv[0]                + C3*(_tDr-_tGr);
      _Frv(s,1) = C1*(_Grv[0]*_tGv[0]+_tGp) + C2*(_Drv[0]*_tDv[0]+_tDp) + C3*(_Drv[0]-_Grv[0]);
      _Frv(s,2) = C1*(_Grv[0]*_tGv[1])      + C2*(_Drv[0]*_tDv[1])      + C3*(_Drv[1]-_Grv[1]);
      _Frv(s,3) = C1*(_Grv[0]*_tGv[2])      + C2*(_Drv[0]*_tDv[2])      + C3*(_Drv[2]-_Grv[2]);
      _FE(s)    = C1*(_tGv[0]*(_GE+_tGp))   + C2*(_tDv[0]*(_DE+_tDp))   + C3*(_DE-_GE);
   }
   else { // Flux = FD
      _Fr(s)    = _Drv[0] ;
      _Frv(s,1) = _Drv[0]*_tDv[0] + _tDp;
      _Frv(s,2) = _Drv[0]*_tDv[1];
      _Frv(s,3) = _Drv[0]*_tDv[2];
      _FE(s)    = _tDv[0]*(_DE+_tDp);
   }
   return Splus;
}


real_t ICPG3DT::RiemannSolver_HLLC(int s)
{
   real_t one_over_rho_G=1./_Gr(s), one_over_rho_D=1./_Dr(s);
   real_t SL, SR, Splus, Sstar;
   real_t rhoGstar, rhoDstar, uGstar, uDstar, vGstar, vDstar, wGstar, wDstar, EGstar, EDstar;
   real_t AL, AR;

   if (_Gr(s)<_min_density) {
      if (_verbose)
         cerr << "WARNING: Left minimum density rho = " << _Gr(s) << endl;
      _Gr(s) = _min_density;
   }
   if (_Dr(s)<_min_density) {
      if (_verbose)
         cerr << "WARNING: Right minimum density rho = " << _Dr(s) << endl;
      _Dr(s) = _min_density;
   }
   _tGr = _Gr(s); _tDr = _Dr(s);

// Compute velocity
   _tGv[0] = _Gv(s,1); _tGv[1] = _Gv(s,2); _tGv[2] = _Gv(s,3);
   _tDv[0] = _Dv(s,1); _tDv[1] = _Dv(s,2); _tDv[2] = _Dv(s,3);

   _Grv[0] = _tGr*_tGv[0]; _Grv[1] = _tGr*_tGv[1]; _Grv[2] = _tGr*_tGv[2];
   _Drv[0] = _tDr*_tDv[0]; _Drv[1] = _tDr*_tDv[1]; _Drv[2] = _tDr*_tDv[2];

// Compute pressure, internal energy and sound speed
   _tGp = _Gp(s);
   _GE = _tGp/(_Gamma-1.) + 0.5*_tGr*(_tGv[0]*_tGv[0]+_tGv[1]*_tGv[1]+_tGv[2]*_tGv[2]);
   _GH = (_GE+_tGp)*one_over_rho_G; 

   _tDp = _Dp(s);
   _DE = _tDp/(_Gamma-1.) + 0.5*_tDr*(_tDv[0]*_tDv[0]+_tDv[1]*_tDv[1]+_tDv[2]*_tDv[2]);
   _DH = (_DE+_tDp)*one_over_rho_D;

// compute S+
   _Gc = (_Gamma-1.)*(_GH-(_tGv[0]*_tGv[0]+_tGv[1]*_tGv[1]+_tGv[2]*_tGv[2])*0.5); // Here c*c
   _Dc = (_Gamma-1.)*(_DH-(_tDv[0]*_tDv[0]+_tDv[1]*_tDv[1]+_tDv[2]*_tDv[2])*0.5); // Here c*c
   if (_Gc<0. || _Dc<0.) {
      if (_verbose)
         cerr << "WARNING: Sound speed negative: " << _Dc << " " << _Gc << ", at point " << s
              << " on boundary ? (0/1) " << _theMesh->getPtrSide(s)->isOnBoundary()
              << ", x = " << Triang3(_theMesh->getPtrSide(s)).getCenter().x
              << ", y = " << Triang3(_theMesh->getPtrSide(s)).getCenter().y
              << ", z = " << Triang3(_theMesh->getPtrSide(s)).getCenter().z << endl;
      exit(2);
   }
   _Gc = sqrt(_Gc); _Dc = sqrt(_Dc);

   Splus = max(fabs(_tGv[0] - _Gc),fabs(_tGv[0] + _Gc));
   Splus = max(Splus,fabs(_tDv[0] - _Dc));
   Splus = max(Splus,fabs(_tDv[0] + _Dc));
   SL = -Splus; SR = Splus;
   Sstar = (_tDp - _tGp + _Grv[0]*(SL-_tGv[0]) - _Drv[0]*(SR-_tDv[0]))/(_tGr*(SL-_tGv[0]) - _tDr*(SR-_tDv[0])) ;
// See other choices of SR, SL

    AL = _tGr*(SL - _tGv[0])/(SL - Sstar);
    AR = _tDr*(SR - _tDv[0])/(SR - Sstar);

// compute star region
   rhoGstar = AL; rhoDstar = AR;
   uGstar = AL*Sstar; uDstar = AR*Sstar;
   vGstar = AL*_tGv[1]; vDstar = AR*_tDv[1];
   wGstar = AL*_tGv[2]; wDstar = AR*_tDv[2];
   EGstar = AL*(_GE*one_over_rho_G + (Sstar-_tGv[0])*(Sstar + _tGp/(_tGr*(SL-_tGv[0]))));
   EDstar = AR*(_DE*one_over_rho_D + (Sstar-_tDv[0])*(Sstar + _tDp/(_tDr*(SR-_tDv[0]))));

// Compute HLLC flux
   if (SL>=0.) { // Flux = FG
      _Fr(s)    = _Grv[0];
      _Frv(s,1) = _Grv[0]*_tGv[0] + _tGp;
      _Frv(s,2) = _Grv[0]*_tGv[1];
      _Frv(s,3) = _Grv[0]*_tGv[2] ;
      _FE(s)    = _tGv[0]*(_GE+_tGp);
   }
   else if (Sstar >=0.) { // Flux = FHLLCG
      _Fr(s)    = _Grv[0] + SL*(rhoGstar-_Gr(s)) ;
      _Frv(s,1) = _Grv[0]*_tGv[0] + _tGp + SL*(uGstar-_Grv[0]);
      _Frv(s,2) = _Grv[0]*_tGv[1] + SL*(vGstar-_Grv[1]);
      _Frv(s,3) = _Grv[0]*_tGv[2] + SL*(wGstar-_Grv[2]);
      _FE(s)    = _tGv[0]*(_GE+_tGp) + SL*(EGstar-_GE);
   }
   else if (SR >=0.) { // Flux = FHLLCD
      _Fr(s)    = _Drv[0]+ SR*(rhoDstar-_Dr(s));
      _Frv(s,1) = _Drv[0]*_tDv[0]+ _tDp + SR*(uDstar-_Drv[0]);
      _Frv(s,2) = _Drv[0]*_tDv[1] + SR*(vDstar-_Drv[1]);
      _Frv(s,3) = _Drv[0]*_tDv[2] + SR*(wDstar-_Drv[2]);
      _FE(s)    = _tDv[0]*(_DE+_tDp) + SR*(EDstar-_DE);
    }
   else { // Flux = FD
      _Fr(s)    = _Drv[0];
      _Frv(s,1) = _Drv[0]*_tDv[0] + _tDp;
      _Frv(s,2) = _Drv[0]*_tDv[1];
      _Frv(s,3) = _Drv[0]*_tDv[2];
      _FE(s)    = _tDv[0]*(_DE+_tDp);
   }
   return Splus;
}


real_t ICPG3DT::getFlux()
{
    /////////////////////////////////////////////////////////////
    // !!! the solver assumes that boundary condition     !!!
    // !!! i.e. DU at boundary, has already been defined  !!!
    /////////////////////////////////////////////////////////////

   real_t lambda=1.e-10;
   real_t Gu, Gv, Gw, Du, Dv, Dw, Fru, Frv, Frw, det;
   Point<real_t> A, M, t, r, n, B, Pt;
   LocalMatrix<real_t,3,3> b;

   MESH_SD {
      size_t s = side_label;
      A = The_side(1)->getCoord();
      B = The_side(2)->getCoord();
      t = A - B; // first tangent
      t.Normalize();
      n = _n(s);
      r = CrossProduct(n,t); // second tangent
      r.Normalize();
      det = 1./(n.x*t.y*r.z-n.x*r.y*t.z-n.y*t.x*r.z+n.y*r.x*t.z+n.z*t.x*r.y-n.z*r.x*t.y);
      b(1,1) =  (t.y*r.z-r.y*t.z)*det;
      b(2,1) = -(t.x*r.z-r.x*t.z)*det;
      b(3,1) =  (t.x*r.y-r.x*t.y)*det;
      b(1,2) = -(n.y*r.z-r.y*n.z)*det;
      b(2,2) =  (n.x*r.z-r.x*n.z)*det;
      b(3,2) = -(n.x*r.y-r.x*n.y)*det;
      b(1,3) =  (n.y*t.z-t.y*n.z)*det;
      b(2,3) = -(n.x*t.z-t.x*n.z)*det;
      b(3,3) =  (n.x*t.y-t.x*n.y)*det;

//    Save the value
      Gu = _Gv(s,1); Gv = _Gv(s,2); Gw = _Gv(s,3);
      Du = _Dv(s,1); Dv = _Dv(s,2); Dw = _Dv(s,3);

//    Change the orientation
      _Gv(s,1) = b(1,1)*Gu + b(2,1)*Gv + b(3,1)*Gw;
      _Gv(s,2) = b(1,2)*Gu + b(2,2)*Gv + b(3,2)*Gw;
      _Gv(s,3) = b(1,3)*Gu + b(2,3)*Gv + b(3,3)*Gw;
      _Dv(s,1) = b(1,1)*Du + b(2,1)*Dv + b(3,1)*Dw;
      _Dv(s,2) = b(1,2)*Du + b(2,2)*Dv + b(3,2)*Dw;
      _Dv(s,3) = b(1,3)*Du + b(2,3)*Dv + b(3,3)*Dw;

//    Compute the flux with the new orientation
      assert(mySolver);       // just in case
      lambda = max((this->*mySolver)(int(s)),lambda);

//    come back to the initial orientation
      Fru = _Frv(s,1);
      Frv = _Frv(s,2);
      Frw = _Frv(s,3);
      _Frv(s,1) = n.x*Fru + t.x*Frv + r.x*Frw;
      _Frv(s,2) = n.y*Fru + t.y*Frv + r.y*Frw;
      _Frv(s,3) = n.z*Fru + t.z*Frv + r.z*Frw;

//    RESTORE THE VALUES // really usefull????
      _Gv(s,1) = Gu; _Gv(s,2) = Gv; _Gv(s,3) = Gw;
      _Dv(s,1) = Du; _Dv(s,2) = Dv; _Dv(s,3) = Dw;
   }
   return lambda;
}


void ICPG3DT::setBC(const Side& sd,
                    real_t      u)
{
   size_t s = sd.n();
   _Dr(s) = _Gr(s);
   real_t tmp  = _Gv(s,1)*_n(s).x + _Gv(s,2)*_n(s).y + _Gv(s,3)*_n(s).z;
   tmp *= u - 1;
   _Dv(s,1) = _Gv(s,1) + tmp*_n(s).x;
   _Dv(s,2) = _Gv(s,2) + tmp*_n(s).y;
   _Dv(s,3) = _Gv(s,3) + tmp*_n(s).z;
   _Dp(s) = _Gp(s);
}


void ICPG3DT::setBC(int    code,
                    real_t u)
{
   MESH_SD
      if (The_side.getCode()==code)
         setBC(The_side,u);
}


void ICPG3DT::setBC(real_t u)
{
   MESH_SD
      setBC(The_side,u);
}


void ICPG3DT::setBC(const Side&                sd,
                    const LocalVect<real_t,5>& u)
{
   size_t s = sd.n();
   _Dr(s)   = u[0];
   _Dv(s,1) = u[1];
   _Dv(s,2) = u[2];
   _Dv(s,3) = u[3];
   _Dp(s)   = u[4];
}


void ICPG3DT::setBC(int                        code,
                    const LocalVect<real_t,5>& u)
{
   MESH_SD
      if (The_side.getCode()==code)
         setBC(The_side,u);
}


void ICPG3DT::setBC(const LocalVect<real_t,5>& u)
{
   MESH_SD
      setBC(The_side,u);
}


void ICPG3DT::setInOutFlowBC(const Side&                sd,
                             const LocalVect<real_t,5>& u)
{
   size_t s = sd.n();
   real_t Gc = sqrt(_Gamma*_Gp(s)/_Gr(s));
   real_t GUn = _Gv(s,1)*_n(s).x + _Gv(s,2)*_n(s).y + _Gv(s,3)*_n(s).z;
   int choice = (GUn>=0.) + 2*(Gc>=1.);

   switch (choice) {

      case 0:                          // GU.n < 0 and Gc < 1  : subsonic inflow
         _Dr(s )  = _Gr(s);             // getting from interior
         _Dv(s,1) = u[1];              // because exterior normal
         _Dv(s,2) = u[2];
         _Dv(s,3) = u[3];
         _Dp(s)   = u[4];
         break;

      case 1:                          // GU.n >=0 and Gc < 1  : subsonic outflow
         _Dr(s)   = _Gr(s);
         _Dv(s,1) = _Gv(s,1);
         _Dv(s,2) = _Gv(s,2);
         _Dv(s,3) = _Gv(s,3);
         _Dp(s)   = u[4];              // imposing pressure
         break;

      case 2:                          // GU.n < 0 and Gc > 1  : supersonic inflow
         _Dr(s)   = u[0];              // imposing density
         _Dv(s,1) = u[1];              // because exterior normal
         _Dv(s,2) = u[2];
         _Dv(s,3) = u[3];
         _Dp(s)   = u[4];
         break;

      case 3:                          // GU.n > 0 and Gc > 1 : supersonic outflow
         _Dr(s)   = _Gr(s);
         _Dv(s,1) = _Gv(s,1);
         _Dv(s,2) = _Gv(s,2);
         _Dv(s,3) = _Gv(s,3);
         _Dp(s)   = _Gp(s);              // getting from interior
         break;

      default :
         cout << "ERROR in setInOutFlowBC" << endl;
         exit(1);
   }
}


void ICPG3DT::setInOutFlowBC(int                        code,
                             const LocalVect<real_t,5>& u)
{
   MESH_SD
      if (The_side.getCode()==code)
         setInOutFlowBC(The_side,u);
}


void ICPG3DT::setInOutFlowBC(const LocalVect<real_t,5>& u)
{
   MESH_SD
      setInOutFlowBC(The_side,u);
}

/*=====================================================
  SECTION ONE STEP METHOD TO EVALUATE VECTOR A TIME N+1
  Computation of the primal variable    n->n+1
  =====================================================*/
/*----------------------------------------------------------------------------
   The variable flux contains all fluxes issued from the Riemann solver with
   as left element getNeighborElement(1) and right element getNeighborElement(2).
   If getNeighborElement(2) is NULL, we are on the boundary and we impose
   the symmetry condition
  ----------------------------------------------------------------------------*/

void ICPG3DT::fromPrimalToConservative()
{
   real_t G = 1./(_Gamma-1.);
   MESH_EL {
      size_t n = element_label;
      real_t rho = (*_r)(n);
      (*_p)(n) = (*_p)(n)*G + 0.5*rho*((*_v)(n,1)*(*_v)(n,1) +
                                       (*_v)(n,2)*(*_v)(n,2) +
                                       (*_v)(n,3)*(*_v)(n,3));
      (*_v)(n,1) *= rho;
      (*_v)(n,2) *= rho;
      (*_v)(n,3) *= rho;
   }
}


void ICPG3DT::fromConservativeToPrimal()
{
   real_t G = _Gamma-1.;
   MESH_EL {
      size_t n = element_label;
      real_t invrho = 1./(*_r)(n);
      (*_v)(n,1) *= invrho;
      (*_v)(n,2) *= invrho;
      (*_v)(n,3) *= invrho;
      (*_p)(n) = ((*_p)(n) - 0.5*(*_r)(n)*((*_v)(n,1)*(*_v)(n,1) +
                                           (*_v)(n,2)*(*_v)(n,2) +
                                           (*_v)(n,3)*(*_v)(n,3)))*G;
   }
}


void ICPG3DT::forward()
{
   Element *elg, *eld;
   size_t ns, ng, nd;
   real_t rate;

   fromPrimalToConservative();
// Be careful, v is now r*v and p is now E
   MESH_SD {
      ns = side_label;
      elg = The_side.getNeighborElement(1);
      ng = elg->n();
      rate = _Lrate(ns)*_TimeStep;
      (*_r)(ng)   -= _Fr(ns)*rate;
      (*_v)(ng,1) -= _Frv(ns,1)*rate;
      (*_v)(ng,2) -= _Frv(ns,2)*rate;
      (*_v)(ng,3) -= _Frv(ns,3)*rate;
      (*_p)(ng)   -= _FE(ns) *rate;
      if (!The_side.isOnBoundary()) {
         eld = theSide->getNeighborElement(2);
         nd = eld->n();
         rate = _Rrate(ns)*_TimeStep;
         (*_r)(nd)   += _Fr(ns)*rate;
         (*_v)(nd,1) += _Frv(ns,1)*rate;
         (*_v)(nd,2) += _Frv(ns,2)*rate;
         (*_v)(nd,3) += _Frv(ns,3)*rate;
         (*_p)(nd)   += _FE(ns,1) *rate;
      }
   }

// Back now to primal variables
   fromConservativeToPrimal();
}


void ICPG3DT::Forward(const Vect<real_t>& flux,
                      Vect<real_t>&       field)
{
   Element *elg, *eld;
   size_t ns, ng, nd;
   real_t rate;
   MESH_SD {
      ns = side_label;
      elg = The_side.getNeighborElement(1);
      ng = elg->n();
      rate = _Lrate(ns)*_TimeStep;
      field(ng,1) -= flux(ns)*rate;
      if (!The_side.isOnBoundary()) {
         eld = The_side.getNeighborElement(2);
         nd = eld->n();
         rate = _Rrate(ns)*_TimeStep;
         field(nd) += flux(ns)*rate;
      }
   }
}


void ICPG3DT::getMomentum(Vect<real_t>& m) const
{
   m.setMesh(*_theMesh,ELEMENT_DOF,3);
   MESH_EL {
      size_t n = element_label;
      m(n,1) = (*_r)(n)*(*_v)(n,1);
      m(n,2) = (*_r)(n)*(*_v)(n,2);
      m(n,3) = (*_r)(n)*(*_v)(n,3);
   }
}


void ICPG3DT::getInternalEnergy(Vect<real_t>& e) const
{
   e.setMesh(*_theMesh,ELEMENT_DOF,1);
   MESH_EL
      e(element_label) = (*_p)(element_label)*(_Gamma -1.);
}


void ICPG3DT::getTotalEnergy(Vect<real_t>& e) const
{
   e.setMesh(*_theMesh,ELEMENT_DOF,1);
   Vect<real_t> tmp;
   getInternalEnergy(tmp);
   MESH_EL {
      size_t n = element_label;
      e(n) = tmp(n) + 0.5*(*_r)(n)*(*_v)(n,1)*(*_v)(n,1) +
                                   (*_v)(n,2)*(*_v)(n,2) +
                                   (*_v)(n,3)*(*_v)(n,3);
   }
}


void ICPG3DT::getSoundSpeed(Vect<real_t>& s) const
{
   s.setMesh(*_theMesh,ELEMENT_DOF,1);
   MESH_EL {
      size_t n = element_label;
      s(n) = sqrt(_Gamma*(*_p)(n)/(*_r)(n));
   }
}


void ICPG3DT::getMach(Vect<real_t>& m) const
{
   m.setMesh(*_theMesh,ELEMENT_DOF,3);
   Vect<real_t> tmp;
   getSoundSpeed(tmp);
   MESH_EL {
      size_t n = element_label;
      real_t inv = 1./tmp(n);
      m(n,1) = (*_v)(n,1)*inv;
      m(n,2) = (*_v)(n,2)*inv;
      m(n,3) = (*_v)(n,3)*inv;
   }
}


real_t ICPG3DT::runOneTimeStep()
{
   _TimeStep = _ReferenceLength/getFlux()*_CFL;
   forward();
   return _TimeStep;
}

} /* namespace OFELI */
