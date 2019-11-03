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

                            Implementation of Class ICPG1D

       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/

#include <iostream>
#include <algorithm>

using std::min;
using std::max;
using std::cout;

#include "equations/cl/ICPG1D.h"
#include "mesh/Mesh.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "linear_algebra/LocalVect_impl.h"
#include "linear_algebra/Vect_impl.h"
#include "shape_functions/Triang3.h"

namespace OFELI {

ICPG1D::e_fptr ICPG1D::fsolver [] = {
        &ICPG1D::RiemannSolver_ROE,
        &ICPG1D::RiemannSolver_VFROE,
        &ICPG1D::RiemannSolver_LF,
        &ICPG1D::RiemannSolver_RUSANOV,
        &ICPG1D::RiemannSolver_HLL,
        &ICPG1D::RiemannSolver_HLLC
       };


ICPG1D::ICPG1D(Mesh &ms) : Muscl1D(ms), _min_density(1.e-9)
{
   init();
   _init_alloc = true;
   _r = new Vect<real_t>(*_theMesh,1,ELEMENT_DOF);
   _v = new Vect<real_t>(*_theMesh,1,ELEMENT_DOF);
   _p = new Vect<real_t>(*_theMesh,1,ELEMENT_DOF);
   _min_density = 1.e-9;
}


ICPG1D::ICPG1D(Mesh&         ms,
               Vect<real_t>& r,
               Vect<real_t>& v,
               Vect<real_t>& p)
      : Muscl1D(ms) 
{
   _r = &r;
   _v = &v;
   _p = &p;
   _init_alloc = false;
   _min_density = 1.e-9;
}


void ICPG1D::init()
{
   _Gr.setSize(_nb_sides);
   _Gv.setSize(_nb_sides);
   _Gp.setSize(_nb_sides);
   _Dr.setSize(_nb_sides);
   _Dv.setSize(_nb_sides);
   _Dp.setSize(_nb_sides);
   _Fr.setSize(_nb_sides);
   _FrU.setSize(_nb_sides);
   _FE.setSize(_nb_sides);
   setGamma(1.4);
   _ReferenceLength = _MinimumLength;
   setCFL(0.2);
   setTimeStep(0.);
   mySolver = fsolver[0];
}


ICPG1D::~ICPG1D()
{
   if (_init_alloc) {
      delete _r;
      delete _v;
      delete _p;
   }
}


void ICPG1D::setInitialCondition_shock_tube(const LocalVect<real_t,3>& BcG,
                                            const LocalVect<real_t,3>& BcD,
                                                  real_t               x0)
{
// Two differents zones (for shock tube initial solution)
   MESH_EL {
      size_t n = theElementLabel;
      Triang3 tria(theElement);
      if (tria.getCenter().x<x0) {
         (*_r)(n) = BcG[0];
         (*_v)(n) = BcG[1];
         (*_p)(n) = BcG[2];
      }
      else {
         (*_r)(n) = BcD[0];
         (*_v)(n) = BcD[1];
         (*_p)(n) = BcD[2];
      }
   }
}


void ICPG1D::setInitialCondition(const LocalVect<real_t,3>& U)
{
   MESH_EL {
      size_t n = theElementLabel;
      (*_r)(n) = U[0];
      (*_v)(n) = U[1];
      (*_p)(n) = U[2];
   }
}


void ICPG1D::setReconstruction()
{
   if (Muscl1D::setReconstruction(*_r,_Gr,_Dr,1)) {
      cout << "ERROR: reconstruction of rho failed" << endl;
      exit(3);
   }
   if (Muscl1D::setReconstruction(*_v,_Gv,_Dv,1)) {
      cout << "ERROR: reconstruction of u failed" << endl;
      exit(3);
   }
   if (Muscl1D::setReconstruction(*_p,_Gp,_Dp,1)) {
      cout << "ERROR: reconstruction of p failed" << endl;
      exit(3);
   }
}


real_t ICPG1D::RiemannSolver_ROE(int s)
{
   real_t lambda1, lambda2, lambda4;                            // eigenvalues
   real_t lambda_entropy_left, lambda_entropy_right;            // entropy correction
   real_t rho_star, u_star, p_star, c_star;	                // entropy correction state
   real_t alpha1, alpha2, alpha4;                               // amplitude of the waves
   real_t Roe_u, Roe_H, Roe_c;                                  // Roe averaged values
   real_t sqrt_rho_G=sqrt(_Gr(s)), sqrt_rho_D=sqrt(_Dr(s));     // square roots for rho_G and rho_D
   real_t one_over_sqrt = 1./(sqrt_rho_G+sqrt_rho_D);           // Normalization factor of the Roe average
   real_t one_over_rho_G=1./_Gr(s), one_over_rho_D=1./_Dr(s);   // Coefficient for Roe average
   sqrt_rho_G *= one_over_sqrt;
   sqrt_rho_D *= one_over_sqrt;

// Averages are normalized
   if (_Gr(s)<_min_density) {
      cerr << "Error: left minimum density rho = " << _Gr(s) << endl;
      exit(2);
   }
   if (_Dr(s)<_min_density) {
      cerr << "Error: right  minimum density rho = " << _Dr(s) << endl;
      exit(2);
   }

// get data
   real_t Gr=_Gr(s), Dr=_Dr(s), GU=_Gv(s), DU=_Dv(s);
   real_t GrU=Gr*GU, DrU=Dr*DU;
// compute pressure, internal energy and sound speed
   real_t GP = _Gp(s);
   real_t GE = GP/(_Gamma-1.) + 0.5*Gr*(GU*GU);
   if (GP<0.)
      cerr << "\nWARNING, GP <0" << endl;
   real_t GH = (GE+GP)*one_over_rho_G; 
   real_t DP = _Dp(s);
   real_t DE = DP/(_Gamma-1.) + 0.5*Dr*(DU*DU);
   if (DP<0.)
      cerr << "\nWARNING: DP <0" << endl;
   real_t DH = (DE+DP)*one_over_rho_D; 
// compute Roe average
   Roe_u = sqrt_rho_G*GU + sqrt_rho_D*DU;
   Roe_H = sqrt_rho_G*GH + sqrt_rho_D*DH;
   Roe_c = (_Gamma-1.)*(Roe_H-(Roe_u*Roe_u)*0.5);   //  c*c
   if (Roe_c<0) {
      cerr << "Error: Roe sound speed negative: " << Roe_c << endl
           << "in point " << s << " on boundary ? (0/1):  " << _theMesh->getPtrSide(s)->isOnBoundary()
           << "x = " << Triang3(_theMesh->getPtrSide(s)).getCenter().x
           << "y = " << Triang3(_theMesh->getPtrSide(s)).getCenter().y << endl;
      exit(2);
   }
   Roe_c = sqrt(Roe_c);

// Compute eigenvalues
   lambda1 = Roe_u - Roe_c;
   lambda2 = Roe_u;
   lambda4 = Roe_u + Roe_c;

// Compute flux
   if (0.<lambda1) {
      _Fr(s) = GrU;
      _FrU(s) = GrU*GU + GP;
      _FE(s) = GU*(GE+GP);
      return (Roe_u+Roe_c);
   }

   alpha2 = (Dr-Gr)*(Roe_H-Roe_u*Roe_u) - (DE-GE);
   alpha2 *= (_Gamma-1)/(Roe_c*Roe_c);
   alpha1 = (Dr-Gr)*(Roe_c+Roe_u) - (DrU-GrU) - alpha2*Roe_c;
   alpha1 *= 0.5/Roe_c;

   if (0.<=lambda2) { 
//    entropy fix if necessary (Toro  1997 p. 338)
      lambda_entropy_left = GU - _Gc;
      rho_star = Gr + alpha1;
      u_star = ( GrU+alpha1*(Roe_u-Roe_c))/rho_star;
      p_star = (_Gamma-1)*(GE+alpha1*(Roe_H-Roe_u*Roe_c) - 0.5*rho_star*u_star*u_star);
      if (p_star<0) {
         cerr << "Error: negative pressure in left entropy fixed method" << endl;
         exit(2);
      }
      c_star = sqrt(_Gamma*p_star/rho_star);
      lambda_entropy_right = u_star - c_star;
      if (lambda_entropy_left<0. && 0.<lambda_entropy_right) //we have to fix the entropy
         lambda1 = lambda_entropy_left*(lambda_entropy_right-lambda1)/(lambda_entropy_right-lambda_entropy_left);
      _Fr(s) = GrU + alpha1*lambda1;
      _FrU(s) = GrU*GU + GP + alpha1*lambda1*(Roe_u-Roe_c);
      _FE(s) = GU*(GE+GP) + alpha1*lambda1*(Roe_H-Roe_u*Roe_c);
      return (Roe_u+Roe_c);
   }

   alpha4 = (Dr-Gr) - alpha2 - alpha1;
   if (0.<lambda4) {
//    entropy fix if necessary (Toro  1997 p. 338)
      lambda_entropy_right = DU + _Dc;
      rho_star = Dr - alpha4;
      u_star = (DrU-alpha4*(Roe_u+Roe_c))/rho_star;
      p_star = (_Gamma-1)*(DE-alpha4*(Roe_H+Roe_u*Roe_c)-0.5*rho_star*u_star*u_star);
      if (p_star<0) {
         cerr << "Error: negative pressure in right entropy fix method" << endl;
         exit(2);
      }
      c_star = sqrt(_Gamma*p_star/rho_star);
      lambda_entropy_left=u_star+c_star;
      if (lambda_entropy_left<0. && 0.<lambda_entropy_right) //we have to fix the entropy
         lambda4 = lambda_entropy_right*(lambda4-lambda_entropy_left)/(lambda_entropy_right-lambda_entropy_left);
      _Fr(s) = DrU - alpha4*lambda4;
      _FrU(s) = DrU*DU + DP - alpha4*lambda4*(Roe_u+Roe_c);
      _FE(s) = DU*(DE+DP) - alpha4*lambda4*(Roe_H+Roe_u*Roe_c);
      return (-Roe_u+Roe_c);
   }

   _Fr(s) = DrU;
   _FrU(s) = DrU*DU + DP;
   _FE(s) = DU*(DE+DP);
   return (-Roe_u+Roe_c);
}


real_t ICPG1D::RiemannSolver_VFROE(int s)
{
   real_t lambda1, lambda2, lambda4;                            // eigenvalues
   real_t lambda_entropy_left, lambda_entropy_right;            // entropy correction
   real_t rho_star, u_star, p_star, c_star;	                // entropy correction state
   real_t alpha1, alpha2, alpha4;                               // amplitude of the waves
   real_t Roe_u, Roe_H, Roe_c;                                  // Roe averaged values
   real_t sqrt_rho_G=sqrt(_Gr(s)), sqrt_rho_D=sqrt(_Dr(s));     // square roots for rho_G and rho_D
   real_t one_over_sqrt = 1./(sqrt_rho_G+sqrt_rho_D);           // normalization factor of the Roe average
   real_t one_over_rho_G=1./_Gr(s), one_over_rho_D=1./_Dr(s);   // coefficient for Roe average
   sqrt_rho_G *= one_over_sqrt;
   sqrt_rho_D *= one_over_sqrt;

   if (_Gr(s)<_min_density) {
      cerr << "Error: left minimum density rho = " << _Gr(s,1) << endl;
      exit(2);
   }
   if (_Dr(s)<_min_density) {
      cerr << "Error: right  minimum density rho = " << _Dr(s,1) << endl;
      exit(2);
   }
   real_t Gr=_Gr(s), Dr=_Dr(s);

// compute velocities
   real_t GU=_Gv(s), DU=_Dv(s);
   real_t GrU=Gr*GU, DrU=Dr*DU;

// compute pressure, internal energy and sound speed
   real_t GP = _Gp(s);
   real_t GE = GP/(_Gamma-1.) + 0.5*Gr*(GU*GU);
   if (GP<0.)
      cerr << "WARNING, GP < 0" << endl;
   real_t GH = (GE+GP)*one_over_rho_G; 
   real_t DP = _Dp(s);
   real_t DE = DP/(_Gamma-1.) + 0.5*Dr*(DU*DU);
   if (DP<0.)
      cerr << "WARNING, DP <0" << endl;
   real_t DH = (DE+DP)*one_over_rho_D;

// compute Roe average
   Roe_u = 0.5*(GU+DU);
   Roe_H = 0.5*(GH+DH);
   Roe_c = (_Gamma-1.)*(Roe_H-(Roe_u*Roe_u)*0.5);   //  c*c !!!!

   if (Roe_c<0) {
      cerr << "Error: Roe sound speed negative : " << Roe_c << endl;
      exit(2);
   }
   Roe_c = sqrt(Roe_c);

// compute eigenvalues
   lambda1 = Roe_u - Roe_c;
   lambda2 = Roe_u;
   lambda4 = Roe_u + Roe_c;
    
// compute flux
   if (0.<lambda1) {
      _Fr(s) = GrU;
      _FrU(s) = GrU*GU + GP;
      _FE(s) = GU*(GE+GP);
      return (Roe_u+Roe_c);
   }

   alpha2  = (Dr-Gr)*(Roe_H-Roe_u*Roe_u) - (DE-GE);
   alpha2 *= (_Gamma-1)/(Roe_c*Roe_c);
   alpha1  = (Dr-Gr)*(Roe_c+Roe_u) - (DrU-GrU) - alpha2*Roe_c;
   alpha1 *= 0.5/Roe_c;

   if (0.<=lambda2)	{ 
//    entropy fix if necessary (Toro 1997 p. 338)
      lambda_entropy_left = GU - _Gc;  //$$$$$ _GC n'a jamais ete initialise
      rho_star = Gr + alpha1;
      u_star = (GrU+alpha1*(Roe_u-Roe_c))/rho_star;
      p_star = (_Gamma-1)*(GE+alpha1*(Roe_H-Roe_u*Roe_c) - 0.5*rho_star*u_star*u_star);
      if (p_star<0) {
         cerr << "Error: negative pressure in left entropy fix method" << endl;
         exit(2);
      }
      c_star = sqrt(_Gamma*p_star/rho_star);
      lambda_entropy_right = u_star - c_star;
      if (lambda_entropy_left<0. && 0.<lambda_entropy_right) //we have to fix the entropy
         lambda1 = lambda_entropy_left*(lambda_entropy_right-lambda1)/(lambda_entropy_right-lambda_entropy_left);
      _Fr(s) = GrU + alpha1*lambda1;
      _FrU(s) = GrU*GU + GP + alpha1*lambda1*(Roe_u-Roe_c);
      _FE(s) = GU*(GE+GP) + alpha1*lambda1*(Roe_H-Roe_u*Roe_c);
      return (Roe_u+Roe_c);
   }

   alpha4 = (Dr-Gr) - alpha2 - alpha1;
   if (0.<lambda4) {
//    entropy fix if necessary (Toro  1997 p. 338)
      lambda_entropy_right = DU + _Dc;
      rho_star = Dr - alpha4;
      u_star = (DrU-alpha4*(Roe_u+Roe_c))/rho_star;
      p_star = (_Gamma-1)*(DE-alpha4*(Roe_H+Roe_u*Roe_c)-0.5*rho_star*u_star*u_star);
      if (p_star<0) {
         cerr << "Error: negative pressure in right entropy fix method" << endl;
         exit(2);
      }
      c_star = sqrt(_Gamma*p_star/rho_star);
      lambda_entropy_left = u_star + c_star;
      if (lambda_entropy_left<0. && 0.<lambda_entropy_right) //we have to fix the entropy
         lambda4 = lambda_entropy_right*(lambda4-lambda_entropy_left)/(lambda_entropy_right-lambda_entropy_left);
      _Fr(s) = DrU - alpha4*lambda4;
      _FrU(s) = DrU*DU + DP - alpha4*lambda4*(Roe_u+Roe_c);
      _FE(s) = DU*(DE+DP) - alpha4*lambda4*(Roe_H+Roe_u*Roe_c);
      return (-Roe_u+Roe_c);
   }

   _Fr(s) = DrU;
   _FrU(s) = DrU*DU + DP;
   _FE(s) = DU*(DE+DP);
   return (-Roe_u+Roe_c);
}


real_t ICPG1D::RiemannSolver_LF(int s)
{
   real_t Roe_u, Roe_H, Roe_c;
   real_t one_over_rho_G=1./_Gr(s,1), one_over_rho_D=1./_Dr(s,1);

   if (_Gr(s)<_min_density)	{
      cerr << "Warning: left minimum density rho = " << _Gr(s) << endl;
      _Gr(s) = _min_density;
   }
   if (_Dr(s)<_min_density)	{
      cerr << "Warning: right minimum density rho = " << _Dr(s) << endl;
      _Dr(s) = _min_density;
   }

// get data
   real_t Gr=_Gr(s), Dr=_Dr(s), GU=_Gv(s), DU=_Dv(s);
   real_t GrU=Gr*GU, DrU=Dr*DU;

// compute pressure, internal energy and sound speed
   real_t GP = _Gp(s);
   real_t GE = GP/(_Gamma-1.) + 0.5*Gr*(GU*GU);
   if (GP<0.)
      cerr << "WARNING: GP <0" << endl;
   real_t GH = (GE+GP)*one_over_rho_G;

   real_t DP = _Dp(s);
   real_t DE = DP/(_Gamma-1.) + 0.5*Dr*(DU*DU);
   if (DP<0.)
      cerr << "WARNING: DP <0" << endl;
   real_t DH = (DE+DP)*one_over_rho_D; 

// use VFRoe method to compute eigenvalues
// compute VFRoe average
   Roe_u = 0.5*(GU+DU);
   Roe_H = 0.5*(GH+DH);
   Roe_c = (_Gamma-1.)*(Roe_H-(Roe_u*Roe_u)*0.5); // Here c*c

   if (Roe_c<0) {
      cerr << "Warning: Roe sound speed negative: " << Roe_c << endl;
      exit(2);
   }
   Roe_c = sqrt(Roe_c);

// compute LF flux
   _Fr(s) = (GrU+DrU)*0.5;
   _FrU(s) = (GrU*GU+GP+DrU*DU+DP)*0.5;
   _FE(s) = (GU*(GE+GP)+DU*(DE+DP))*0.5;
   if (_TimeStep) {
      _Fr(s) -= 0.5*_ReferenceLength/_TimeStep*(Dr-Gr);
      _FrU(s) -= 0.5*_ReferenceLength/_TimeStep*(DrU-GrU);
      _FE(s) -= 0.5*_ReferenceLength/_TimeStep*(DE-GE);
   }
   return (_ReferenceLength/_TimeStep);
}


real_t ICPG1D::RiemannSolver_RUSANOV(int s)
{
   real_t one_over_rho_G=1./_Gr(s), one_over_rho_D=1./_Dr(s);
   real_t Splus;

   if (_Gr(s)<_min_density) {
      cerr << "Warning: left minimum density rho = " << _Gr(s,1) << endl;
      _Gr(s) = _min_density;
   }
   if (_Dr(s)<_min_density)	{
      cerr << "Warning: right  minimum density rho = " << _Dr(s,1) << endl;
      _Dr(s) = _min_density;
   }

// get data
   real_t Gr=_Gr(s), Dr=_Dr(s);
   real_t GU=_Gv(s), DU=_Dv(s);
   real_t GrU=Gr*GU, DrU=Dr*DU;

// compute pressure, internal energy and sound speed
   real_t GP = _Gp(s);
   real_t GE = GP/(_Gamma-1.) + 0.5*Gr*(GU*GU);
   if (GP<0.)
      cerr << "WARNING: GP <0" << endl;
   real_t GH = (GE+GP)*one_over_rho_G; 
   real_t DP = _Dp(s);
   real_t DE = DP/(_Gamma-1.) + 0.5*Dr*(DU*DU);
   if (DP<0.)
      cerr << "WARNING: DP <0" << endl;
   real_t DH = (DE+DP)*one_over_rho_D;
    
// compute S+
   real_t Gc = (_Gamma-1.)*(GH-(GU*GU)*0.5); // Here c*c
   real_t Dc = (_Gamma-1.)*(DH-(DU*DU)*0.5); // Here c*c

   if (Gc<0. || Dc<0.) {
      cerr << "Error: sound speed negative: " << Gc << ", " << Dc << endl;
      exit(2);
   }
   Gc = sqrt(Gc);
   Dc = sqrt(Dc);
   Splus = max(fabs(GU-Gc),fabs(GU+Gc));
   Splus = max(Splus,fabs(DU-Dc));
   Splus = max(Splus,fabs(DU+Dc));

// compute Rusanov flux
   _Fr(s) = (GrU+DrU)*0.5;
   _FrU(s) = (GrU*GU+GP+DrU*DU+DP)*0.5;
   _FE(s) = (GU*(GE+GP)+DU*(DE+DP))*0.5;
   _Fr(s) -= 0.5*Splus*(Dr-Gr);
   _FrU(s) -= 0.5*Splus*(DrU-GrU);
   _FE(s) -= 0.5*Splus*(DE-GE);
   return Splus;
}


real_t ICPG1D::RiemannSolver_HLL(int s)
{
   real_t one_over_rho_G=1./_Gr(s), one_over_rho_D=1./_Dr(s);
   real_t SL, SR, Splus, C1, C2, C3;

   if (_Gr(s)<_min_density)	{
      cerr << "Warning: left minimum density rho = " << _Gr(s) << endl;
      _Gr(s) = _min_density;
   }
   if (_Dr(s)<_min_density)	{
      cerr << "Warning: right minimum density rho = " << _Dr(s) << endl;
      _Dr(s) = _min_density;
   }

// get data
   real_t Gr=_Gr(s), Dr=_Dr(s);
   real_t GU=_Gv(s), DU=_Dv(s);
   real_t GrU=Gr*GU, DrU=Dr*DU;

// compute pressure, internal energy and sound speed
   real_t GP = _Gp(s);
   real_t GE = GP/(_Gamma-1.) + 0.5*Gr*(GU*GU);
   if (GP<0.)
      cerr << "WARNING: GP <0" << endl;
   real_t GH = (GE+GP)*one_over_rho_G; 
   real_t DP = _Dp(s);
   real_t DE = DP/(_Gamma-1.) + 0.5*Dr*(DU*DU);
   if (DP<0.)
      cerr << "WARNING: DP < 0" << endl;
   real_t DH = (DE+DP)*one_over_rho_D;

// compute S+
   real_t Gc = (_Gamma-1.)*(GH-(GU*GU)*0.5); // Here c*c
   real_t Dc = (_Gamma-1.)*(DH-(DU*DU)*0.5); // Here c*c
   if ((Gc<0.) || (Dc<0.)) {
      cerr << "Warning: sound speed negative: " << Gc << ", " << Dc << endl;
      exit(2);
   }
   Gc = sqrt(Gc);
   Dc = sqrt(Dc);

   Splus = max(fabs(GU - Gc),fabs(GU + Gc));
   Splus = max(Splus,fabs(DU - Dc));
   Splus = max(Splus,fabs(DU + Dc));
   SL = -Splus;
   SR = Splus;
   C1 = SR/(SR-SL);
   C2 = -SL/(SR-SL);
   C3 = SL*SR/(SR-SL);

//See other choices of SR, SL

// compute HLL flux
   if (SL>=0.) {  // Flux = FG
      _Fr(s) = GrU ;
      _FrU(s) = GrU*GU + GP;
      _FE(s) = GU*(GE+GP);
      return(Splus);
   }
   if (SR<=0.) { // Flux = FD
      _Fr(s) = DrU ;
      _FrU(s) = DrU*DU + DP;
      _FE(s) = DU*(DE+DP);
   }

// (SL<=0.) && (SR >=0.)) => Flux = FHLL
   _Fr(s)  = C1*GrU + C2*DrU + C3*(Dr-Gr);
   _FrU(s) = C1*(GrU*GU+GP) + C2*(DrU*DU+DP) + C3*(DrU-GrU);
   _FE(s)  = C1*(GU*(GE+GP)) + C2*(DU*(DE+DP)) + C3*(DE-GE);
   return Splus;
}


real_t ICPG1D::RiemannSolver_HLLC(int s)
{
   real_t one_over_rho_G=1./_Gr(s), one_over_rho_D=1./_Dr(s);
   real_t SL, SR, Splus, Sstar;
   real_t rhoGstar, rhoDstar, uGstar, uDstar, EGstar, EDstar;
   real_t AL, AR;

   if (_Gr(s)<_min_density) {
      cerr << "WARNING: left minimum density rho = " << _Gr(s) << endl;
      _Gr(s) = _min_density;
   }
   if (_Dr(s)<_min_density) {
      cerr << "WARNING: right  minimum density rho = " << _Dr(s) << endl;
      _Dr(s) = _min_density;
   }

// get data
   real_t Gr=_Gr(s), Dr=_Dr(s);
   real_t GU=_Gv(s), DU=_Dv(s);
   real_t GrU=Gr*GU, DrU=Dr*DU;

// compute pressure, internal energy and sound speed
   real_t GP = _Gp(s);
   real_t GE = GP/(_Gamma-1.) + 0.5*Gr*(GU*GU);
   real_t GH = (GE+GP)*one_over_rho_G; 
   real_t DP = _Dp(s);
   real_t DE = DP/(_Gamma-1.) + 0.5*Dr*(DU*DU);
   real_t DH = (DE+DP)*one_over_rho_D;

// compute S+
   real_t Gc = (_Gamma-1.)*(GH-(GU*GU)*0.5); // Here c*c
   real_t Dc = (_Gamma-1.)*(DH-(DU*DU)*0.5); // Here c*c

   if (Gc<0. || Dc<0.) {
      cerr << "Warning: sound speed negative " << Gc << ", " << Dc << endl;
      exit(2);
   }
   Gc = sqrt(Gc);
   Dc = sqrt(Dc);

   Splus = max(fabs(GU - Gc),fabs(GU+ Gc));
   Splus = max(Splus,fabs(DU - Dc));
   Splus = max(Splus,fabs(DU + Dc));
   SL = -Splus;
   SR = Splus;
   Sstar = (DP - GP +GrU*(SL-GU) - DrU* (SR-DU))/(Gr *(SL-GU) - Dr *(SR-DU)) ;
// see other choices of SR, SL
   AL = Gr* (SL - GU)/(SL - Sstar);
   AR = Dr* (SR - DU)/(SR - Sstar);
// compute star region
   rhoGstar = AL;
   rhoDstar = AR;
   uGstar = AL * Sstar;
   uDstar = AR * Sstar;
   EGstar = AL*(GE*one_over_rho_G + (Sstar-GU)*(Sstar + GP/(Gr*(SL-GU))));
   EDstar = AR*(DE*one_over_rho_D + (Sstar-DU)*(Sstar + DP/(Dr*(SR-DU))));
// compute HLLC flux
   if (SL>=0.) { // Flux = FG
      _Fr(s) = GrU;
      _FrU(s) = GrU*GU + GP;
      _FE(s) = GU*(GE+GP);
      return Splus;
   }
   if (SR<=0.) {  // Flux = FD
      _Fr(s) = DrU;
      _FrU(s) = DrU*DU + DP;
      _FE(s) = DU*(DE+DP);
      return Splus;
   }
   if ( Sstar >=0.) { // Flux = FHLLCG
      _Fr(s) = GrU + SL*(rhoGstar-_Gr(s)) ;
      _FrU(s) = GrU*GU + GP + SL*(uGstar-GrU);
      _FE(s) = GU*(GE+GP) + SL*(EGstar-GE);
   }
   else { // Flux = FHLLCD
      _Fr(s) = DrU + SR*(rhoDstar-_Dr(s));
      _FrU(s) = DrU*DU + DP + SR*(uDstar-DrU);
      _FE(s) = DU*(DE+DP) + SR*(EDstar-DE);
   }
   return Splus;
}


real_t ICPG1D::getFlux()
{
   real_t lambda=1.e-10;
   MESH_SD {
      assert(mySolver);
      lambda = max(lambda,(this->*mySolver)(int(theSideLabel)));
   }
   return lambda;
}


void ICPG1D::setBC(const Side&  sd,
                         real_t u)
{
   size_t s = sd.n();
   _Dr(s) = _Gr(s);
   _Dv(s) = u*_Gv(s);
   _Dp(s) = _Gp(s);
}


void ICPG1D::setBC(int    code,
                   real_t u)
{
   MESH_SD
      if (theSide->getCode(1)==code)
         setBC(*theSide,u);
}


void ICPG1D::setBC(real_t u)
{
   MESH_SD
      setBC(*theSide,u);
}


void ICPG1D::setBC(const Side&                sd,
                   const LocalVect<real_t,3>& u)
{
   size_t s = sd.n();
   _Dr(s) = u(1);
   _Dv(s) = u(2);
   _Dp(s) = u(3);
}


void ICPG1D::setBC(      int                  code,
                   const LocalVect<real_t,3>& u)
{
   MESH_SD
      if (theSide->getCode()==code)
         setBC(*theSide,u);
}


void ICPG1D::setBC(const LocalVect<real_t,3>& u)
{
   MESH_SD
      setBC(*theSide,u);
}


void ICPG1D::setInOutflowBC(const Side&                sd,
                            const LocalVect<real_t,3>& u)
{
   size_t s = sd.n();
   real_t Gc = sqrt(_Gamma*_Gp(s)/_Gr(s)), GUn = _Gv(s);
   int choice = (GUn>=0.) + 2*(Gc>=1.);

   switch (choice) {

      case 0:                        // GU.n < 0 and Gc < 1  : subsonic inflow
        _Dr(s) = _Gr(s);             // getting from interior
        _Dv(s) = u[1];               // because exterior normal
        _Dp(s) = u[2]; 
        break;

      case 1:                       // GU.n >=0 and Gc < 1  : subsonic outflow
        _Dr(s) = _Gr(s);
        _Dv(s) = _Gv(s);
        _Dp(s) = u[2];              // imposing pressure
        break;

     case 2:                        // GU.n < 0 and Gc > 1  : supersonic inflow
        _Dr(s) = u[0];	            // imposing density
        _Dv(s) = u[1];              // because exterior normal
        _Dp(s) = u[2];
        break;

     case 3:                        //  GU.n > 0 and Gc > 1 : supersonic outflow
        _Dr(s) = _Gr(s);
        _Dv(s) = _Gv(s);
        _Dp(s) = _Gp(s);	    // getting from interior
        break;

     default:
        cerr << "ERROR in setBCInOutflow" << endl;
        exit(1);
   }
}


void ICPG1D::setInOutflowBC(      int                  code,
                            const LocalVect<real_t,3>& u)
{
   MESH_SD
      if (theSide->getCode(1)==code)
         setInOutflowBC(TheSide,u);
}


void ICPG1D::setInOutflowBC(const LocalVect<real_t,3>& u)
{
   MESH_SD
      setInOutflowBC(TheSide,u);
}

/*=====================================================
  SECTION ONE STEP METHOD TO EVALUATE VECTOR A TIME N+1
  Computation of the primal variable    n->n+1
  ===================================================*/

/*----------------------------------------------------------------------------
   The variable flux contains all fluxes issued from the Riemann solver with
   as left element getNeighborElement(1) and right element getNeighborElement(2).
   If getNeighborElement(2) is NULL, we are on the boundary and we impose
   the symmetry condition
  ----------------------------------------------------------------------------*/

void ICPG1D::fromPrimalToConservative()
{
   MESH_EL {
      size_t n = theElementLabel;
      real_t rho = (*_r)(n);
      (*_p)(n) = (*_p)(n)/(_Gamma-1.) + 0.5*rho*(*_v)(n)*(*_v)(n);
      (*_v)(n) *= rho;
   }
}


void ICPG1D::fromConservativeToPrimal()
{
   MESH_EL {
      size_t n = theElementLabel;
      (*_v)(n) /= (*_r)(n);
      (*_p)(n) = ((*_p)(n) - 0.5*(*_r)(n)*(*_v)(n)*(*_v)(n))*(_Gamma-1.) ;
   }
}


void ICPG1D::forward()
{
   Element *elg, *eld;
   size_t ns, ng, nd;
   real_t rate;
   fromPrimalToConservative();

// Warning: v is now rU and p is now E
   MESH_SD {
      ns = theSide->n();
      elg = theSide->getNeighborElement(1);
      ng = elg->n();
      rate = _Lrate(ns) * _TimeStep;
      (*_r)(ng) -=_Fr(ns) * rate;
      if ((*_r)(ng)<0.)
         cerr << "Warning: Negative density, Fr = " << _Fr(ns,1)*rate << endl;
      (*_v)(ng) -= _FrU(ns)*rate;
      (*_p)(ng) -= _FE(ns) *rate;
      if (!theSide->isOnBoundary()) {
         eld = theSide->getNeighborElement(2);
         nd = eld->n();
         rate = _Rrate(ns)*_TimeStep;
         (*_r)(nd) += _Fr(ns)*rate;
         (*_v)(nd) += _FrU(ns)*rate;
         (*_p)(nd) += _FE(ns)*rate;
      }
   }

// WARNING: we are now back to primal variables
   fromConservativeToPrimal();
}


void ICPG1D::Forward(const Vect<real_t>& flux,
                           Vect<real_t>& field)
{
   Element *lel, *rel;
   size_t ns, ln, rn;
   real_t rate;
   MESH_SD {
      ns = theSideLabel;
      lel = TheSide.getNeighborElement(1);
      ln = lel->n();
      rate = _Lrate(ns)*_TimeStep;
      field.add(ln,-flux(ns)*rate);
      if (!theSide->isOnBoundary()) {
         rel = theSide->getNeighborElement(2);
         rn = rel->n();
         rate = _Rrate(ns)*_TimeStep;
         field.add(rn,flux(ns)*rate);
      }
   }
}


void ICPG1D::getMomentum(Vect<real_t>& m) const
{
   m.setMesh(*_theMesh,1,ELEMENT_FIELD);
   MeshElements(*_theMesh) {
      size_t n = theElementLabel;
      m.set(n,(*_r)(n)*(*_v)(n));
   }
}


void ICPG1D::getInternalEnergy(Vect<real_t>& ie) const
{
   ie.setMesh(*_theMesh,1,ELEMENT_FIELD);
   MeshElements(*_theMesh) {
      size_t n = theElementLabel;
      ie.set(n,(*_p)(n)*(_Gamma -1.));
   }
}


void ICPG1D::getTotalEnergy(Vect<real_t>& te) const
{
   te.setMesh(*_theMesh,1,ELEMENT_FIELD);
   Vect<real_t> tmp_e;
   getInternalEnergy(tmp_e);
   MeshElements(*_theMesh) {
      size_t n = theElementLabel;
      te.set(n,tmp_e(n) + 0.5*(*_r)(n)*(*_v)(n)*(*_v)(n));
   }
}


void ICPG1D::getSoundSpeed(Vect<real_t>& s) const
{
   s.setMesh(*_theMesh,1,ELEMENT_FIELD);
   MeshElements(*_theMesh) {
      size_t n = theElementLabel;
      s.set(n,sqrt(_Gamma*(*_p)(n)/(*_r)(n)));
   }
}


void ICPG1D::getMach(Vect<real_t>& m) const
{
   m.setMesh(*_theMesh,1,ELEMENT_FIELD);
   Vect<real_t> temp(*_theMesh);
   getSoundSpeed(temp);
   MeshElements(*_theMesh) {
      size_t n = theElementLabel;
      m.set(n,(*_v)(n)/temp(n));
   }
}


real_t ICPG1D::runOneTimeStep()
{
   real_t V = getFlux();
   _TimeStep = _ReferenceLength/V*_CFL;
   forward();
   return _TimeStep;
}

} /* namespace OFELI */
