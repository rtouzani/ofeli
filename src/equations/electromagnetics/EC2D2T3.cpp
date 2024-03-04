/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2024 Rachid Touzani

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

            Class EC2D2T3 : Eddy Current Problems in Two-Dimensions
            with a vector magnetic field using the 3-Node triangle

  ==============================================================================*/


#include "equations/electromagnetics/EC2D2T3.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"
#include "linear_algebra/Vect_impl.h"


namespace OFELI {


EC2D2T3::EC2D2T3()
        : Equation<3,6,2,4>()
{ }


EC2D2T3::EC2D2T3(Mesh& ms)
        : Equation<3,6,2,4>(ms)
{ 
   if (Equa::_nb_dof!=2)
      throw OFELIException("In EC2D2T3::EC2D2T3(..): Nodes must have "
                            "2 degrees of freedom because of complex-valued formulation.");

}


EC2D2T3::EC2D2T3(Mesh&         ms,
                 Vect<real_t>& u)
        : Equation<3,6,2,4>(ms,u)
{
   if (Equa::_nb_dof!=2)
      throw OFELIException("In EC2D2T3::EC2D2T3(..): Nodes must have "
                            "2 degrees of freedom because of complex-valued formulation.");

}


EC2D2T3::EC2D2T3(const Side* sd1,
                 const Side* sd2)
{
   _ns = sd1->n(); _nt = sd2->n();
   _theSide = sd1, _theElement = nullptr;
   Node *nd1 = sd1->getPtrNode(1), *nd2 = sd1->getPtrNode(2);
   _i1 = nd1->n(); _j1 = nd2->n();
   _N1 = nd1->getCoord(); _N2 = nd2->getCoord();
   _ll1 = sqrt((_N1.x-_N2.x)*(_N1.x-_N2.x) + (_N1.y-_N2.y)*(_N1.y-_N2.y));
   nd1 = sd2->getPtrNode(1); nd2 = sd2->getPtrNode(2);
   _i2 = nd1->n(); _j2 = nd2->n();
   _M1 = nd1->getCoord(); _M2 = nd2->getCoord();
   _ll2 = sqrt((_M1.x-_M2.x)*(_M1.x-_M2.x) + (_M1.y-_M2.y)*(_M1.y-_M2.y));
}


void EC2D2T3::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   _label = el->n();
   setMaterial();
   Triang3 tr(el);
   _el_geo.area = tr.getArea();
   _el_geo.center = tr.getCenter();
   _x[0] = (*_theElement)(1)->getCoord();
   _x[1] = (*_theElement)(2)->getCoord();
   _x[2] = (*_theElement)(3)->getCoord();
   _dSh = tr.DSh();
   if (_omega_set)
      _omega = _omega_fct(_el_geo.center,0.);
   if (_Mu_set)
      _Mu = _Mu_fct(_el_geo.center,0.);
   if (_sigma_set)
      _sigma = _sigma_fct(_el_geo.center,0.);
   eMat = 0;
   eRHS = 0;
}


void EC2D2T3::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Node *nd1 = (*_theSide)(1), *nd2 = (*_theSide)(2);
   _i1 = nd1->n(); _j1 = nd2->n();
   _N1 = nd1->getCoord(); _N2 = nd2->getCoord();
   _ll1 = sqrt((_N1.x-_N2.x)*(_N1.x-_N2.x) + (_N1.y-_N2.y)*(_N1.y-_N2.y));
   _x[0] = (*_theSide)(1)->getCoord();
   _x[1] = (*_theSide)(2)->getCoord();
   Line2 ln(sd);
   _dSh = ln.DSh();
   _ns = sd->n();
   sMat = 0;
   sRHS = 0;
}


void EC2D2T3::RHS(real_t coef)
{
   eRHS(1) = eRHS(3) = eRHS(5) = coef*OFELI_THIRD*_el_geo.area*_sigma;
   eRHS(2) = eRHS(4) = eRHS(6) = 0;
}


void EC2D2T3::FEBlock()
{
   real_t c = 0.25*_Mu*_omega*OFELI_THIRD*_el_geo.area*_sigma;
   for (size_t i=1; i<=3; i++) {
      eMat(2*i-1,2*i  ) -= c;
      eMat(2*i  ,2*i-1) += c;
      for (size_t j=1; j<=3; j++) {
         real_t e = OFELI_THIRD*_el_geo.area*(_dSh[i-1],_dSh[j-1]);
         eMat(2*i-1,2*j-1) += e; eMat(2*i  ,2*j  ) += e;
         eMat(2*i-1,2*j  ) -= c; eMat(2*i  ,2*j-1) += c;
      }
   }
}


void EC2D2T3::BEBlocks(size_t            n1,
                       size_t            n2,
                       SpMatrix<real_t>& L,
                       SpMatrix<real_t>& U,
                       SpMatrix<real_t>& D)
{
   real_t aa = -0.5*_ll1;
   if (_ns==_nt) {
      U.add(2*_i1-1,2*_ns-1,aa);
      U.add(2*_i1  ,2*_ns  ,aa);
      U.add(2*_j1-1,2*_ns-1,aa);
      U.add(2*_j1  ,2*_ns  ,aa);
      L.add(2*_ns-1,2*_i1-1,aa);
      L.add(2*_ns  ,2*_i1  ,aa);
      L.add(2*_ns-1,2*_j1-1,aa);
      L.add(2*_ns  ,2*_j1  ,aa);
   }
   real_t d1, d2;
   _Dkl(_N1, _N2, _M1, _M2, d1);
   D.add(2*(_ns+n1)-1,2*(_nt+n2)-1,d1);
   D.add(2*(_ns+n1)  ,2*(_nt+n2)  ,d1);
   if (_ns != _nt) {
      _Lkl(_N1, _N2, _M1, _M2, d1, d2);
      L.add(2*_ns-1,2*_i2-1,d1);
      L.add(2*_ns  ,2*_i2  ,d1);
      L.add(2*_ns-1,2*_j2-1,d2);
      L.add(2*_ns  ,2*_j2  ,d2);
   }
}


complex_t EC2D2T3::Constant(const Vect<real_t>& u,
                            complex_t&          I)
{
   real_t ur = OFELI_THIRD*_el_geo.area*(u(1)+u(3)+u(5));
   real_t ui = OFELI_THIRD*_el_geo.area*(u(2)+u(4)+u(6));
   real_t c1 = 1./_sigma*(I.real() - _sigma*_omega*ui)/_el_geo.area;
   real_t c2 = 1./_sigma*(I.imag() + _sigma*_omega*ur)/_el_geo.area;
   return std::complex<real_t> (c1,c2);
}


real_t EC2D2T3::MagneticPressure(const Vect<real_t>& u)
{
   real_t c = 0.5*_Mu*_el_geo.area;
   real_t dxr=0, dxi=0, dyr=0, dyi=0;
   for (size_t i=1; i<=3; i++) {
      dxr += _dSh[i-1].x*u(2*i-1);
      dxi += _dSh[i-1].x*u(2*i  );
      dyr += _dSh[i-1].y*u(2*i-1);
      dyi += _dSh[i-1].y*u(2*i  );
   }
   return (c*(dxr*dxr+dxi*dxi+dyr*dyr+dyi*dyi));
}


size_t EC2D2T3::_log_det(const Point<real_t>& ck,
                         const Point<real_t>& cl)
/*-----------------------------------------------------------------------
      Choice of the log determination

      ck,cl: centers of lines [ak,bk] and [al,bl]

      det = 1 : cut on the negative real half axis
      det = 2 : cut on the negative imaginary half axis
      det = 3 : cut on the positive real half axis
      det = 4 : cut on the positive imaginary half axis
  -----------------------------------------------------------------------*/
{
   size_t ret = 0;
   Point<real_t> u = cl - ck;
   real_t t = Max(u.x,u.y,-u.x,-u.y);
   if (t == u.x)
      ret = 1;
   else if (t == u.y)
      ret = 2;
   else if (t == -u.x)
      ret = 3;
   else if (t == -u.y)
      ret = 4;
   return ret;
}


complex<real_t> EC2D2T3::_ablog(size_t          det,
                                complex<real_t> a,
                                complex<real_t> b,
                                real_t          t)
/*-----------------------------------------------------------------------
      ablog = a*b*(log(a)-t)
  -----------------------------------------------------------------------*/
{
   if (a==complex<real_t>(0))
      return complex<real_t>(0);
   real_t dd = 0.;
   if (det==3 && a.imag()<0)
      dd = 2*OFELI_PI;
   if (det==4 && a.imag()>=0. && a.real()<0)
      dd = -2*OFELI_PI;
   if (det==2 && a.imag()<0 && a.real()<0)
      dd = 2*OFELI_PI;
   return a*b*(Log(a)+complex<real_t>(0,dd)-complex<real_t>(t,0));
}


void EC2D2T3::_Lkl(const Point<real_t>& uk,
                   const Point<real_t>& vk,
                   const Point<real_t>& ul,
                   const Point<real_t>& vl,
                   real_t&              d1,
                   real_t&              d2)
/*-----------------------------------------------------------------------
     Coefficients of matrix L

     uk,vk : Vertex coordinates for line [ak,bk]
     ul,vl : Vertex coordinates for line [al,bl]
     d1,d2 : coefficients of basis functions for nodes al and bl

  -----------------------------------------------------------------------*/
{
   Point<real_t> ak = uk, bk = vk, al = ul, bl = vl;
   Point<real_t> ck = 0.5*(ak+bk), cl = 0.5*(al+bl);
   real_t gk = sqrt((bk.x-ak.x)*(bk.x-ak.x)+(bk.y-ak.y)*(bk.y-ak.y));

// If [ak,bk] and [al,bl] are colinear, then d1 = d2 = 0.

   Point<real_t> v1 = al - ak, v2 = bl - ak, v = (bk - ak)/gk;
   real_t u1 = v.x*v1.y - v.y*v1.x;
   real_t u2 = v.x*v2.y - v.y*v2.x;
   if (abs(u1)<=OFELI_EPSMCH && abs(u2)<=OFELI_EPSMCH) {
     d1 = d2 = 0;
     return;
   }
   size_t det = _log_det(ck,cl);
   complex_t alpha(al.x-ak.x,al.y-ak.y);
   complex_t beta(bl.x-al.x,bl.y-al.y);
   complex_t gamma(bk.x-ak.x,bk.y-ak.y);

// Compute second coefficient (Node bl)

   complex_t t1 = _ablog(det,alpha+beta-gamma,0.5*(alpha+beta-gamma),0.5);
   complex_t t2 = _ablog(det,alpha-gamma,0.5*(alpha-gamma),1.5);
   complex_t t3 = _ablog(det,alpha+beta-gamma,alpha-gamma,1.);
   complex_t t4 = _ablog(det,alpha+beta,0.5*(alpha+beta),0.5);
   complex_t t5 = _ablog(det,alpha,0.5*alpha,1.5);
   complex_t t6 = _ablog(det,alpha+beta,alpha,1.);
   t1 = (-t1-t2+t3+t4+t5-t6)/(beta*beta*gamma);
   complex_t tt = t1;
   d2 = -gk/OFELI_PI*(beta.imag()*t1.real() + beta.real()*t1.imag());

// Compute first coefficient (node al)

   complex_t one(1,0);
   t1 = _ablog(det,alpha+beta-gamma,one,1.);
   t2 = _ablog(det,alpha-gamma,one,1.);
   t3 = _ablog(det,alpha+beta,one,1.);
   t4 = _ablog(det,alpha,one,1.);
   t1 = (-t1+t2+t3-t4)/(beta*gamma) - tt;
   d1 = -gk/OFELI_PI*(beta.imag()*t1.real() + beta.real()*t1.imag());
}

void EC2D2T3::_Dkl(const Point<real_t>& uk,
                   const Point<real_t>& vk,
                   const Point<real_t>& ul,
                   const Point<real_t>& vl,
                   real_t&              d)
/*-----------------------------------------------------------------------
     Coefficient of matrix D

     uk,vk : Vertex coordinates for line [ak,bk]
     ul,vl : Vertex coordinates for line [al,bl]

  -----------------------------------------------------------------------*/
{
   Point<real_t> ak = uk, bk = vk, al = ul, bl = vl;
   Point<real_t> ck = 0.5*(ak + bk), cl = 0.5*(al + bl);
   real_t gk = sqrt((bk.x-ak.x)*(bk.x-ak.x)+(bk.y-ak.y)*(bk.y-ak.y));
   real_t gl = sqrt((bl(1)-al(1))*(bl(1)-al(1))+(bl(2)-al(2))*(bl(2)-al(2)));
   size_t det = _log_det(ck,cl);
   complex_t alpha(al.x-ak.x,al.y-ak.y);
   complex_t beta(bl.x-al.x,bl.y-al.y);
   complex_t gamma(bk.x-ak.x,bk.y-ak.y);
   complex_t t1 = _ablog(det,alpha+beta-gamma,0.5*(alpha+beta-gamma),1.5);
   complex_t t2 = _ablog(det,alpha-gamma,0.5*(alpha-gamma),1.5);
   complex_t t3 = _ablog(det,alpha+beta,0.5*(alpha+beta),1.5);
   complex_t t4 = _ablog(det,alpha,0.5*alpha,1.5);
   t1 = (-t1+t2+t3-t4)/(beta*gamma);
   d = gk*gl/OFELI_PI*t1.real();
}

} /* namespace OFELI */
