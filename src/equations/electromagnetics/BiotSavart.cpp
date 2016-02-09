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

                     Implementation of class 'BiotSavart'

  ==============================================================================*/

#include "equations/electromagnetics/BiotSavart.h"
#include "shape_functions/Tetra4.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"

namespace OFELI {


BiotSavart::BiotSavart()
{
   _mu = 4*OFELI_PI*1.e-7;
   _C = false;
   _bound = false;
   _code = 0;
   _theEdgeList = NULL;
   _theSideList = NULL;
}


BiotSavart::BiotSavart(Mesh &ms)
{
   _mu = 4*OFELI_PI*1.e-7;
   _C = false;
   _bound = false;
   _code = 0;
   _theEdgeList = NULL;
   _theSideList = NULL;
}


BiotSavart::BiotSavart(Mesh &ms, const Vect<real_t> &J, Vect<real_t> &B, int code)
{
   _mu = 4*OFELI_PI*1.e-7;
   _theMesh = &ms;
   _code = code;
   setSurface();
   _B = &B;
   _C = false;
   _bound = false;
   _type = J.getDOFType();
   _theEdgeList = NULL;
   _theSideList = NULL;
   _J = &J;
   if (_type==SIDE_FIELD)
      _theSideList = new SideList(*_theMesh);
   if (_type==EDGE_FIELD)
      _theEdgeList = new EdgeList(*_theMesh);
}


BiotSavart::BiotSavart(Mesh &ms, const Vect<complex_t> &J, Vect<complex_t> &B, int code)
{
   _mu = 4*OFELI_PI*1.e-7;
   _theMesh = &ms;
   _code = code;
   setSurface();
   _BC = &B;
   _C = true;
   _bound = false;
   _type = J.getDOFType();
   _theEdgeList = NULL;
   _theSideList = NULL;
   _JC = &J;
   if (_type==SIDE_FIELD)
      _theSideList = new SideList(*_theMesh);
   if (_type==EDGE_FIELD)
      _theEdgeList = new EdgeList(*_theMesh);
}
   

BiotSavart::~BiotSavart()
{
   if (_theEdgeList)
      delete _theEdgeList;
   if (_theSideList)
      delete _theSideList;
}


void BiotSavart::setVolume()
{
   _meas.setSize(_theMesh->getNbElements());
   _center.setSize(_theMesh->getNbElements());
   size_t k=1;
   MeshElements(*_theMesh) {
      Tetra4 t(theElement);
      _meas(k) = t.getVolume();
      _center(k++) = t.getCenter();
   }
}


void BiotSavart::setSurface()
{
   _meas.setSize(_theMesh->getNbSides());
   _center.setSize(_theMesh->getNbSides());
   size_t k=1;
   MeshSides(*_theMesh) {
      Triang3 t(theSide);
      if (theSide->getCode()==_code) {
         _meas(k) = t.getArea();
         _center(k++) = t.getCenter();
      }
   }
}
   
   
void BiotSavart::setLine()
{
   _meas.setSize(_theMesh->getNbEdges());
   _center.setSize(_theMesh->getNbEdges());
   size_t k=1;
   MeshEdges(*_theMesh) {
      Line2 l(theEdge);
      if (theEdge->getCode()==_code) {
         _meas(k) = l.getLength();
         _center(k++) = l.getCenter();
      }
   }
}
   

void BiotSavart::setPermeability(real_t mu)
{
   _mu = mu;
}

   
void BiotSavart::setCurrentDensity(const Vect<real_t> &J)
{
   _J = &J;
   _C = false;
   _type = J.getDOFType();
}


void BiotSavart::setCurrentDensity(const Vect<complex_t> &J)
{
   _JC = &J;
   _C = true;
   _type = J.getDOFType();
}
   

void BiotSavart::setMagneticInduction(Vect<real_t> &B)
{
   _B = &B;
   _C = false;
}
      

void BiotSavart::setMagneticInduction(Vect<complex_t> &B)
{
   _BC = &B;
   _C = true;
}
   

int BiotSavart::run()
{
   Point<real_t> B;
   Point<complex_t> BC;
   size_t k=1;

   switch (_type) {

      case 3: 
         if (_bound==false) {
            if (_C) {
               MeshNodes(*_theMesh) {
                  BC = getBC3(theNode->getCoord());
                  (*_BC)(k,1) = BC.x; (*_BC)(k,2) = BC.y; (*_BC)(k,3) = BC.z;
                  k++;
               }
            }
            else {
               MeshNodes(*_theMesh) {
                  B = getB3(theNode->getCoord());
                  (*_B)(k,1) = B.x; (*_B)(k,2) = B.y; (*_B)(k,3) = B.z;
                  k++;
               }
            }
         }
         else {
            _theMesh->getBoundaryNodes();
            if (_C) {
               MeshBoundaryNodes(*_theMesh) {
                  BC = getBC3(theNode->getCoord());
                  (*_BC)(k,1) = BC.x; (*_BC)(k,2) = BC.y; (*_BC)(k,3) = BC.z;
                  k++;
               }
            }
            else {
               MeshBoundaryNodes(*_theMesh) {
                  B = getB3(theNode->getCoord());
                  (*_B)(k,1) = B.x; (*_B)(k,2) = B.y; (*_B)(k,3) = B.z;
                  k++;
               }
            }
         }
         break;

   case 2:
      _theSideList = new SideList(*_theMesh);
      _theSideList->selectCode(_code);
      if (_C) {
         MeshNodes(*_theMesh) {
            BC = getBC2(theNode->getCoord());
            (*_BC)(k,1) = BC.x; (*_BC)(k,2) = BC.y; (*_BC)(k,3) = BC.z;
            k++;
         }
      }
      else {
         MeshNodes(*_theMesh) {
            B = getB2(theNode->getCoord());
            (*_B)(k,1) = B.x; (*_B)(k,2) = B.y; (*_B)(k,3) = B.z;
            k++;
         }
      }
      break;

         
   case 1:
      _theEdgeList = new EdgeList(*_theMesh);
      _theEdgeList->selectCode(_code);
      if (_C) {
         MeshNodes(*_theMesh) {
            BC = getBC1(theNode->getCoord());
            (*_BC)(k,1) = BC.x; (*_BC)(k,2) = BC.y; (*_BC)(k,3) = BC.z;
            k++;
         }
      }
      else {
         MeshNodes(*_theMesh) {
            B = getB1(theNode->getCoord());
            (*_B)(k,1) = B.x; (*_B)(k,2) = B.y; (*_B)(k,3) = B.z;
            k++;
         }
      }
      break;
   }
   return 0;
}


Point<real_t> BiotSavart::getB3(Point<real_t> x)
{
   real_t xy, c=0.25*_mu/OFELI_PI;
   Point<real_t> y, J, B=0.;
   size_t k=1;
   MeshElements(*_theMesh) {
      J.x = (*_J)(k,1);
      J.y = (*_J)(k,2);
      J.z = (*_J)(k,3);
      y = x - _center(k);
      xy = y.Norm();
      B += c*_meas(k++)/(xy*xy*xy)*CrossProduct(J,y);
   }
   return B;
}


Point<complex_t> BiotSavart::getBC3(Point<real_t> x)
{
   real_t xy, c=0.25*_mu/OFELI_PI;
   Point<real_t> y, Jr, Ji, Br=0., Bi=0.;
   Point<complex_t> B;
   size_t k=1;
   MeshElements(*_theMesh) {
      Jr.x = (*_JC)(k,1).real();
      Jr.y = (*_JC)(k,2).real();
      Jr.z = (*_JC)(k,3).real();
      Ji.x = (*_JC)(k,1).imag();
      Ji.y = (*_JC)(k,2).imag();
      Ji.z = (*_JC)(k,3).imag();
      y = x - _center(k);
      xy = y.Norm();
      Br += c*_meas(k  )/(xy*xy*xy)*CrossProduct(Jr,y);
      Bi += c*_meas(k++)/(xy*xy*xy)*CrossProduct(Ji,y);
   }
   B.x = complex_t(Br.x,Bi.x);
   B.y = complex_t(Br.y,Bi.y);
   B.z = complex_t(Br.z,Bi.z);
   return B;
}


Point<real_t> BiotSavart::getB2(Point<real_t> x)
{
   real_t xy, c=0.25*_mu/OFELI_PI;
   Point<real_t> y, J, B=0.;
   size_t k=1;
   for (_theSideList->top(); (theSide=_theSideList->get());) {
      J.x = (*_J)(k,1);
      J.y = (*_J)(k,2);
      J.z = (*_J)(k,3);
      y = x - _center(k);
      xy = y.Norm();
      B += c*_meas(k++)/(xy*xy*xy)*CrossProduct(J,y);
   }
   return B;
}
   
   
Point<complex_t> BiotSavart::getBC2(Point<real_t> x)
{
   real_t xy, c=0.25*_mu/OFELI_PI;
   Point<real_t> y, Jr, Ji, Br=0., Bi=0.;
   Point<complex_t> B;
   size_t k=1;
   Side *sd;
   for (_theSideList->top(); (sd=_theSideList->get());) {
      Jr.x = (*_JC)(k,1).real();
      Jr.y = (*_JC)(k,2).real();
      Jr.z = (*_JC)(k,3).real();
      Ji.x = (*_JC)(k,1).imag();
      Ji.y = (*_JC)(k,2).imag();
      Ji.z = (*_JC)(k,3).imag();
      y = x - _center(k);
      xy = y.Norm();
      Br += c*_meas(k  )/(xy*xy*xy)*CrossProduct(Jr,y);
      Bi += c*_meas(k++)/(xy*xy*xy)*CrossProduct(Ji,y);
   }
   B.x = complex_t(Br.x,Bi.x);
   B.y = complex_t(Br.y,Bi.y);
   B.z = complex_t(Br.z,Bi.z);
   return B;
}


Point<real_t> BiotSavart::getB1(Point<real_t> x)
{
   real_t xy, c=0.25*_mu/OFELI_PI;
   Point<real_t> y, J, B=0.;
   size_t k=1;
   for (_theEdgeList->top(); (theEdge=_theEdgeList->get());) {
      J.x = (*_J)(k,1);
      J.y = (*_J)(k,2);
      J.z = (*_J)(k,3);
      y = x - _center(k);
      xy = y.Norm();
      B += c*_meas(k++)/(xy*xy*xy)*CrossProduct(J,y);
   }
   return B;
}


Point<complex_t> BiotSavart::getBC1(Point<real_t> x)
{
   real_t xy, c=0.25*_mu/OFELI_PI;
   Point<real_t> y, Jr, Ji, Br=0., Bi=0.;
   Point<complex_t> B;
   size_t k=1;
   for (_theEdgeList->top(); (theEdge=_theEdgeList->get());) {
      cout<<theEdge->n()<<endl;
      Jr.x = (*_JC)(k,1).real();
      Jr.y = (*_JC)(k,2).real();
      Jr.z = (*_JC)(k,3).real();
      Ji.x = (*_JC)(k,1).imag();
      Ji.y = (*_JC)(k,2).imag();
      Ji.z = (*_JC)(k,3).imag();
      y = x - _center(k);
      xy = y.Norm();
      Br += c*_meas(k  )/(xy*xy*xy)*CrossProduct(Jr,y);
      Bi += c*_meas(k++)/(xy*xy*xy)*CrossProduct(Ji,y);
   }
   B.x = complex_t(Br.x,Bi.x);
   B.y = complex_t(Br.y,Bi.y);
   B.z = complex_t(Br.z,Bi.z);
   return B;
}

} /* namespace OFELI */
