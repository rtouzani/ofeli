/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2026 Rachid Touzani

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

                         Implementation of class Muscl3DT

       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/


#include "equations/cl/Muscl3DT.h"
#include "mesh/Mesh.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/LocalVect_impl.h"
#include "linear_algebra/LocalMatrix_impl.h"
#include "linear_algebra/Point.h"
#include "shape_functions/Tetra4.h"
#include "shape_functions/Triang3.h"
#include <fstream>

using std::min;
using std::max;
using std::ofstream;

namespace OFELI {

Muscl3DT::Muscl3DT(Mesh& m) : Muscl(m)
{
   size_t nb_el = _theMesh->getNbElements(), nb_sd = _theMesh->getNbSides();
   _n.setSize(nb_sd);
   _Lrate.setSize(nb_sd);
   _Rrate.setSize(nb_sd);
   _BBv.setSize(nb_el,4);               // = B_j.x - B_i.x      4 neighbors
   _BBl.setSize(nb_el,4);               // = ||B_j - B_i||
   _BQv.setSize(nb_el,4);               // just like BQ
   _BQl.setSize(nb_el,4);               // just like BQ
   _beta12.setSize(nb_el);              // 3 coeff == all the matrix
   _beta13.setSize(nb_el);
   _beta14.setSize(nb_el);
   _order.setSize(nb_el);
   _BMv.setSize(nb_el,4);
   _BMl.setSize(nb_el,4);
   _betabis12.setSize(nb_el);
   _betabis13.setSize(nb_el);
   _betabis14.setSize(nb_el);
   _mu1.setSize(nb_el,4);               // decomposition r_1 with t_j
   _mu2.setSize(nb_el,4);               // decomposition r_2 with t_j
   _mu3.setSize(nb_el,4);               // decomposition r_3 with t_j
   _mu4.setSize(nb_el,4);               // decomposition r_4 with t_j

   _MinimumFaceArea = getMinSideMeasure(*_theMesh);
   _MinimumElementVolume = getMinElementMeasure(*_theMesh);
   _MaximumFaceArea = getMaxSideMeasure(*_theMesh);
   _MaximumElementVolume = getMaxElementMeasure(*_theMesh);
   _MeanFaceArea = getMeanSideMeasure(*_theMesh);
   _MeanElementVolume = getMinElementMeasure(*_theMesh);
   _MinimumEdgeLength = getMinSize(*_theMesh);
   _MinimumVolumebyArea = getMinSize32();
   _MaximumEdgeLength = getMaxSize(*_theMesh);

   Initialize();
   _method = FIRST_ORDER_METHOD;
   _limiter = MINMOD_LIMITER;
   _betalim = 0.05;
}


Muscl3DT::~Muscl3DT() { }


Point<real_t> Muscl3DT::getDetCrammer_33(const LocalMatrix<real_t,3,3>& A,
                                         const Point<real_t>&           B)
{
   real_t detA = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(2,1)*A(3,2)*A(1,3) -
                 A(1,3)*A(2,2)*A(3,1) - A(1,1)*A(2,3)*A(3,2) - A(3,3)*A(1,2)*A(2,1);
   Point<real_t> x;
   x.x = B.x*A(2,2)*A(3,3) + B.y*A(3,2)*A(1,3) + A(1,2)*A(2,3)*B.z -
         B.z*A(2,2)*A(1,3) - B.x*A(2,3)*A(3,2) - A(3,3)*B.y*A(1,2);
   x.y = A(1,1)*B.y*A(3,3) + A(2,1)*B.z*A(1,3) + B.x*A(2,3)*A(3,1) -
         A(1,3)*B.y*A(3,1) - A(1,1)*B.z*A(2,3) - A(3,3)*A(2,1)*B.x;
   x.z = A(1,1)*A(2,2)*B.z + A(2,1)*A(3,2)*B.x + A(1,2)*B.y*A(3,1) -
         A(3,1)*A(2,2)*B.x - A(1,1)*A(3,2)*B.y - B.z*A(1,2)*A(2,1);
   return x/detA;
}


void Muscl3DT::Initialize()
{
   Point<real_t> normal, center, AB, AC;
   size_t i;
   int bad_tetra=0;
   MESH_SD {
      _n(side_label) = the_side->getUnitNormal();
      the_element = The_side.getNeighborElement(1);

//    Compute the surface rate
      _Lrate(the_side->n()) = Triang3(the_side).getArea()/Tetra4(the_element).getVolume();
      if (!the_side->isOnBoundary()) {
         the_element = the_side->getNeighborElement(2);
         _Rrate(side_label) = Triang3(the_side).getArea()/Tetra4(the_element).getVolume();
      }
   }

// Precalculate the B_iB_j and B_iQ_j abd B_iM_j vectors
   MESH_EL {
      size_t n = element_label;
      if (the_element->isOnBoundary()==false) {
         Side *sd1=the_element->getPtrSide(1), *sd2=the_element->getPtrSide(2),
              *sd3=the_element->getPtrSide(3), *sd4=the_element->getPtrSide(4);
         Element *el1=sd1->getNeighborElement(1), *el2=sd2->getNeighborElement(1),
                 *el3=sd3->getNeighborElement(1), *el4=sd4->getNeighborElement(1);
         if (el1==the_element)
            el1 = sd1->getNeighborElement(2);
         if (el2==the_element)
            el2 = sd2->getNeighborElement(2);
         if (el3==the_element)
            el3 = sd3->getNeighborElement(2);
         if (el4==the_element)
            el4 = sd4->getNeighborElement(2);
         Tetra4 tetra(the_element), tetra1(el1), tetra2(el2), tetra3(el3), tetra4(el4);
         Triang3 tria1(sd1), tria2(sd2), tria3(sd3), tria4(sd4);

//       compute BiBj
         _BBv(n,1) = tetra1.getCenter() - tetra.getCenter();
         _BBv(n,2) = tetra2.getCenter() - tetra.getCenter();
         _BBv(n,3) = tetra3.getCenter() - tetra.getCenter();
         _BBv(n,4) = tetra4.getCenter() - tetra.getCenter();

//       compute BiMj : distance between face centers
         _BMv(n,1) = tria1.getCenter() - tetra.getCenter();
         _BMv(n,2) = tria2.getCenter() - tetra.getCenter();
         _BMv(n,3) = tria3.getCenter() - tetra.getCenter();
         _BMv(n,4) = tria4.getCenter() - tetra.getCenter();

//       compute BiQj
//       intersection plan (3pts) / droite (2pts)
         Point<real_t> a, A, B, C, D;
         real_t lambda;
         D = tetra.getCenter();
         for (i=1; i<=4; ++i) {
            A = the_element->getPtrSide(i)->getPtrNode(1)->getCoord();
            B = the_element->getPtrSide(i)->getPtrNode(2)->getCoord();
            C = the_element->getPtrSide(i)->getPtrNode(3)->getCoord();
            a = CrossProduct(B-A,C-A);
            lambda = -(a*(D-A)) / (a*_BBv(n,i));
            _BQv(n,i) = lambda*_BBv(n,i);
         }
      }

//    compute all the norms
      for (i=1; i<=4; ++i) {
         _BBl(n,i) = _BBv(n,i).Norm();
         _BQl(n,i) = _BQv(n,i).Norm();
         _BMl(n,i) = _BMv(n,i).Norm();
      }
   }

// compute geometric coefficients
   Point<real_t> x, B;
   LocalMatrix<real_t,3,3> A;
   LocalVect<Point<real_t>,4> t, r, t_ortho;
   LocalVect<real_t,4> alpha, alpha_ortho;
   LocalMatrix<real_t,4,4> PI, beta;
   MESH_EL {
      size_t n = element_label;

//    we first compute _beta and _betabis
      if (the_element->isOnBoundary()) {
         _beta12(n) = 0.; _beta13(n) = 0.; _beta14(n) = 0.;
         _betabis12(n) = 0.; _betabis13(n) = 0.; _betabis14(n) = 0.;
         _mu1(n,1) = 0.; _mu2(n,1) = 0.; _mu3(n,1) = 0.; _mu4(n,1) = 0.;
         _mu1(n,2) = 0.; _mu2(n,2) = 0.; _mu3(n,2) = 0.; _mu4(n,2) = 0.;
         _mu1(n,3) = 0.; _mu2(n,3) = 0.; _mu3(n,3) = 0.; _mu4(n,3) = 0.;
         _mu1(n,4) = 0.; _mu2(n,4) = 0.; _mu3(n,4) = 0.; _mu4(n,4) = 0.;
      }
      else {

//       compute the normalized vector
         for (i=1; i<=4; ++i) {
            t(i) = _BBv(n,i)/_BBl(n,i);
            r(i) = _BMv(n,i)/_BMl(n,i);
         }

//       From p- to p+
         for (i=2; i<=4; ++i) {
            A(1,i-1) = t(i).x;
            A(2,i-1) = t(i).y;
            A(3,i-1) = t(i).z;
         }
         B(1) = t(1).x;
         B(2) = t(1).y;
         B(3) = t(1).z;
         x = getDetCrammer_33(A,B);
         _beta12(n) = x.x; _beta13(n) = x.y; _beta14(n) = x.z;

//       From q- to q+
         for (i=2; i<=4; ++i) {
            A(1,i-1) = r(i).x;
            A(2,i-1) = r(i).y;
            A(3,i-1) = r(i).z;
         }
         B(1) = r(1).x; B(2) = r(1).y; B(3) = r(1).z;
         x = getDetCrammer_33(A,B);
         _betabis12(n) = x.x;
         _betabis13(n) = x.y;
         _betabis14(n) = x.z;

//       Compute coefficient _mu

//       first: compute t_ortho
         for (size_t k=1; k<=4; ++k) {
//          alpha = dot product r_k and t_k
            alpha(k) = r(k)*t(k);
//          build t_ortho
            t_ortho(k) = r(k) - alpha(k)*t(k);
//          Normalize. Warning t_ortho should be small even zero so we do not normalize
            if (t_ortho(k).Norm()<OFELI_EPSMCH) { // use a home made orthogonal vector but it does not interfere since m=Q
               t_ortho(k).x = -t(k).y;
               t_ortho(k).y = -t(k).x;
               t_ortho(k).z = 0.;
            }
            t_ortho(k) /= t_ortho(k).Norm();      // and we normalize
            alpha_ortho(k) = r(k)*t_ortho(k);
         }
//       second: compute PI  t_ortho with respect to t
         for (size_t k=1; k<=4; ++k) {
            for (i=1; i<=4; ++i) {
               size_t j = i;
               if (i==k)
                  continue;
               if (i>k)
                  j = i-1;
               A(1,j) = t(i).x; A(2,j) = t(i).y; A(3,j) = t(i).z;
            }
            B(1) = t_ortho(k).x; B(2) = t_ortho(k).y; B(3) = t_ortho(k).z;
            x = getDetCrammer_33(A,B);
            for (i=1; i<=4; ++i) {
               size_t j=i;
               if (i==k)
                  continue;
               if (i>k)
                  j = i-1;
               PI(k,i) = x(j);
            }
         }
//       third: compute _mu : r with respect t.
         _mu1(1) = alpha(1);
         for (i=1; i<=4; ++i) {
            if (i!=1)
               _mu1(i) = alpha_ortho(1)*PI(1,i);
         }
         _mu2(2) = alpha(2);
         for (i=1; i<=4; ++i) {
            if (i!=2)
               _mu2(i) = alpha_ortho(2)*PI(2,i);
         }
         _mu3(3) = alpha(3);
         for (i=1; i<=4; ++i) {
            if (i!=3)
               _mu3(i) = alpha_ortho(3)*PI(3,i);
         }
         _mu4(4) = alpha(4);
         for (i=1; i<=4; ++i)
            if (i!=4)
               _mu4(i) = alpha_ortho(4)*PI(4,i);

//       check the coefficient for H-compatibility
         beta(1,1) = -1;
         beta(1,2) = _beta12(n);
         beta(1,3) = _beta13(n);
         beta(1,4) = _beta14(n);
         beta(2,1) = 1./_beta12(n);
         beta(2,2) = -1;
         beta(2,3) = -beta(2,1)*_beta13(n);
         beta(2,4) = -beta(2,1)*_beta14(n);
         beta(3,1) = 1./_beta13(n,1);
         beta(3,2) = 1./beta(2,3);
         beta(3,3) = -1;
         beta(3,4) = -beta(3,1)*_beta14(n);
         beta(4,1) = 1./_beta14(n);
         beta(4,2) = 1./beta(2,4);
         beta(4,3) = 1./beta(3,4);
         beta(4,4) = -1;

         bool verif = false;
         for (size_t j=1; j<=4; j++) {
            for (size_t k=1; k<=4; k++)
               if (beta(j,k) > -_betalim)
                  verif = true;
         }
         if (verif) {
            bad_tetra++;
            _order(n) = 0;
         }
         else
            _order(n) = 1;
      }
   }

   cerr << "WARNING: There are " << bad_tetra << " bad elements among " << int(_theMesh->getNbElements()) << endl;
   _taulim = gettaulim();
   getgraphtau();
   _Comega = getComega();
   getgraphComega();
}


bool Muscl3DT::setReconstruction(const Vect<real_t>& U,
                                       Vect<real_t>& LU,
                                       Vect<real_t>& RU,
                                       size_t        dof)
{
   switch (_method) {

      case FIRST_ORDER_METHOD:
         FirstOrder(U,LU,RU,dof);
         break;

      case MULTI_SLOPE_M_METHOD:
         MultiSlopeM(U,LU,RU,dof);
         break;

      case MULTI_SLOPE_Q_METHOD:
         MultiSlopeQ(U,LU,RU,dof);
         break;

      default:
         return true;
   };
   return false;
}


void Muscl3DT::FirstOrder(const Vect<real_t>& U,
                                Vect<real_t>& LU,
                                Vect<real_t>& RU,
                                size_t        dof)
{
   Element *ell, *elr;
   MESH_SD {
      ell = the_side->getNeighborElement(1);
      if (the_side->isOnBoundary())
         elr = ell;
      else
         elr = the_side->getNeighborElement(2);
      LU(side_label,dof) = U(ell->n(),dof);
      RU(side_label,dof) = U(elr->n(),dof);
   }
}


void Muscl3DT::MultiSlopeQ(const Vect<real_t>& U,
                                 Vect<real_t>& LU,
                                 Vect<real_t>& RU,
                                 size_t        dof)
{
   MESH_EL {
      size_t n = element_label;
      Side *sd1 = the_element->getPtrSide(1); size_t ns1 = sd1->n();
      Side *sd2 = the_element->getPtrSide(2); size_t ns2 = sd2->n();
      Side *sd3 = the_element->getPtrSide(3); size_t ns3 = sd3->n();
      Side *sd4 = the_element->getPtrSide(4); size_t ns4 = sd4->n();
      real_t ue = U(n,dof);

//    The element is on the boundary, we keep first order scheme
      real_t val1=ue, val2=ue, val3=ue, val4=ue;
      if (the_element->isOnBoundary() || !(_order(n))) {
         if (sd1->isOnBoundary())
            RU(ns1,dof) = ue;
         if (sd2->isOnBoundary())
            RU(ns2,dof) = ue;
         if (sd3->isOnBoundary())
            RU(ns3,dof) = ue;
         if (sd4->isOnBoundary())
            RU(ns4,dof) = ue;
      }

//    We have an inner element, second order method
      else {
         Element *el1=sd1->getNeighborElement(1), *el2=sd2->getNeighborElement(1),
                 *el3=sd3->getNeighborElement(1), *el4=sd4->getNeighborElement(1);
         if (el1==the_element)
            el1 = sd1->getNeighborElement(2);
         if (el2==the_element)
            el2 = sd2->getNeighborElement(2);
         if (el3==the_element)
            el3 = sd3->getNeighborElement(2);
         if (el4==the_element)
            el4 = sd4->getNeighborElement(2);

         LocalVect<real_t,4> dU;
         dU[0] = U(el1->n(),dof) - ue; dU[1] = U(el2->n(),dof) - ue;
         dU[2] = U(el3->n(),dof) - ue; dU[3] = U(el4->n(),dof) - ue;

//       Compute slopeQ_plus
         LocalVect<real_t,4> SlopeQForward;
         SlopeQForward[0] = dU[0]/_BBl(n,1);
         SlopeQForward[1] = dU[1]/_BBl(n,2);
         SlopeQForward[2] = dU[2]/_BBl(n,3);
         SlopeQForward[3] = dU[3]/_BBl(n,4);

//       Compute slope_minus
         LocalVect<real_t,4> SlopeQBackward;
         SlopeQBackward[0] = SlopeQForward[1]*_beta12(n) + SlopeQForward[2]*_beta13(n) + SlopeQForward[3]*_beta14(n);
         SlopeQBackward[1] = (SlopeQForward[0] - SlopeQForward[2]*_beta13(n) - SlopeQForward[3]*_beta14(n))/_beta12(n);
         SlopeQBackward[2] = (SlopeQForward[0] - SlopeQForward[1]*_beta12(n) - SlopeQForward[3]*_beta14(n))/_beta13(n);
         SlopeQBackward[3] = (SlopeQForward[0] - SlopeQForward[1]*_beta12(n) - SlopeQForward[2]*_beta13(n))/_beta14(n);

//       Limiting slope
         LocalVect<real_t,4> SlopeQ;
         SlopeQ[0] = minmod(SlopeQBackward[0],SlopeQForward[0])*SlopeQForward[0];
         SlopeQ[1] = minmod(SlopeQBackward[1],SlopeQForward[1])*SlopeQForward[1];
         SlopeQ[2] = minmod(SlopeQBackward[2],SlopeQForward[2])*SlopeQForward[2];
         SlopeQ[3] = minmod(SlopeQBackward[3],SlopeQForward[3])*SlopeQForward[3];

         val1 = ue + SlopeQ[0]*_BQl(n,1); val2 = ue + SlopeQ[1]*_BQl(n,2);
         val3 = ue + SlopeQ[2]*_BQl(n,3); val4 = ue + SlopeQ[3]*_BQl(n,4);
      }
      (sd1->getNeighborElement(1)==the_element) ? LU(ns1,dof) = val1 : RU(ns1,dof) = val1;
      (sd2->getNeighborElement(1)==the_element) ? LU(ns2,dof) = val2 : RU(ns2,dof) = val2;
      (sd3->getNeighborElement(1)==the_element) ? LU(ns3,dof) = val3 : RU(ns3,dof) = val3;
      (sd4->getNeighborElement(1)==the_element) ? LU(ns4,dof) = val4 : RU(ns4,dof) = val4;
   }
}


void Muscl3DT::MultiSlopeM(const Vect<real_t>& U,
                                 Vect<real_t>& LU,
                                 Vect<real_t>& RU,
                                 size_t        dof)
{
   MESH_EL {
      size_t n = element_label;
      real_t ue = U(n,dof);
      Side *sd1 = the_element->getPtrSide(1); size_t ns1 = sd1->n();
      Side *sd2 = the_element->getPtrSide(2); size_t ns2 = sd2->n();
      Side *sd3 = the_element->getPtrSide(3); size_t ns3 = sd3->n();
      Side *sd4 = the_element->getPtrSide(4); size_t ns4 = sd4->n();

//    The element is on the boundary, we keep first order scheme
      real_t val1=ue, val2=ue, val3=ue, val4=ue;
      if (the_element->isOnBoundary() || !(_order(n))) {
         if (sd1->isOnBoundary())
            RU(ns1,dof) = ue;
         if (sd2->isOnBoundary())
            RU(ns2,dof) = ue;
         if (sd3->isOnBoundary())
            RU(ns3,dof) = ue;
         if (sd4->isOnBoundary())
            RU(ns4,dof) = ue;
      }

//    We have an inner element, second order method
      else {
         Element *el1=sd1->getNeighborElement(1), *el2=sd2->getNeighborElement(1),
                 *el3=sd3->getNeighborElement(1), *el4=sd4->getNeighborElement(1);
         if (el1==the_element)
            el1 = sd1->getNeighborElement(2);
         if (el2==the_element)
            el2 = sd2->getNeighborElement(2);
         if (el3==the_element)
            el3 = sd3->getNeighborElement(2);
         if (el4==the_element)
            el4 = sd4->getNeighborElement(2);

         LocalVect<real_t,4> dU;
         dU[0] = U(el1->n(),dof) - ue; dU[1] = U(el2->n(),dof) - ue;
         dU[2] = U(el3->n(),dof) - ue; dU[3] = U(el4->n(),dof) - ue;

//       Compute slopeQ_plus
         LocalVect<real_t,4> SlopeQForward;
         SlopeQForward[0] = dU[0]/_BBl(n,1); SlopeQForward[1] = dU[1]/_BBl(n,2);
         SlopeQForward[2] = dU[2]/_BBl(n,3); SlopeQForward[3] = dU[3]/_BBl(n,4);

         LocalVect<real_t,4> SlopeMForward;
         SlopeMForward[0] = _mu1(n,1)*SlopeQForward[0] + _mu1(n,2)*SlopeQForward[1] + _mu1(n,3)*SlopeQForward[2] + _mu1(n,4)*SlopeQForward[3];
         SlopeMForward[1] = _mu2(n,1)*SlopeQForward[0] + _mu2(n,2)*SlopeQForward[1] + _mu2(n,3)*SlopeQForward[2] + _mu2(n,4)*SlopeQForward[3];
         SlopeMForward[2] = _mu3(n,1)*SlopeQForward[0] + _mu3(n,2)*SlopeQForward[1] + _mu3(n,3)*SlopeQForward[2] + _mu3(n,4)*SlopeQForward[3];
         SlopeMForward[3] = _mu4(n,1)*SlopeQForward[0] + _mu4(n,2)*SlopeQForward[1] + _mu4(n,3)*SlopeQForward[2] + _mu4(n,4)*SlopeQForward[3];

//       Compute slopeM_moins
         LocalVect<real_t,4> SlopeMBackward;
         SlopeMBackward[0] = SlopeMForward[1]*_betabis12(n) + SlopeMForward[2]*_betabis13(n) + SlopeMForward[3]*_betabis14(n);
         SlopeMBackward[1] = (SlopeMForward[0] - SlopeMForward[2]*_betabis13(n) - SlopeMForward[3]*_betabis14(n))/_betabis12(n);
         SlopeMBackward[2] = (SlopeMForward[0] - SlopeMForward[1]*_betabis12(n) - SlopeMForward[3]*_betabis14(n))/_betabis13(n);
         SlopeMBackward[3] = (SlopeMForward[0] - SlopeMForward[1]*_betabis12(n) - SlopeMForward[2]*_betabis13(n))/_betabis14(n);

//       Compute limited slope
         LocalVect<real_t,4> SlopeM;
         SlopeM[0] = minmod(SlopeMBackward[0],SlopeMForward[0])*SlopeMForward[0];
         SlopeM[1] = minmod(SlopeMBackward[1],SlopeMForward[1])*SlopeMForward[1];
         SlopeM[2] = minmod(SlopeMBackward[2],SlopeMForward[2])*SlopeMForward[2];
         SlopeM[3] = minmod(SlopeMBackward[3],SlopeMForward[3])*SlopeMForward[3];

         val1 = ue + _order(n)*SlopeM[0]*_BMl(n,1);
         val2 = ue + _order(n)*SlopeM[1]*_BMl(n,2);
         val3 = ue + _order(n)*SlopeM[2]*_BMl(n,3);
         val4 = ue + _order(n)*SlopeM[3]*_BMl(n,4);
         real_t min1 = min(ue,U(el1->n(),dof)), min2 = min(ue,U(el2->n(),dof));
         real_t min3 = min(ue,U(el3->n(),dof)), min4 = min(ue,U(el4->n(),dof));
         real_t max1 = max(ue,U(el1->n(),dof)), max2 = max(ue,U(el2->n(),dof));
         real_t max3 = max(ue,U(el3->n(),dof)), max4 = max(ue,U(el4->n(),dof));
         val1 = min(val1,min1); val2 = min(val2,min2);
         val3 = min(val3,min3); val4 = min(val4,min4);
         val1 = max(val1,max1); val2 = max(val2,max2);
         val3 = max(val3,max3); val4 = max(val4,max4);
      }
      (sd1->getNeighborElement(1)==the_element) ? LU(ns1,dof) = val1 : RU(ns1,dof) = val1;
      (sd2->getNeighborElement(1)==the_element) ? LU(ns2,dof) = val2 : RU(ns2,dof) = val2;
      (sd3->getNeighborElement(1)==the_element) ? LU(ns3,dof) = val3 : RU(ns3,dof) = val3;
      (sd4->getNeighborElement(1)==the_element) ? LU(ns4,dof) = val4 : RU(ns4,dof) = val4;
   }
}


void Muscl3DT::GradientQ(const Vect<real_t>& U,
                               Vect<real_t>& LU,
                               Vect<real_t>& RU,
                               size_t        dof)
{

   MESH_EL {
      size_t n = the_element->n();
      real_t ue = U(n,dof);
      Side *sd1 = the_element->getPtrSide(1); size_t ns1 = sd1->n();
      Side *sd2 = the_element->getPtrSide(2); size_t ns2 = sd2->n();
      Side *sd3 = the_element->getPtrSide(3); size_t ns3 = sd3->n();
      Side *sd4 = the_element->getPtrSide(4); size_t ns4 = sd4->n();

//    The element is on the boundary: First order scheme
      real_t val1=ue, val2=ue, val3=ue, val4=ue;
      if (the_element->isOnBoundary()) {
         if (sd1->isOnBoundary())
            RU(ns1,dof) = ue;
         if (sd2->isOnBoundary())
            RU(ns2,dof) = ue;
         if (sd3->isOnBoundary())
            RU(ns3,dof) = ue;
         if (sd4->isOnBoundary())
            RU(ns4,dof) = ue;
      }

//    The element is on the boundary: First order scheme
      else {
         Element *el1=sd1->getNeighborElement(1), *el2=sd2->getNeighborElement(1),
                 *el3=sd3->getNeighborElement(1), *el4=sd4->getNeighborElement(1);
         if (el1==the_element)
            sd1->getNeighborElement(2);
         if (el2==the_element)
            sd2->getNeighborElement(2);
         if (el3==the_element)
            sd3->getNeighborElement(2);
         if (el4==the_element)
            sd4->getNeighborElement(2);

         LocalVect<real_t,4> qU, dU, ar;
         qU(1) = U(el1->n(),dof); qU(2) = U(el2->n(),dof);
         qU(3) = U(el3->n(),dof); qU(4) = U(el4->n(),dof);
         dU(1) = qU(1) - U(the_element->n(),dof); dU(2) = qU(2) - U(the_element->n(),dof);
         dU(3) = qU(3) - U(the_element->n(),dof); dU(4) = qU(4) - U(the_element->n(),dof);

//       MP Limiter
         Point<real_t> slope = getGrad(qU,the_element);
         ar(1) = slope*_BQv(n,1); ar(2) = slope*_BQv(n,2);
         ar(3) = slope*_BQv(n,3); ar(4) = slope*_BQv(n,4);
         real_t phi = 1.;  // 1 for a second order scheme, 0 for a first order scheme

         (ar(1)*dU(1)<=0.0) ? phi = 0.0 : phi = min(phi,dU(1)/ar(1));
         (ar(2)*dU(2)<=0.0) ? phi = 0.0 : phi = min(phi,dU(2)/ar(2));
         (ar(3)*dU(3)<=0.0) ? phi = 0.0 : phi = min(phi,dU(3)/ar(3));
         (ar(4)*dU(4)<=0.0) ? phi = 0.0 : phi = min(phi,dU(4)/ar(4));

         val1 = ue + phi*slope*_BQv(n,1); val2 = ue + phi*slope*_BQv(n,2);
         val3 = ue + phi*slope*_BQv(n,3); val4 = ue + phi*slope*_BQv(n,4);
      }
      (sd1->getNeighborElement(1)==the_element) ? LU(ns1,dof) = val1 : RU(ns1,dof) = val1;
      (sd2->getNeighborElement(1)==the_element) ? LU(ns2,dof) = val2 : RU(ns2,dof) = val2;
      (sd3->getNeighborElement(1)==the_element) ? LU(ns3,dof) = val3 : RU(ns3,dof) = val3;
      (sd4->getNeighborElement(1)==the_element) ? LU(ns4,dof) = val3 : RU(ns4,dof) = val4;
   }
}


void Muscl3DT::GradientM(const Vect<real_t>& U,
                               Vect<real_t>& LU,
                               Vect<real_t>& RU,
                               size_t        dof)
{
   MESH_EL {
      size_t n = element_label;
      real_t ue = U(n,dof);
      Side *sd1 = the_element->getPtrSide(1); size_t ns1 = sd1->n();
      Side *sd2 = the_element->getPtrSide(2); size_t ns2 = sd2->n();
      Side *sd3 = the_element->getPtrSide(3); size_t ns3 = sd3->n();
      Side *sd4 = the_element->getPtrSide(4); size_t ns4 = sd4->n();

//    The element is on the boundary: First order scheme
      real_t val1=ue, val2=ue, val3=ue, val4=ue;
      if (the_element->isOnBoundary()) {
         if (sd1->isOnBoundary())
            RU(ns1,dof) = ue;
         if (sd2->isOnBoundary())
            RU(ns2,dof) = ue;
         if (sd3->isOnBoundary())
            RU(ns3,dof) = ue;
         if (sd4->isOnBoundary())
            RU(ns4,dof) = ue;
      }
//    The element is inner: Second order method
      else {
         Element *el1=sd1->getNeighborElement(1), *el2=sd2->getNeighborElement(1),
                 *el3=sd3->getNeighborElement(1), *el4=sd4->getNeighborElement(1);
         if (el1==the_element)
            el1 = sd1->getNeighborElement(2);
         if (el2==the_element)
            el2 = sd2->getNeighborElement(2);
         if (el3==the_element)
            el3 = sd3->getNeighborElement(2);
         if (el4==the_element)
            el4 = sd4->getNeighborElement(2);

         LocalVect<real_t,4> qU, dU, ar;
         qU(1) = U(el1->n(),dof); qU(2) = U(el2->n(),dof);
         qU(3) = U(el3->n(),dof); qU(4) = U(el4->n(),dof);
         dU(1) = qU(1) - ue; dU(2) = qU(2) - ue;
         dU(3) = qU(3) - ue; dU(4) = qU(4) - ue;

//       MP LIMITER
         Point<real_t> slope = getGrad(qU,the_element);
         ar(1) = slope*_BMv(n,1); ar(2) = slope*_BMv(n,2);
         ar(3) = slope*_BMv(n,3); ar(4) = slope*_BMv(n,4);

         real_t phi = 1.;     //  phi=1 for a second order scheme phi=0 for a first order scheme
         (ar(1)*dU(1)<=0.0) ? phi = 0.0 : phi = min(phi,dU(1)/ar(1));
         (ar(2)*dU(2)<=0.0) ? phi = 0.0 : phi = min(phi,dU(2)/ar(2));
         (ar(3)*dU(3)<=0.0) ? phi = 0.0 : phi = min(phi,dU(3)/ar(3));
         (ar(4)*dU(4)<=0.0) ? phi = 0.0 : phi = min(phi,dU(4)/ar(4));

         val1 = ue + phi*slope*_BMv(n,1);
         val2 = ue + phi*slope*_BMv(n,2);
         val3 = ue + phi*slope*_BMv(n,3);
         val4 = ue + phi*slope*_BMv(n,4);
      }

//    set new left and right values
      (sd1->getNeighborElement(1)==the_element) ? LU(ns1,dof) = val1 : RU(ns1,dof) = val1;
      (sd2->getNeighborElement(1)==the_element) ? LU(ns2,dof) = val2 : RU(ns2,dof) = val2;
      (sd3->getNeighborElement(1)==the_element) ? LU(ns3,dof) = val3 : RU(ns3,dof) = val3;
      (sd4->getNeighborElement(1)==the_element) ? LU(ns4,dof) = val3 : RU(ns4,dof) = val4;
   }
}


Point<real_t> Muscl3DT::getGrad(const LocalVect<real_t,4>& qU,
                                      Element*             el)
{
/* qU contains U1, U2, U3 of the three adjacent triangles                            */
/* nel is the label of the central triangle to be treated                            */
/* On return, slope contains the gradient                                            */
   Point<real_t> x, B;
   LocalMatrix<real_t,3,3> A;
   Side *sd1 = el->getPtrSide(1), *sd2 = el->getPtrSide(2), *sd3 = el->getPtrSide(3), *sd4= el->getPtrSide(4);
   Element *el1=sd1->getNeighborElement(1), *el2=sd2->getNeighborElement(1);
   Element *el3=sd3->getNeighborElement(1), *el4=sd4->getNeighborElement(1);
   if (el1==el)
      el1 = sd1->getNeighborElement(2);
   if (el2==el)
      el2 = sd2->getNeighborElement(2);
   if (el3==el)
      el3 = sd3->getNeighborElement(2);
   if (el4==el)
      el4 = sd4->getNeighborElement(2);

// compute geometric stuff
   Tetra4 t1(el1), t2(el2), t3(el3), t4(el4);
   A(1,1) = t2.getCenter().x - t1.getCenter().x;
   A(1,2) = t2.getCenter().y - t1.getCenter().y;
   A(1,3) = t2.getCenter().z - t1.getCenter().z;
   A(2,1) = t3.getCenter().x - t1.getCenter().x;
   A(2,2) = t3.getCenter().y - t1.getCenter().y;
   A(2,3) = t3.getCenter().z - t1.getCenter().z;
   A(3,1) = t4.getCenter().x - t1.getCenter().x;
   A(3,2) = t4.getCenter().y - t1.getCenter().y;
   A(3,3) = t4.getCenter().z - t1.getCenter().z;
   B.x = qU(2) - qU(1);
   B.y = qU(3) - qU(1);
   B.z = qU(4) - qU(1);
   x = getDetCrammer_33(A,B);
   return x;
}


real_t Muscl3DT::getMinSize32()
{
   real_t maxS=OFELI_EPSMCH, minVsS=1./OFELI_EPSMCH;
   MESH_EL {
      Tetra4 tt(the_element);
      for (size_t i=1; i<=4; ++i) {
         Triang3 t(the_element->getPtrSide(i));
         if (t.getArea()>maxS)
            maxS = t.getArea();
      }
      if (tt.getVolume()/maxS<minVsS)
         minVsS = tt.getVolume()/maxS;
   }
   return minVsS;
}


real_t Muscl3DT::gettaulim()
{
   real_t taulim = 2.;

   MESH_EL {
      if (the_element->isOnBoundary()==false) {
         size_t n = the_element->n();
         for (size_t i=1; i<=4; ++i)
            taulim = min(taulim,_BBl(n,i)/_BQl(n,i));
      }
   }
   return taulim;
}


void Muscl3DT::getgraphtau()
{
   const size_t nb_tau_max = 10000;
   Vect<int> nbtau(nb_tau_max);
   MESH_EL {
      if (The_element.isOnBoundary()==false) {
         size_t n = element_label;
         for (int i=1; i<=4; ++i)
            nbtau[min(int((_BBl(n,i)/_BQl(n,i)-1)/4.*(int(nb_tau_max)-1)),int(nb_tau_max)-1)]++;
      }
   }
   ofstream F_graph( "tau.graph",std::ios_base::out );
   for (size_t i=0; i<nb_tau_max ; i++)
      F_graph << 1+ 4*real_t(i)/(int(nb_tau_max)-1) << "\t\t" << nbtau[i] << endl;
   F_graph.close();
}


real_t Muscl3DT::getComega()
{
   real_t Comega = 0.;
   LocalMatrix<real_t,4,4> beta;
   MESH_EL {
      size_t n = element_label;
      if (the_element->isOnBoundary()==false && _order(n)) {
         beta(1,1) = -1;
         beta(1,2) = _beta12(n);
         beta(1,3) = _beta13(n);
         beta(1,4) = _beta14(n);
         beta(2,1) = 1./_beta12(n);
         beta(2,2) = -1;
         beta(2,3) = -beta(2,1)*_beta13(n);
         beta(2,4) = -beta(2,1)*_beta14(n);
         beta(3,1) = 1./_beta13(n);
         beta(3,2) = 1./beta(2,3);
         beta(3,3) = -1;
         beta(3,4) = -beta(3,1)*_beta14(n);
         beta(4,1) = 1./_beta14(n);
         beta(4,2) = 1./beta(2,4),
         beta(4,3) = 1./beta(3,4);
         beta(4,4) = -1;
         for (size_t j=1; j<=4; ++j) {
            real_t tmp = 0;
            for (size_t k=1; k<=4; ++k) {
               if (k!=j)
                  tmp += beta(j,k)/_BBl(n,k);
            }
            tmp *= -_BBl(n,j);
            Comega = max(Comega,tmp);
        }
     }
   }
   return Comega;
}


void Muscl3DT::getgraphComega()
{
   const size_t nb_comega_max = 10000;
   LocalMatrix<real_t,4,4> beta;
   Vect<int> nbcomega(nb_comega_max);
   MESH_EL {
      size_t n = the_element->n();
      if (the_element->isOnBoundary()==false && _order(n)) {
         beta(1,1) = -1;
         beta(1,2) = _beta12(n);
         beta(1,3) = _beta13(n);
         beta(1,4) = _beta14(n);
         beta(2,1) = 1./_beta12(n);
         beta(2,2) = -1;
         beta(2,3) = -beta(2,1)*_beta13(n);
         beta(2,4) = -beta(2,1)*_beta14(n);
         beta(3,1) = 1./_beta13(n);
         beta(3,2) = 1./beta(2,3);
         beta(3,3) = -1;
         beta(3,4) = -beta(3,1)*_beta14(n);
         beta(4,1) = 1./_beta14(n);
         beta(4,2) = 1./beta(2,4);
         beta(4,3) = 1./beta(3,4);
         beta(4,4) = -1;
         for (int j=1; j<=4; ++j){
            real_t tmp = 0;
            for (int k=1; k<=4; ++k)
               if (k!=j)
                  tmp += beta(j,k) / _BBl(n,k);
            tmp *= -_BBl(n,j);
            nbcomega[min(max(int(tmp*(int(nb_comega_max)-1)/_Comega),0),int(nb_comega_max)-1)]++;
         }
      }
   }
   ofstream F_graph("comega.graph",std::ios_base::out);
   for (size_t i=0; i<nb_comega_max ; i++)
      F_graph << real_t(i)*_Comega/(int(nb_comega_max)-1) << "\t\t" << nbcomega[i] << endl;
   F_graph.close();
}


ostream& operator<<(      ostream&  s,
                    const Muscl3DT& m)
{
   s.setf(ios::right|ios::scientific);
   s << "Mesh statistics as calculated in Muscl3DT:" << endl;
   s << "* Smin\t\t" << m.getMinimumFaceArea() << endl;
   s << "* Smax\t\t" << m.getMaximumFaceArea() << endl;
   s << "* Smean\t\t" << m.getMeanFaceArea() << endl;
   s << "* Vmin\t\t" << m.getMinimumElementVolume() << endl;
   s << "* Vmax\t\t" << m.getMaximumElementVolume() << endl;
   s << "* Vmean\t\t" << m.getMeanElementVolume() << endl;
   s << "* hmin\t\t" << m.getMinimumEdgeLength() << endl;
   s << "* hmax\t\t" << m.getMaximumEdgeLength() << endl;
   s << "* V/S m\t\t" << m.getMinimumVolumebyArea() << endl << endl;
   s << "* hmin\t\t" << m.getMinimumEdgeLength() << "\t\n* sqrt(S)\t" << sqrt(m.getMinimumFaceArea()) << "\t" << endl;
   return s;
}

} /* namespace OFELI */
