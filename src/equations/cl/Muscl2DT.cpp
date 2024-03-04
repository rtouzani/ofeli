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

                         Implementation of class Muscl2DT

       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/


#include "equations/cl/Muscl2DT.h"
#include "mesh/Mesh.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/LocalVect_impl.h"
#include "linear_algebra/LocalMatrix_impl.h"
#include "shape_functions/Line2.h"
#include "shape_functions/Triang3.h"
#include <fstream>
#include <algorithm>

using std::min;
using std::max;
using std::ofstream;
using std::cout;
using std::to_string;

namespace OFELI {

Muscl2DT::Muscl2DT(Mesh& m) : Muscl(m)
{
   _n.setSize(_nb_sides);
   _Lrate.setSize(_nb_sides);         // ratio surface/left length of side
   _Rrate.setSize(_nb_sides);         // ratio surface/right length of side
   _BBv.setSize(_nb_elements,3);      //  B_j - B_i
   _BBl.setSize(_nb_elements,3);      // ||B_j - B_i||
   _BQv.setSize(_nb_elements,3);      // Q_ij - B_i
   _BQl.setSize(_nb_elements,3);      // ||Q_ij - B_i||
   _BMv.setSize(_nb_elements,3);      //  M_ij - B_i
   _BMl.setSize(_nb_elements,3);      // ||M_ij - B_i||
   _beta12.setSize(_nb_elements);     // Decomposition t_1 with t_2 and t_3
   _beta13.setSize(_nb_elements);
   _betabis12.setSize(_nb_elements);  // Decomposition r_1 with r_2 and r_3
   _betabis13.setSize(_nb_elements);
   _mu1.setSize(_nb_elements,3);      // Decomposition r_1 with t_j
   _mu2.setSize(_nb_elements,3);      // Decomposition r_2 with t_j
   _mu3.setSize(_nb_elements,3);      // Decomposition r_3 with t_j
   _order.setSize(_nb_elements);
   _MinimumFaceArea = getMinSideMeasure(*_theMesh);
   _MinimumElementVolume = getMinElementMeasure(*_theMesh);
   _MaximumFaceArea = getMaxSideMeasure(*_theMesh);
   _MaximumElementVolume = getMaxElementMeasure(*_theMesh);
   _MeanFaceArea = getMeanSideMeasure(*_theMesh);
   _MeanElementVolume = getMeanElementMeasure(*_theMesh);
   _MinimumEdgeLength = getMinSize(*_theMesh);
   _MinimumVolumebyArea	= getMinSize32();
   _MaximumEdgeLength = getMaxSize(*_theMesh);
   Initialize();
   _method = FIRST_ORDER_METHOD;
   _limiter = MINMOD_LIMITER;
   _betalim = 0.05;
   if (_verbose) {
      std::cout << "   Automatically set betalim to " << _betalim << endl;
      std::cout << "   Automatically set limiter to Minmod" << endl;
      std::cout << "   Automatically set method to First Order" << endl;
   }
}


Muscl2DT::~Muscl2DT() { }


void Muscl2DT::Initialize()
{
   MESH_SD {
      size_t s = side_label;
      _n.set(s,The_side.getUnitNormal());
      real_t L = Line2(the_side).getLength();

//    compute surface rate
      _Lrate.set(s,L/Triang3(The_side.getNeighborElement(1)).getArea());
      _Rrate.set(s,0.);
      if (The_side.isOnBoundary()==false)
         _Rrate.set(s,L/Triang3(The_side.getNeighborElement(2)).getArea());
   }

// Precalculate the B_iB_j and B_iQ_j abd B_iM_j vectors
// (for the second order method)
   size_t bad_triang=0;
   MESH_EL {
      size_t n = element_label;
      Side *sd1 = The_element.getPtrSide(1),
           *sd2 = The_element.getPtrSide(2),
           *sd3 = The_element.getPtrSide(3);
      if (The_element.isOnBoundary()==true) {
         _BBv.set(n,1,0.); _BBv.set(n,2,0.); _BBv.set(n,3,0.);
         _BQv.set(n,1,0.); _BQv.set(n,2,0.); _BQv.set(n,3,0.);
         _BMv.set(n,1,0.); _BMv.set(n,2,0.); _BMv.set(n,3,0.);
      }
      else {
         Element *el1 = sd1->getNeighborElement(1),
                 *el2 = sd2->getNeighborElement(1),
                 *el3 = sd3->getNeighborElement(1);
         if (el1==the_element)
            el1 = sd1->getNeighborElement(2);
         if (el2==the_element)
            el2 = sd2->getNeighborElement(2);
         if (el3==the_element)
            el3 = sd3->getNeighborElement(2);

//       compute BiBj: distance center to center
         Triang3 tr(the_element), tr1(el1), tr2(el2), tr3(el3);
         Point<real_t> c(tr.getCenter());
         _BBv.set(n,1,tr1.getCenter() - c);
         _BBv.set(n,2,tr2.getCenter() - c);
         _BBv.set(n,3,tr3.getCenter() - c);

//       compute BiMj : distance edge center to edge center
         _BMv.set(n,1,Line2(sd1).getCenter() - c);
         _BMv.set(n,2,Line2(sd2).getCenter() - c);
         _BMv.set(n,3,Line2(sd3).getCenter() - c);

//       compute BiQj
//       1
         Point<real_t> P1=(*sd1)(1)->getCoord(), P2=(*sd1)(2)->getCoord();
         real_t det1 = _determ(P1-c,P2-P1),
                det2 = _determ(tr1.getCenter()-c,P2-P1);
         if (det1/det2<0.0)
            throw OFELIException("Muscl2DT::Initialize(): Element " + to_string(element_label) +
                                 " is inconsistent with Muscl method.");
         _BQv.set(n,1,_BBv(n,1)*det1/det2);
//       2
         P1=(*sd2)(1)->getCoord(); P2=(*sd2)(2)->getCoord();
         det1 = _determ(P1-c,P2-P1);
         det2 = _determ(tr2.getCenter()-c,P2-P1);
         if (det1/det2<0.0)
            throw OFELIException("Muscl2DT::Initialize(): Element " + to_string(element_label) +
                                 " is inconsistent with Muscl method.");
         _BQv.set(n,2,_BBv(n,2)*det1/det2);
//       3
         P1=(*sd3)(1)->getCoord(), P2=(*sd3)(2)->getCoord();
         det1 = _determ(P1-c,P2-P1);
         det2 = _determ(tr3.getCenter()-c,P2-P1);
         if (det1/det2<0.0)
            throw OFELIException("Muscl2DT::Initialize(): Element "
                                 + to_string(element_label) +
                                 " is inconsistent with Muscl method.");
         _BQv.set(n,3,_BBv(n,3)*det1/det2);
      }

//    Compute norms
      for (int i=1; i<=3; ++i) {
         _BBl.set(n,i,_BBv(n,i).Norm());
         _BQl.set(n,i,_BQv(n,i).Norm());
         _BMl.set(n,i,_BMv(n,i).Norm());
      }
   }

// Q and M Methods
   MESH_EL {
      size_t n = element_label;
      if (The_element.isOnBoundary()==true) {
         _beta12.set(n,0.); _beta13.set(n,0.);
         _betabis12.set(n,0.); _betabis13.set(n,0.);
         _mu1.set(n,1,0.); _mu1.set(n,2,0.); _mu1.set(n,3,0.);
         _mu2.set(n,1,0.); _mu2.set(n,2,0.); _mu2.set(n,3,0.);
         _mu3.set(n,1,0.); _mu3.set(n,2,0.); _mu3.set(n,3,0.);
      }
      else {

//       compute normalized directions t_i
         Point<real_t> t1=_BBv(n,1)/_BBl(n,1), t2=_BBv(n,2)/_BBl(n,2), t3=_BBv(n,3)/_BBl(n,3);
         Point<real_t> s1=_BMv(n,1)/_BMl(n,1), s2=_BMv(n,2)/_BMl(n,2), s3=_BMv(n,3)/_BMl(n,3);

//       compute t_ortho an orthogonal vector to t
         Point<real_t> t1_ortho(t1.y,-t1.x), t2_ortho(t2.y,-t2.x), t3_ortho(t3.y,-t3.x);

//       compute det
         real_t det12=_determ(t1,t2), det13=_determ(t1,t3), det23=_determ(t2,t3);
         real_t dets12=_determ(s1,s2), dets13=_determ(s1,s3), dets23=_determ(s2,s3);

//       Compute dot products
         real_t ps11=t1*s1, ps22=t2*s2, ps33=t3*s3;
         real_t pso11=t1_ortho*s1, pso22=t2_ortho*s2, pso33=t3_ortho*s3;

//       compute beta: the decomposition of t_i related to the two other directions t_i
         _beta12.set(n, det13/det23);
         _beta13.set(n,-det12/det23);

//       compute _betabis12 and _betabis_13:
//       the decomposition of s_i related to the two other directions s_j
         _betabis12.set(n, dets13/dets23);
         _betabis13.set(n,-dets12/dets23);

//       Compute mu: decompose s_i in t_i and t_ortho_i
         _mu1.set(n,1,ps11);
         _mu1.set(n,2,pso11*_determ(t1_ortho,t3)/_determ(t2,t3));
         _mu1.set(n,3,pso11*_determ(t1_ortho,t2)/_determ(t3,t2));
         _mu2.set(n,1,pso22*_determ(t2_ortho,t3)/_determ(t1,t3));
         _mu2.set(n,2,ps22);
         _mu2.set(n,3,pso22*_determ(t2_ortho,t1)/_determ(t3,t1));
         _mu3.set(n,1,pso33*_determ(t3_ortho,t2)/_determ(t1,t2));
         _mu3.set(n,2,pso33*_determ(t3_ortho,t1)/_determ(t2,t1));
         _mu3.set(n,3,ps33);

//       Order test of negativity condition on beta's
         LocalMatrix<real_t,3,3> beta;
         int verif = 0;
         beta(1,1) = -1; beta(1,2) = _beta12(n); beta(1,3) = _beta13(n);
         beta(2,1) = 1./_beta12(n); beta(2,2) = -1.; beta(2,3) = -beta(2,1)*_beta13(n);
         beta(3,1) = 1./_beta13(n); beta(3,2) = 1./beta(2,3); beta(3,3) = -1.;
         for (size_t j=1; j<=3; j++) {
            for (size_t k=1; k<=3; k++)
               if (beta(j,k) > -_betalim)
                  verif = 1;
         }
         _order.set(n,1);
         if (verif) {
            bad_triang++;
            _order.set(n,0);
         }
      }
   }

   if (bad_triang)
      cerr << "WARNING: There are " << bad_triang << " bad triangles among " 
           << _nb_elements << endl;
   _taulim = gettaulim();
   _Comega = getComega();
}


bool Muscl2DT::setReconstruction(const Vect<real_t>& U,
                                 Vect<real_t>&       LU,
                                 Vect<real_t>&       RU,
                                 size_t              dof)
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


void Muscl2DT::FirstOrder(const Vect<real_t>& U,
                          Vect<real_t>&       LU,
                          Vect<real_t>&       RU,
                          size_t              dof)
{
   MESH_SD {
      Element *Lel = The_side.getNeighborElement(1), *Rel = nullptr;
      (The_side.isOnBoundary()) ? Rel = Lel : Rel = The_side.getNeighborElement(2);
      LU.set(side_label,dof,U(Lel->n(),dof));
      RU.set(side_label,dof,U(Rel->n(),dof));
   }
}


void Muscl2DT::MultiSlopeQ(const Vect<real_t>& U,
                           Vect<real_t>&       LU,
                           Vect<real_t>&       RU,
                           size_t              dof)
{
   MESH_EL {
      size_t n = element_label;
      real_t ue = U(n,dof);
      Side *sd1 = The_element.getPtrSide(1),
           *sd2 = The_element.getPtrSide(2),
           *sd3 = The_element.getPtrSide(3);
      size_t ns1=sd1->n(), ns2=sd2->n(), ns3=sd3->n();
      real_t val1=ue, val2=ue, val3=ue;
      if (The_element.isOnBoundary() || _order(n)==0) {
         if (sd1->isOnBoundary())
            RU.set(ns1,dof,ue);
         if (sd2->isOnBoundary())
            RU.set(ns2,dof,ue);
         if (sd3->isOnBoundary())
            RU.set(ns3,dof,ue);
      }
      else {
         Element *el1 = sd1->getNeighborElement(1),
                 *el2 = sd2->getNeighborElement(1),
                 *el3 = sd3->getNeighborElement(1);
         if (el1==the_element)
            el1 = sd1->getNeighborElement(2);
         if (el2==the_element)
            el2 = sd2->getNeighborElement(2);
         if (el3==the_element)
            el3 = sd3->getNeighborElement(2);

//       compute slope_plus
         real_t dU1 = U(el1->n(),dof) - ue,
                dU2 = U(el2->n(),dof) - ue,
                dU3 = U(el3->n(),dof) - ue;
         real_t SlopeQForward1 = dU1/_BBl(n,1),
                SlopeQForward2 = dU2/_BBl(n,2),
                SlopeQForward3 = dU3/_BBl(n,3);

//       compute slope_moins
         real_t b12 = _beta12(n,1), b13 = _beta13(n,1), b23 = -b12/b13;
         real_t SlopeQBackward1 = SlopeQForward2*b12 + SlopeQForward3*b13;
         real_t SlopeQBackward2 = SlopeQForward1/b12 + SlopeQForward3*b23;
         real_t SlopeQBackward3 = SlopeQForward1/b13 + SlopeQForward2/b23;

//       compute limited slope
         real_t SlopeQ1 = minmod(SlopeQBackward1,SlopeQForward1)*SlopeQForward1;
         real_t SlopeQ2 = minmod(SlopeQBackward2,SlopeQForward2)*SlopeQForward2;
         real_t SlopeQ3 = minmod(SlopeQBackward3,SlopeQForward3)*SlopeQForward3;
         val1 += SlopeQ1*_BQl(n,1);
         val2 += SlopeQ2*_BQl(n,2);
         val3 += SlopeQ3*_BQl(n,3);
      }
      (sd1->getNeighborElement(1)==the_element) ? LU.set(ns1,dof,val1) : RU.set(ns1,dof,val1);
      (sd2->getNeighborElement(1)==the_element) ? LU.set(ns2,dof,val2) : RU.set(ns2,dof,val2);
      (sd3->getNeighborElement(1)==the_element) ? LU.set(ns3,dof,val3) : RU.set(ns3,dof,val3);
   }
}


void Muscl2DT::MultiSlopeM(const Vect<real_t>& U,
                           Vect<real_t>&       LU,
                           Vect<real_t>&       RU,
                           size_t              dof)
{
   MESH_EL {
      size_t n = element_label;
      Side *sd1 = The_element.getPtrSide(1),
           *sd2 = The_element.getPtrSide(2),
           *sd3 = The_element.getPtrSide(3);
      size_t ns1=sd1->n(), ns2=sd2->n(), ns3=sd3->n();
      real_t ue = U(n,dof);
      real_t val1=ue, val2=ue, val3=ue;
      if (The_element.isOnBoundary() || _order(n)==0) {
         if (sd1->isOnBoundary())
            RU.set(ns1,dof,ue);
         if (sd2->isOnBoundary())
            RU.set(ns2,dof,ue);
         if (sd3->isOnBoundary())
            RU.set(ns3,dof,ue);
      }
      else {
         Element *el1 = sd1->getNeighborElement(1),
                 *el2 = sd2->getNeighborElement(1),
                 *el3 = sd3->getNeighborElement(1);
         if (el1==the_element)
            el1 = sd1->getNeighborElement(2);
         if (el2==the_element)
            el2 = sd2->getNeighborElement(2);
         if (el3==the_element)
            el3 = sd3->getNeighborElement(2);

//       compute slope_plus
         real_t SlopeQF1 = (U(el1->n(),dof) - ue)/_BBl(n,1),
                SlopeQF2 = (U(el2->n(),dof) - ue)/_BBl(n,2),
                SlopeQF3 = (U(el3->n(),dof) - ue)/_BBl(n,3);

//       compute slopeM_plus using _pi decomposition
         real_t SlopeMF1 = _mu1(n,1)*SlopeQF1 + _mu1(n,2)*SlopeQF2 + _mu1(n,3)*SlopeQF3,
                SlopeMF2 = _mu2(n,1)*SlopeQF1 + _mu2(n,2)*SlopeQF2 + _mu2(n,3)*SlopeQF3,
                SlopeMF3 = _mu3(n,1)*SlopeQF1 + _mu3(n,2)*SlopeQF2 + _mu3(n,3)*SlopeQF3;

//       compute slopeM_minus using betabis decomposition
         real_t b12bis = _betabis12(n), b13bis = _betabis13(n);
         real_t b23bis = -b12bis/b13bis;
         real_t SlopeMB1 = SlopeQF2*b12bis + SlopeMF3*b13bis,
                SlopeMB2 = SlopeQF1/b12bis + SlopeMF3*b23bis,
                SlopeMB3 = SlopeQF1/b13bis + SlopeMF2/b23bis;

//       compute limited slope using limiter
         val1 += minmod(SlopeMB1,SlopeMF1)*SlopeMF1*_BMl(n,1);
         val2 += minmod(SlopeMB2,SlopeMF2)*SlopeMF2*_BMl(n,2);
         val3 += minmod(SlopeMB3,SlopeMF3)*SlopeMF3*_BMl(n,3);

//       stability limiter
         real_t min1 = min(ue,U(el1->n(),dof)), 
                min2 = min(ue,U(el2->n(),dof)),
                min3 = min(ue,U(el3->n(),dof));
         real_t max1 = max(ue,U(el1->n(),dof)),
                max2 = max(ue,U(el2->n(),dof)),
                max3 = max(ue,U(el3->n(),dof));
         if (val1<min1)
            val1 = min1;
         if (val2<min2)
            val2 = min2;
         if (val3<min3)
            val3 = min3;
         if (val1>max1)
            val1 = max1;
         if (val2>max2)
            val2 = max2;
         if (val3>max3)
            val3 = max3;
      }
      (sd1->getNeighborElement(1)==the_element) ? LU.set(ns1,dof,val1) : RU.set(ns1,dof,val1);
      (sd2->getNeighborElement(1)==the_element) ? LU.set(ns2,dof,val2) : RU.set(ns2,dof,val2);
      (sd3->getNeighborElement(1)==the_element) ? LU.set(ns3,dof,val3) : RU.set(ns3,dof,val3);
   }
}


void Muscl2DT::Grad(const LocalVect<real_t,3>& U,
                    Point<real_t>&             slope,
                    size_t                     n)
{
   Point<real_t> u(U(2)-U(1),U(3)-U(1));
   Point<real_t> A=_BBv(n,2)-_BBv(n,1), B=_BBv(n,3)-_BBv(n,1);
   real_t det = _determ(A,B);
   slope = Point<real_t>(_determ(u,B)/det,_determ(A,u)/det);
}


void Muscl2DT::LeastSquare(const LocalVect<real_t,3>& U,
                           Point<real_t>&             slope,
                           size_t                     n)
{
   real_t a11 = _BBv(n,1).x*_BBv(n,1).x + _BBv(n,2).x*_BBv(n,2).x + _BBv(n,3).x*_BBv(n,3).x;
   real_t a12 = _BBv(n,1).x*_BBv(n,1).y + _BBv(n,2).x*_BBv(n,2).y + _BBv(n,3).x*_BBv(n,3).y;
   real_t a22 = _BBv(n,1).y*_BBv(n,1).y + _BBv(n,2).y*_BBv(n,2).y + _BBv(n,3).y*_BBv(n,3).y;
   real_t b1 = _BBv(n,1).x*U(1) + _BBv(n,2).x*U(2) + _BBv(n,3).x*U(3);
   real_t b2 = _BBv(n,1).y*U(1) + _BBv(n,2).y*U(2) + _BBv(n,3).y*U(3);
   real_t a = a11*a22 - a12*a12;
   slope = Point<real_t>((b1*a22-b2*a12)/a,(a11*b2-a12*b1)/a);
}


void Muscl2DT::GradientQ(const Vect<real_t>& U,
                         Vect<real_t>&       LU,
                         Vect<real_t>&       RU,
                         size_t              dof)
{
   LocalVect<real_t,3> qU, dU, ar;
   MESH_EL {
      size_t n = element_label;
      real_t ue = U(n,dof);
      Side *sd1 = The_element.getPtrSide(1); size_t ns1 = sd1->n();
      Side *sd2 = The_element.getPtrSide(2); size_t ns2 = sd2->n();
      Side *sd3 = The_element.getPtrSide(3); size_t ns3 = sd3->n();
      real_t val1=ue, val2=ue, val3=ue;
      if (The_element.isOnBoundary()) {
         if (sd1->isOnBoundary())
            RU.set(ns1,dof,ue);
         if (sd2->isOnBoundary())
            RU.set(ns2,dof,ue);
         if (sd3->isOnBoundary())
            RU.set(ns3,dof,ue);
      }
      else {
         Element *el1=sd1->getNeighborElement(1),
                 *el2=sd2->getNeighborElement(1),
                 *el3=sd3->getNeighborElement(1);
         if (el1==the_element)
            el1 = sd1->getNeighborElement(2);
         if (el2==the_element)
            el2 = sd2->getNeighborElement(2);
         if (el3==the_element)
            el3 = sd3->getNeighborElement(2);
         qU(1) = U(el1->n(),dof), qU(2) = U(el2->n(),dof), qU(3) = U(el3->n(),dof);
         dU(1) = qU(1) - ue, dU(2) = qU(2) - ue, dU(3) = qU(3) - ue;

//       compute gradient
         Point<real_t> slope = 0;
         Grad(qU,slope,n);

//       Limiter MP
         ar(1) = slope*_BQv(n,1), ar(2) = slope*_BQv(n,2), ar(3) = slope*_BQv(n,3);
         real_t phi = 1.; //  phi=1 for a second order scheme phi=0 for a first order scheme

//       Very rough limiter!!!! could be improved
         (ar(1)*dU(1)<=0.0) ? phi = 0.0 : phi = min(phi,dU(1)/ar(1));
         (ar(2)*dU(2)<=0.0) ? phi = 0.0 : phi = min(phi,dU(2)/ar(2));
         (ar(3)*dU(3)<=0.0) ? phi = 0.0 : phi = min(phi,dU(3)/ar(3));

//       compute new values
         val1 += phi*slope*_BQv(n,1);
         val2 += phi*slope*_BQv(n,2);
         val3 += phi*slope*_BQv(n,3);
      }
      (sd1->getNeighborElement(1)==the_element) ? LU.set(ns1,dof,val1) : RU.set(ns1,dof,val1);
      (sd2->getNeighborElement(1)==the_element) ? LU.set(ns2,dof,val2) : RU.set(ns2,dof,val2);
      (sd3->getNeighborElement(1)==the_element) ? LU.set(ns3,dof,val3) : RU.set(ns3,dof,val3);
   }
}


void Muscl2DT::GradientM(const Vect<real_t>& U,
                         Vect<real_t>&       LU,
                         Vect<real_t>&       RU,
                         size_t              dof)
{
   LocalVect<real_t,3> qU, dU, aU, ar;
   MESH_EL {
      size_t n = element_label;
      real_t ue = U(n,dof);
      Side *sd1 = The_element.getPtrSide(1); size_t ns1 = sd1->n();
      Side *sd2 = The_element.getPtrSide(2); size_t ns2 = sd2->n();
      Side *sd3 = The_element.getPtrSide(3); size_t ns3 = sd3->n();
      real_t val1=ue, val2=ue, val3=ue;
      if (The_element.isOnBoundary()) {
         if (sd1->isOnBoundary())
            RU.set(ns1,dof,ue);
         if (sd2->isOnBoundary())
            RU.set(ns2,dof,ue);
         if (sd3->isOnBoundary())
            RU.set(ns3,dof,ue);
      }
      else {
         Element *el1=sd1->getNeighborElement(1), 
                 *el2=sd2->getNeighborElement(1), 
                 *el3=sd3->getNeighborElement(1);
         if (el1==the_element)
            el1 = sd1->getNeighborElement(2);
         if (el2==the_element)
            el2 = sd2->getNeighborElement(2);
         if (el3==the_element)
            el3 = sd3->getNeighborElement(2);
         qU(1) = U(el1->n(),dof), qU(2) = U(el2->n(),dof), qU(3) = U(el3->n(),dof);
         dU(1) = qU(1) - ue, dU(2) = qU(2) - ue, dU(3) = qU(3) - ue;

//       compute gradient
         Point<real_t> slope = 0;
         Grad(qU,slope,n);

//       MP limiter
         ar(1) = slope*_BMv(n,1), ar(2) = slope*_BMv(n,2), ar(3) = slope*_BMv(n,3);

//       phi=1 for a second order scheme phi=0 for a first order scheme
         real_t phi = 1.;
         (ar(1)*dU(1)<=0) ? phi = 0 : phi = min(phi,dU(1)/ar(1));
         (ar(2)*dU(2)<=0) ? phi = 0 : phi = min(phi,dU(2)/ar(2));
         (ar(3)*dU(3)<=0) ? phi = 0 : phi = min(phi,dU(3)/ar(3));

//       set the new left and right value
         val1 += phi*slope*_BMv(n,1);
         val2 += phi*slope*_BMv(n,2);
         val3 += phi*slope*_BMv(n,3);
      }
      (sd1->getNeighborElement(1)==the_element) ? LU.set(ns1,dof,val1) : RU.set(ns1,dof,val1);
      (sd2->getNeighborElement(1)==the_element) ? LU.set(ns2,dof,val2) : RU.set(ns2,dof,val2);
      (sd3->getNeighborElement(1)==the_element) ? LU.set(ns3,dof,val3) : RU.set(ns3,dof,val3);
   }
}


real_t Muscl2DT::getMinSize32()
{
   real_t maxS=1e-10, minVsS=1e10;
   MESH_EL {
      for (size_t i=1; i<=3; ++i) {
         Line2 l(The_element.getPtrSide(i));
         if (l.getLength()>maxS)
            maxS = l.getLength();
      }
      Triang3 t(the_element);
      if (t.getArea()/maxS<minVsS)
         minVsS = t.getArea()/maxS;
   }
   return minVsS;
}


real_t Muscl2DT::gettaulim()
{
   real_t taulim = 2.;
   MESH_EL {
      if (The_element.isOnBoundary() == false) {
         size_t n = element_label;
         for (size_t i=1; i<=3; ++i)
            taulim = min(taulim,_BBl(n,i)/_BQl(n,i));
      }
   }
   return taulim;
}


real_t Muscl2DT::getComega()
{
   real_t Comega = 0.;
   LocalMatrix<real_t,3,3> beta;
   MESH_EL {
      size_t n = element_label;
      if (The_element.isOnBoundary()==false && _order(n)) {
         beta(1,1) = -1; beta(1,2) = _beta12(n); beta(1,3) = _beta13(n);
         beta(2,1) = 1./_beta12(n); beta(2,2) = -1; beta(2,3) = -beta(2,1)*_beta13(n);
         beta(3,1) = 1./_beta13(n); beta(3,2) = 1./beta(2,3); beta(3,3) = -1;
         for (size_t j=1; j<=3; ++j) {
            real_t tmp = 0;
            for (size_t k=1; k<=3; ++k)
               if (k!=j)
                  tmp += beta(j,k)/_BBl(n,k);
            tmp *= -_BBl(n,j);
            Comega = max(Comega,tmp);
         }
      }
   }
   return Comega;
}


void Muscl2DT::getgraphComega()
{
   const size_t nb_comega_max = 10000;
   LocalMatrix<real_t,3,3> beta;
   Vect<int> nbcomega(nb_comega_max);
   MESH_EL {
      size_t n = element_label;
      if (The_element.isOnBoundary()==false && _order(n)) {
         beta(1,1) = -1;
         beta(1,2) = _beta12(n);
         beta(1,3) = _beta13(n);
         beta(2,1) = 1./_beta12(n);
         beta(2,2) = -1;
         beta(2,3) = -beta(2,1)*_beta13(n);
         beta(3,1) = 1./_beta13(n);
         beta(3,2) = 1./beta(2,3);
         beta(3,3) = -1;
         for (size_t j=1; j<=3; ++j) {
            real_t tmp = 0;
            for (size_t k=1; k<=3; ++k) {
               if (k!=j)
                  tmp += beta(j,k) / _BBl(n,k);
            }
            tmp *= -_BBl(n,j);
            nbcomega.add(min(max(int(tmp*(nb_comega_max-1)/_Comega),0),int(nb_comega_max)-1)+1,1);
         }
      }
   }
}


void Muscl2DT::print_mesh_stat()
{
   cout << "****************************************************************************" << endl;
   cout << "Mesh statistics :" << endl;
   cout << "* Smin\t\t" << getMinimumFaceArea() << endl;
   cout << "* Smax\t\t" << getMaximumFaceArea() << endl;
   cout << "* Smean\t\t" << getMeanFaceArea() << endl;
   cout << "* Vmin\t\t" << getMinimumElementVolume() << endl;
   cout << "* Vmax\t\t" << getMaximumElementVolume() << endl;
   cout << "* Vmean\t\t" << getMeanElementVolume() << endl;
   cout << "* hmin\t\t" << getMinimumEdgeLength() << endl;
   cout << "* hmax\t\t" << getMaximumEdgeLength() << endl;
   cout << "* V/S m\t\t" << getMinimumVolumebyArea() << endl;
   cout << "*----------------------------------\n" << endl;
   cout << "* hmin\t\t" << getMinimumEdgeLength() << "\t\n sqrt(S)\t" <<  sqrt(getMinimumFaceArea())
        << "\t\n* sqrt3(V) \t" <<  sqrt(getMinimumElementVolume()) << "\t" << endl;
   cout << "****************************************************************************" << endl;
}

} /* namespace OFELI */
