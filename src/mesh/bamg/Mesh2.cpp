// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY: Bamg: Bidimensional Anisotrope Mesh Generator
// RELEASE: 0 
// USAGE  : You may copy freely these files and use it for    
//          teaching or research. These or part of these may   
//          not be sold or used for a commercial purpose with- 
//          out our consent: fax (33) 1 39 63 55 14       
//
// AUTHOR:  F. Hecht,    
// ORG   :  INRIA
// E-MAIL:  Frederic.Hecht@Inria.fr
//
// ORIG-DATE:     Dec 97

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <iostream>
using namespace std;

#include "mesh/bamg/Mesh2.h"
#include "mesh/bamg/QuadTree.h"
#include "mesh/bamg/SetOfE4.h"

namespace bamg {

int Triangles::counter = 0;
Triangles *CurrentTh = 0;
int hinterpole = 1;
long NbUnSwap = 0;
int  ForDebugging = 0;
const Direction NoDirOfSearch = Direction();

inline void MyAssert(int   i,
                     char* ex,
                     char* file,
                     long  line) 
{
   if (i) {
      cerr << "Error Assert:" << ex << " in " << file << " line: " << line << endl;
      exit(1);
   }
}


long AGoodNumberPrimeWith(long n)
{
   const long BigPrimeNumber[] = {567890359L, 567890431L, 567890437L, 567890461L, 567890471L,
                                  567890483L, 567890489L, 567890497L, 567890507L, 567890591L,
                                  567890599L, 567890621L, 567890629L, 0};
  
   long o=0, pi=BigPrimeNumber[1];
   for (int i=0; BigPrimeNumber[i]; i++) {
      long r=BigPrimeNumber[i] % n;
      long oo=Min(Min(r,n-r),Min(Abs(n-2*r),Abs(n-3*r)));
      if (o<oo) 
         o = oo, pi = BigPrimeNumber[i];
   }
   return pi;
}


void MeshError(int Err)
{
   cerr << " Fatal error in the mesh generator " << Err << endl;
   exit(1);
}


ostream& operator <<(ostream&        f,
                     const Triangle& ta)
{
   if (CurrentTh)
      f << "[" << CurrentTh->Number(ta) << "::"
        <<  CurrentTh->Number(ta.ns[0]) << ","
        <<  CurrentTh->Number(ta.ns[1]) << ","
        <<  CurrentTh->Number(ta.ns[2]) << ","
        << "{" <<  CurrentTh->Number(ta.at[0]) << " " << ta.aa[0] << "} "
        << "{" <<  CurrentTh->Number(ta.at[1]) << " " << ta.aa[1] << "} "
        << "{" <<  CurrentTh->Number(ta.at[2]) << " " << ta.aa[2] << "} ]";
   else
      f << "[" 
        << ta.ns[0] << "," << ta.ns[1] << "," << ta.ns[2] << "," 
        << "{" << ta.at[0] << " " << ta.aa[0] << "} " 
        << "{" << ta.at[1] << " " << ta.aa[1] << "} " 
        << "{" << ta.at[2] << " " << ta.aa[2] << "} ]";
   return f;
}


void swap(Triangle* t1,
          short     a1,
          Triangle* t2,
          short     a2,
          Vertex*   s1,
          Vertex*   s2,
          Icoor2    det1,
          Icoor2    det2)
{
/*--------------------------------------------------------------
     short a2=aa[a];// les 2 numero de l arete dans les 2 triangles

                   sb                     sb    
                 / | \                   /   \     
             as1/  |  \                 /a2   \    
               /   |   \               /    t2 \   
           s1 /t1  | t2 \s2  -->   s1 /___as2___\s2
              \  a1|a2  /             \   as1   /  
               \   |   /               \ t1    /   
                \  |  / as2             \   a1/    
                 \ | /                   \   /     
                  sa                       sa   
  -------------------------------------------------------------*/

   int as1 = NextEdge[a1], as2 = NextEdge[a2];
   int ap1 = PreviousEdge[a1], ap2 = PreviousEdge[a2];
   (*t1)(VerticesOfTriangularEdge[a1][1]) = s2;  // before sb
   (*t2)(VerticesOfTriangularEdge[a2][1]) = s1;  // before sa
// update both external adjacencies 
   TriangleAdjacent taas1=t1->Adj(as1),
   taas2 = t2->Adj(as2),
   tas1(t1,as1), tas2(t2,as2),
   ta1(t1,a1),ta2(t2,a2);
// external top left
   taas1.SetAdj2(ta2,taas1.GetAllFlag_UnSwap());
// externe bas droite
   taas2.SetAdj2(ta1,taas2.GetAllFlag_UnSwap());
// remove the Mark UnMarkSwap 
   t1->SetUnMarkUnSwap(ap1);
   t2->SetUnMarkUnSwap(ap2);
// internal
   tas1.SetAdj2(tas2);
   t1->det = det1;
   t2->det = det2;
   t1->SetTriangleContainingTheVertex();
   t2->SetTriangleContainingTheVertex();
}


long FindTriangle(Triangles& Th,
                  double     x,
                  double     y,
                  double*    a,
                  int&       inside)
{
   CurrentTh = &Th;
//   assert(&Th);
   I2 I = Th.toI2(R2(Min(Max(Th.pmin.x,x),Th.pmax.x),Min(Max(Th.pmin.y,y),Th.pmax.y))); 
   Icoor2 dete[3];
   Triangle &tb = *Th.FindTriangleContaining(I,dete);
   
   if (tb.link) { // internal point in a true triangle
      a[0] = double(dete[0]) / tb.det;
      a[1] = double(dete[1]) / tb.det;
      a[2] = double(dete[2]) / tb.det;
      inside = 1;
      return Th.Number(tb);
   }
   else  {
      inside = 0; 
      double aa, bb;
      TriangleAdjacent ta=CloseBoundaryEdgeV2(I,&tb,aa,bb);	 
      int k = ta;
      Triangle *tc = ta;
      if (!tc->link) {
         ta = ta.Adj();
         tc = ta;
         k = ta;
         Exchange(aa,bb);
         assert(tc->link);
      }
      a[VerticesOfTriangularEdge[k][0]] = aa;
      a[VerticesOfTriangularEdge[k][1]] = bb;
      a[OppositeVertex[k]] = 1 - aa - bb;
      return Th.Number(tc);
   }
}


TriangleAdjacent CloseBoundaryEdge(I2        A,
                                   Triangle* t,
                                   double&   a,
                                   double&   b)
{
   int k=(*t)(0) ? (((*t)(1) ? ((*t)(2) ? -1 : 2) : 1)) : 0;
   int dir=0;
   assert(k>=0);
   Icoor2 IJ_IA, IJ_AJ;
   TriangleAdjacent edge(t,OppositeEdge[k]);          
   for (;;edge = dir >0 ? Next(Adj(Next(edge))) : Previous(Adj(Previous(edge)))) {
      Vertex &vI=*edge.EdgeVertex(0);
      Vertex &vJ=*edge.EdgeVertex(1);
      I2 I=vI, J=vJ, IJ= J-I;
      IJ_IA = (IJ,(A-I));
      if (IJ_IA<0) {
         if (dir>0) {
            a = 1;
            b = 0;
            return edge;
         }
         else {
            dir = -1;
            continue;
         }
      }
      IJ_AJ = (IJ,(J-A));
      if (IJ_AJ<0) {
         if (dir<0) {
            a = 0;
            b = 1;
            return edge;
         }
         else {
            dir = 1;
            continue;
         }
      }
      double IJ2 = IJ_IA + IJ_AJ;
      assert(IJ2);
      a = IJ_AJ/IJ2;
      b = IJ_IA/IJ2;
      return edge;
   }
}


TriangleAdjacent Triangle::FindBoundaryEdge(int i) const
{
   Triangle *t=(Triangle *)this, *ttc;
   int k=0, j=EdgesVertexTriangle[i][0], jc, ext=!link;
   do {
      int exterieurp=ext;
      k++; 
      ttc = t->at[j];
      ext = !ttc->link;
      if (ext+exterieurp==1) 
         return TriangleAdjacent(t,j);
      jc = NextEdge[t->aa[j]&3];
      t = ttc;
      j = NextEdge[jc];
      assert(k<2000);
   } while ((this!=t)); 
   return TriangleAdjacent(0,0);
}


TriangleAdjacent CloseBoundaryEdgeV2(I2        C,
                                     Triangle* t,
                                     double&   a,
                                     double&   b) 
{
   assert(t->link==0);
   Vertex *s=0, *s1=0, *s0=0;
   Icoor2 imax=MaxICoor22, l0=imax, l1=imax;
   double dd2=imax;
   TriangleAdjacent er;
   for (int j=0; j<3; j++) { 
      TriangleAdjacent ta=t->FindBoundaryEdge(j);
      if (!(Triangle *) ta)
         continue;
      s0 = ta.EdgeVertex(0);
      s1 = ta.EdgeVertex(1);
      I2 A = *s0, B = *ta.EdgeVertex(1);
      I2 AB = B-A, AC=C-A, BC=B-C;
      Icoor2 ACAC=(AC,AC), BCBC=(BC,BC);
      Icoor2 AB2 = Norme2_2(AB); //  ||AB||^2
      Icoor2 ABAC = (AB,AC);         //  AB.AC|
      double d2;
      if (ABAC<0) {
         if ((d2=(double)ACAC) < dd2) {
            er = ta;
            l0 = ACAC;
            l1 = BCBC;
            s = s0;
         }
      }
      else if (ABAC>AB2) {
         if ((d2=double(BCBC)) < dd2) {
            dd2 = d2;
            er = Adj(ta);
            l0 = BCBC;
            l1 = ACAC;
            s = s1;
         }
      }
      else {
         double det_2=double(Det(AB,AC)); 
         det_2 *= det_2; // square of area*2 of triangle ABC
         d2 = det_2/double(AB2); // hauteur^2 in C of of triangle ABC
         if (d2<dd2) {
            dd2 = d2;
            er = ta;
            l0 = (AC,AC);
            l1 = (BC,BC);
            s = 0;
            b = (double(ABAC)/double(AB2));
            a = 1 - b;
         }
      }
   }
   if (s) {
      t = er;
      TriangleAdjacent edge(er);
      int kkk=0;
      int linkp = t->link == 0;
      Triangle *tt=t=edge=Adj(Previous(edge));
      do {
         assert(edge.EdgeVertex(0)==s && kkk<10000);
         kkk++;
         int link = tt->link == 0;
         if ((link + linkp) == 1) {
            Vertex *st = edge.EdgeVertex(1);
            I2 I=*st;
            Icoor2 ll = Norme2_2 (C-I);
            if (ll < l1) {  // the other vertex is neart 
               s1 = st;
               l1 = ll;
               er = edge;
               if (ll<l0) { // change of direction --
                  s1 = s;
                  l1 = l0;
                  s = st;
                  l0 = ll;
                  t = tt;
                  edge = Adj(edge);
                  link = linkp;
                  er = edge;
               }
            }
         }
         linkp = link;
         edge = Adj(Previous(edge));
         tt = edge;
      } while (t!=tt);

      assert((Triangle *) er);
      I2 A((I2)*er.EdgeVertex(0));
      I2 B((I2)*er.EdgeVertex(1));
      I2 AB=B-A, AC=C-A, CB=B-C;
      double aa=double((AB,AC)), bb=double((AB,CB));
      if (aa<0)
         a = 1, b = 0;
      else if (bb<0)
         a = 0, b = 1;
      else {
         a = bb/(aa+bb);
         b = aa/(aa+bb);
      }
   }
   return er;
}


Metric Triangles::MetricAt(const R2& A) const
{
   I2 a = toI2(A);
   Icoor2 deta[3];
   Triangle *t=FindTriangleContaining(a,deta);
   if (t->det <0) {
      double ba, bb;
      TriangleAdjacent edge = CloseBoundaryEdge(a,t,ba,bb) ;
      return Metric(ba,*edge.EdgeVertex(0),bb,*edge.EdgeVertex(1));
   }
   else {
      double aa[3];
      double s = deta[0]+deta[1]+deta[2];
      aa[0] = deta[0]/s;
      aa[1] = deta[1]/s;
      aa[2] = deta[2]/s;
      return Metric(aa,(*t)[0],(*t)[1],(*t)[2]);
   }
}


void ListofIntersectionTriangles::SplitEdge(const Triangles& Bh,
                                            const R2&        A,
                                            const R2&        B,
                                            int              nbegin)
{
   Triangle *tbegin, *t;
   Icoor2 deta[3], deti=0, detj=0;
   double ba[3];
   int nbt=0, ifirst=-1, ilast;
   int i0, i1, i2;
   int ocut=0, i, j, k=-1;
   Icoor2 dt[3];
   I2 a=Bh.toI2(A), b=Bh.toI2(B), vi, vj;  
   int iedge=-1;

   if (nbegin) {
      t = tbegin = lIntTria[ilast=(Size-1)].t;
      if (tbegin->det>=0) 
         ifirst = ilast;
   }
   else {// not optimisation 
      init();
      t = tbegin = Bh.FindTriangleContaining(a,deta);
      if (t->det>=0)
         ilast = NewItem(t,double(deta[0])/t->det,double(deta[1])/t->det,double(deta[2])/t->det);
      else {
//       find the nearest boundary edge  of the vertex A
//       find a edge or such normal projection a the edge IJ is on the edge
//   <=> IJ.IA >=0 && IJ.AJ >=0
         ilast = ifirst;
         double ba, bb;
         TriangleAdjacent edge=CloseBoundaryEdge(a,t,ba,bb);
         Vertex &v0 = *edge.EdgeVertex(0), &v1 = *edge.EdgeVertex(1);
         NewItem(A,Metric(ba,v0,bb,v1));
         t = edge;
//       test if the point b is in the same side
         if (det(v0.i,v1.i,b)>=0) {
            TriangleAdjacent edge=CloseBoundaryEdge(a,t,ba,bb);
	    Vertex &v0 = *edge.EdgeVertex(0), &v1 = *edge.EdgeVertex(1);
            NewItem(A,Metric(ba,v0,bb,v1));
            return;
         }
      }
   }
   if (t->det<0) {  // outside departure
      while (t->det <0) { // intersection boundary edge and a,b,
         k = (*t)(0) ? (((*t)(1) ? ((*t)(2) ? -1 : 2) : 1)) : 0;
         assert(k>=0);
         ocut = OppositeEdge[k];
         i = VerticesOfTriangularEdge[ocut][0];
         j = VerticesOfTriangularEdge[ocut][1];
         vi = (*t)[i];
         vj = (*t)[j];
         deti = bamg::det(a,b,vi);
         detj = bamg::det(a,b,vj);
         if (deti>0) // go to  i direction on gamma
            ocut = PreviousEdge[ocut];      
         else if (detj<=0) // go to j direction on gamma
            ocut = NextEdge[ocut];         
         TriangleAdjacent tadj=t->Adj(ocut);
         t = tadj;
         iedge = tadj;
         if (t == tbegin) { 
            double ba, bb;
            if (verbosity>7)
               cout << "       SplitEdge: All the edge " << A << B << nbegin << det(vi,vj,b) 
                    << " deti= " << deti <<  " detj=" << detj << endl;
            TriangleAdjacent edge=CloseBoundaryEdge(a,t,ba,bb);
            Vertex &v0 = *edge.EdgeVertex(0), &v1 = *edge.EdgeVertex(1);
            NewItem(A,Metric(ba,v0,bb,v1));
            return;
         }
      }
      if (det(vi,vj,b)>=0) {
         if (verbosity>7)
            cout << "       SplitEdge: all AB outside " << A << B << endl;
         t = tbegin;
         double ba, bb;
         TriangleAdjacent edge=CloseBoundaryEdge(b,t,ba,bb);
         NewItem(B,Metric(ba,*edge.EdgeVertex(0),bb,*edge.EdgeVertex(1)));
         return;
      }
      else {
         k = OppositeVertex[iedge];
         i = VerticesOfTriangularEdge[iedge][0];
         j = VerticesOfTriangularEdge[iedge][1];
         double dij = detj-deti;
         assert(i+j+k == 0+1+2);
         ba[j] =  detj/dij;
         ba[i] = -deti/dij;
         ba[k] = 0;
         ilast = NewItem(t,ba[0],ba[1],ba[2]);
      }
   }

   for (;;) {
      if (iedge < 0) {
         i0 = 0; i1 = 1; i2 = 2;
         dt[0] = bamg::det(a,b,(*t)[0]);
         dt[1] = bamg::det(a,b,(*t)[1]);
         dt[2] = bamg::det(a,b,(*t)[2]);
      }
      else {
         i2 = iedge;
         i0 = NextEdge[i2];
         i1 = NextEdge[i0];
//       We revert i, j because we take the triangle by the other side
         dt[VerticesOfTriangularEdge[iedge][0]] = detj;
         dt[VerticesOfTriangularEdge[iedge][1]] = deti;
         dt[iedge] = det(a,b,(*t)[OppositeVertex[iedge]]);
      }
      if ((dt[i=VerticesOfTriangularEdge[i0][0]] < 0) &&
          (dt[j=VerticesOfTriangularEdge[i0][1]] > 0))
         ocut = i0;
      else if ((dt[i=VerticesOfTriangularEdge[i1][0]] < 0) &&
               (dt[j=VerticesOfTriangularEdge[i1][1]] > 0))
         ocut = i1;
      else if ((dt[i=VerticesOfTriangularEdge[i2][0]] < 0) && 
               (dt[j=VerticesOfTriangularEdge[i2][1]] > 0))
         ocut = i2;
      else if ((dt[i=VerticesOfTriangularEdge[i0][0]] == 0) &&
               (dt[j=VerticesOfTriangularEdge[i0][1]] >  0))
         ocut = i0;
      else if ((dt[i=VerticesOfTriangularEdge[i1][0]] == 0) &&
               (dt[j=VerticesOfTriangularEdge[i1][1]] >  0))
         ocut = i1;
      else if ((dt[i=VerticesOfTriangularEdge[i2][0]] == 0) && 
               (dt[j=VerticesOfTriangularEdge[i2][1]] >  0))
         ocut = i2;
      else if ((dt[i=VerticesOfTriangularEdge[i0][0]] <  0) &&
             (dt[j=VerticesOfTriangularEdge[i0][1]] == 0))
         ocut = i0;
      else if ((dt[i=VerticesOfTriangularEdge[i1][0]] <  0) &&
               (dt[j=VerticesOfTriangularEdge[i1][1]] == 0))
         ocut = i1;
      else if ((dt[i=VerticesOfTriangularEdge[i2][0]] <  0) && 
               (dt[j=VerticesOfTriangularEdge[i2][1]] == 0))
         ocut = i2;
      else {
         k = 0;
         if (dt[0])
            ocut = 0, k++; 
         if (dt[1])
            ocut = 1, k++; 
         if (dt[2])
            ocut = 2, k++;
         if (k==1) {
            if (dt[ocut]>0) // triangle upper AB
               ocut = NextEdge[ocut];
            i = VerticesOfTriangularEdge[ocut][0];
            j = VerticesOfTriangularEdge[ocut][1];
         }
         else {
            cerr << " Bug Split Edge " << endl;
            cerr << " dt[0] = " << dt[0] 
                 << " dt[1] = " << dt[1] 
                 << " dt[2] = "<< dt[2] << endl;
            cerr << i0 << " " << i1 << " " << i2 << endl;
            cerr << " A = " << A << " B= " << B << endl;
            cerr << " Triangle t = " <<  *t << endl;
            cerr << (*t)[0] << (*t)[1] << (*t)[0] << endl;
            cerr << " nbt = " << nbt << endl;
            MeshError(100);
         }
      }
      k = OppositeVertex[ocut];
      Icoor2 detbij = bamg::det((*t)[i],(*t)[j],b);
      if (detbij >= 0) { //we find the triangle containing b
         dt[0] = bamg::det((*t)[1],(*t)[2],b);
         dt[1] = bamg::det((*t)[2],(*t)[0],b);
         dt[2] = bamg::det((*t)[0],(*t)[1],b);
         double dd = t->det;
         NewItem(t,dt[0]/dd,dt[1]/dd,dt[2]/dd);      
         return;
      }
      else { // next triangle by  adjacent by edge ocut 
         deti = dt[i];
         detj = dt[j];
         double dij=detj-deti;
         ba[i] = detj/dij;
         ba[j] = -deti/dij;
         ba[3-i-j] = 0;
         ilast = NewItem(t,ba[0],ba[1],ba[2]);
         TriangleAdjacent ta = t->Adj(ocut);
         t = ta;
         iedge = ta; 
         if (t->det <= 0)  {
            double ba, bb;
            TriangleAdjacent edge=CloseBoundaryEdge(b,t,ba,bb);
            NewItem(B,Metric(ba,*edge.EdgeVertex(0),bb,*edge.EdgeVertex(1)));
            return;
         }
      } 
   }
}


int ListofIntersectionTriangles::NewItem(Triangle* tt,
                                         double    d0,
                                         double    d1,
                                         double    d2)
{
   int n;
   R2 x(0,0);
   if (d0)
      x = (*tt)[0].r * d0;
   if (d1)
      x = x + (*tt)[1].r * d1;
   if (d2)
      x = x + (*tt)[2].r * d2;

// newer add same point 
   if (!Size || Norme2_2(lIntTria[Size-1].x-x)) {
      if (Size==MaxSize)
         ReShape();
      lIntTria[Size].t = tt;
      lIntTria[Size].bary[0] = d0;
      lIntTria[Size].bary[1] = d1;
      lIntTria[Size].bary[2] = d2;
      lIntTria[Size].x = x;
      Metric m0, m1, m2;
      Vertex* v;
      if ((v=(*tt)(0)))
         m0 = v->m;
      if ((v=(*tt)(1)))
         m1 = v->m;
      if ((v=(*tt)(2)))
         m2 = v->m;
      lIntTria[Size].m = Metric(lIntTria[Size].bary,m0,m1,m2);
      n = Size++;
   }
   else
      n = Size - 1;
   return n;
}


int ListofIntersectionTriangles::NewItem(R2            A,
                                         const Metric& mm)
{
   int n;
   if (!Size || Norme2_2(lIntTria[Size-1].x-A)) {
      if (Size==MaxSize)
         ReShape();
      lIntTria[Size].t = 0;
      lIntTria[Size].x = A;
      lIntTria[Size].m = mm;
      n = Size++;
   }
   else
      n = Size - 1;
   return n;
}


double ListofIntersectionTriangles::Length()
{
   assert(Size>0);
   Metric Mx,My;
   int ii, jj;
   R2 C, x, y, xy;
   SegInterpolation *SegI=lSegsI;
   SegI = lSegsI;
   lSegsI[NbSeg].last = Size;// improvement
   int EndSeg=Size;
   y = lIntTria[0].x;
   double sxy, s=0;
   lIntTria[0].s = 0;
   SegI->lBegin = s;
   for (jj=0, ii=1; ii<Size; jj=ii++) {  
      x=y;
      y = lIntTria[ii].x;
      xy = y-x;
      Mx = lIntTria[ii].m;
      My = lIntTria[jj].m;
      sxy =  LengthInterpole(Mx,My,xy);
      s += sxy;
      lIntTria[ii].s = s;
      if (ii == EndSeg) 
         SegI->lEnd=s,
         SegI++,
         EndSeg = SegI->last,
         SegI->lBegin = s;
   }
   len = s;
   SegI->lEnd = s;
   return s;
}


long ListofIntersectionTriangles::NewPoints(Vertex* vertices,
                                            long&   nbv,
                                            long    nbvx)
{
   const long nbvold=nbv;
   double s=Length();
   if (s<1.5)
      return 0;
   int ii=1;
   R2 y, x;
   Metric My, Mx;
   double sx=0, sy;
   int nbi=Max(2,int(s+0.5));
   double sint=s/nbi;
   double si=sint;
   int EndSeg=Size;
   SegInterpolation *SegI=0;
   if (NbSeg)
      SegI=lSegsI, EndSeg=SegI->last;
   for (int k=1; k<nbi; k++) {
      while ((ii<Size) && (lIntTria[ii].s<=si)) 
         if (ii++ == EndSeg) 
            SegI++, EndSeg=SegI->last;
      int ii1=ii-1;
      x = lIntTria[ii1].x;
      sx = lIntTria[ii1].s;
      Metric Mx=lIntTria[ii1].m;
      y = lIntTria[ii].x;
      sy = lIntTria[ii].s;
      Metric My=lIntTria[ii].m;
      double lxy=sy-sx;
      double cy=abscisseInterpole(Mx,My,y-x,(si-sx)/lxy);
      R2 C;
      double cx = 1-cy;
      C = SegI ? SegI->F(si): x * cx + y *cy;
      si += sint;
      if (nbv<nbvx) {
         vertices[nbv].r = C;
         vertices[nbv++].m = Metric(cx,lIntTria[ii-1].m,cy,lIntTria[ii].m);
         if ((verbosity/100%10)==2)
            cout << "   -- Add point " << nbv-1 << " " << vertices[nbv-1] << " "
                 << vertices[nbv-1].m << endl;
      }
      else
         return nbv-nbvold;
   }
   return nbv-nbvold;
}


int SwapForForcingEdge(Vertex* &         pva,
                       Vertex* &         pvb,
                       TriangleAdjacent& tt1,
                       Icoor2&           dets1,
                       Icoor2&           detsa,
                       Icoor2&           detsb,
                       int&              NbSwap)
{ // l'arete ta coupe l'arete pva pvb
  // de cas apres le swap sa coupe toujours
  // on cherche l'arete suivante 
  // on suppose que detsa >0 et detsb <0
  // attention la routine echange pva et pvb 

   if (tt1.Locked())
      return 0; // frontiere croise 
   TriangleAdjacent tt2 = Adj(tt1);
   Triangle *t1=tt1, *t2=tt2;// les 2 triangles adjacent
   short a1=tt1, a2=tt2;// les 2 numero de l arete dans les 2 triangles
   assert(a1>=0 && a1<3);

   Vertex& sa = (*t1)[VerticesOfTriangularEdge[a1][0]];
   Vertex& s1 = (*t1)[OppositeVertex[a1]];
   Vertex& s2 = (*t2)[OppositeVertex[a2]];
   Icoor2 dets2 = det(*pva,*pvb,s2);
   Icoor2 det1 = t1->det, det2 = t2->det;
   Icoor2 detT = det1 + det2;
   assert(det1>0 && det2>0);
   assert(detsa<0 && detsb>0); // [a,b] cut infinite line va,bb
   Icoor2 ndet1=bamg::det(s1,sa,s2);
   Icoor2 ndet2=detT-ndet1;
   int ToSwap=0;             // no swap
   if (ndet1>0 && ndet2>0) { // we can swap  
      if ((dets1<=0 && dets2<=0) || (dets2>=0 && detsb>=0))
         ToSwap = 1;
      else 
         if (BinaryRand()) 
            ToSwap = 2; 
   }
   if (ToSwap)
      NbSwap++;
   bamg::swap(t1,a1,t2,a2,&s1,&s2,ndet1,ndet2);

   int ret=1;
   if (dets2<0) {
      dets1 = ToSwap ? dets1 : detsa;
      detsa = dets2;
      tt1 =  Previous(tt2);
   }
   else if (dets2>0) {
      dets1 = ToSwap ? dets1 : detsb;
      detsb = dets2;
      if (!ToSwap)
         tt1 = Next(tt2);
   }
   else { 
      if (ForDebugging)
         cout << "changement de sens" << endl;
      ret = -1;
      Exchange(pva,pvb);
      Exchange(detsa,detsb);
      Exchange(dets1,dets2);
      Exchange(tt1,tt2);
      dets1 = -dets1;
      dets2 = -dets2;
      detsa = -detsa;
      detsb = -detsb;
      if (ToSwap) { 
         if (dets2<0) {
            dets1 = (ToSwap ? dets1 : detsa);
            detsa = dets2;
            tt1 = Previous(tt2);
         }
         else if (dets2>0){
            dets1 = (ToSwap ? dets1 : detsb);
            detsb = dets2;
            if (!ToSwap)
               tt1 = Next(tt2);
         }
         else {
            tt1 = Next(tt2);
            ret = 0;
         }
      }
   }
   return ret;
}


int ForceEdge(Vertex&           a,
              Vertex&           b,
              TriangleAdjacent& taret) 
{
   int NbSwap=0;
   assert(a.t && b.t);
   int k=0;
   taret = TriangleAdjacent(0,0);
   TriangleAdjacent tta(a.t,EdgesVertexTriangle[a.vint][0]);
   Vertex *v2=tta.EdgeVertex(0), *vbegin =v2;
   Icoor2 det2 = v2 ? det(*v2,a,b) : -1, det1;
   if (v2) 
      det2 = det(*v2,a,b);
   else {
      tta = Previous(Adj(tta));
      v2 = tta.EdgeVertex(0);
      vbegin = v2;
      assert(v2);
      det2 = det(*v2,a,b);
   }

   while (v2 != &b) {
      TriangleAdjacent tc = Previous(Adj(tta));    
      v2 = tc.EdgeVertex(0);
      det1 = det2;
      det2 =  v2 ? det(*v2,a,b) : det2;
      if (det1<0 && det2>0) {
         Vertex *va=&a, *vb=&b;
         tc = Previous(tc);
         Icoor2 detss=0, l=0, ks;
         while ((ks=SwapForForcingEdge(va,vb,tc,detss,det1,det2,NbSwap))) {
            if (l++ > 1000000) {
               cerr << " Loop in forcing Egde AB" 
                    << "\n vertex A " << a
                    << "\n vertex B " <<  b
                    << "\n nb de swap " << NbSwap
                    << "\n nb of try swap too large = " <<  l << " greater than " <<  100000 << endl;
               if (CurrentTh)
                  cerr << " vertex number " << CurrentTh->Number(a) << " " << CurrentTh->Number(b) << endl;
               MeshError(990);
            }
         }
         Vertex *aa=tc.EdgeVertex(0), *bb=tc.EdgeVertex(1);
         if (((aa==&a) && (bb==&b)) || ((bb==&a) && (aa==&b))) {
            tc.SetLock();
            a.Optim(1,0);
            b.Optim(1,0);
            taret = tc;
            return NbSwap;
         }
      }
      tta = tc;
      assert(k<2000);
      k++;
      if (vbegin==v2)
         return -1;
   }
   tta.SetLock();
   taret = tta;
   a.Optim(1,0);
   b.Optim(1,0);
   return NbSwap; 
}


int Triangle::swap(short a,
                   int   koption)
{
   if (a/4 !=0)
      return 0;

   Triangle *t1=this, *t2=at[a];// les 2 triangles adjacent
   short a1=a, a2=aa[a];         // les 2 numero de l arete dans les 2 triangles
   if (a2/4!=0)
      return 0; // arete lock or MarkUnSwap
   Vertex *sa=t1->ns[VerticesOfTriangularEdge[a1][0]];
   Vertex *sb=t1->ns[VerticesOfTriangularEdge[a1][1]];
   Vertex *s1=t1->ns[OppositeVertex[a1]];
   Vertex *s2=t2->ns[OppositeVertex[a2]];

   Icoor2 det1=t1->det, det2=t2->det;
   Icoor2 detT=det1+det2;
   Icoor2 detA=Abs(det1) + Abs(det2);
   Icoor2 detMin=Min(det1,det2);

   int OnSwap=0;
// si 2 triangle infini (bord) => detT = -2;
   if (sa==0) {
      det2 = bamg::det(s2->i,sb->i,s1->i);
      OnSwap = det2>0;
   }
   else if (sb==0) {
      det1 = bamg::det(s1->i,sa->i,s2->i);
      OnSwap = det1>0;
   }
   else if ((s1!=0) && (s2!=0)) {
      det1 = bamg::det(s1->i,sa->i,s2->i);
      det2 = detT - det1;
      OnSwap = (Abs(det1) + Abs(det2)) < detA;
      Icoor2 detMinNew=Min(det1,det2);
      if (!OnSwap && detMinNew>0) {
         OnSwap = detMin==0;
         if (!OnSwap) {
	    int kopt = koption;
            while (1)
               if (kopt) {
//                critere de Delaunay pure isotrope
                  Icoor2 xb1 = sb->i.x - s1->i.x,
                  x21 = s2->i.x - s1->i.x,
                  yb1 = sb->i.y - s1->i.y,
                  y21 = s2->i.y - s1->i.y,
                  xba = sb->i.x - sa->i.x, 
                  x2a = s2->i.x - sa->i.x,
                  yba = sb->i.y - sa->i.y,
                  y2a = s2->i.y - sa->i.y;
                  double cosb12 = double(xb1*x21 + yb1*y21),
                                  cosba2 = double(xba*x2a + yba*y2a) ,
                                  sinb12 = double(det2),
                                  sinba2 = double(t2->det);
                  OnSwap = (double(cosb12)*double(sinba2)) < (double(cosba2)*double(sinb12));
                  break;
               }
               else {
//                critere de Delaunay anisotrope 
                  double som;
                  I2 AB=(I2)*sb-(I2)*sa, MAB2=((I2)*sb+(I2)*sa);
                  R2 MAB(MAB2.x*0.5,MAB2.y*0.5);
                  I2 A1=(I2)*s1-(I2)*sa, D=(I2)*s1-(I2)*sb;
                  R2 S2(s2->i.x,s2->i.y), S1(s1->i.x,s1->i.y);
                  {
                     Metric M=s1->m;
                     R2 ABo=M.Orthogonal(AB), A1o=M.Orthogonal(A1);
                     double dd=Abs(ABo.x*A1o.y)+Abs(ABo.y*A1o.x);
                     double d=(ABo.x*A1o.y - ABo.y*A1o.x)*2; // because D/2
                     if (Abs(d) > dd*1.e-3) {
                        R2 C(MAB+ABo*((D.x*A1o.y - D.y*A1o.x)/d));
                        som = M(C-S2)/M(C-S1);
                     }
                     else {
                        kopt = 1;
                        continue;
                     }	
	          }
                  {
                     Metric M=s2->m;
                     R2 ABo = M.Orthogonal(AB), A1o = M.Orthogonal(A1);
                     double dd = Abs(ABo.x*A1o.y)+Abs(ABo.y*A1o.x);
                     double d = (ABo.x*A1o.y - ABo.y*A1o.x)*2; // because D/2
                     if (Abs(d) > dd*1.e-3) {
                        R2 C(MAB+ABo*((D.x*A1o.y-D.y*A1o.x)/d));
                        som += M(C-S2)/M(C-S1);
                     }
                     else {
                        kopt = 1;
                        continue;
                     }
                  }
                  OnSwap = som < 2;
                  break;
               }
            } 
         }
      }
      if (OnSwap) 
         bamg::swap(t1,a1,t2,a2,s1,s2,det1,det2);
      else {
         NbUnSwap++;
         t1->SetMarkUnSwap(a1);     
      }
      return OnSwap;
}


double Vertex::Smoothing(Triangles&       Th,
                         const Triangles& BTh,
                         Triangle*        &tstart,
                         double           omega)
{
   Vertex *s=this;
   Vertex &vP=*s, vPsave=vP;
   Triangle *tbegin=t, *tria=t, *ttc;
   int k=0, kk=0, j=EdgesVertexTriangle[vint][0], jc;
   R2 P(s->r), PNew(0,0);
   do {
      k++;
      if (!tria->Hidden(j)) {
         Vertex &vQ=(*tria)[VerticesOfTriangularEdge[j][0]];
         R2 Q=vQ, QP(P-Q);
         double lQP = LengthInterpole(vP,vQ,QP);
         PNew += Q + QP/Max(lQP,1e-20);
         kk++;
      }
      ttc =  tria->TriangleAdj(j);
      jc = NextEdge[tria->NuEdgeTriangleAdj(j)];
      tria = ttc;
      j = NextEdge[jc];
      assert(k<2000);
   } while (tbegin!=tria); 
   if (kk<4)
      return 0;
   PNew = PNew/double(kk);
   R2 Xmove((PNew-P)*omega);
   PNew = P + Xmove;
   double delta=Norme2_2(Xmove); 
   Icoor2 deta[3];
   I2 IBTh = BTh.toI2(PNew);
   tstart = BTh.FindTriangleContaining(IBTh,deta,tstart);  
   if (tstart->det<0) { // outside
      double ba, bb;
      TriangleAdjacent edge=CloseBoundaryEdge(IBTh,tstart,ba,bb);
      tstart = edge;
      vP.m = Metric(ba,*edge.EdgeVertex(0),bb,*edge.EdgeVertex(1));
   }
   else { // inside
      double aa[3];
      double s=deta[0]+deta[1]+deta[2];
      aa[0] = deta[0]/s;
      aa[1] = deta[1]/s;
      aa[2] = deta[2]/s;
      vP.m = Metric(aa,(*tstart)[0],(*tstart)[1],(*tstart)[2]);
   }

// recompute the det of the triangle
   vP.r = PNew;
   vP.i = Th.toI2(PNew);
   Vertex vPnew=vP;
   int ok=1, loop=1;
   k = 0;
   while (ok) {
      ok = 0;
      do {
         k++;
         double detold=tria->det;
         tria->det = bamg::det((*tria)[0],(*tria)[1],(*tria)[2]);
         if (loop) {
            Vertex *v0, *v1, *v2, *v3;
            if (tria->det<0)
               ok = 1;			       
            else if (tria->Quadrangle(v0,v1,v2,v3)) {
               vP = vPsave;
               double qold=QuadQuality(*v0,*v1,*v2,*v3);
               vP = vPnew;
               double qnew=QuadQuality(*v0,*v1,*v2,*v3);
               if (qnew<qold)
                  ok = 1;
            }
            else if ((double)tria->det < detold/2)
               ok = 1;
         }
         tria->SetUnMarkUnSwap(0);
         tria->SetUnMarkUnSwap(1);
         tria->SetUnMarkUnSwap(2);
         ttc = tria->TriangleAdj(j);
         jc = NextEdge[tria->NuEdgeTriangleAdj(j)];
         tria = ttc;
         j = NextEdge[jc];
         assert(k<2000);
      } while (tbegin != tria); 
      if (ok && loop)
         vP = vPsave; // no move 
      loop = 0;
   }
   return delta;
}


void Triangles::Add(Vertex&   s,
                    Triangle* t,
                    Icoor2*   det3) 
{
/*-------------------------------------------
                 s2

                 /|\
                / | \
               /  |  \
        tt1   /   |   \ tt0
             /    |s   \
            /     .     \
           /  .      `   \
          / .           ` \
          ----------------
       s0       tt2       s1
  --------------------------------------------*/

   Vertex *s0 = &(*t)[0], *s1 = &(*t)[1], *s2 = &(*t)[2];
   Icoor2 det3local[3] = { Icoor2(0.), Icoor2(0.), Icoor2(0.) };

// The following lines are added to avoid a bug in the Mac Os compiler when
// an optimization option is used
   static int qq=0;
   qq++;
   if (qq==0)
      s0 = NULL;

   int infv = s0 ? ((s1?(s2?-1:2):1)) : 0;
   int nbd0 = 0, izerodet = -1, iedge = 0;
   Icoor2 detOld = t->det;

   if ((infv<0 && detOld<0) || (infv>=0 && detOld>0)) {
      cout << "infv = " << infv << ", det = " << detOld << endl;
      cout << Number(s) << " " << Number(*s0) << " "  
           << Number(*s1) << " " << Number(*s2) << endl;
      MeshError(3);
   }

   if (!det3) {
      det3 = det3local;
      if (infv<0) {
         det3[0] = bamg::det(s,*s1,*s2);
         det3[1] = bamg::det(*s0,s,*s2);
         det3[2] = bamg::det(*s0,*s1,s);
      }
      else { 
//       one of &s1, &s2, &s0 is NULL so (&si || &sj) <=> !&sk
         det3[0] = s0 ? -1 : bamg::det(s,*s1,*s2);
         det3[1] = s1 ? -1 : bamg::det(*s0,s,*s2);
         det3[2] = s2 ? -1 : bamg::det(*s0,*s1,s);
      }
   }

   if (!det3[0])
      izerodet = 0, nbd0++;
   if (!det3[1])
      izerodet = 1, nbd0++;
   if (!det3[2])
      izerodet = 2, nbd0++;

   if (nbd0>0) { // point s on an egde or on a vertex 
      if (nbd0==1) {
         iedge = OppositeEdge[izerodet];
         TriangleAdjacent ta=t->Adj(iedge);
//       the point is on the edge if the point is on the boundary 
//       add the point in outside part
         if (t->det>=0) { // inside triangle
            if (((Triangle *)ta)->det<0) {
//             add in outside triangle 
               Add(s,(Triangle *)ta);
               return;
            }
         }
      }
      else {
         cerr << " bug " << nbd0 <<endl;
         cerr << " Bug double points in " << endl ;
         cerr << " s = " << Number(s) << " " << s << endl;
         cerr << " s0 = " << Number(*s0) << " " << *s0 << endl;
         cerr << " s1 = " << Number(*s1) << " " << *s1 << endl;
         cerr << " s2 = " << Number(*s2) << " " << *s2 << endl;
         MeshError(5);
      }
   }

// remove of MarkUnSwap edge
   t->SetUnMarkUnSwap(0);
   t->SetUnMarkUnSwap(1);
   t->SetUnMarkUnSwap(2);
   Triangle *tt[3]; // the 3 new triangles
   tt[0] = t;
   tt[1] = &triangles[nbt++];
   tt[2] = &triangles[nbt++];
   if (nbt>nbtx) {
      cerr << " Not enough triangles." << endl;
      MeshError(999);
   }

   *tt[1] = *tt[2] = *t;
   tt[0]->link = tt[1];
   tt[1]->link = tt[2];
   (*tt[0])(OppositeVertex[0]) = &s;
   (*tt[1])(OppositeVertex[1]) = &s;
   (*tt[2])(OppositeVertex[2]) = &s;
   tt[0]->det = det3[0];
   tt[1]->det = det3[1];
   tt[2]->det = det3[2];

// update adj of external triangles
   tt[0]->SetAdjAdj(0);
   tt[1]->SetAdjAdj(1);
   tt[2]->SetAdjAdj(2);
// update adj of the 3 internal triangles
   const int i0=0;
   const int i1=NextEdge[i0], i2=PreviousEdge[i0];

   tt[i0]->SetAdj2(i2,tt[i2],i0);
   tt[i1]->SetAdj2(i0,tt[i0],i1);
   tt[i2]->SetAdj2(i1,tt[i1],i2);

   tt[0]->SetTriangleContainingTheVertex();
   tt[1]->SetTriangleContainingTheVertex();
   tt[2]->SetTriangleContainingTheVertex();

// swap if the point s is on an edge
   if (izerodet>=0) {
      int rswap=tt[izerodet]->swap(iedge);
      if (!rswap)
         cout << " Pb swap the point s is on an edge =>swap " << iedge << " "  << *tt[izerodet] << endl;
      assert(rswap);
   }
}


long Triangles::SplitInternalEdgeWithBorderVertices()
{
   long NbSplitEdge=0;
   SetVertexFieldOn();
   long nbvold = nbv;
   for (long it=0; it<nbt; it++) {
      Triangle &t=triangles[it];
      if (t.link) {
         for (int j=0; j<3; j++) {
            if (!t.Locked(j) && !t.Hidden(j)) {
               Triangle *tt=t.TriangleAdj(j);
               if (tt && tt->link && it<Number(*tt)) { // an internal edge 
                  Vertex &v0 = t[VerticesOfTriangularEdge[j][0]];
                  Vertex &v1 = t[VerticesOfTriangularEdge[j][1]];
                  if (v0.on && v1.on) {
                     R2 P = (R2(v0) + R2(v1))*0.5;
                     if (nbv<nbvx) {
                        _vertices[nbv].r = P;
                        _vertices[nbv++].m = Metric(0.5,v0.m,0.5,v1.m);
                        _vertices[nbv].ReferenceNumber = 0;
                        _vertices[nbv].DirOfSearch = NoDirOfSearch;
                     }
                     NbSplitEdge++;
                     if (verbosity>7)
                        cout << " Internal edge with two vertices on boundary" 
                             << Number(v0) << " " << Number(v1) << " by " << endl;
                  }
               }
            }
         }
      }
   }
   ReMakeTriangleContainingTheVertex();    
   if (nbvold!=nbv) {
      long iv = nbvold, NbSwap = 0;
      Icoor2 dete[3];
      for (long i=nbvold; i<nbv; i++) {
         Vertex &vi=_vertices[i];
         vi.i = toI2(vi.r);
         vi.r = toR2(vi.i);
         vi.ReferenceNumber = 0; 
         vi.DirOfSearch = NoDirOfSearch;
         Triangle *tcvi = FindTriangleContaining(vi.i,dete);
         if (tcvi && !tcvi->link) {
            cout << i <<  " PB insert point " << Number(vi) << vi << Number(vi) 
                 << " tcvi = " << tcvi << " " << tcvi->link << endl;
            cout << (*tcvi)[1] << (*tcvi)[2] << endl;
            tcvi = FindTriangleContaining(vi.i,dete);
            cout << (*tcvi)[1] << (*tcvi)[2] << endl;
            MeshError(1001);
         }
         _quadtree->Add(vi);
         assert(tcvi && tcvi->det>=0);
         Add(vi,tcvi,dete);
         NbSwap += vi.Optim(1);
         iv++;
      }
      if (verbosity>3) {
         cout << "    Nb Of new points " << iv;
         cout << " Nb swap = " << NbSwap << " to  split internal edges with border vertices";
      }
      nbv = iv;
   }
   if (NbSplitEdge>nbv-nbvold)
      cout << " Warning not enough vertices to split all internal edges " << endl
           << "    we lost " << NbSplitEdge-(nbv-nbvold) << " edges !!" << endl;
   if (verbosity>2)
      cout << "SplitInternalEdgeWithBorderVertices: Number of splitted edges " << NbSplitEdge << endl;
   return NbSplitEdge;
}


long Triangles::InsertNewPoints(long  nbvold,
                                long& NbTSwap)
{
   double seuil = 0.707;
   long i;

// insertion part ---
   const long nbvnew=nbv-nbvold;
   if (verbosity>5) 
      cout << "    Try to insert the " << nbvnew << " new points " << endl;  
   long NbSwap=0;
   Icoor2 dete[3];

// construction of random order 
   if (!nbvnew)
      return 0;
   if (nbvnew) {
      const long PrimeNumber=AGoodNumberPrimeWith(nbv);
      long k3=rand()%nbvnew;
      for (long is3=0; is3<nbvnew; is3++) {
         long j = nbvold +(k3 = (k3+PrimeNumber)%nbvnew);
         long i = nbvold+is3; 
         ordre[i] = _vertices + j;
         ordre[i]->ReferenceNumber = i;
      }

//    Be careful
      long iv=nbvold;
      for (i=nbvold; i<nbv; i++) {// for all the new point
         Vertex &vi=*ordre[i];
         vi.i = toI2(vi.r);
         vi.r = toR2(vi.i);
         double hx, hy;
         vi.m.Box(hx,hy);
         Icoor1 hi=(Icoor1)(hx*coefIcoor), hj=(Icoor1)(hy*coefIcoor);
         if (!_quadtree->ToClose(vi,seuil,hi,hj)) {
            Vertex &vj=_vertices[iv];
            long j=vj.ReferenceNumber; 
            assert(&vj== ordre[j]);
            Exchange(vi,vj);
            Exchange(ordre[j],ordre[i]);
            vj.ReferenceNumber = 0;
            Triangle *tcvj=FindTriangleContaining(vj.i,dete);
            if (tcvj && !tcvj->link) {
               cerr << i <<  " PB insert point " << Number(vj) << vj << Number(vi) 
                    << " tcvj = " << tcvj << " " << tcvj->link << endl;
               cerr << (*tcvj)[1] << (*tcvj)[2] << endl;
               tcvj = FindTriangleContaining(vj.i,dete);
               MeshError(1001);
            }
            _quadtree->Add(vj);
            assert(tcvj && tcvj->det>=0);
            Add(vj,tcvj,dete);
            NbSwap += vj.Optim(1);          
            iv++;
         }
      }
      if (verbosity>3) {
         cout << "    Nb Of new points " << iv << ", Nb of too close points " << nbv-iv;
         cout << " Nb swap = " << NbSwap << " after ";
      }
      nbv = iv;
   }
   for (i=nbvold; i<nbv; i++) 
      NbSwap += _vertices[i].Optim(1);  
   if (verbosity>3) 
      cout << " NbSwap = " << NbSwap << endl;
   NbTSwap +=  NbSwap ;
   return nbv-nbvold;
}


void  Triangles::NewPoints(Triangles& Bh,
                           int        KeepBackVertex)
{
   long nbtold(nbt), nbvold(nbv);
   if (verbosity>2) 
      cout << " -- Triangles::NewPoints ";
   if (verbosity>3)
      cout <<  " nbv (in)  on Boundary  = " << nbv << endl;
   long i, k;
   int j;
   long *first_np_or_next_t = new long[nbtx];
   long NbTSwap=0;
   nbtold = nbt;
   if (KeepBackVertex && (&Bh != this) && (nbv+Bh.nbv<nbvx)) {
      for (i=0; i<Bh.nbv; i++) {
         Vertex &bv = Bh[i];
         if (!bv.on) {
            _vertices[nbv].r = bv.r;
            _vertices[nbv++].m = bv.m;
         }
      }
      int nbv1=nbv;
      Bh.ReMakeTriangleContainingTheVertex();
      InsertNewPoints(nbvold,NbTSwap);
      if (verbosity>2)
         cout <<  "      (Nb of Points from background mesh  = " 
              << nbv-nbvold  << " / " << nbv1-nbvold << ")" << endl;
   }
   else
      Bh.ReMakeTriangleContainingTheVertex();     

   Triangle *t;
// generation of the list of next Triangle 
// at 1 time we test all the triangles
   long Headt=0, next_t;
   for (i=0; i<nbt; i++)
      first_np_or_next_t[i] = -(i+1);
// end list i >= nbt 
// the list of test triangle is 
// the next traingle on i is  -first_np_or_next_t[i]
   int iter=0;
// Big loop
   do {
      iter++;
      nbtold = nbt;
      nbvold = nbv;
//    default size of  IntersectionTriangle

      i = Headt;
      next_t = -first_np_or_next_t[i];
      for (t=&triangles[i]; i<nbt; t=&triangles[i=next_t], next_t=-first_np_or_next_t[i]) {
         assert(i>=0 && i < nbt );
         first_np_or_next_t[i] = iter; 
         for (j=0; j<3; j++) { // for each edge
            TriangleAdjacent tj(t,j);
            Vertex& vA = * tj.EdgeVertex(0);
            Vertex& vB = * tj.EdgeVertex(1);
            if (!t->link)
               continue;// boundary
            if (t->det<0)
               continue;
            if (t->Locked(j))
               continue;

            TriangleAdjacent tadjj=t->Adj(j);	  
            Triangle *ta=tadjj;
            if (ta->det <0)
               continue;
            R2 A=vA, B=vB;
            k = Number(ta);
            if (first_np_or_next_t[k]==iter)  // this edge is done before 
               continue; // next edge of the triangle 
            lIntTria.SplitEdge(Bh,A,B);
            lIntTria.NewPoints(_vertices,nbv,nbvx);
         }
      }

      if (!InsertNewPoints(nbvold,NbTSwap)) 
         break;

      for (i=nbtold;i<nbt;i++)
         first_np_or_next_t[i] = iter;

      Headt = nbt; // empty list
      for (i=nbvold; i<nbv; i++) { // for all triangles containing the vertex i
         Vertex *s = _vertices + i;
         TriangleAdjacent ta(s->t, EdgesVertexTriangle[s->vint][1]);
         Triangle *tbegin = (Triangle*)ta;
         long kt;
         do { 
            kt = Number((Triangle*)ta);
            if (first_np_or_next_t[kt]>0)
               first_np_or_next_t[kt] = -Headt, Headt = kt;
            assert(ta.EdgeVertex(0)==s);
            ta = Next(Adj(ta));
         } while ((tbegin != (Triangle*) ta)); 
      }
   } while (nbv!=nbvold);
   delete [] first_np_or_next_t;

   long NbSwapf=0;
   for (i=0; i<nbv; i++)
      NbSwapf += _vertices[i].Optim(0);
   NbTSwap +=  NbSwapf;
   if (verbosity>3)
      cout << "   ";
   if (verbosity>2)
      cout << " Nb Of Vertices = " << nbv << " Nb of triangles = " << nbt-NbOutT 
           << " NbSwap final = " << NbSwapf << " Nb Total Of Swap = " << NbTSwap << endl;
}


void  Triangles::NewPointsOld(Triangles &Bh)
{
   double seuil= 0.7 ;// for two nearest points 
   if (verbosity>1) 
      cout << " begin: Triangles::NewPointsOld " << endl;
   long i,k;
   int j;
   long BeginNewPoint[3], EndNewPoint[3];
   int step[3];
   long *first_np_or_next_t = new long[nbtx];
   long color=-1;
   Triangle *t;
// generation of the list of next Triangle 
// at 1 time we test all the triangles
   long Headt =0,next_t;
   for (i=0; i<nbt; i++)
      first_np_or_next_t[i] = -(i+1);
// end list i >= nbt 
// the list of test triangle is 
// the next Triangle on i is  -first_np_or_next_t[i]
   long nbtold, nbvold;

// Big loop 
   do {
      nbtold = nbt;
      nbvold = nbv;
//    default size of  IntersectionTriangle

      i = Headt;
      next_t = -first_np_or_next_t[i];
      for (t=&triangles[i]; i<nbt; t=&triangles[i=next_t], next_t=-first_np_or_next_t[i]) {
#ifdef TRACETRIANGLE
         trace =  TRACETRIANGLE <0 ? 1 : i == TRACETRIANGLE;
#endif
         assert(i>=0 && i < nbt );
         first_np_or_next_t[i] = nbv; // to save the fist new point of triangle
         for (j=0; j<3; j++) { // for each edge
            TriangleAdjacent tj(t,j);
            color = 3*i + j;
            BeginNewPoint[j] = nbv;
            EndNewPoint[j] = nbv - 1;
            step[j] = 1; // right sense
            Vertex &vA = * tj.EdgeVertex(0);
            Vertex &vB = * tj.EdgeVertex(1);
#ifdef TRACETRIANGLE
            if (trace)
               cout << " i " << Number(vA) << " j " << Number(vB) << " " << t->Locked(j);
#endif
            if (!t->link)
               continue;// boundary
            if (t->det<0)
               continue;
            if (t->Locked(j))
               continue;
            TriangleAdjacent tadjj = t->Adj(j);
            Triangle *ta=tadjj;
            if (ta->det<0)
               continue;
            R2 A=vA, B=vB;
            k = Number(ta);
//          the 2 opposite vertices
            const Vertex &vC1=*tj.OppositeVertex();
            const Vertex &vC2=*tadjj.OppositeVertex();
	  
#ifdef TRACETRIANGLE
            trace = trace || k == TRACETRIANGLE;
            if (trace) {
               cout << "Test Arete " << i << " AB = " << A << B 
                    << "i "  << Number(vA) << "j" << Number(vB); 
               cout << " link" <<(int)t->link << " ta=" << Number(ta) 
                    << " det " << ta->det;
               cout << " hA = " << vA.m.h << " hB = " << vB.m.h;
               cout << " loc " << ta->Locked(j) << endl;
            }
#endif
            if (first_np_or_next_t[k]>0) { // this edge is done before 
//             find the color of the edge and begin, end of newpoint
               int kk = t->NuEdgeTriangleAdj(j);
               assert((*t)(VerticesOfTriangularEdge[j][0]) == (*ta)(VerticesOfTriangularEdge[kk][1]));
               assert((*t)(VerticesOfTriangularEdge[j][1]) == (*ta)(VerticesOfTriangularEdge[kk][0]));
               long kolor=3*k+kk, kkk=1;
               step[j] = -1;// other sens 
               BeginNewPoint[j] = 0;
               EndNewPoint[j] = -1; // empty list          
               for (long iv=first_np_or_next_t[k]; iv<nbv; iv++) 
               if (_vertices[iv].color > kolor)
                  break; // the color is passed            
               else if (_vertices[iv].color == kolor) {
                  EndNewPoint[j] = iv; 
                  if (kkk) // one time test
                     kkk=0, BeginNewPoint[j] = iv;
               }
               continue; // next edge of the triangle 
            } // end  if( k < i) 
            const long NbvOld=nbv;
            lIntTria.SplitEdge(Bh,A,B);
            long nbvNew=nbv;
            nbv = NbvOld;
            for (long iv=NbvOld; iv<nbvNew; iv++) {
               _vertices[nbv].color = color;
               _vertices[nbv].ReferenceNumber = nbv;  // circular link
               R2 C = _vertices[iv].r;
               _vertices[nbv].r = C;
               _vertices[nbv].m = _vertices[iv].m ;
//             test if the new point is not to close to the 2 opposite vertex
               R2 CC1=C-vC1, CC2=C-vC2;
               if (((vC1.m(CC1)+_vertices[nbv].m(CC1))>seuil) && ((vC2.m(CC2)+_vertices[nbv].m(CC2))>seuil))
                  nbv++;
            }
            EndNewPoint[j] = nbv-1;
         } // end loop for each edge 

#ifdef TRACETRIANGLE
         if (trace) {
            cout << "\n ------------ " << t->link << " " << t->det 
                 << " b " << BeginNewPoint[0] << " " << BeginNewPoint[1]
                 << " " << BeginNewPoint[2] << " " 
                 << " e " << EndNewPoint[0] << " " << EndNewPoint[1] 
                 << " " << EndNewPoint[2] << " " 
                 << " s " << step[0] << " " << step[1] << " " << step[2] << " " << endl;
         }
#endif
         if (!t->link)
            continue;  // boundary
         if (t->det<=0)
            continue;// outside 

         for (int j0=0; j0<3; j0++) {
            for (long i0=BeginNewPoint[j0]; i0<=EndNewPoint[j0]; i0++) {
//             find the nearest point one the opposite edge to compute i1
               Vertex &vi0=_vertices[i0];
               int kstack=0, j1=j0;
               long stack[10];
               while (j0 != (j1=NextEdge[j1])) { // loop on the 2 other edge
//                computation of the intersection of edge j1 and DOrto
//                take the good sens
                  if (BeginNewPoint[j1] > EndNewPoint[j1])
                     continue;
                  else if (EndNewPoint[j1] - BeginNewPoint[j1] < 1) {
                     for (long ii1=BeginNewPoint[j1]; ii1<=EndNewPoint[j1]; ii1++)
                        stack[kstack++] = ii1;
                     continue;
                  }
                  int k0=0, k1=0;
                  if (step[j1]<0)
                     k0 = 1; // reverse
                  else
                     k1 = 1; 
                  R2 V10 = R2((*t)[VerticesOfTriangularEdge[j1][k0]]);
                  R2 V11 = R2((*t)[VerticesOfTriangularEdge[j1][k1]]);
                  R2 D=V11-V10;
                  double c0=vi0.m(D,R2(vi0)), c10=vi0.m(D,V10), c11=vi0.m(D,V11);
                  int sss = (c11-c10)>0 ? 1 : -1;
                  long ii0=BeginNewPoint[j1], ii1=EndNewPoint[j1];
                  double ciii=-1, cii0=-1, cii1=-1;
                  if (sss * ((cii0=vi0.m(D,R2(_vertices[ii0])))-c0)>0)
                     stack[kstack++] = ii0;
                  else if (sss * ((cii1=vi0.m(D,R2(_vertices[ii1])))-c0)<0)
                     stack[kstack++] = ii1;
                  else {
                     while ((ii1-ii0)> 1) {
                        long iii = (ii0+ii1)/2;
                        ciii = vi0.m(D,R2(_vertices[iii]));
                        if (sss*(ciii-c0)<0)
                           ii0 = iii;
                        else
                           ii1 = iii;
                     }
                     stack[kstack++] = ii0;
                     if (ii1 != ii0)
                        stack[kstack++] = ii1;
                  }
                  if (kstack >5) // bug ?
                     cout << "NewPoints: bug????? " << kstack << " stack  " << stack[kstack] << endl;
               }
               stack[kstack++] = -1; // to stop
               long i1;
               kstack = 0; 
               while ((i1=stack[kstack++]) >= 0) { // the two parameters are i0 and i1 
                  assert(i1<nbv && i1>=0);
                  assert(i0<nbv && i0>=0);
                  assert(i1!=i0);
                  R2 v01 = R2(_vertices[i1])-R2(_vertices[i0]);
                  double d01 = (_vertices[i0].m(v01) + _vertices[i1].m(v01));
#ifdef TRACETRIANGLE
                  if (trace) {
                     cout << "\n test j" << j <<" "  << i0 
                          << " " << i1 << " d01=" << d01 <<endl;}
#endif
                  assert(i0>=nbvold);
                  assert(i1>=nbvold);
                  assert(i0!=i1);
                  if (d01==0) 
                     break;
                  if (d01<seuil) {
                     if (i1<nbvold) {
//                      remove all the points i0;
                        long ip, ipp;
                        for (ip=i0; i0!=(ipp=_vertices[ip].ReferenceNumber); ip=ipp)
                           _vertices[ip].ReferenceNumber = -1;// mark remove
                        _vertices[ip].ReferenceNumber = -1;// mark remove
                     }
                     else {
//                      remove on of two points
                        long ip0, ip1, ipp0, ipp1;
                        int kk0=1, kk1=1;
//                      count the number of common points to compute weight w0,w1
                        for (ip0=i0; i0!=(ipp0=_vertices[ip0].ReferenceNumber); ip0=ipp0)
                           kk0++;
                        for (ip1=i1; i1!=(ipp1=_vertices[ip1].ReferenceNumber); ip1=ipp1)
                           kk1++;
                        double w0=double(kk0)/(kk0+kk1), w1=double(kk1)/(kk0+kk1);
//                      make a circular link
                        Exchange(_vertices[i0].ReferenceNumber,_vertices[i1].ReferenceNumber);
//                      the new coordinate 
                        R2 C = _vertices[i0].r*w0 + _vertices[i1].r*w1;
#ifdef TRACETRIANGLE
                        if (trace)
                           cout << "\n ref = " << _vertices[i0].ref << " " << _vertices[i1].ref << endl;
#endif
//                      update the new point points of the list 
                        for (ip0=i0; i0!=(ipp0=_vertices[ip0].ReferenceNumber); ip0=ipp0)
                           _vertices[ip0].r = C;
                        _vertices[ip0].r = C;
                     }
                  }
               }
            }
         }
      }

//    remove all double points

      long ip, ipp, kkk=nbvold;
      for (i=nbvold; i<nbv; i++) {
         if (_vertices[i].ReferenceNumber>=0) {
            for (ip=i; i!=(ipp=_vertices[ip].ReferenceNumber); ip=ipp)
               _vertices[ip].ReferenceNumber = -1;
            _vertices[ip].ReferenceNumber = -1;
            _vertices[kkk] = _vertices[i];
            _vertices[kkk].i = toI2(_vertices[kkk].r);
            _vertices[kkk++].ReferenceNumber = 0; 
         }
      }

//    insertion part --- 
      const long nbvnew=kkk-nbvold;
      cout << "    Remove " << nbv - kkk  << " to close  vertex " ;
      cout << " and insert the " << nbvnew << " new points " << endl;  
      nbv = kkk;
      long NbSwap=0;
      Icoor2 dete[3];

//    Construction of a random order 
      if (!nbvnew)
         break;
      if (nbvnew) {
         const long PrimeNumber=AGoodNumberPrimeWith(nbv);
         long k3=rand()%nbvnew;
         for (long is3=0; is3<nbvnew; is3++)
            ordre[nbvold+is3]= &_vertices[nbvold+(k3=(k3+PrimeNumber)%nbvnew)];

         for (i=nbvold; i<nbv; i++) {
            Vertex *vi=ordre[i];
            Triangle *tcvi=FindTriangleContaining(vi->i,dete);
            _quadtree->Add(*vi);
            assert (tcvi->det>=0);
            Add(*vi,tcvi,dete);
            NbSwap += vi->Optim(1);          
         }
      }
      cout << " Nb swap = " << NbSwap << " after ";

      for (i=nbvold;i<nbv;i++) 
         NbSwap += _vertices[i].Optim(1);  
 
      for (i=nbtold; i<nbt; i++)
         first_np_or_next_t[i] = 1;
      Headt = nbt; // empty list 
      for (i=nbvold; i<nbv; i++) { // for all triangles containing the vertex i
         Vertex *s = _vertices + i;
         TriangleAdjacent ta(s->t, EdgesVertexTriangle[s->vint][1]);
         Triangle *tbegin=(Triangle*)ta;
         long kt;
         do {
            kt = Number((Triangle*)ta);
            if (first_np_or_next_t[kt]>0)
               first_np_or_next_t[kt] = -Headt, Headt = kt;
            assert(ta.EdgeVertex(0) == s);
            ta = Next(Adj(ta));
         } while ((tbegin != (Triangle*) ta)); 
      }
   } while (nbv!=nbvold);
   delete [] first_np_or_next_t;
   cout << " end: Triangles::NewPoints old  nbv=" << nbv << endl;
}


void Triangles::Insert() 
{
   if (verbosity>2)
      cout << " -- Initial insert " << nbv << " vertices " << endl ;
   Triangles *OldCurrentTh=CurrentTh;
   CurrentTh = this;
   double time0=CPUtime(), time1, time2, time3;
   SetIntCoor();
   for (long i=0; i<nbv; i++) 
      ordre[i] = &_vertices[i];

// construction d'un ordre aleatoire 
   const long PrimeNumber = AGoodNumberPrimeWith(nbv);
   long k3 = rand()%nbv;
   for (int is3=0; is3<nbv; is3++) 
      ordre[is3] = & _vertices[k3 = (k3 + PrimeNumber)% nbv];
   long i;
   for (i=2; det(ordre[0]->i,ordre[1]->i,ordre[i]->i)==0;) {
      if (++i>=nbv) {
         cerr << " All vertices are aligned." << endl;
         MeshError(998);
      }
   }

// echange i et 2 dans ordre afin 
// que les 3 premiers ne soit pas aligne
   Exchange(ordre[2],ordre[i]);

// on ajoute un point a l'infini pour construire le maillage
// afin d'avoir une definition simple des aretes frontieres
   nbt = 2;

// on construit un maillage trivale forme
// d'une arete et de 2 triangles
// construit avec le 2 aretes orientes et 
   Vertex *v0=ordre[0], *v1=ordre[1];

   triangles[0](0) = 0; // sommet pour infini 
   triangles[0](1) = v0;
   triangles[0](2) = v1;

   triangles[1](0) = 0; // sommet pour infini 
   triangles[1](2) = v0;
   triangles[1](1) = v1;
   const int e0=OppositeEdge[0], e1=NextEdge[e0], e2=PreviousEdge[e0];
   triangles[0].SetAdj2(e0,&triangles[1],e0);
   triangles[0].SetAdj2(e1,&triangles[1],e2);
   triangles[0].SetAdj2(e2,&triangles[1],e1);
   triangles[0].det = -1;  // fake triangles
   triangles[1].det = -1;  // fake triangles
   triangles[0].SetTriangleContainingTheVertex();
   triangles[1].SetTriangleContainingTheVertex();
   triangles[0].link = &triangles[1];
   triangles[1].link = &triangles[0];
   if (!_quadtree)
      _quadtree = new QuadTree(this,0);
   _quadtree->Add(*v0);
   _quadtree->Add(*v1);
   long NbSwap=0;
   time1 = CPUtime();
   if (verbosity>3)
      cout << " -- Begin of insertion process " << endl;

   for (long icount=2; icount<nbv; icount++) {
      Vertex *vi=ordre[icount];
      Icoor2 dete[3];
      Triangle *tcvi=FindTriangleContaining(vi->i,dete);
      _quadtree->Add(*vi);
      Add(*vi,tcvi,dete);
      NbSwap += vi->Optim(1,0);
   }
   time2 = CPUtime();
   if (verbosity>3) 
      cout << " NbSwap of insertion " << NbSwap 
           << " NbSwap/Nbv " << float(NbSwap)/float(nbv) 
           << " NbUnSwap " << NbUnSwap << ", Nb UnSwap/Nbv " 
           << float(NbUnSwap)/float(nbv) << endl;
   NbUnSwap = 0;
#ifdef NBLOOPOPTIM
   k3 = rand()%nbv;
   for (int is4=0; is4<nbv; is4++) 
      ordre[is4] = &_vertices[k3 = (k3 + PrimeNumber)% nbv];
   double timeloop=time2;
   for (int Nbloop=0; Nbloop<NBLOOPOPTIM; Nbloop++) {
      double time000 = timeloop;
      long NbSwap=0;
      for (int is1=0; is1<nbv; is1++) 
         NbSwap += ordre[is1]->Optim(0,0);
      timeloop = CPUtime();
      if (verbosity>3)
         cout << "    Optim Loop " << Nbloop << " NbSwap: " << NbSwap 
              << " NbSwap/Nbv " << float(NbSwap)/float(nbv)
              << " CPU = " << timeloop - time000 << "  s, "
              << " NbUnSwap/Nbv " << float(NbUnSwap)/float(nbv) << endl;
      NbUnSwap = 0;
      if (!NbSwap)
         break;
   }
   ReMakeTriangleContainingTheVertex(); 
#endif
   time3 = CPUtime();
   if (verbosity>4) 
      cout << "    init " << time1 - time0 << " initialisation,  " 
           << time2 - time1 << "s, insert point  "
           << time3 -time2 << "s, optim " << endl
           << "     Init Total Cpu Time = " << time3-time0 << "s " << endl;
   CurrentTh = OldCurrentTh;
}


void Triangles::ForceBoundary()
{
   if (verbosity>2)
      cout << " -- ForceBoundary nb of edge " << nbe << endl;
   int k=0;
   long nbfe=0, nbswp=0, Nbswap=0;
   for (long t=0; t<nbt; t++)  
      if (!triangles[t].det)
         k++, cerr << " det T" << t << " = 0" << endl;
   if (k!=0) {
      cerr << " the is  " << k << " 0 triangles " << endl;
      MeshError(11);
   }
   TriangleAdjacent ta(0,0);
   for (long i=0; i<nbe; i++) {
      nbswp = ForceEdge(edges[i][0],edges[i][1],ta);
      if (nbswp<0)
         k++;
      else
         Nbswap += nbswp;
      if (nbswp)
         nbfe++;
      if (nbswp>=0 && edges[i].on->Cracked())
         ta.SetCracked();
   }
   if (k!=0) {
      cerr << " there are " << k << " lost edges " << endl;
      cerr << " The boundary is maybe crossing!" << endl;
      MeshError(10);
   }
   for (long j=0; j<nbv; j++)
      Nbswap += _vertices[j].Optim(1,0);
   if (verbosity>3)
      cout << "     Nb of inforced edge = " << nbfe << " Nb of Swap "
           << Nbswap << endl;
}


void Triangles::FindSubDomain(int OutSide=0)
{
   if (verbosity >2) {
      if (OutSide)
         cout << " -- Find all external sub-domains ";	
      else
         cout << " -- Find all internal sub-domains ";
   }
   short *HeapArete = new short[nbt];
   Triangle **HeapTriangle = new Triangle* [nbt];
   Triangle *t, *t1;
   long k, it;
   for (long itt=0; itt<nbt; itt++) 
      triangles[itt].link = 0;
   long NbSubDomTot=0;
   for (it=0; it<nbt; it++)  { 
      if (!triangles[it].link) {
         t = triangles + it;
         NbSubDomTot++; // new composante connexe
         long i=0; // niveau de la pile 
         t->link = t; // sd forme d'un triangle cicular link
         HeapTriangle[i] = t; 
         HeapArete[i] = 3;
         while (i>=0) { // boucle sur la pile
            while (HeapArete[i]--) { // boucle sur les 3 aretes 
	       int na=HeapArete[i];
	       Triangle *tc = HeapTriangle[i]; // triangle courant
	       if (!tc->Locked(na)) {  // arete non frontiere
                  Triangle *ta = tc->TriangleAdj(na);
                  if (ta->link == 0 ) {
                     i++;
                     ta->link = t->link;  // on chaine les triangles
                     t->link = ta;        // d'un meme sous domaine
                     HeapArete[i] = 3;     // pour les 3 triangles adjacents
                     HeapTriangle[i] = ta;
                  }
               }
            }
            i--;
         }
      }
   }

// Supression de tous les sous domaine infini <=> contient le sommet NULL
   it = 0;
   NbOutT = 0;
   while (it<nbt) {
      if (triangles[it].link) {
         if (!(triangles[it](0) && triangles[it](1) && triangles[it](2))) {
            NbSubDomTot --;
            t = &triangles[it];
            while (t) {
               NbOutT++;
               t1 = t;
               t = t->link;
               t1->link = 0;
            }
         }
      }
      it++;
   }
   if (nbt==NbOutT || !NbSubDomTot) {
      cerr << "Error: The boundary is not closed => All triangles are outside " << endl;
      MeshError(888);
   }
   delete [] HeapArete;
   delete [] HeapTriangle;
   if (OutSide|| !Gh.subdomains || !Gh.NbSubDomains) {
      long i;
      if (subdomains)
         delete [] subdomains;
      subdomains = new SubDomain [NbSubDomTot];
      NbSubDomains = NbSubDomTot;
      for (i=0; i<NbSubDomains; i++) {
         subdomains[i].head = NULL;
         subdomains[i].ref = i + 1;
      }
      long *mark = new long [nbt];
      for (it=0; it<nbt; it++)
         mark[it] = triangles[it].link ? -1 : -2;
      it = 0;
      k = 0;
      while (it<nbt) {
         if (mark[it] == -1) {
            t1 = &triangles[it];
            t = t1->link;
            mark[it] = k;
            subdomains[k].head = t1;
            do {
               mark[Number(t)] = k;
               t = t->link;
            } while (t!=t1);
            mark[it] = k++;
         }
         it++;
      }
      assert(k==NbSubDomains);
      if (OutSide) {
         long nbk = NbSubDomains;
         while (nbk) {
            for (it=0; it<nbt && nbk; it++) {
               for (int na=0; na<3 && nbk; na++) {
                  Triangle *ta = triangles[it].TriangleAdj(na);
                  long kl = ta ? mark[Number(ta)] : -2;
                  long kr = mark[it];
                  if (kr != kl) {
                     if (kl>=0 && subdomains[kl].ref<0 && kr>=0 && subdomains[kr].ref>=0)
                        nbk--, subdomains[kr].ref = subdomains[kl].ref - 1;
                     if (kr>=0 && subdomains[kr].ref<0 && kl>=0 && subdomains[kl].ref>=0)
                        nbk--, subdomains[kl].ref = subdomains[kr].ref - 1;
                     if (kr<0 && kl>=0 && subdomains[kl].ref>=0)
                        nbk--, subdomains[kl].ref=-1;
                     if (kl<0 && kr>=0 && subdomains[kr].ref>=0)
                        nbk--, subdomains[kr].ref = -1;
                  }
               }
            }
         }
         long j=0;
         for (i=0; i<NbSubDomains; i++) {
            if ((-subdomains[i].ref)%2) {
               if (i != j) 
                  Exchange(subdomains[i],subdomains[j]);
               j++;
            }
            else {
               t = subdomains[i].head;
               while (t) {
                  NbOutT++;
                  t1 = t;
                  t = t->link;
                  t1->link = 0;
               }
            }
         }
         if (verbosity>4)
            cout << " Number of remove sub domain (OutSideMesh) =" << NbSubDomains-j << endl;
         NbSubDomains=j;
      }
      delete [] mark;
   }
   else { // find the head for all sub domaine
      if (Gh.NbSubDomains != NbSubDomains && subdomains)
         delete [] subdomains, subdomains=NULL;
      if (!subdomains)
         subdomains = new SubDomain [Gh.NbSubDomains];
      NbSubDomains = Gh.NbSubDomains;
      if (verbosity>4)
         cout << "     find the " << NbSubDomains << " sub domain " << endl;
      long err=0;
      ReMakeTriangleContainingTheVertex();
      long *mark = new long[nbt];
      Edge **GeometricalEdgetoEdge = MakeGeometricalEdgeToEdge();
      for (it=0; it<nbt; it++)
         mark[it] = triangles[it].link ? -1 : -2;
      long inew =0;
      for (long i=0; i<NbSubDomains; i++) {
         GeometricalEdge &eg = *Gh.subdomains[i].edge;
         subdomains[i].ref = Gh.subdomains[i].ref;
         Edge &e = *GeometricalEdgetoEdge[Gh.Number(eg)];
//         assert(&e);
         Vertex *v0=e(0), *v1 = e(1);
         Triangle *t = v0->t;
         int sens = Gh.subdomains[i].sens;
         if (((eg[0].r-eg[1].r),(e[0].r-e[1].r))<0)
            sens = -sens;
         subdomains[i].sens = sens;
         subdomains[i].edge = &e;
         assert(t && sens);
         TriangleAdjacent ta(t,EdgesVertexTriangle[v0->vint][0]);// previous edges
         while (1) {
            assert(v0==ta.EdgeVertex(1));
            if (ta.EdgeVertex(0)==v1) { // ok we find the edge
               if (sens>0)  
                  subdomains[i].head = t = Adj(ta);
               else
                  subdomains[i].head = t = ta;
               if (t<triangles || t >= triangles+nbt || t->det < 0) {
                  cerr << " Error in the def of subdomain: Wrong orientation " << i << " " << "Edge "
                       << Gh.Number(eg) << " " << sens << endl;
                  err++;
                  break;
               }
               long it=Number(t);
               if (mark[it]>=0) {
                  if (verbosity>10)
                     cerr << "     Warning: the sub domain " << i << " ref = " << subdomains[i].ref 
                          << " is previouly defined with "  << mark[it] << " ref = " << subdomains[mark[it]].ref
                          << " skip this def " << endl;
                  break;
               }
               if (i!=inew) 
                  Exchange(subdomains[i],subdomains[inew]);
               inew++;
               Triangle *tt=t;
	       long kkk=0;
               do {
                  kkk++;
                  assert(mark[Number(tt)]<0);
                  mark[Number(tt)] = i;
                  tt = tt->link;
               } while (tt!=t);
               if (verbosity>7)
                  cout << "     Nb of triangles in the subdomain " << i << " with ref " << subdomains[i].ref
                       << " = " << kkk << endl;
               break;
            }
            ta = Previous(Adj(ta));         
            if (t == (Triangle *) ta) {
               err++;
      	       cerr << " Error in the def of subdomain " << i
                    << " edge=" << Gh.Number(eg) << " " << sens << endl;
               break;
            }
         }
      }
      if (err)
         MeshError(777);
      if (inew<NbSubDomains) {
         if (verbosity>5) 
            cout << "     Warning: We remove " << NbSubDomains-inew << " subdomains " << endl;
         NbSubDomains = inew;
      }
      for (it=0; it<nbt; it++)
         if (mark[it]==-1) 
            NbOutT++, triangles[it].link = 0;
      delete [] GeometricalEdgetoEdge;
      delete [] mark;
   }
   NbOutT = 0;
   for (it=0; it<nbt; it++)
      if (!triangles[it].link)
         NbOutT++;
   if (verbosity> 4)
      cout << "    ";
   if (verbosity> 2)
      cout << " Nb of Sub bounded Domain  = " << NbSubDomTot << " NbOutTriangles = " << NbOutT << endl;
}


void Triangles::ReNumberingVertex(long *renu)
{
// warning be careful because pointer
// from on mesh to over mesh 
//  --  so do ReNumbering a the beginning
   Vertex *ve = _vertices + nbv;
   long it, ie, i;
   for (it=0; it<nbt; it++)
      triangles[it].ReNumbering(_vertices,ve,renu);
   for (ie=0; ie<nbe; ie++)
      edges[ie].ReNumbering(_vertices,ve,renu);
   for (i=0; i<NbVerticesOnGeomVertex; i++) {
      Vertex *v=VerticesOnGeomVertex[i].mv;
      if (v>=_vertices && v<ve)
         VerticesOnGeomVertex[i].mv = _vertices + renu[Number(v)];
   }
   for (i=0; i<NbVerticesOnGeomEdge; i++) {
      Vertex *v=VerticesOnGeomEdge[i].mv;
      if (v>=_vertices && v<ve)
	 VerticesOnGeomEdge[i].mv = _vertices + renu[Number(v)];
   }
   for (i=0; i<NbVertexOnBThVertex; i++) {
      Vertex *v=VertexOnBThVertex[i].v;
      if (v>=_vertices && v<ve)
         VertexOnBThVertex[i].v = _vertices + renu[Number(v)];
   }
   for (i=0; i<NbVertexOnBThEdge; i++) {
      Vertex *v=VertexOnBThEdge[i].v;
      if (v>=_vertices && v<ve)
         VertexOnBThEdge[i].v = _vertices + renu[Number(v)];
   }

// move the vertices without a copy of the array 
// be careful: non trivial code 
   long j;
   for (it=0; it<nbv; it++) { // for all sub cycles of the permutation renu
      if (renu[it]>=0) { // a new sub cycle
         i = it;
         Vertex ti=_vertices[i],tj;
         while ( (j=renu[i]) >= 0) { // i is old, and j is new 
            renu[i] = -1 - renu[i]; // mark 
            tj = _vertices[j]; // save new
            _vertices[j] = ti; // new <- old
            i = j;     // next 
            ti = tj;
         }
      }
   }
   if (_quadtree) {
      delete _quadtree;
      _quadtree = new QuadTree(this);
   }
   for (it=0; it<nbv; it++)
      renu[i] = -renu[i] - 1;  
}


void Triangles::ReNumberingTheTriangleBySubDomain(bool justcompress)
{
   long *renu = new long[nbt];
   Triangle *t0, *t, *te=triangles+nbt;
   long k=0, it, i, j;
   for (it=0; it<nbt; it++) 
      renu[it] = -1;
   for (i=0; i<NbSubDomains; i++) {
      t = t0 = subdomains[i].head;
      assert(t0);
      do { 
         long kt = Number(t);
         assert(kt>=0 && kt<nbt);
         assert(renu[kt]==-1);
         renu[kt] = k++;
      }
      while (t0 != (t=t->link));
   }
   if (verbosity>9)
      cout << " number of inside triangles " << k << " nbt = " << nbt << endl;
// take is same numbering if possible    
   if (justcompress)
      for (k=0, it=0; it<nbt; it++)
         if (renu[it]>=0) 
            renu[it] = k++;

// put the outside triangles at the end
   for (it=0; it<nbt; it++) 
      if (renu[it]==-1) 
         renu[it] = k++; 
    assert(k == nbt);

// do the change on all the pointer 
   for (it=0; it<nbt; it++)
      triangles[it].ReNumbering(triangles,te,renu);     
   for (i=0; i<NbSubDomains; i++)
      subdomains[i].head = triangles + renu[Number(subdomains[i].head)];

// move the Triangles  without a copy of the array 
// be careful not trivial code 
   for (it=0; it<nbt; it++) { // for all sub cycles of the permutation renu
      if (renu[it] >= 0) {  // a new sub cycle
         i = it;
         Triangle ti=triangles[i], tj;
         while ((j=renu[i]) >= 0) { // i is old, and j is new 
            renu[i] = -1; // mark 
            tj = triangles[j]; // save new
            triangles[j] = ti; // new <- old
            i = j;
            ti = tj;
         }
      }
   }
   delete [] renu;
   nt = nbt - NbOutT;
}


long Triangles::ConsRefTriangle(long *reft) const
{
   assert(reft);
   Triangle *t0, *t;
   long k=0, num;
   for (long it=0; it<nbt; it++) 
      reft[it] = -1;
   for (long i=0; i<NbSubDomains; i++) {
      t = t0 = subdomains[i].head;
      assert(t0);
      do {
         k++;
         num = Number(t);
         assert(num>=0 &&num < nbt);
         reft[num] = i;
      }
      while (t0 != (t=t->link));
   }
   if (verbosity>5) 
      cout << " Nb of Sub Domain =" << NbSubDomains  << " Nb of In Triangles " << k 
           << " Nbt = " << nbt << " Out Triangles = " << nbt-k << endl;
   return k;
}


Vertex *Triangles::NearestVertex(Icoor1 i,
                                 Icoor1 j) 
{
   return _quadtree->NearestVertex(i,j);
}


void Triangles::PreInit(long  inbvx,
                        char* fname)
{
   srand(19999999);
   OnDisk = 0;
   NbRef = 0;
   identity = 0;
   NbOfTriangleSearchFind = 0;
   NbOfSwapTriangle = 0;
   nbiv = 0;
   nbv = 0;
   nbvx = inbvx;
   nbt = 0;
   NbOfQuad = 0;
   nbtx = 2*inbvx - 2;
   NbSubDomains = 0;
   NbVertexOnBThVertex = 0;
   NbVertexOnBThEdge = 0;
   VertexOnBThVertex = 0;
   VertexOnBThEdge = 0;
   NbCrackedVertices = NbCrackedEdges = 0;
   CrackedEdges = 0;  
   nbe = 0; 
   name = fname;
   if (inbvx) {
      _vertices = new Vertex[nbvx];
      assert(_vertices);
      ordre = new Vertex *[nbvx];
      assert(ordre);
      triangles = new Triangle[nbtx];
      assert(triangles);
   }
   else {
      _vertices = 0;
      ordre = 0;
      triangles = 0;
      nbtx = 0;
   }
   if (name || inbvx) {
      time_t timer = time(0);
      char buf[70];     
      strftime(buf ,70,", Date: %y/%m/%d %H:%M %Ss",localtime(&timer));
      counter++; 
      char countbuf[30];   
      snprintf(countbuf,30,"%d",counter);
      int lg = 0;
      if (&BTh != this && BTh.name)
         lg = strlen(BTh.name) + 4;
      identity = new char[lg+strlen(buf)+strlen(countbuf)+2+10+(Gh.name?strlen(Gh.name)+4:0)];
      identity[0] = 0;
      if (lg)
         strcat(strcat(strcat(identity,"B="),BTh.name),", ");
      if (Gh.name)
         strcat(strcat(identity,"G="),Gh.name);
      strcat(strcat(identity,";"),countbuf);
      strcat(identity,buf);
   }
   _quadtree = 0;
   edges = 0;
   VerticesOnGeomVertex = VerticesOnGeomEdge = 0;
   NbVerticesOnGeomVertex = NbVerticesOnGeomEdge = NbSubDomains = 0;
   subdomains = 0;
   if (verbosity>98) 
      cout << "Triangles::PreInit() " << nbvx << " " << nbtx 
           << " " << _vertices << " " << ordre << " " <<  triangles << endl;
}


void Triangles::GeomToTriangles1(long inbvx,
                                 int  KeepBackVertices) 
{
   Gh.NbRef++;

/************************************************************************* 
  Method in 2 step:
   1 - compute the number of new edge to allocate
   2 - construct the edge

   Remark: 
   In this part we suppose to have a background mesh with the same geometry 

   To construct the discretization of the new mesh we have to 
   rediscretize the boundary of background Mesh 
   because we have only the pointer from the background mesh to the geometry.
   We need the abcisse of the background mesh vertices on geometry
   so a vertex is 
      0 on GeometricalVertex
      1 on GeometricalEdge + abcisse
      2 internal
  *************************************************************************/

   assert(&BTh.Gh == &Gh);
   BTh.NbRef++;
   PreInit(inbvx);
   BTh.SetVertexFieldOn();
   int *bcurve = new int[Gh.NbOfCurves];
  
// We have 2 ways to make the loop 
// 1) on the geometry 
// 2) on the background mesh
//    If you do the loop on geometry, we don't have the pointeur on background,
//    and if you do the loop in background we have the pointeur on geometry
//    so do the walk on  background

   NbVerticesOnGeomVertex = NbVerticesOnGeomEdge = 0;
   int i; 
   for (i=0; i<Gh.nbv; i++)
      if (Gh[i].Required())
         NbVerticesOnGeomVertex++;

   VerticesOnGeomVertex = new VertexOnGeom [NbVerticesOnGeomVertex];
   VertexOnBThVertex = new VertexOnVertex [NbVerticesOnGeomVertex];
   if (NbVerticesOnGeomVertex>=nbvx) {
      cerr << " Too much vertices on geometry " << NbVerticesOnGeomVertex << " >= " << nbvx << endl; 
      MeshError(1);
   }
   assert(_vertices);
   for (i=0; i<Gh.nbv; i++) {
      if (Gh[i].Required()) {
         _vertices[nbv] = Gh[i];
         _vertices[nbv].i = I2(0,0);
         Gh[i].to = _vertices + nbv;
         VerticesOnGeomVertex[nbv] = VertexOnGeom(_vertices[nbv],Gh[i]);
         nbv++;
      }
      else
         Gh[i].to = 0;
   }

   for (i=0; i<BTh.NbVerticesOnGeomVertex; i++) { 
      VertexOnGeom &vog = BTh.VerticesOnGeomVertex[i];
      if (vog.IsRequiredVertex()) {
         GeometricalVertex *gv=vog;
         Vertex *bv=vog;
         assert(gv->to);
         VertexOnBThVertex[NbVertexOnBThVertex++] = VertexOnVertex(gv->to,bv);
         gv->to->m = bv->m;
      }
   }
   assert(NbVertexOnBThVertex == NbVerticesOnGeomVertex);

   {
      Gh.UnMarkEdges();	
      int bfind=0;
      for (int i=0; i<Gh.NbOfCurves; i++)
         bcurve[i] = -1;
      for (int iedge=0; iedge<BTh.nbe; iedge++) {
         Edge &ei = BTh.edges[iedge];
         for (int je=0; je<2; je++) {
            if (!ei.on->Mark() && ei[je].on->IsRequiredVertex()) {
               int nc = ei.on->CurveNumber;
               if (ei.on==Gh.curves[nc].be && 
                   (GeometricalVertex *) *ei[je].on == &(*Gh.curves[nc].be)[Gh.curves[nc].kb]) { 
                  bcurve[nc] = iedge*2 + je;
                  bfind++;
               }
            }
         }
      }
      assert(bfind==Gh.NbOfCurves);
   }

// Method in 2 + 1 step 
//  0.0) compute the length and the number of vertex to do allocation
//  1.0) recompute the length
//  1.1) compute the vertex
   for (int step=0; step<2; step++) {
      long NbOfNewPoints=0, NbOfNewEdge=0, iedge;
      Gh.UnMarkEdges();	
      double L=0;
      for (int icurve=0; icurve<Gh.NbOfCurves; icurve++) {
         iedge = bcurve[icurve]/2;
         int jedge = bcurve[icurve]%2;
         if (!Gh.curves[icurve].master)
            continue;
         Edge &ei = BTh.edges[iedge];
         double Lstep=0;// step between two points   (phase==1)
         long NbCreatePointOnCurve=0;// Nb of new points on curve     (phase==1) 
         for (int phase=0; phase<=step; phase++) {
            for (Curve *curve=Gh.curves+icurve; curve; curve=curve->next) {
               int icurveequi=Gh.Number(curve);
               if (phase==0 && icurveequi!=icurve)
                  continue;
               int k0=jedge, k1;
               Edge *pe=BTh.edges+iedge;
//             GeometricalEdge *ong = ei.on;
               int iedgeequi=bcurve[icurveequi]/2, jedgeequi=bcurve[icurveequi]%2;
               int k0equi=jedgeequi, k1equi;
               Edge *peequi=BTh.edges+iedgeequi;
               GeometricalEdge *ongequi = peequi->on;
               double sNew=Lstep; 
               L = 0;
               long i=0;
               GeometricalVertex *GA0 = *(*peequi)[k0equi].on;
               Vertex *A0;
               A0 = GA0->to;  // the vertex in new mesh
               Vertex *A1;
               VertexOnGeom *GA1;
               Edge *PreviousNewEdge=0;
               assert(A0-_vertices>=0 && A0-_vertices <nbv);
               if (ongequi->Required()) {
                  GeometricalVertex *GA1 = *(*peequi)[1-k0equi].on;
                  A1 = GA1->to;
               }
               else {
                  for (;;) {
                     Edge &ee=*pe; 
                     Edge &eeequi=*peequi; 
                     k1 = 1 - k0; // next vertex of the edge 
                     k1equi = 1 - k0equi;
                     assert(pe && ee.on);
                     ee.on->SetMark();
                     Vertex &v0=ee[0], &v1=ee[1];
                     R2 AB = R2(v1) - R2(v0);
                     double L0=L, LAB;
                     LAB = LengthInterpole(v0.m,v1.m,AB);
                     L += LAB;
                     if (phase) {// computation of the new points
                        while ((i!=NbCreatePointOnCurve) && sNew <= L) { 
                           assert (sNew >= L0);
                           assert(LAB);
                           assert(_vertices && nbv<nbvx);
                           A1 = _vertices + nbv++;
                           GA1 = VerticesOnGeomEdge + NbVerticesOnGeomEdge;
                           Edge *e = edges + nbe++;
                           double se=(sNew-L0)/LAB;
                           assert(se>=0 && se<1.000000001);
                           se = abscisseInterpole(v0.m,v1.m,AB,se,1);
                           assert(se>=0 && se<=1);
                           se = k1 ? se : 1. - se;
                           se = k1==k1equi ? se : 1. - se;
                           VertexOnBThEdge[NbVerticesOnGeomEdge++] = VertexOnEdge(A1,&eeequi,se);
                           ongequi = Gh.ProjectOnCurve(eeequi,se,*A1,*GA1); 
                           A1->ReferenceNumber = eeequi.ref;
                           A1->DirOfSearch = NoDirOfSearch;
                           e->on = ongequi;
                           e->v[0] = A0;
                           e->v[1] = A1;
                           e->ref = eeequi.ref;
                           e->adj[0] = PreviousNewEdge;
                           if (PreviousNewEdge)
                              PreviousNewEdge->adj[1] = e;
                           PreviousNewEdge = e;
                           A0 = A1;
                           sNew += Lstep;
                           if (++i== NbCreatePointOnCurve)
                              break;
                        }
                     }
                     assert(ee.on->CurveNumber==ei.on->CurveNumber);
                     ei.on->CurveNumber = 0;
                     if (ee[k1].on->IsRequiredVertex()) {
                        assert(eeequi[k1equi].on->IsRequiredVertex());
                        GeometricalVertex *GA1 = *eeequi[k1equi].on;
                        A1 = GA1->to;
                        assert(A1-_vertices>=0 && A1-_vertices <nbv);
                        break;
                     }
                     if (!ee.adj[k1]) {
                        cerr << "Error adj edge " << BTh.Number(ee) << ", nbe = " << nbe 
                             << " Gh.vertices " << Gh._vertices 
                             << " k1 = " << k1 << " on=" << *ee[k1].on << endl;
                        cerr << ee[k1].on->gv-Gh._vertices << endl;
                     }
                     pe = ee.adj[k1]; // next edge
                     k0 = pe->Intersection(ee);
                     peequi = eeequi.adj[k1equi];  // next edge
                     k0equi = peequi->Intersection(eeequi);
                  }  
               }

               if (phase) { // construction of the last edge
                  Edge *e = edges + nbe++;
                  if (verbosity>10) 
                     cout << " Fin curve A1" << *A1 << " " << icurve << " <=> " << icurveequi <<"-----" <<
                  NbCreatePointOnCurve << " == " << i << endl;
                  e->on = ongequi;
                  e->v[0] = A0;
                  e->v[1] = A1;
                  e->ref = peequi->ref;
                  e->adj[0] = PreviousNewEdge;
                  e->adj[1] = 0;
                  if (PreviousNewEdge)
                     PreviousNewEdge->adj[1] = e;
                  PreviousNewEdge = e;
                  assert(i==NbCreatePointOnCurve);
               }
            }
            if (!phase) {
               long NbSegOnCurve=Max(long(L+0.5),long(1));
               Lstep = L/NbSegOnCurve; 
               NbCreatePointOnCurve = NbSegOnCurve - 1;
               for (Curve *curve=Gh.curves+icurve; curve; curve=curve->next) {
                  NbOfNewEdge += NbSegOnCurve;
                  NbOfNewPoints += NbCreatePointOnCurve;
               }
               if (verbosity>5)
                  cout << " NbSegOnCurve = " <<  NbSegOnCurve << " Lstep = " << Lstep
                       << " " << NbOfNewPoints << " NBPC= " << NbCreatePointOnCurve << endl;
            }
         }
      }
      if (step==0) {
         if (nbv+NbOfNewPoints > nbvx) {
            cerr << " Too many vertices on geometry " << nbv+NbOfNewPoints  << " >= " << nbvx << endl;
            MeshError(3);
         }
         edges = new Edge [NbOfNewEdge];
         if (NbOfNewPoints) {
            VerticesOnGeomEdge = new VertexOnGeom[NbOfNewPoints];
            NbVertexOnBThEdge = NbOfNewPoints;
            VertexOnBThEdge = new VertexOnEdge[NbOfNewPoints];
         }
         NbOfNewPoints = 0;
         NbOfNewEdge = 0;
      }
   }
   assert(nbe);
   delete [] bcurve;
   Insert();
   ForceBoundary();
   FindSubDomain();
   NewPoints(BTh,KeepBackVertices);
   CurrentTh = 0;
}


void Triangles::GeomToTriangles0(long inbvx) 
{
   Gh.NbRef++;
   long i, NbOfCurves=0, NbNewPoints, NbEdgeCurve=0;
   double lcurve, lstep, s;
   R2 AB;
   GeometricalVertex *a, *b;
   Vertex *va=NULL, *vb=NULL;
   GeometricalEdge *e=NULL;
   PreInit(inbvx);
   int background = &BTh != this;
   nbv = 0;
   NbVerticesOnGeomVertex = NbVerticesOnGeomEdge = 0;
   for (i=0; i<Gh.nbv; i++)
      if (Gh[i].Required() && Gh[i].IsThe())
         NbVerticesOnGeomVertex++;
   VerticesOnGeomVertex = new VertexOnGeom[NbVerticesOnGeomVertex];
   if (NbVerticesOnGeomVertex >= nbvx) {
      cerr << " Too many vertices on geometry " << NbVerticesOnGeomVertex << " >= "
           << nbvx << endl;
      MeshError(1);
   }
   for (i=0; i<Gh.nbv; i++)
      if (Gh[i].Required() && Gh[i].IsThe()) {
         if (nbv<nbvx)
            _vertices[nbv] = Gh[i];
         Gh[i].to = _vertices + nbv;
         VerticesOnGeomVertex[nbv]= VertexOnGeom(*Gh[i].to,Gh[i]);
         nbv++;
   }
// Method in 2 step:  0 and 1 
// 1) compute de nb of edge 
// 2) construct the edge
// generation of the curves
   assert(!edges);
// 2 step 
// --step=0 to compute the number of edges + alloc at end
// --step=1 to construct the edges
   for (int step=0; step<2; step++) { 
      long nbex=0;
      nbe = 0;
      long NbVerticesOnGeomEdge0=NbVerticesOnGeomEdge;
      Gh.UnMarkEdges();	
      NbOfCurves = 0;
      for (i=0; i<Gh.nbe; i++) {
         GeometricalEdge &ei = Gh.edges[i];
         if (!ei.Dup()) {
            for (int j=0; j<2; j++) {
               if (!ei.Mark() && ei[j].Required()) {
                  long nbvend=0;
                  Edge *PreviousNewEdge=0;
                  lstep = -1;   // to do not create points
                  if (ei.Required()) {
                     if (j==0) {
                        if (step==0)
                           nbe++;
                        else {
                           e = &ei;
                           a = ei(0)->The();
                           b = ei(1)->The();
                           assert(edges);
                           edges[nbe].v[0] = a->to;
                           edges[nbe].v[1] = b->to;;
                           edges[nbe].ref = e->ref;
                           edges[nbe].on = e;
                           edges[nbe].adj[0] = 0;
                           edges[nbe].adj[1] = 0;
                           nbe++;
                        }
                     }
                  }
                  else {
                     for (int kstep=0; kstep<=step; kstep++) {
                        PreviousNewEdge = 0;
                        NbNewPoints = NbEdgeCurve = 0;
                        assert(nbvend < nbvx); 
                        lcurve = 0;
                        s = lstep;
                        int k=j;
                        e = &ei;
                        a = ei(k)->The();
                        va = a->to;
                        e->SetMark();
                        for (;;) { 
                           k = 1 - k;
                           b = (*e)(k)->The();
                           AB = b->r - a->r;
                           Metric MA = background ? BTh.MetricAt(a->r) : a->m;
                           Metric MB =  background ? BTh.MetricAt(b->r) : b->m;
                           double ledge = (MA(AB) + MB(AB))/2;
                           const int MaxSubEdge = 10;
                           int NbSubEdge = 1;
                           double lSubEdge[MaxSubEdge];
                           R2 A, B;
                           if (ledge<1.5) 
                              lSubEdge[0] = ledge;
                           else {
                              NbSubEdge = Min(MaxSubEdge, (int)(ledge+0.5));
                              A = a->r;
                              Metric MAs=MA, MBs;
                              ledge = 0; 
                              double x=0, xstep=1./NbSubEdge;
                              for (int kk=0; kk<NbSubEdge; kk++, A=B, MAs=MBs) {
                                 x += xstep;
                                 B = e->F(k ? x : 1-x);
                                 MBs = background ? BTh.MetricAt(B) : Metric(1-x,MA,x,MB);
                                 AB = A - B;
                                 lSubEdge[kk] = (ledge += (MAs(AB)+MBs(AB))/2);
                              }
                           }
                           double lcurveb=lcurve+ledge;
                           while (lcurve<=s && s<=lcurveb && nbv<nbvend) {
                              double ss=s-lcurve;
//                            1) find the SubEdge containing ss by dichotomie
                              int kk0=-1, kk1=NbSubEdge-1, kkk;
                              double ll0=0, ll1=ledge, llk;
                              while (kk1-kk0>1) {
                                 if (ss < (llk=lSubEdge[kkk=(kk0+kk1)/2]))
                                    kk1 = kkk, ll1 = llk;
                                 else
                                    kk0 = kkk,ll0=llk;
                              }
                              assert(kk1!=kk0);
                              double sbb=(ss-ll0)/(ll1-ll0);
                              double bb=(kk1+sbb)/NbSubEdge, aa=1-bb;
                              vb = &_vertices[nbv++];
                              vb->m = Metric(aa,a->m,bb,b->m);
                              vb->ReferenceNumber = e->ref;
                              vb->DirOfSearch = NoDirOfSearch;
                              double abcisse = k ? bb : aa;
                              vb->r = e->F(abcisse);
                              VerticesOnGeomEdge[NbVerticesOnGeomEdge++] = VertexOnGeom(*vb,*e,abcisse);
                              s += lstep;
                              edges[nbe].v[0] = va;
                              edges[nbe].v[1] = vb;
                              edges[nbe].ref = e->ref;
                              if (e->ref>0)
                                 edges[nbe].ref = 0;
                              edges[nbe].on = e;
                              edges[nbe].adj[0] = PreviousNewEdge;
                              if (PreviousNewEdge)
                                 PreviousNewEdge->adj[1] = &edges[nbe];
                              PreviousNewEdge = edges + nbe;
                              nbe++;
                              va = vb;
                           }
                           lcurve = lcurveb;
                           e->SetMark();
                           a = b;
                           if (b->Required())
                              break;
                           int kprev=k;
                           k = e->SensAdj[kprev];
                           e = e->Adj[kprev];
                           assert(e);
                        }
                        vb = b->to;
                        NbEdgeCurve = Max((long)(lcurve +0.5),(long)1);
                        NbNewPoints = NbEdgeCurve - 1;
                        if (!kstep) {
                           NbVerticesOnGeomEdge0 += NbNewPoints;
                           NbOfCurves++;
                        } 
                        nbvend = nbv + NbNewPoints;  
                        lstep = lcurve / NbEdgeCurve;
                     }
                     if (edges) { 
                        edges[nbe].v[0] = va;
                        edges[nbe].v[1] = vb;
                        edges[nbe].ref = e->ref;
                        edges[nbe].on = e;
                        edges[nbe].adj[0] = PreviousNewEdge;
                        edges[nbe].adj[1] = 0;
                        if (PreviousNewEdge)
                           PreviousNewEdge->adj[1] = &edges[nbe];
                        nbe++;
                     }
                     else
                        nbe += NbEdgeCurve;
                  }
               }
            }
         }
      }
      if (!step) {
         assert(!edges);
         assert(!VerticesOnGeomEdge);
         edges = new Edge[nbex=nbe];
         if (NbVerticesOnGeomEdge0)
            VerticesOnGeomEdge = new VertexOnGeom[NbVerticesOnGeomEdge0];
         assert(edges);
         assert(VerticesOnGeomEdge || NbVerticesOnGeomEdge0==0);
         NbVerticesOnGeomEdge0 = NbVerticesOnGeomEdge;       
      }
      else
         assert(NbVerticesOnGeomEdge==NbVerticesOnGeomEdge0);
      assert((nbex=nbe));
   }

   Insert();
   ForceBoundary();
   FindSubDomain();
   NewPoints(*this,0);
   CurrentTh = 0;
}


Edge** Triangles::MakeGeometricalEdgeToEdge()
{
   assert(Gh.nbe);
   Edge **e = new Edge* [Gh.nbe];
   long i;
   for (i=0; i<Gh.nbe ;i++)
      e[i] = NULL;
   for (i=0; i<nbe; i++) {
      Edge *ei = edges+i;
      GeometricalEdge *on = ei->on; 
      e[Gh.Number(on)] = ei;
   }
   for (i=0; i<nbe ;i++) {
      for (int ii=0; ii<2; ii++) { 
         Edge *ei = edges + i;
         GeometricalEdge *on = ei->on;
         int j=ii;
         while (!(*on)[j].Required()) {
            Adj(on,j);
            j = 1 - j;
            if (e[Gh.Number(on)])
               break;
            e[Gh.Number(on)] = ei; 
         }
      }
   }
   int kk=0;
   for (i=0; i<Gh.nbe; i++) {
      if (!e[i])
         if (kk++<10)
            cerr << "Error: The geometrical edge " << i << " is on no edge curve = "
                 << Gh.edges[i].CurveNumber << " s0 " << Gh.Number(Gh.edges[i][0])
                 << " s1 " << Gh.Number(Gh.edges[i][1]) << endl;
   }
   if (kk)
      MeshError(997);
   return e;
}


Triangles::~Triangles() 
{
   assert(NbRef<=0);
   if (CurrentTh==this)
      CurrentTh = 0;
   if (verbosity>10)
      cout << " ~Triangles "<< this <<" "<< identity << endl;
   if (_vertices)
      delete [] _vertices;
   if (edges)
      delete [] edges;
   if (triangles)
      delete [] triangles;
   if (_quadtree)
      delete _quadtree;
   if (ordre)
      delete [] ordre;
   if (subdomains)
      delete [] subdomains;
   if (VerticesOnGeomEdge)
      delete [] VerticesOnGeomEdge;
   if (VerticesOnGeomVertex)
      delete [] VerticesOnGeomVertex;
   if (name)
      delete [] name;
   if (identity)
      delete [] identity;
   if (VertexOnBThVertex)
      delete [] VertexOnBThVertex;
   if (VertexOnBThEdge)
      delete [] VertexOnBThEdge;
   Geometry *gh=&Gh;
   if (gh) {
      if (gh->NbRef>0)
         gh->NbRef--;
      else if (gh->NbRef==0)
         delete gh;
   }
   Triangles *bth=&BTh;
   if (bth && (bth!=this)) {
      if (bth->NbRef>0)
         bth->NbRef--;
      else if (bth->NbRef==0)
         delete bth;
   }
   PreInit(0);
}


void Triangles::SetIntCoor(const char* strfrom)
{
   pmin = pmax = _vertices[0].r;

// Compute extrema for vertices
   for (long i=0; i<nbv; i++) {
      pmin.x = Min(pmin.x,_vertices[i].r.x);
      pmin.y = Min(pmin.y,_vertices[i].r.y);
      pmax.x = Max(pmax.x,_vertices[i].r.x);
      pmax.y = Max(pmax.y,_vertices[i].r.y);
   }
   R2 DD=(pmax-pmin)*0.05;
   pmin = pmin - DD;
   pmax = pmax + DD; 
   coefIcoor = (MaxICoor)/(Max(pmax.x-pmin.x,pmax.y-pmin.y));
   assert(coefIcoor >0);

// Generation of integer coord  
   for (long i=0; i<nbv; i++)
      _vertices[i].i = toI2(_vertices[i].r);

// computation of the det 
   int Nberr=0;
   for (long i=0; i<nbt; i++) {
      Vertex *v0=&triangles[i][0], *v1=&triangles[i][1], *v2=&triangles[i][2];
      if (v0 && v1 && v2) {
         triangles[i].det = det(*v0,*v1,*v2);
         if (triangles[i].det<=0 && Nberr++<10) {
            if (Nberr==1) {
               if (strfrom)
                  cerr << "+++ Fatal Error " << strfrom << "(SetInCoor) Error: area of Triangle < 0" << endl; 
               else 
                  cerr << "+++  Fatal Error Triangle (in SetInCoor) area of Triangle < 0" << endl;
            }
            cerr << " Triangle " << i << "  det  (I2) = " << triangles[i].det ;
            cerr << " (R2) " << Det(v1->r-v0->r,v2->r-v0->r);
            cerr << "; The 3 vertices " << endl;
            cerr << Number(*v0) << " " << Number(*v1) << " " << Number(*v2) << ": ";
            cerr << v0->r << v1->r << v2->r << " ; ";
            cerr << v0->i << v1->i << v2->i << endl;
         }
      }
      else
         triangles[i].det = -1; // Null triangle; 
   }
   if (Nberr)
      MeshError(899);
}


void Triangles::FillHoleInMesh() 
{
   Triangles *OldCurrentTh=CurrentTh;
   CurrentTh = this;
   {
      long i;
      if (verbosity>2)
         cout << " -- FillHoleInMesh: Nb of vertices =" << nbv 
              << " Pmin = "<< pmin << " Pmax = "<< pmax << endl;
      assert(ordre);
      for (i=0; i<nbv; i++) 
         ordre[i]= 0;
      NbSubDomains = 0;

//    Generation of the adjacence of the triangles
      SetOfEdges4 *edge4= new SetOfEdges4(nbt*3,nbv);
      long *st = new long[nbt*3];
      for (i=0; i<nbt*3; i++)
         st[i] = -1;
      long kk=0;
      for (i=0; i<nbe; i++)
         kk += (i==edge4->addtrie(Number(edges[i][0]),Number(edges[i][1])));
      if (kk != nbe) { 
         cerr << " Some Double edge in the mesh, the number is " << kk-nbe << endl;
         MeshError(1002);
      }
      for (i=0; i<nbt; i++) {
         for (int j=0; j<3; j++) {
            long k=edge4->addtrie(Number(triangles[i][VerticesOfTriangularEdge[j][0]]),
                                  Number(triangles[i][VerticesOfTriangularEdge[j][1]]));
            long invisible = triangles[i].Hidden(j);
            if (st[k]==-1)
               st[k] = 3*i+j;
            else if (st[k]>=0) {
               assert(!triangles[i].TriangleAdj(j) && !triangles[st[k]/3].TriangleAdj(int(st[k]%3)));
               triangles[i].SetAdj2(j,triangles + st[k]/3,int(st[k]%3));
               if (invisible)
                  triangles[i].SetHidden(j);
               if (k<nbe)
                  triangles[i].SetLocked(j);
               st[k] = -2 - st[k];
            }
            else {
               cerr << " The edge (" 
                    << Number(triangles[i][VerticesOfTriangularEdge[j][0]])
                    << " , " << Number(triangles[i][VerticesOfTriangularEdge[j][1]])
                    << " ) is in more than 2 triangles " << k <<endl;
               cerr << " Edge " << j << " Of Triangle " << i << endl;
               cerr << " Edge " << (-st[k]+2)%3 << " Of Triangle " << (-st[k]+2)/3 << endl;
               cerr << " Edge " << triangles[(-st[k]+2)/3].NuEdgeTriangleAdj(int((-st[k]+2)%3))
                    << " Of Triangle " << Number(triangles[(-st[k]+2)/3].TriangleAdj(int((-st[k]+2)%3))) << endl;
               MeshError(9999);
            }
         }
      }
      if (verbosity>5) {
         cout << "    On Mesh " << name << endl;
         cout << "    - The number of Vertices   = " << nbv << endl;
         cout << "    - The number of Triangles  = " << nbt << endl;
         cout << "    - The number of given edges = " << nbe << endl;
         cout << "    - The number of all edges = " << edge4->nb() << endl;
         cout << "    - The Euler number = 1-Nb Of Hole = " << nbt-edge4->nb()+nbv << endl;
      }

//    Check consistency of edge[].adj and the geometrical required vertex
      long k=0;
      for (i=0; i<edge4->nb(); i++) {
         if (st[i] >=0) {
            if (i<nbe) {
               long i0=edge4->i(i);
               ordre[i0] = _vertices + i0;
               long i1=edge4->j(i);
               ordre[i1] = _vertices + i1;
            }
            else {
               k++;
               if (verbosity>20 && k<20) {
                  long i0=edge4->i(i), i1=edge4->j(i);
                  cerr << " Lose boundary edges " << i << " : " << i0 << " " << i1 << endl;
               }
            }
         }
      }
      if (k != 0) {
         if (verbosity>20) {
            cout << " The given edges are " << endl;
            for (int i=0; i<nbe; i++)
               cout << " Edge " << i << " : " <<  Number(edges[i][0]) << " " <<  Number(edges[i][1]) 
                    << " " << edges[i].ref << endl; 
         }
         cerr << k << " boundary edges are not defined as edges" << endl;
         MeshError(9998);
      }

//    Generation of the mesh with boundary points   
      long nbvb = 0;
      for (i=0; i<nbv; i++) {
         _vertices[i].t = 0;
         _vertices[i].vint = 0;
         if (ordre[i]) 
            ordre[nbvb++] = ordre[i];
      }
      Triangle *savetriangles=triangles;
      long savenbt=nbt, savenbtx=nbtx;
      SubDomain *savesubdomains = subdomains;
      subdomains = 0;
      long Nbtriafillhole = 2*nbvb;
      Triangle *triafillhole = new Triangle[Nbtriafillhole];
      if (verbosity>9)
         cout << " Nbtriafillhole triafillhole*" << triafillhole << endl; 
      triangles =  triafillhole;
      nbt = 2;
      nbtx = Nbtriafillhole;
      for (i=2; det(ordre[0]->i,ordre[1]->i,ordre[i]->i) == 0;) {
         if (++i >= nbvb) {
            cerr << "FillHoleInMesh: All vertices are aligned! " << nbvb << endl;
            MeshError(998);
         }
      }
      Exchange(ordre[2],ordre[i]);

      Vertex *v0=ordre[0], *v1=ordre[1];
      triangles[0](0) = 0; triangles[0](1) = v0; triangles[0](2) = v1;
      triangles[1](0) = 0; triangles[1](2) = v0; triangles[1](1) = v1;
      const int e0=OppositeEdge[0], e1=NextEdge[e0], e2=PreviousEdge[e0];
      triangles[0].SetAdj2(e0,&triangles[1],e0);
      triangles[0].SetAdj2(e1,&triangles[1],e2);
      triangles[0].SetAdj2(e2,&triangles[1],e1);
      triangles[0].det = -1;  // fake triangles
      triangles[1].det = -1;  // fake triangles
      triangles[0].SetTriangleContainingTheVertex();
      triangles[1].SetTriangleContainingTheVertex();
      triangles[0].link = &triangles[1];
      triangles[1].link = &triangles[0];
      if (!_quadtree)
         delete _quadtree; // ->ReInitialise();
      _quadtree = new QuadTree(this,0);
      _quadtree->Add(*v0);
      _quadtree->Add(*v1);
	
//    Add vertices one by one
      for (long icount=2; icount<nbvb; icount++) {
         Vertex *vi = ordre[icount];
         Icoor2 dete[3];
         Triangle *tcvi = FindTriangleContaining(vi->i,dete);
         _quadtree->Add(*vi);
         Add(*vi,tcvi,dete);
      }

//    Inforce the boundary 
      TriangleAdjacent ta(0,0);
      long nbloss=0, knbe=0;
      for (i=0; i<nbe; i++) {
         Vertex &a=edges[i][0],&b=edges[i][1];
         if (a.t && b.t) {
            knbe++;
            if (ForceEdge(a,b,ta)<0)
               nbloss++;
         }
      }
      if (nbloss) {
         cerr << " We lost some  " << nbloss << " edges other " << knbe << endl;
         MeshError(1100);
      }
      FindSubDomain(1);

//    Remove all the hole 
//    remove all the good sub domain
      long krm=0;
      for (i=0; i<nbt; i++) {
         if (triangles[i].link) {
            krm++;
            for (int j=0; j<3; j++) {
               TriangleAdjacent ta =  triangles[i].Adj(j);
               Triangle &tta = *(Triangle *) ta;
               if (!tta.link) {
                  int ja = ta;
                  Vertex *v0=ta.EdgeVertex(0), *v1=ta.EdgeVertex(1);
                  long k = edge4->addtrie(v0?Number(v0):nbv,v1? Number(v1):nbv);
                  assert(st[k]>=0);
                  tta.SetAdj2(ja,savetriangles+st[k]/3,int(st[k]%3));
                  ta.SetLock();
                  st[k] = -2 - st[k]; 
               }
            }
         }
      }
      long NbTfillHoll=0;
      for (i=0; i<nbt; i++) {
         if (triangles[i].link) {
            triangles[i] = Triangle((Vertex *) NULL,(Vertex *) NULL,(Vertex *) NULL);
            triangles[i].color = -1;
         }
         else
            triangles[i].color = savenbt + NbTfillHoll++;
      }
      assert(savenbt+NbTfillHoll <= savenbtx );
      for (i=0;i <nbt; i++) {
         if (triangles[i].color>=0) {
            savetriangles[savenbt] = triangles[i];
            savetriangles[savenbt].link = 0;
            savenbt++;
         }
      }
      k = 0;
      Triangle *tmax = triangles + nbt;
      for (i=0; i<savenbt; i++) { 
         Triangle &ti = savetriangles[i];
         for (int j=0; j<3; j++) {
            Triangle *ta=ti.TriangleAdj(j);
            int aa=ti.NuEdgeTriangleAdj(j);
            int lck=ti.Locked(j);
            if (!ta)
               k++; // bug 
            else if (ta>=triangles && ta<tmax) {
               ta = savetriangles + ta->color;
               ti.SetAdj2(j,ta,aa);
               if (lck)
                  ti.SetLocked(j);
            }
         }
      }
 
//    Restore triangles;
      nbt = savenbt;
      nbtx = savenbtx;
      delete [] triangles;
      delete [] subdomains;
      triangles = savetriangles;
      subdomains = savesubdomains;
      if (k) {
         cerr << "Error Nb of triangles edge alone = " << k << endl;
         MeshError(9997);
      }
      FindSubDomain();
      delete edge4;
      delete [] st;
      for (i=0; i<nbv; i++)
         _quadtree->Add(_vertices[i]);
      SetVertexFieldOn();
      for (i=0; i<nbe; i++) {
         if (edges[i].on) { 
            for (int j=0; j<2; j++) {
               if (!edges[i].adj[j]) {
                  if (!edges[i][j].on->IsRequiredVertex()) {
                     cerr << " Erreur adj et sommet requis edges [" << i <<  "][ " << j << "]= "
                          <<  Number(edges[i][j]) << ": "  << " on = " << Gh.Number(edges[i].on) ;
                     if (edges[i][j].on->OnGeomVertex())
                        cerr << " vertex " << Gh.Number(edges[i][j].on->gv);
                     else if (edges[i][j].on->OnGeomEdge())
                        cerr << "Edges " << Gh.Number(edges[i][j].on->ge);
                     else
                        cerr << " = " << edges[i][j].on;
                     cerr << endl;
                  }
               }
            }
         }
      }
   }
   CurrentTh = OldCurrentTh;
}


Triangles::Triangles(Triangles& Th,
                     Geometry*  pGh,
                     Triangles* pBth,
                     long       nbvxx)
           : Gh(*(pGh?pGh:&Th.Gh)), BTh(*(pBth?pBth:this))
{
   Gh.NbRef++;
   nbvxx = Max(nbvxx,Th.nbv); 
// do all the allocation to be sure all the pointer exists

   char *cname = 0;
   if (Th.name) {
      cname = new char[strlen(Th.name)+1];
      strcpy(cname,Th.name);
   }
   PreInit(nbvxx,cname);// to make the allocation

// copy of triangles
   nt = Th.nt;
   nbv = Th.nbv;
   nbt = Th.nbt;
   nbiv = Th.nbiv;
   nbe = Th.nbe;
   NbSubDomains = Th.NbSubDomains;
   NbOutT = Th.NbOutT;
   NbOfQuad = Th.NbOfQuad;
   NbOfSwapTriangle = 0;
   NbVerticesOnGeomVertex = Th.NbVerticesOnGeomVertex;
   if (NbVerticesOnGeomVertex)
      VerticesOnGeomVertex = new VertexOnGeom[NbVerticesOnGeomVertex];
   NbVerticesOnGeomEdge = Th.NbVerticesOnGeomEdge;
   if (NbVerticesOnGeomEdge)
      VerticesOnGeomEdge = new VertexOnGeom[NbVerticesOnGeomEdge] ;
   if (&BTh == &Th.BTh) {
      BTh.NbRef++;
      NbVertexOnBThVertex = Th.NbVertexOnBThVertex;
      if (NbVertexOnBThVertex)
         VertexOnBThVertex = new VertexOnVertex[NbVertexOnBThVertex];
      NbVertexOnBThEdge = Th.NbVertexOnBThEdge;
      if (NbVertexOnBThEdge)
         VertexOnBThEdge = new VertexOnEdge[NbVertexOnBThEdge];
   }
   else {
      BTh.NbRef++;
      NbVertexOnBThVertex = 0;
      VertexOnBThVertex = 0;
      NbVertexOnBThEdge = 0;
      VertexOnBThEdge = 0;
   }
   if (nbe)
      edges = new Edge[nbe];
   if (NbSubDomains)
      subdomains = new SubDomain[NbSubDomains];
   pmin = Th.pmin;
   pmax = Th.pmax;
   coefIcoor = Th.coefIcoor;
   for (long i=0; i<nbt; i++)
      triangles[i].Set(Th.triangles[i],Th,*this);
   for (long i=0; i<nbe; i++)
      edges[i].Set(Th,i,*this);
   for (long i=0; i<nbv; i++)
      _vertices[i].Set(Th._vertices[i],Th,*this);
   for (long i=0; i<NbSubDomains; i++)
      subdomains[i].Set(Th,i,*this);
   for (long i=0; i<NbVerticesOnGeomVertex; i++)
      VerticesOnGeomVertex[i].Set(Th.VerticesOnGeomVertex[i],Th,*this);
   for (long i=0; i<NbVerticesOnGeomEdge; i++)
      VerticesOnGeomEdge[i].Set(Th.VerticesOnGeomEdge[i],Th,*this);
   _quadtree = 0;
}


long Triangle::Optim(short i,
                     int   koption)
{
   long NbSwap=0;
   Triangle *t=this;
   int k=0, j=OppositeEdge[i], jp=PreviousEdge[j];
   Triangle *tp=at[jp];
   jp = aa[jp]&3;
   do {
      while (t->swap(j,koption)) {
         NbSwap++;
         assert(k<20000);
         k++;
         t = tp->at[jp];
         j = NextEdge[tp->aa[jp]&3];
      }
      tp = t;
      jp = NextEdge[j];
      t = tp->at[jp];      // set unchange t qnd j for previous triangles
      j = NextEdge[tp->aa[jp]&3];
   } while(t != this);
   return NbSwap;
}


void Triangles::SmoothingVertex(int    nbiter,
                                double omega)
{
// If quadtree exists remove it end reconstruct
   if (_quadtree)
      delete _quadtree;
   _quadtree = NULL;
   ReMakeTriangleContainingTheVertex();
   Triangle vide; // a triangle to mark the boundary vertex
   Triangle ** tstart= new Triangle* [nbv];
   long i, j, k;
// Note: If Background == Triangle, we cannot use fast research
   if (this==&BTh) {
      for (i=0; i<nbv; i++)
         tstart[i] = _vertices[i].t;
   }     
   else {
      for (i=0; i<nbv; i++)
         tstart[i] = 0;
   }
   for (j=0; j<NbVerticesOnGeomVertex; j++) 
      tstart[ Number(VerticesOnGeomVertex[j].mv)] = &vide;
   for (k=0; k<NbVerticesOnGeomEdge; k++) 
      tstart[Number(VerticesOnGeomEdge[k].mv)] = &vide;
   if (verbosity>2) 
      cout << " -- SmoothingVertex: nb Iteration = " << nbiter << " Omega = " << omega << endl;
   for (k=0; k<nbiter; k++) {
      long i,NbSwap =0;
      double delta=0;
      for (i=0; i<nbv; i++)
         if (tstart[i] != &vide)    // not a boundary vertex 
	    delta = Max(delta,_vertices[i].Smoothing(*this,BTh,tstart[i],omega));
      if (!NbOfQuad) {
         for (i=0; i<nbv; i++)
            if (tstart[i] != &vide) // not a boundary vertex 
               NbSwap += _vertices[i].Optim(1);
      }
      if (verbosity>3)
         cout << "    Move max = " <<  sqrt(delta) << " iteration = " 
              << k << " Nb of Swap = " << NbSwap << endl;
   }
   delete [] tstart;
   if (_quadtree)
      _quadtree = new QuadTree(this);
}


void Triangles::MakeQuadTree()
{  
   if (verbosity>8)
      cout << "      MakeQuadTree" << endl;
   if (!_quadtree)
      _quadtree = new QuadTree(this);
}


void Triangles::ShowHistogram() const
{
   const long kmax=10;
   const double llmin=0.5, llmax=2;
   const double lmin=log(llmin), lmax=log(llmax), delta=kmax/(lmax-lmin);
   long histo[kmax+1];
   long i, it, k, nbedges=0;
   for (i=0; i<=kmax; i++)
      histo[i] = 0;
   for (it=0; it<nbt; it++) {
      if (triangles[it].link) {
         for (int j=0; j<3; j++) {
            Triangle *ta=triangles[it].TriangleAdj(j);
            if (!ta || !ta->link || Number(ta) >= it) {
               Vertex *vP=&triangles[it][VerticesOfTriangularEdge[j][0]];
               Vertex *vQ=&triangles[it][VerticesOfTriangularEdge[j][1]];
               if (!vP || !vQ)
                  continue;
               R2 PQ = vQ->r - vP->r;
               double l=log(LengthInterpole(*vP,*vQ,PQ));
               nbedges++;
               k = int((l-lmin)*delta);
               k = Min(Max(k,0L),kmax);
               histo[k]++;
            }
         }
      }
   }
   cout << " -- Histogram of the unit mesh,  nb of edges" << nbedges << endl <<endl;
   cout << "        length of edge in   | % of edge  | Nb of edges " << endl;
   cout << "        ------------------- | ---------- | ----------- " << endl;
   for (i=0; i<=kmax; i++) {
      cout << "    ";
      cout.width(10);
      if (i==0)
         cout << " 0 ";
      else
         cout << exp(lmin+i/delta);
      cout.width();
      cout << ",";
      cout.width(10);
      if (i==kmax)
         cout << " +infty ";
      else
         cout << exp(lmin+(i+1)/delta);
      cout.width();
      cout << "   |   ";
      cout.precision(4);
      cout.width(6);
      cout <<  ((long)((10000.0 * histo[i])/nbedges))/100.0;
      cout.width();
      cout.precision();
      cout << "   |   " << histo[i] << endl;
   }
   cout << "        ------------------- | ---------- | ----------- " << endl << endl;
}


int  Triangles::Crack()
{
   assert(NbCrackedEdges==0 || NbCrackedVertices>0); 
   for (int i=0; i<NbCrackedEdges; i++)
      CrackedEdges[i].Crack();
   return NbCrackedEdges;
}


int Triangles::UnCrack() 
{ 
   assert(NbCrackedEdges==0 || NbCrackedVertices>0);
   for (int i=0; i<NbCrackedEdges; i++)
      CrackedEdges[i].UnCrack();
   return NbCrackedEdges;
}


int Triangles::CrackMesh()
{
   Triangles *CurrentThOld = CurrentTh;
// compute the number of cracked edges
   int i, k;
   for (k=i=0; i<nbe; i++)
      if (edges[i].on->Cracked())
         k++;
   if (k==0)
      return 0;
   CurrentTh = this;
   cout << " Nb of Cracked Edges = " << k << endl;
   NbCrackedEdges = k;
   CrackedEdges = new CrackedEdge[k];
   Edge *e = new Edge[nbe+k];
   for (i=0; i<nbe; i++) 
      e[i] = edges[i];
   delete edges;
   edges = e;

   const int nbe0=nbe;
   for (k=i=0; i<nbe0; i++) {
      if (edges[i].on->Cracked()) {
         e[nbe] = e[i];
         e[nbe].v[0] = e[i].v[1];
         e[nbe].v[1] = e[i].v[0];
         e[nbe].on = e[i].on->link;
         CrackedEdges[k++] = CrackedEdge(edges,i,nbe);
         nbe++;
      }
   }
   ReMakeTriangleContainingTheVertex();
   int nbcrakev=0;
   Vertex *vlast = _vertices + nbv;
   Vertex *vend = _vertices + nbvx;
   for (int iv=0; iv<nbv; iv++) {
      Vertex &v = _vertices[iv];
      Vertex *vv = &v;  
      int kk=0, kc=0, kkk=0;
      Triangle *tbegin=v.t;
      int i=v.vint;       
      assert(tbegin && (i>=0) && (i<3));
      TriangleAdjacent ta(tbegin,EdgesVertexTriangle[i][0]);
      int k=0;
      do {
         int kv = VerticesOfTriangularEdge[ta][1];
         k++; 
         Triangle *tt(ta);
         if (ta.Cracked()) {   
            if (kk== 0)
               tbegin = ta,kkk = 0;  //  begin by a cracked edge  => restart                
            if (kkk) {
               kc = 1;
               vv = vlast++;
               kkk = 0;
            } // new vertex if use 
	    kk++;// number of cracked edge view                 
	 }
         if (tt->link) { // if good triangles store the value 
            assert(Number(tt) < nt);
            (*tt)(kv) = vv; //   Change the vertex of triangle 
            if (vv<vend) {
               *vv = v;
               vv->ReferenceNumber = iv;
            } // copy the vertex value + store the old vertex number in ref 
            kkk++;
         }
         else if (kk) { // crack + boundary 
            if (kkk) {
               kc = 1;
               vv = vlast++;
               kkk = 0;
            } // new vertex if use 
         }
         ta = Next(ta).Adj(); 
      } while ( (tbegin != ta)); 
      assert(k);
      if (kc)
         nbcrakev++;
   }

   if (nbcrakev) 
      for (int iec=0; iec<NbCrackedEdges; iec++)
         CrackedEdges[iec].Set();

// set the ref
   cout << " set the ref " << endl;
   NbCrackedVertices = nbcrakev;
   nbv = vlast - _vertices;
   int nbnewv =  nbv - nbv; // nb of new vrtices 
   if (nbcrakev && verbosity>1)
      cout << " Nb of cracked vertices: " << nbcrakev 
           << " Nb of created vertices: " << nbnewv << endl;
   if (nbnewv) {
      long n = nbnewv + NbVerticesOnGeomVertex;
      long i, j, k;
      VertexOnGeom *vog = new VertexOnGeom[n];
      for (i=0; i<NbVerticesOnGeomVertex; i++) 
         vog[i] = VerticesOnGeomVertex[i];
      delete [] VerticesOnGeomVertex;
      VerticesOnGeomVertex = vog;
      Vertex *LastOld = _vertices + nbv - nbnewv;
      for (int iec=0; iec< NbCrackedEdges; ++iec) {
         for (k=0; k<2; k++) {
            Edge &e = *(k ? CrackedEdges[iec].a.edge : CrackedEdges[iec].b.edge);
            for (j=0; j<2; j++) {
               Vertex *v = e(j);
               if (v>=LastOld) {
                  long old = v->ReferenceNumber;
                  long i = (v - LastOld);
                  vog[i] = vog[old];
               }
            }
         }
      }
      NbVerticesOnGeomVertex = n;
   }
   SetVertexFieldOn();
   if (vlast >= vend) {
      cerr << " Not enough vertices to crack the mesh we need " << nbv << " vertices " << endl;
      MeshError(555);
   }
   cout << "  NbCrackedVertices " <<  NbCrackedVertices << endl;
   CurrentTh = CurrentThOld;
   return NbCrackedVertices;
}


Triangles::Triangles(const Triangles& Tho,
                     const int*       flag,
                     const int*       bb)
          : Gh(*(new Geometry())), BTh(*this)
{
   char cname[] = "trunc";
   int i, k, itadj, kt=0;
   int *kk = new int [Tho.nbv];
   long *reft = new long[Tho.nbt];
   long nbInT = Tho.ConsRefTriangle(reft);
   long *refv = new long[Tho.nbv];
   for (i=0; i<Tho.nbv; i++) {
      kk[i] = -1;
      refv[i] = 0;
   }
   int nbNewBedge=0;
   for (i=0; i<Tho.nbt; i++) {
      if (reft[i] >=0 && flag[i]) {
         const Triangle &t=Tho.triangles[i];
         kt++;
         kk[Tho.Number(t[0])] = kk[Tho.Number(t[1])] = kk[Tho.Number(t[2])] = 1;
         itadj = Tho.Number(t.TriangleAdj(0));
         if (reft[itadj]>=0 && !flag[itadj]) {
            nbNewBedge++;
            refv[Tho.Number(t[VerticesOfTriangularEdge[0][0]])] = bb[i];
            refv[Tho.Number(t[VerticesOfTriangularEdge[0][1]])] = bb[i];
         }
         itadj = Tho.Number(t.TriangleAdj(1));
         if (reft[itadj]>=0 && !flag[itadj]) {
            nbNewBedge++;
            refv[Tho.Number(t[VerticesOfTriangularEdge[1][0]])] = bb[i];
            refv[Tho.Number(t[VerticesOfTriangularEdge[1][1]])] = bb[i];
         }
         itadj = Tho.Number(t.TriangleAdj(2));
         if (reft[itadj]>=0 && !flag[itadj]) {
            nbNewBedge++;
            refv[Tho.Number(t[VerticesOfTriangularEdge[2][0]])] = bb[i];
            refv[Tho.Number(t[VerticesOfTriangularEdge[2][1]])] = bb[i];
         }
      }
   }
   k = 0;
   for (i=0; i<Tho.nbv; i++)
      if (kk[i]>=0) 
         kk[i] = k++;
   cout << " number of vertices " << k << " remove = " << Tho.nbv - k << endl;
   cout << " number of triangles " << kt << " remove = " << nbInT-kt << endl;
   cout << " number of new boundary edges " << nbNewBedge << endl;
   long inbvx=k;
   PreInit(inbvx,cname);
   for (i=0; i<Tho.nbv; i++) {
      if (kk[i]>=0) {
         _vertices[nbv] = Tho._vertices[i];
         if (!_vertices[nbv].ref())
            _vertices[nbv].ReferenceNumber = refv[i];
         nbv++;
      }
   }
   assert(inbvx==nbv);
   for (i=0; i<Tho.nbt; i++) {
      if (reft[i]>=0 && flag[i]) {
         const Triangle &t = Tho.triangles[i];
         int i0=Tho.Number(t[0]), i1=Tho.Number(t[1]), i2=Tho.Number(t[2]);
         assert(i0>=0 && i1>=0 && i2>=0);
         assert(i0<Tho.nbv && i1<Tho.nbv && i2<Tho.nbv);
         triangles[nbt] = Triangle(this,kk[i0],kk[i1],kk[i2]);
         triangles[nbt].color = Tho.subdomains[reft[i]].ref; 
         nbt++;           
      }
   }
   assert(kt==nbt);
   if (nbt==0 && nbv==0) {
      cout << "Error: all triangles were remove." << endl;
      MeshError(999);
   }
   delete [] kk;
   delete [] reft;
   delete [] refv;
   double cutoffradian = 10.0/180.0*Pi;
   ConsGeometry(cutoffradian);
   Gh.AfterRead(); 
   SetIntCoor();
   FillHoleInMesh();
   assert(NbSubDomains);
   assert(subdomains[0].head && subdomains[0].head->link);       
}


Triangle* Triangles::FindTriangleContaining(const I2& B,
                                            Icoor2    dete[3],
                                            Triangle* tstart) const
{ // in: B
  // out: t
  // out: dete[3]
  // t the triangle with vertices s0,s1,s2
  // in dete[3] = { det(B,s1,s2), det(s0,B,s2), det(s0,s1,B) }
  // with det(a,b,c) = -1 if one of 3 vertices a,b,c is NULL 
   Triangle *t = NULL;
   if (tstart)
      t = tstart;
   else {
      assert(_quadtree);
      Vertex *a = _quadtree->NearestVertex(B.x,B.y);
      if (!a || !a->t) {
         if (a) {
            cerr << "Attention PB TriangleContainingTheVertex vertex number = " 
                 << Number(a) << endl;
            cerr << "We forget a call to ReMakeTriangleContainingTheVertex" << endl;
         }
         cerr << "Pb with " << B << toR2(B) << endl;
         MeshError(7777);
      }
      assert(a>=_vertices && a<_vertices+nbv);
      t = a->t;
      assert(t>=triangles && t<triangles+nbt);
   }
   Icoor2 detop;
   int kk=0;            // number of test triangles
   while (t->det<0) {   // the initial triangle is outside
      int k0=(*t)(0) ? (((*t)(1) ? ((*t)(2) ? -1 : 2) : 1)) : 0;
      assert(k0>=0);    // k0 the NULL vertex
      int k1=NextVertex[k0], k2=PreviousVertex[k0];
      dete[k0] = det(B,(*t)[k1],(*t)[k2]);
      dete[k1] = dete[k2] = -1;
      if (dete[k0]>0)   // outside B
         return t;
      t = t->TriangleAdj(OppositeEdge[k0]);
      assert(kk<2);
      kk++;
   }
   int jj = 0;
   detop = det(*(*t)(VerticesOfTriangularEdge[jj][0]),
               *(*t)(VerticesOfTriangularEdge[jj][1]),B); 
   while (t->det>0) {
      assert(kk++<2000);
      int j = OppositeVertex[jj];
      dete[j] = detop;  // det(*b,*s1,*s2);
      int jn=NextVertex[j], jp=PreviousVertex[j];
      dete[jp] = det(*(*t)(j),*(*t)(jn),B);
      dete[jn] = t->det - dete[j] - dete[jp];
      int k=0, ii[3];
      if (dete[0]<0)
         ii[k++] = 0; 
      if (dete[1]<0)
         ii[k++] = 1;
      if (dete[2]<0)
         ii[k++] = 2;
      if (k==0)
         break;
      if (k==2 && BinaryRand())
         Exchange(ii[0],ii[1]);
      assert(k<3);
      TriangleAdjacent t1=t->Adj(jj=ii[0]);
      if (t1.det()<0 && k==2)
         t1 = t->Adj(jj=ii[1]);
      t = t1;
//    for optimization we know the -det[OppositeVertex[j]]
      j = t1; 
      detop = -dete[OppositeVertex[jj]];
      jj = j;
   }
   if (t->det<0) {
//    outside triangle 
      dete[0] = dete[1] = dete[2] = -1;
      dete[OppositeVertex[jj]] = detop;
   }
   return t;
}

} /* namespace bamg */
