// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY: Bamg: Bidimensional Anisotrope Mesh Generator
// RELEASE: 0 
// USAGE  : You may copy freely these files and use it for    
//          teaching or research. These or part of these may   
//          not be sold or used for a commercial purpose with- 
//          out our consent : fax (33) 1 39 63 55 14       
// --------------------------------------------------------
// AUTHOR:   F. Hecht,    
// ORG    :  INRIA
// E-MAIL :   Frederic.Hecht@Inria.fr   
//
// ORIG-DATE:     frev 98
//---------------------------------------------------------
//  to make quad 
// -------------------

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "mesh/bamg/Meshio.h"
#include "mesh/bamg/Mesh2.h"
#include "mesh/bamg/QuadTree.h"
#include "mesh/bamg/SetOfE4.h"

namespace bamg {

static const Direction NoDirOfSearch=Direction();

class DoubleAndlong
{
 public:
    double q;
    long i3j;
    int operator<(DoubleAndlong a) { return q > a.q; }
};

template<class T>
inline void HeapSort(T*   c,
                     long n)
{
   long l, j, r, i;
   T crit;
   c--;
   if (n <= 1)
      return;
   l = n/2 + 1;
   r = n;
   while (1) {
      if (l <= 1) {
         crit = c[r];
         c[r--] = c[1];
         if (r == 1) {
            c[1] = crit;
            return;
         }
      }
      else
         crit = c[--l]; 
      j = l;
      while (1) {
         i = j;
         j = 2*j;
         if (j>r) {
            c[i] = crit;
            break;
         }
         if ((j<r) && (c[j] < c[j+1]))
            j++;
         if (crit < c[j])
            c[i] = c[j];
         else {
            c[i] = crit;
            break;
         }
      }
   }
}


Triangle* swapTest(Triangle* t1,
                   short     a);

double QuadQuality(const Vertex& a,
                   const Vertex& b,
                   const Vertex& c,
                   const Vertex& d)
{
   R2 A((R2)a), B((R2)b), C((R2)c), D((R2)d);
   R2 AB(B-A), BC(C-B), CD(D-C), DA(A-D);
   const Metric &Ma = a, &Mb = b, &Mc = c, &Md = d;

   double lAB = Norme2(AB), lBC = Norme2(BC), lCD = Norme2(CD), lDA = Norme2(DA);
   AB /= lAB;  BC /= lBC;  CD /= lCD;  DA /= lDA;

// version aniso 
   double cosDAB = Ma(DA,AB)/(Ma(DA)*Ma(AB)), sinDAB = Det(DA,AB);
   double cosABC = Mb(AB,BC)/(Mb(AB)*Mb(BC)), sinABC = Det(AB,BC);
   double cosBCD = Mc(BC,CD)/(Mc(BC)*Mc(CD)), sinBCD = Det(BC,CD);
   double cosCDA = Md(CD,DA)/(Md(CD)*Md(DA)), sinCDA = Det(CD,DA);
   double sinmin =Min(Min(sinDAB,sinABC),Min(sinBCD,sinCDA));
   if (sinmin<=0)
      return sinmin;
   return 1.0-Max(Max(Abs(cosDAB),Abs(cosABC)),Max(Abs(cosBCD),Abs(cosCDA)));
}


GeometricalEdge* Triangles::ProjectOnCurve(Edge&         BhAB,
                                           Vertex&       vA,
                                           Vertex&       vB,
					   double        theta,
					   Vertex&       R,
                                           VertexOnEdge& BR,
                                           VertexOnGeom& GR) 
{
   void *pA=0, *pB=0;
   double tA=0, tB=0;
   R2 A=vA, B=vB;
   Vertex *pvA=&vA, *pvB=&vB;
   if (vA.vint == IsVertexOnVertex)
      pA = vA.onbv;
   else if (vA.vint == IsVertexOnEdge) {
      pA = vA.onbe->be;
      tA = vA.onbe->abcisse;
   }
   else {
      cerr << "ProjectOnCurve On Vertex " << BTh.Number(vA) << " " << endl;
      cerr << " forget call to SetVertexFieldOnBTh" << endl;
      MeshError(-1);
   }
   if (vB.vint == IsVertexOnVertex)
      pB = vB.onbv;
   else if (vB.vint == IsVertexOnEdge) {
      pB = vB.onbe->be;
      tB = vB.onbe->abcisse;
   }
   else {
      cerr << "ProjectOnCurve On Vertex " << BTh.Number(vB) << " " << endl;
      cerr << " forget call to SetVertexFieldOnBTh" << endl;
      MeshError(-1);
   } 
   Edge *e = &BhAB;
   assert (pA && pB && e);
// Be careful the background edge e is on same geom edge 
// of the initial edge def by the 2 vertex A B;
   assert(e>=BTh.edges && e<BTh.edges+BTh.nbe);// Is a background Mesh;   
// walk on BTh edge 
// 1 first find a back ground edge contening the vertex A
// 2 walk n back gound boundary to find the final vertex B
   if (vA.vint == IsVertexOnEdge) // find the starting edge 
      e = vA.onbe->be;
   else if (vB.vint == IsVertexOnEdge) {
      theta = 1-theta;
      Exchange(tA,tB);
      Exchange(pA,pB);
      Exchange(pvA,pvB);
      Exchange(A,B);
      e =  vB.onbe->be;
   } 
   else
      assert(0);

// find the direction of walking with sens of edge and pA,PB;
   R2 AB=B-A;
   int kkk=0;
   double cosE01AB = (((R2)(*e)[1] - (R2)(*e)[0]),AB);
   int sens = (cosE01AB>0) ? 1 : 0;
   double abscisse = -1;
   for (int cas=0; cas<2; cas++) {
      int iii;
      Vertex *v0=pvA, *v1=NULL; 
      Edge *neee, *eee;
      double lg =0, te0;
      for (eee=e, iii=sens, te0=tA;
             eee && (((void*) eee) != pB) && (( void*) (v1=&((*eee)[iii]))) != pB ;
             neee = eee->adj[iii],iii = 1-neee->Intersection(*eee), eee=neee, v0=v1, te0=1-iii) { 
         assert(kkk<100);
         kkk++;
         assert(eee);
         double lg0 = lg;
         double dp = LengthInterpole(v0->m,v1->m,(R2) *v1 - (R2) *v0);
         lg += dp;
         if (cas && abscisse <= lg) {
            double sss = (abscisse-lg0)/dp;
            double thetab = te0*(1-sss)+ sss*iii;
            assert(thetab>=0 && thetab<=1);
            BR = VertexOnEdge(&R,eee,thetab);
            return Gh.ProjectOnCurve(*eee,thetab,R,GR);
         }
      }
      if (v1 != pvB) {
         if ((void*)v1 == pB)
            tB = iii;
         double lg0 = lg;
         assert(eee);
         v1 = pvB;
         double dp = LengthInterpole(v0->m,v1->m,(R2) *v1 - (R2) *v0);
         lg += dp;
         abscisse = lg*theta;
         if (abscisse<=lg && abscisse>=lg0) {
            double sss = (abscisse-lg0)/dp;
            double thetab = te0*(1-sss)+ sss*tB;
            assert(thetab>=0 && thetab<=1);
            BR = VertexOnEdge(&R,eee,thetab);
            return Gh.ProjectOnCurve(*eee,thetab,R,GR);
         }
      }
      abscisse = lg*theta;    
   }
   cerr << " Big Bug" << endl;
   MeshError(678);
   return 0;   
}


void Triangles::MakeQuadrangles(double costheta)
{
   if (verbosity>2) 
      cout << " -- MakeQuadrangles costheta = " << costheta << endl;
   if (verbosity>5)  
      cout << "    (in)  Nb of Quadrilaterals = " << NbOfQuad 
           << " Nb Of Triangles = " << nbt-NbOutT- NbOfQuad*2 
           << " Nb of outside triangles = " << NbOutT << endl;
   if (costheta >1) {
      if (verbosity>5)
         cout << "     do nothing costheta >1 "<< endl;
      return;
   }
   long nbqq = (nbt*3)/2;
   DoubleAndlong  *qq = new DoubleAndlong[nbqq];
   long i, ij, k=0;
   int j;
   for (i=0; i<nbt; i++)
      for (j=0; j<3; j++)
         if ((qq[k].q=triangles[i].QualityQuad(j))>=costheta)
            qq[k++].i3j=i*3+j;
   HeapSort(qq,k);
   long kk=0;
   for (ij=0; ij<k; ij++) {
      i = qq[ij].i3j/3;
      j = int(qq[ij].i3j%3);
      if (triangles[i].QualityQuad(j,0) >=costheta) 
         triangles[i].SetHidden(j),kk++;
   }
   NbOfQuad = kk;
   if (verbosity>2) {
      cout << "    (out)  Nb of Quadrilaterals = " << NbOfQuad 
           << " Nb Of Triangles = " << nbt-NbOutT- NbOfQuad*2 
           << " Nb of outside triangles = " << NbOutT << endl;
   }
   delete [] qq;
}


int Triangles::SplitElement(int choice)
{
   Direction NoDirOfSearch;
   Triangles *bth=&BTh;
   const int withBackground = bth != this && bth;
   if (verbosity>2) 
      cout << " -- SplitElement " << (choice? " Q->4Q and T->4T " : " Q->4Q or T->3Q " ) << endl;
   if (verbosity>5)
      cout << endl << "    (in)  Nb of Quadrilaterals = " << NbOfQuad 
           << " Nb Of Triangles = " << nbt-NbOutT- NbOfQuad*2 
           << " Nb of outside triangles = " << NbOutT << endl;

   ReNumberingTheTriangleBySubDomain();
   const long nfortria( choice ? 4 : 6);
   if (withBackground) {
      BTh.SetVertexFieldOn();
      SetVertexFieldOnBTh();
    }
    else
       BTh.SetVertexFieldOn();

   long newnbt=0, newnbv=0, kkk=0, *kedge=0;
   long newNbOfQuad=NbOfQuad, *ksplit=0, *ksplitarray=0;
   int ret=0;
   if (nbvx<nbv+nbe)
      return 1;
   Triangles *OCurrentTh = CurrentTh;
   CurrentTh = this;

// 1) create  the new points by spliting the internal edges 
   long nbvold=nbv, nbtold=nbt, NbOutTold=NbOutT, NbEdgeOnGeom=0, i;
   nbt -= NbOutT; // remove all the  the ouside triangles 
   long nbtsave = nbt;
   Triangle *lastT = triangles + nbt;
   for (i=0; i<nbe; i++)
      if (edges[i].on)
         NbEdgeOnGeom++;
   long newnbe = nbe + nbe;
   long newNbVerticesOnGeomEdge = NbVerticesOnGeomEdge + NbEdgeOnGeom;
   long newNbVertexOnBThEdge = withBackground ? NbVertexOnBThEdge + NbEdgeOnGeom : 0;

// do allocation for pointer to the geometry and background
   VertexOnGeom *newVerticesOnGeomEdge = new VertexOnGeom[newNbVerticesOnGeomEdge];
   VertexOnEdge *newVertexOnBThEdge = newNbVertexOnBThEdge ? new VertexOnEdge[newNbVertexOnBThEdge] : 0;
   if (NbVerticesOnGeomEdge)
      memcpy(newVerticesOnGeomEdge,VerticesOnGeomEdge,sizeof(VertexOnGeom)*NbVerticesOnGeomEdge);
   if (NbVertexOnBThEdge)
      memcpy(newVertexOnBThEdge,VertexOnBThEdge,sizeof(VertexOnEdge)*NbVertexOnBThEdge);
   Edge *newedges = new Edge [newnbe];
   SetOfEdges4 *edge4 = new SetOfEdges4(nbe,nbv);

   long k=nbv, kk=0, kvb=NbVertexOnBThEdge, kvg=NbVerticesOnGeomEdge, ie=0;
   Edge **edgesGtoB=0;
   if (withBackground)
      edgesGtoB = BTh.MakeGeometricalEdgeToEdge();
   long ferr=0;
   for (i=0; i<nbe; i++)
      newedges[ie].on = 0;

   for (i=0; i<nbe; i++) {
      GeometricalEdge *ong =  edges[i].on;
      newedges[ie] = edges[i];
      newedges[ie].adj[0] = newedges + (edges[i].adj[0]-edges) ;
      newedges[ie].adj[1] = newedges + ie +1;
      R2 A = edges[i][0],B = edges[i][1];
      kk += (i == edge4->addtrie(Number(edges[i][0]),Number(edges[i][1])));
      if (ong) { 
         if (withBackground) {
            assert(edgesGtoB); 
            ong = ProjectOnCurve(*edgesGtoB[Gh.Number(edges[i].on)],
                                 edges[i][0],edges[i][1],0.5,vertices[k],
                                 newVertexOnBThEdge[kvb],
                                 newVerticesOnGeomEdge[kvg++]);
            vertices[k].ReferenceNumber = edges[i].ref;
            vertices[k].DirOfSearch = NoDirOfSearch;        

//          get the Info on background mesh 
            double s = newVertexOnBThEdge[kvb];
            Vertex &bv0 = newVertexOnBThEdge[kvb][0];
            Vertex &bv1 = newVertexOnBThEdge[kvb][1];
//          compute the metrix of the new points 
            vertices[k].m =  Metric(1-s,bv0,s,bv1); 
            kvb++;
         }
         else {
            ong = Gh.ProjectOnCurve(edges[i],0.5,vertices[k],newVerticesOnGeomEdge[kvg++]);
//          vertices[k].i = toI2( vertices[k].r);
            vertices[k].ReferenceNumber = edges[i].ref;
            vertices[k].DirOfSearch = NoDirOfSearch;
            vertices[k].m = Metric(0.5,edges[i][0],0.5,edges[i][1]);	      
         }
      }
      else {
         vertices[k].r = ((R2) edges[i][0] + (R2)  edges[i][1] )*0.5;
         vertices[k].m = Metric(0.5,edges[i][0],0.5,edges[i][1]);
         vertices[k].on = 0;
      }
      R2 AB = vertices[k].r;
      R2 AA = (A+AB)*0.5;
      R2 BB = (AB+B)*0.5;
      vertices[k].ReferenceNumber = edges[i].ref;
      vertices[k].DirOfSearch = NoDirOfSearch;	    
      newedges[ie].on = Gh.Containing(AA,ong);
      newedges[ie++].v[1] = vertices + k;
      newedges[ie] = edges[i];
      newedges[ie].adj[0] = newedges + ie -1;
      newedges[ie].adj[1] = newedges + (edges[i].adj[1]-edges) ;
      newedges[ie].on =  Gh.Containing(BB,ong);
      newedges[ie++].v[0] = vertices + k;
      k++;
   }
   if (edgesGtoB)
      delete [] edgesGtoB;
   edgesGtoB = 0;
   newnbv = k;
   newNbVerticesOnGeomEdge=kvg;
   if (newnbv>nbvx)
      goto Error;
   nbv = k;

   kedge = new long[3*nbt+1];
   ksplitarray = new long[nbt+1];
   ksplit = ksplitarray + 1; // because ksplit[-1] == ksplitarray[0]

   for (i=0; i<3*nbt; i++)
      kedge[i] = -1;
   for (i=0; i<nbt; i++) {
      Triangle &t=triangles[i];
      assert(t.link);
      for (int j=0; j<3; j++) {
         const TriangleAdjacent ta=t.Adj(j);
         const Triangle &tt=ta;
         if (&tt>=lastT)
            t.SetAdj2(j,0,0);// unset adj
         const Vertex &v0=t[VerticesOfTriangularEdge[j][0]];
         const Vertex &v1=t[VerticesOfTriangularEdge[j][1]];
         long ke = edge4->findtrie(Number(v0),Number(v1));
         if (ke>0) {
            long ii=Number(tt);
            int jj=ta;
            long ks = ke + nbvold;
            kedge[3*i+j] = ks;
            if (ii<nbt) // good triangle
               kedge[3*ii+jj] = ks;
            Vertex &A=vertices[ks];
            double aa=0, bb=0, cc=0, dd;
            if ((dd=Area2(v0.r,v1.r,A.r)) >= 0) { // warning PB roundoff error
               if (t.link && ((aa=Area2(A.r,t[1].r,t[2].r)) < 0.0 
                          ||  (bb=Area2(t[0].r,A.r,t[2].r)) < 0.0  
                          ||  (cc=Area2(t[0].r,t[1].r,A.r)) < 0.0))
                  ferr++, cerr << " Error : " <<  ke + nbvold << " not in triangle " 
                               << i << " In=" << !!t.link
                               << " " <<  aa  << " " << bb << " " << cc << " " << dd << endl;
            }
            else {
              if (tt.link && ( (aa=Area2(A.r    , tt[1].r, tt[2].r)) < 0 
                           ||   (bb=Area2(tt[0].r, A.r    , tt[2].r)) < 0 
                           ||   (cc=Area2(tt[0].r,tt[1].r, A.r    )) < 0)) 
                  ferr++, cerr << " Warning : " << ke+nbvold << " not in triangle " << ii 
                               << " In=" << !!tt.link 
                               << " " << aa << " " << bb << " " << cc << " " << dd << endl;	
            }
         }
      }
   }
   if (ferr) {
      cerr << " Number of triangles with P2 interpolation Problem " << ferr << endl;;
      MeshError(9);
   }

   for (i=0; i<nbt; i++) {
      ksplit[i] = 1; // no split by default
      const Triangle &t=triangles[ i];
      int nbsplitedge=0, nbinvisible=0, invisibleedge=0;
      int kkk[3];
      for (int j=0; j<3; j++) {
         if (t.Hidden(j))
            invisibleedge = j, nbinvisible++;
         const TriangleAdjacent ta=t.Adj(j);
         const Triangle &tt=ta;
         const Triangle *ttt=&tt;
         const Vertex &v0=t[VerticesOfTriangularEdge[j][0]];
         const Vertex &v1=t[VerticesOfTriangularEdge[j][1]];
         if (kedge[3*i+j]<0) {
            long ke = edge4->findtrie(Number(v0),Number(v1));
            if (ke<0) {
               if (ttt) {
                  long ii=Number(tt);
                  int jj=ta;
                  kedge[3*i+j] = k; 
                  kedge[3*ii+jj] = k;
                  if (k<nbvx) {
                     vertices[k].r = ((R2)v0+(R2)v1)/2;
                     vertices[k].ReferenceNumber = 0;
                     vertices[k].DirOfSearch = NoDirOfSearch;
                     vertices[k].m =  Metric(0.5,v0,0.5,v1);
                  }
                  k++;
                  kkk[nbsplitedge++] = j;	      
               }
               else
                  cerr << endl << " Bug " << i << " " << j << " t = " << t << endl;
            }
	    else {
               kedge[3*i+j] = nbvold + ke;
               kkk[nbsplitedge++] = j;// previously splited
            }
         }
         else 
            kkk[nbsplitedge++] = j;// previously splited 
      }
      assert(nbinvisible<2);
      switch (nbsplitedge) {

         case 0: ksplit[i] = 10;
                 newnbt++;
                 break;

         case 1: ksplit[i] = 20 + kkk[0];
                 newnbt += 2;
                 break;
 
         case 2: ksplit[i] = 30+3-kkk[0]-kkk[1];
                 newnbt += 3;
                 break;

         case 3: if (nbinvisible)
                    ksplit[i] = 40 + invisibleedge, newnbt += 4;
                 else
                    ksplit[i] = 10*nfortria, newnbt+=nfortria;
                  break;
      } 
      assert(ksplit[i]>=40);
   }

// now do the element split
   newNbOfQuad = 4*NbOfQuad;
   nbv = k;
   kkk = nbt;
   ksplit[-1] = nbt;

// look on old true  triangles
   for (i=0; i<nbtsave; i++) {
      long kk=ksplit[i]/10;
      int ke=int(ksplit[i]%10);
      assert(kk<7 && kk >0);

//    def the numbering   k (edge) i vertex 
      int k0=ke, k1=NextEdge[k0], k2=PreviousEdge[k0];
      int i0=OppositeVertex[k0], i1=OppositeVertex[k1], i2=OppositeVertex[k2];
      Triangle &t0 = triangles[i];
      Vertex *v0=t0(i0), *v1=t0(i1), *v2=t0(i2);

//    save the flag Hidden
      int hid[] = {t0.Hidden(0),t0.Hidden(1),t0.Hidden(2)};
//    unset all adj -- save Hidden flag --
      t0.SetAdj2(0,0,hid[0]);
      t0.SetAdj2(1,0,hid[1]);
      t0.SetAdj2(2,0,hid[2]);

      switch  (kk) {

         case 1: break;

         case 2:
            {
               Triangle &t1=triangles[kkk++];
               t1 = t0;
               assert(kedge[3*i+i0]>=0);
               Vertex *v3 = vertices + kedge[3*i+k0];
               t0(i2) = v3;
               t1(i1) = v3;
               t0.SetAllFlag(k2,0);
               t1.SetAllFlag(k1,0);
            }
            break;
 
         case 3:
            {
               Triangle &t1=triangles[kkk++], &t2=triangles[kkk++];
               t2 = t1 = t0;
               assert(kedge[3*i+k1]>=0);
               assert(kedge[3*i+k2]>=0);
               Vertex *v01 = vertices + kedge[3*i+k2];
               Vertex *v02 = vertices + kedge[3*i+k1]; 
               t0(i1) = v01; 
               t0(i2) = v02; 
               t1(i2) = v02;
               t1(i0) = v01; 
               t2(i0) = v02; 
               t0.SetAllFlag(k0,0);
               t1.SetAllFlag(k1,0);
               t1.SetAllFlag(k0,0);
               t2.SetAllFlag(k2,0);
            }
            break;

         case 4:

         case 6: 
            {
               Triangle &t1=triangles[kkk++], &t2=triangles[kkk++], &t3=triangles[kkk++];
               t3 = t2 = t1 = t0;
               assert(kedge[3*i+k0] >=0 && kedge[3*i+k1] >=0 && kedge[3*i+k2] >=0);
               Vertex *v12 = vertices + kedge[3*i+k0];
               Vertex *v02 = vertices + kedge[3*i+k1]; 
               Vertex *v01 = vertices + kedge[3*i+k2];
               t0(i1) = v01;
               t0(i2) = v02;
               t0.SetAllFlag(k0,hid[k0]);
               t1(i0) = v01;
               t1(i2) = v12;
               t0.SetAllFlag(k1,hid[k1]);
               t2(i0) = v02;
               t2(i1) = v12;
               t2.SetAllFlag(k2,hid[k2]);
               t3(i0) = v12;
               t3(i1) = v02;
               t3(i2) = v01;
               t3.SetAllFlag(0,hid[0]);	   
               t3.SetAllFlag(1,hid[1]);	   
               t3.SetAllFlag(2,hid[2]);
               if (kk == 6) {
                  Triangle &t4 = triangles[kkk++];
                  Triangle &t5 = triangles[kkk++];
                  t4 = t3;
                  t5 = t3;
                  t0.SetHidden(k0);
                  t1.SetHidden(k1);
                  t2.SetHidden(k2);
                  t3.SetHidden(0);
                  t4.SetHidden(1);
                  t5.SetHidden(2);
                  if (nbv < nbvx) {
                     vertices[nbv].r = ((R2)*v01 + (R2)*v12  + (R2)*v02 ) / 3.0;
                     vertices[nbv].ReferenceNumber = 0;
                     vertices[nbv].DirOfSearch = NoDirOfSearch;
                     double a3[] = {1./3.,1./3.,1./3.};
                     vertices[nbv].m = Metric(a3,v0->m,v1->m,v2->m);
                     Vertex *vc = vertices + nbv++;
                     t3(i0) = vc;
                     t4(i1) = vc;
                     t5(i2) = vc;
                  }
               else
                  goto Error;
            }
         }
         break;
      }
      long jj;
      if (t0.link) 
         for (jj=nbt; jj<kkk; jj++) {
            triangles[jj].link = t0.link;
            t0.link = triangles + jj;
         }
      if (kk==6)
         newNbOfQuad += 3;
      nbt = kkk;
      ksplit[i] = nbt; // save last adress of the new triangles
      kkk = nbt;
   }
   for (i=0; i<nbv; i++)
      vertices[i].m = vertices[i].m*2.;
   if (withBackground)
      for (i=0;i<BTh.nbv;i++)
         BTh.vertices[i].m =  BTh.vertices[i].m*2.;
   ret = 2;
   if (nbt>= nbtx)
      goto Error; 
   if (nbv>= nbvx)
      goto Error;

// Generation of the new triangles 
   SetIntCoor("In SplitElement");
   ReMakeTriangleContainingTheVertex();
   if (withBackground)
      BTh.ReMakeTriangleContainingTheVertex();
   delete [] edges;
   edges = newedges;
   nbe = newnbe;
   NbOfQuad = newNbOfQuad;

   for (i=0; i<NbSubDomains; i++) {
      long k = subdomains[i].edge - edges;
      subdomains[i].edge =  edges + 2*k; // split all edge in 2 
   }
    
   if (ksplitarray)
      delete [] ksplitarray;
   if (kedge)
      delete [] kedge;
   if (edge4)
      delete edge4;
   if (VerticesOnGeomEdge)
      delete [] VerticesOnGeomEdge;
   VerticesOnGeomEdge = newVerticesOnGeomEdge;
   if (VertexOnBThEdge)
      delete [] VertexOnBThEdge;
   VertexOnBThEdge = newVertexOnBThEdge;
   NbVerticesOnGeomEdge = newNbVerticesOnGeomEdge;
   NbVertexOnBThEdge=newNbVertexOnBThEdge;
   FillHoleInMesh();

   if (verbosity>2)
      cout << "    (out) Nb of Quadrilaterals = " << NbOfQuad 
           << " Nb Of Triangles = " << nbt-NbOutT- NbOfQuad*2 
           << " Nb of outside triangles = " << NbOutT << endl;

   CurrentTh = OCurrentTh;
   return 0;

 Error:
   nbv = nbvold;
   nbt = nbtold;
   NbOutT = NbOutTold;

   delete [] newedges;
   if (ksplitarray)
      delete [] ksplitarray;
   if (kedge)
      delete [] kedge;
   if (newVerticesOnGeomEdge)
      delete [] newVerticesOnGeomEdge;
   if (edge4)
      delete edge4;
   if (newVertexOnBThEdge)
      delete [] newVertexOnBThEdge;
   CurrentTh = OCurrentTh;
   return ret;
}

} //  end of namespcae bamg 
