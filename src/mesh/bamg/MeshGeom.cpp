// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY: Bamg: Bidimensional Anisotrope Mesh Generator
// RELEASE: 0 
// USAGE  : You may copy freely these files and use it for    
//          teaching or research. These or part of these may   
//          not be sold or used for a commercial purpose with- 
//          out our consent : fax (33) 1 39 63 55 14       
//
// AUTHOR:   F. Hecht,    
// ORG    :  INRIA
// E-MAIL :   Frederic.Hecht@Inria.fr   
//
// ORIG-DATE:     Dec 97

// #define TRACETRIANGLE 3
extern long verbosity;

#include <cstdio>
#include <string.h>
#include <cmath>
#include <time.h>
#include <iostream>
using namespace std;

#include "mesh/bamg/Mesh2.h"
#include "mesh/bamg/QuadTree.h"
#include "mesh/bamg/SetOfE4.h"

namespace bamg {

void Triangles::ConsGeometry(double cutoffradian,
                             int*   equiedges) // construct a geometry if no geo 
{
  //  if equiedges existe taille nbe 
  //   equiedges[i]/2 == i  original
  //   equiedges[i]/2 = j  =>   equivalence entre i et j => meme maillage
  //   equiedges[i]%2   : 0 meme sens , 1 pas meme sens 
  // --------------------------
   if (verbosity>1) 
      cout << " -- construction of the geometry from the 2d mesh " << endl;
   if (nbt<=0 || nbv <=0 ) { MeshError(101);}

// construction of the edges 
   CurrentTh = this;
//  long NbTold = nbt;
// generation of the integer coor
// generation of the adjacence of the triangles
   if (cutoffradian>=0)
      Gh.MaximalAngleOfCorner = cutoffradian;
   SetOfEdges4 * edge4= new SetOfEdges4(nbt*3,nbv);
   long * st = new long[nbt*3];
   long i,k;
   int j; 
   if (Gh.name)
      delete Gh.name;
   Gh.name = new char [ name ? strlen(name) + 15 : 50 ];
   Gh.name[0] = 0;
   strcat(Gh.name,"cons from: ");
   if (name)
      strcat(Gh.name,name);
   else
      strcat(Gh.name," a mesh with no name");
   for (i=0; i<nbt*3; i++)
      st[i] = -1;
   long kk=0, nbeold=nbe;
   for (i=0; i<nbe; i++)
      edge4->addtrie(Number(edges[i][0]),Number(edges[i][1]));
   if (nbe != edge4->nb()) { 
      cerr << " Some Double edge in the mesh, the number is " << nbe 
	   << " nbe4=" << edge4->nb()  << endl;
      MeshError(1002);
   }
   for (i=0; i<nbt; i++) {
      for (j=0; j<3; j++) {
         long k = edge4->addtrie(Number(triangles[i][VerticesOfTriangularEdge[j][0]]),
                                 Number(triangles[i][VerticesOfTriangularEdge[j][1]]));
         long invisible = triangles[i].Hidden(j);
         if (st[k]==-1)
            st[k] = 3*i + j;
         else if (st[k]>=0) {
            assert(!triangles[i].TriangleAdj(j) && !triangles[st[k] / 3].TriangleAdj(int(st[k]%3)));
            triangles[i].SetAdj2(j,triangles + st[k]/3, int(st[k]%3));
            if (invisible)
               triangles[i].SetHidden(j);
            if (k<nbe)
               triangles[i].SetLocked(j);
            st[k] = -2 - st[k];
         }
         else {
            cerr << " The edge (" 
                 << Number(triangles[i][VerticesOfTriangularEdge[j][0]]) << " , "
                 << Number(triangles[i][VerticesOfTriangularEdge[j][1]])
                 << " ) is in more than 2 triangles " << k << endl;
            cerr << " Edge " << j << " of Triangle " << i << endl;
            cerr << " Edge " << (-st[k]+2)%3 << " Of Triangle " << (-st[k]+2)/3 << endl;
            cerr << " Edge " << triangles[(-st[k]+2)/3].NuEdgeTriangleAdj((int)((-st[k]+2)%3))
                 << " of Triangle " <<  Number(triangles[(-st[k]+2)/3].TriangleAdj((int)((-st[k]+2)%3))) 
                 << endl;
            MeshError(9999);
         }
      }
   }
   long nbedges = edge4->nb();
   delete edge4;
   edge4 = NULL;

   if (verbosity>5) {
      if (name)
         cout << "    On Mesh " << name << endl;
      cout << "    - The number of Vertices  = " << nbv << endl;
      cout << "    - The number of Triangles = " << nbt << endl;
      cout << "    - The number of given edge = " << nbe << endl;
      cout << "    - The number of all edges = " << nbedges << endl;
      cout << "    - The Euler number = 1-Nb Of Hole = " << nbt-nbedges+nbv << endl;
   }

// check consistency of edge[].adj and the geometrical required  vertex
   k = kk = 0;
   long it;
   for (i=0; i<nbedges; i++) {
      if (st[i]<-1) { 
         it = (-2-st[i])/3;
         j = (int) ((-2-st[i])%3);
         Triangle &tt = *triangles[it].TriangleAdj(j);
         if (triangles[it].color != tt.color|| i<nbeold)
            k++;
      }
      else if (st[i]>=0)
         kk++;
   }

   if (verbosity>4 && (k+kk))
      cout << "    Nb of ref edge " << kk+k << " (internal " << k << ")"
           << " in file " << nbe << endl;
   k += kk;
   kk = 0;
   if (k) {
      nbe = k;
      Edge *edgessave = edges;
      edges = new Edge[nbe];
      k =0;
      if (verbosity>4)
         cout << "    Construction of the edges  " << nbe << endl;
      for (i=0; i<nbedges; i++) {
         long add = -1;
         if (st[i] < -1) { 
            it = (-2-st[i])/3;
            j = int((-2-st[i])%3);
            Triangle &tt = * triangles[it].TriangleAdj(j);
            if (triangles[it].color !=  tt.color || i<nbeold)
               add = k++;
         }
         else if (st[i]>=0) {
            it = st[i]/3;
            j  = int(st[i]%3);
            add = k++;
         }
         if (add>=0 && add<nbe) {
            edges[add].v[0] = &triangles[it][VerticesOfTriangularEdge[j][0]];
            edges[add].v[1] = &triangles[it][VerticesOfTriangularEdge[j][1]];
            if (i<nbeold) // in file edge // Modif FH 06122055 
               edges[add].ref = edgessave[i].ref; 
            else
               edges[add].ref = Min(edges[add].v[0]->ref(),edges[add].v[1]->ref());
         }
      }
      assert(k==nbe);
      if (edgessave)
         delete [] edgessave;
   }

// construction of edges[].adj 
   for (i=0; i<nbv; i++) 
      vertices[i].color = 0;
   for (i=0; i<nbe; i++)
      for (j=0; j<2; j++) 
	 edges[i].v[j]->color++;

   for (i=0; i<nbv; i++) 
      vertices[i].color = (vertices[i].color==2) ? -1 : -2;
   for (i=0; i<nbe; i++) {
      for (j=0; j<2; j++) {
         Vertex *v=edges[i].v[j];
         long i0=v->color, j0;
         if (i0==-1)
            v->color = i*2 + j;
         else if (i0>=0) {
            j0 = i0%2; i0 /= 2;
            assert(v==edges[i0].v[j0]);
            edges[i].adj[j] = edges + i0;
            edges[i0].adj[j0] = edges + i;
            assert(edges[i0].v[j0] == v);
            v->color = -3;
         }
      }
   }

// now reconstruct the sub domain info 
   assert(!NbSubDomains);
   NbSubDomains = 0;
   { 
      long it;
//    find all subdomains
      long *colorT = new long[nbt];
      Triangle *tt;
      long k;
      for (it=0; it<nbt; it++)
         colorT[it] = -1;
      for (it=0; it<nbt; it++) {
         if (colorT[it]<0) {
            colorT[it] = NbSubDomains;
            long level=1, j, jt, kolor=triangles[it].color;
            st[0] = it; st[1] = 0;
            k = 1;
            while (level>0)
            if ((j=st[level]++)<3) {
              tt = triangles[st[level-1]].TriangleAdj(int(j));
               if (tt && (colorT[jt = Number(tt)] == -1) && (tt->color==kolor)) {
                  colorT[jt] = NbSubDomains;
                  st[++level] = jt;
                  st[++level] = 0;
                  k++;
               }
            }
            else
               level -= 2;
            if (verbosity>5) 
               cout << "   Nb of triangles " << k << " of Subdomain "  
                    << NbSubDomains << " " << kolor << endl;
            NbSubDomains++;
         }
      }
      if (verbosity>3)
         cout << "   Number of sub domains = " << NbSubDomains << endl;

      long isd;
      subdomains = new SubDomain[NbSubDomains];
      for (isd=0; isd<NbSubDomains; isd++)
         subdomains[isd].head = 0;
      k = 0;
      for (it=0; it<nbt; it++) {
         for (int j=0; j<3; j++) {
            tt = triangles[it].TriangleAdj(j);
            if ((!tt || tt->color != triangles[it].color) && !subdomains[isd=colorT[it]].head) {
               subdomains[isd].head = triangles + it;
               subdomains[isd].ref = triangles[it].color;
               subdomains[isd].sens = j; // hack
               subdomains[isd].edge = 0;
               k++;
            }
         }
      }
      assert(k==NbSubDomains);
      delete [] colorT;
   }
   delete [] st;

// now make the geometry
// 1 compress the vertices 
   long *colorV = new long[nbv];
   for (i=0; i<nbv; i++) 
      colorV[i] = -1;
   for (i=0; i<nbe; i++)
      for (j=0; j<2; j++)
        colorV[Number(edges[i][j])] = 0;
   k = 0;
   for (i=0; i<nbv; i++) 
      if (!colorV[i])
         colorV[i] = k++;
   Gh.nbv = k;
   Gh.nbe = nbe;
   Gh.vertices = new GeometricalVertex[k];
   Gh.edges = new GeometricalEdge[nbe];
   Gh.NbSubDomains = NbSubDomains;
   Gh.subdomains = new GeometricalSubDomain[NbSubDomains];
   if (verbosity>3)
      cout << "    Nb of vertices = " << Gh.nbv << " Nb of edges = " << Gh.nbe << endl;
   NbVerticesOnGeomVertex = Gh.nbv;
   VerticesOnGeomVertex = new VertexOnGeom[NbVerticesOnGeomVertex];
   NbVerticesOnGeomEdge = 0;
   VerticesOnGeomEdge = 0;
   {
      long j;
      for (i=0; i<nbv; i++) {
         if ((j=colorV[i])>=0) {
            Vertex &v = Gh.vertices[j];
            v = vertices[i];
            v.color = 0;
            VerticesOnGeomVertex[j] = VertexOnGeom(vertices[i], Gh.vertices[j]);
         }
      }
   }
   edge4 = new SetOfEdges4(nbe,nbv);
   double *len = new double [Gh.nbv];
   for (i=0; i<Gh.nbv; i++)
      len[i] = 0;
   Gh.pmin =  Gh.vertices[0].r;
   Gh.pmax =  Gh.vertices[0].r;

   for (i=0; i<Gh.nbv; i++) {
     Gh.pmin.x = Min(Gh.pmin.x,Gh.vertices[i].r.x);
     Gh.pmin.y = Min(Gh.pmin.y,Gh.vertices[i].r.y);
     Gh.pmax.x = Max(Gh.pmax.x,Gh.vertices[i].r.x);
     Gh.pmax.y = Max(Gh.pmax.y,Gh.vertices[i].r.y);
   }
   R2 DD05 = (Gh.pmax-Gh.pmin)*0.05;
   Gh.pmin -= DD05;
   Gh.pmax += DD05;
   Gh.coefIcoor = (MaxICoor)/(Max(Gh.pmax.x-Gh.pmin.x,Gh.pmax.y-Gh.pmin.y));
   assert(Gh.coefIcoor >0);
   double hmin=HUGE_VAL;
  
   for (i=0; i<nbe; i++) {
      long i0 = Number(edges[i][0]), i1 = Number(edges[i][1]);
      long j0 = colorV[i0], j1 = colorV[i1];
      Gh.edges[i].v[0] = Gh.vertices + j0;
      Gh.edges[i].v[1] = Gh.vertices + j1;
      Gh.edges[i].flag = 0;
      Gh.edges[i].tg[0] = R2();
      Gh.edges[i].tg[1] = R2();
      edges[i].on =  Gh.edges + i;
      if (equiedges && i<nbeold) {
         int j = equiedges[i]/2;
         int sens = equiedges[i]%2;
         if (i!=j) {
            if (verbosity>9)  
               cout << " Edges Equi " << i << " <=> " << j << " sens = " << sens << endl;
            if (sens==0)
               Gh.edges[i].SetEqui();
            else 
               Gh.edges[i].SetReverseEqui();
            Gh.edges[i].link = & Gh.edges[j];
         }
      }

      R2 x12 = Gh.vertices[j0].r-Gh.vertices[j1].r;
      double l12 = Norme2(x12);        
      hmin = Min(hmin,l12);
      Gh.vertices[j1].color++;
      Gh.vertices[j0].color++;
      len[j0] += l12;
      len[j1] += l12;
      hmin = Min(hmin,l12);
      Gh.edges[i].ref = edges[i].ref;      
      k = edge4->addtrie(i0,i1);
      assert(k==i);
   }
   for (i=0; i<Gh.nbv; i++) {
      if (Gh.vertices[i].color > 0) 
         Gh.vertices[i].m = Metric(len[i] /(double)Gh.vertices[i].color);
      else
         Gh.vertices[i].m = Metric(hmin);
   }
   delete [] len;
   for (i=0; i<NbSubDomains; i++) {
      long it=Number(subdomains[i].head);
      int j=subdomains[i].sens;
      long i0=Number(triangles[it][VerticesOfTriangularEdge[j][0]]);
      long i1=Number(triangles[it][VerticesOfTriangularEdge[j][1]]);
      k = edge4->findtrie(i0,i1);
      if (k>=0) {
         subdomains[i].sens = (vertices + i0 == edges[k].v[0]) ? 1 : -1;
         subdomains[i].edge = edges+k;
         Gh.subdomains[i].edge = Gh.edges + k;
         Gh.subdomains[i].sens = subdomains[i].sens;
         Gh.subdomains[i].ref = subdomains[i].ref;
      }
      else
         MeshError(103);
   }
   delete edge4;
   delete [] colorV;

   for (i=0; i<nbt; i++)
      for (j=0; j<3; j++)
         triangles[i].SetAdj2(j,0,triangles[i].GetAllflag(j));
}


void Geometry::EmptyGeometry()
{
   OnDisk = 0;
   NbRef = 0;
   name = 0;
   _quadtree = 0;
   curves = 0;
   triangles = 0;
   edges = 0;
   vertices = 0;
   NbSubDomains = 0;
   nbiv = nbv = nbvx = 0;
   nbe = nbt = nbtx = 0;
   NbOfCurves = 0;
   subdomains = 0;
   MaximalAngleOfCorner = 10*Pi/180;
}


Geometry::Geometry(const Geometry& Gh)
{
   long i;
   *this = Gh;
   NbRef = 0;
   _quadtree = 0;
   name = new char[strlen(Gh.name)+4];
   strcpy(name,"cp:");
   strcat(name,Gh.name);
   vertices = nbv ? new GeometricalVertex[nbv] : NULL;
   triangles = nbt ? new Triangle[nbt]:NULL;
   edges = nbe ? new GeometricalEdge[nbe]:NULL;
   curves = NbOfCurves ? new Curve[NbOfCurves]:NULL;
   subdomains = NbSubDomains ? new GeometricalSubDomain[NbSubDomains]:NULL;
   for (i=0; i<nbv; i++)
      vertices[i].Set(Gh.vertices[i],Gh,*this);
   for (i=0; i<nbe; i++)
      edges[i].Set(Gh.edges[i],Gh,*this);
   for (i=0; i<NbOfCurves; i++)
      curves[i].Set(Gh.curves[i],Gh,*this);
   for (i=0; i<NbSubDomains; i++)
      subdomains[i].Set(Gh.subdomains[i],Gh,*this);
   assert(!nbt);   
}


GeometricalEdge* Geometry::Containing(const R2               P, 
                                            GeometricalEdge* start) const
{
   GeometricalEdge *on=start, *pon=0;
   int k=0;
   while (pon != on) { 
      pon = on;
      assert(k<100);
      k++;
      R2 A = (*on)[0], B = (*on)[1];
      R2 AB = B-A, AP = P-A, BP = P-B;
      if ((AB,AP)<0) 
         on = on->Adj[0];
      else if ((AB,BP)>0) 
         on = on->Adj[1];
      else
         return on;
   }
   return on;
}


GeometricalEdge* Geometry::ProjectOnCurve(const Edge&         e,
                                                double        s,
                                                Vertex&       V,
                                                VertexOnGeom& GV) const 
{
   double save_s=s;
   int NbTry=0;

retry:    
   s = save_s;
   GeometricalEdge *on=e.on;
   assert(on);
   assert(e[0].on && e[1].on);
   const Vertex &v0=e[0], &v1=e[1];
   V.m = Metric(1.0-s,v0,s,v1);
#define MXE__LINE  __LINE__+1
   const int mxe=100;
   GeometricalEdge *ge[mxe+1];
   int sensge[mxe+1];
   double lge[mxe+1];
   int bge=mxe/2, tge=bge;
   ge[bge] = e.on;
   sensge[bge] = 1;
   R2 V0=v0, V1=v1, V01=V1-V0;
   VertexOnGeom vg0=*v0.on, vg1=*v1.on;
   if (NbTry)
      cout << "bug: s==== " << s << " e=" <<  V0 << " " << V1 << endl;

   GeometricalEdge *eg0=on, *eg1=on;
   R2 Ag=(R2)(*on)[0], Bg=(R2)(*on)[1], AB=Bg-Ag; 
   if (NbTry)
      cout <<" G edge= " << Ag << Bg << endl << " v edge" << V01
           << " v geom " << AB << (V01,AB) <<endl; 
   int OppositeSens = (V01,AB) < 0;
   int sens0=0, sens1=1;
   if (OppositeSens)
      s = 1 - s, Exchange(vg0,vg1), Exchange(V0,V1);
   if (NbTry)
      cout << "bug: edge = " << v0.r << " -> " << v1.r << endl 
	   << "sg 0 = " << vg0 << " on = " << Number(on) << ":" << Ag << Bg << "; " 
           << "sg 1 = " << vg1 
           << "--------------------------------------------" << endl;

   while (eg0 != (GeometricalEdge*)vg0 && (*eg0)(sens0) != (GeometricalVertex*)vg0) { 
      if (bge<=0) {
         if (NbTry) {
            cerr << " --Fatal Error: on the class triangles before call Geometry::ProjectOnCurve" << endl; 
            cerr << "   The mesh of the Geometry is too fine: ";
            cerr << "     1)  a mesh edge containing more than "<< mxe/2 << " geometrical edges." << endl;
            cerr << "     2)  code bug: be sure that we call Triangles::SetVertexFieldOn() before " << endl;
            cerr << "   To solve the problem do a coarsening of the geometrical mesh " << endl;
            cerr << " or change the constant value of mxe in " << __FILE__ << " line " << MXE__LINE << "( dangerous way )" << endl;	  
            MeshError(222);
         }
         NbTry++;
         goto retry;
      }
      GeometricalEdge* tmpge = eg0;
      if (NbTry)
         cout << "bug: --Edge @" << Number(tmpge) << " = " << Number(eg0) << ":" << Number(eg0->Adj[0])
              << "," << Number(eg0->Adj[1]) << ",";
      ge[--bge] = eg0 = eg0->Adj[sens0];
      assert(bge>=0 && bge <= mxe);
      sens0 = 1-( sensge[bge] = tmpge->SensAdj[sens0]);
      if (NbTry)
         cout << "bug: Edge " << Number(eg0) << " " << 1-sens0 << " S "
              << Number((*eg0)[1-sens0]) << ":" << Number(eg0->Adj[0]) << "," 
              << Number(eg0->Adj[1]) << "," << endl
              << Number(eg0)<< (*eg0)[sens0].r << "v = " << Number((*eg1)(sens0)) << " e = " << eg0 << endl;
   }
   if (NbTry)
      cout << Number((GeometricalEdge*) vg1) << " " << Number((GeometricalVertex*) vg1) << endl;
   while (eg1 != (GeometricalEdge*) vg1 && (*eg1)(sens1) != (GeometricalVertex*)vg1) {
      if (tge>=mxe) {
         cerr << " --Fatal Error: on the class triangles before call Geometry::ProjectOnCurve" << endl; 
         NbTry++;
         if (NbTry<2)
            goto retry;
         cerr << "   The mesh of the Geometry is to fine:";
         cerr << "     1)  a mesh edge  containing more than "<< mxe/2 << " geometrical edges." << endl;
         cerr << "     2)  code bug : be sure that we call   Triangles::SetVertexFieldOn() before " << endl;
         cerr << "   To solve the problem do a coarsening of the geometrical mesh " << endl;
         cerr << " or change the constant value of mxe in " << __FILE__ << " line " << MXE__LINE << "( dangerous way )" << endl;	
         MeshError(223);
      }

      GeometricalEdge* tmpge=eg1;
      if (NbTry)
         cout << "++Edge @" << tmpge << " = " << Number(eg1) <<"%" << Number(eg1->Adj[0]) << "," 
              << Number(eg1->Adj[1]) << ",";
      ge[++tge] = eg1 = eg1->Adj[sens1];
      sensge[tge] = sens1 = 1-tmpge->SensAdj[sens1];
      assert(tge>=0 && tge<=mxe);
      if (NbTry)
         cout << "  Edge " << Number(eg1) << " " << sens1 << " S "
              << Number((*eg1)[sens1]) << "%"<< Number(eg1->Adj[0]) << "," << Number(eg1->Adj[1]) << ","
              << Number(eg1) << (*eg1)[sens1].r << "v = " << Number((*eg1)(sens1)) << " e = " << Number(eg1) << endl;
   }
   if (NbTry)
      cout << endl;

   if ((*eg0)(sens0) == (GeometricalVertex*) vg0)
      vg0 = VertexOnGeom( *(Vertex *)vg0,*eg0,sens0);
   if ((*eg1)(sens1) == (GeometricalVertex*) vg1)
      vg1 = VertexOnGeom( *(Vertex *) vg1,*eg1,sens1);

   double sg;
   if (eg0==eg1) { 
      double s0=vg0, s1=vg1;
      sg = s0*(1.0-s) + s*s1;
      on = eg0;
   }
   else {
      R2 AA=V0, BB;
      double s0, s1;
      int i=bge;
      double ll=0;
      for (i=bge; i<tge; i++) {
         assert( i>=0 && i <= mxe);
         BB = (*ge[i])[sensge[i]];
         lge[i] = ll += Norme2(AA-BB);
         AA = BB;
      }
      lge[tge] = ll += Norme2(AA-V1); 

//    search the geometrical edge
      assert(s<=1.0);
      double ls=s*ll;
      on = 0;
      s0 = vg0;
      s1= sensge[bge];
      double l0=0, l1;
      i = bge;
      while ((l1=lge[i]) < ls) {
         assert(i>=0 && i<=mxe);
         i++, s0=1-(s1=sensge[i]), l0=l1;
      }
      on = ge[i];
      if (i==tge) 
         s1 = vg1;
      s = (ls-l0)/(l1-l0);
      sg = s0*(1.0-s) +  s*s1;    
   } 
   assert(on);
   V.r = on->F(sg);
   GV = VertexOnGeom(V,*on,sg);
   return on;
}


void Geometry::AfterRead()
{
   if (verbosity>20)
      cout << "Geometry::AfterRead()" << nbv << " " << nbe << endl;
   long i, k=0;
   int jj;
   long *hv = new long [nbv];
   long *ev = new long [2*nbe];
   double *eangle = new double [nbe];
   {
      double eps=1e-20;
      QuadTree quadtree; // to find same vertices
      Vertex *v0=vertices; 
      GeometricalVertex *v0g = (GeometricalVertex *)(void *)v0;
      int k=0;
      for (i=0; i<nbv; i++) 
         vertices[i].link = vertices + i;
      for (i=0; i<nbv; i++) {
         vertices[i].i = toI2(vertices[i].r);
         Vertex *v = quadtree.NearestVertex(vertices[i].i.x,vertices[i].i.y);
         if (v && Norme1(v->r-vertices[i])<eps) { 
            GeometricalVertex *vg = (GeometricalVertex *)(void *)v;
            int j=vg-v0g;
            assert(v== &(Vertex &)vertices[j]);
            vertices[i].link = vertices + j;
            k++;
         }
         else
            quadtree.Add(vertices[i]);
      }
      if (k) {
         cout << " Number of distinct vertices " << nbv - k << " Over " << nbv << endl;
         cout << " The duplicate vertex " << endl;
         for (i=0; i<nbv; i++)
            if (!vertices[i].IsThe())
               cout << " " << i << " and " << Number(vertices[i].The()) << endl;
         MeshError(102);
      }

//    verification of cracked edge
      int err=0;
      for (i=0; i<nbe; i++) {
         if (edges[i].Cracked()) {
//          verification of crack
           GeometricalEdge &e1 = edges[i], &e2 = *e1.link;
            cerr << i << " " << e1[0].The() << " " << e2[0].The() << " " 
                 << e1[1].The() << " " << e2[1].The() << endl;
            if (e1[0].The() == e2[0].The() && e1[1].The() == e2[1].The())
               ;
            else {
               if (e1[0].The()==e2[1].The() && e1[1].The()==e2[0].The())
                  ;
               else {
                  err++;
                  cerr << " Cracked edges with no same vertex " << &e1-edges << " " << &e2-edges << endl;
               }
            }
         }
         else
            ;
      }
      if (err) {
         cerr << " Some vertices are not distint and not on cracked edge " << err<< endl;
         MeshError(222);
      }
   }
   if (verbosity>7)
      for (i=0; i<nbv; i++)
         if (vertices[i].Required())
            cout << "     The geo vertices  " << i << " is required" << endl;
   for (i=0; i<nbv; i++) 
      hv[i] = -1;
   for (i=0; i<nbe; i++) {
      R2 v10 = edges[i].v[1]->r - edges[i].v[0]->r;
      double lv10 = Norme2(v10);
      if (lv10==0) {
         cerr << "The length of " << i << "th Egde is 0." << endl;
         MeshError(1);
      }
      eangle[i] = atan2(v10.y,v10.x);
      if (verbosity>9)
         cout << "     angle edge " << i << " " << eangle[i]*180/Pi << v10 << endl;
      for (jj=0; jj<2; jj++) {
         long v=Number(edges[i].v[jj]);
         ev[k] = hv[v];
         hv[v] = k++;
      }
   }

// bulle sort on the angle of edge
   for (i=0; i<nbv; i++) {
      int exch = 1, ord = 0;
      while (exch) {
         exch = 0;
         long *p = hv+i, *po = p;
         long n = *p;
         double angleold = -1000 ; // angle = - infini 
         ord = 0;
         while (n>=0) {
            ord++;
            long i1 = n/2, j1 = n%2, *pn = ev+n;
            double angle = j1 ? OppositeAngle(eangle[i1]) : eangle[i1];
            n = *pn;
            if (angleold>angle) // exch to have  po -> pn -> p 
               exch = 1, *pn = *po, *po = *p, *p = n, po = pn;
            else // to have: po -> p -> pn 
               angleold = angle, po = p, p = pn;
         }
      }
      if (ord>=1) {
         long n=hv[i];
         while (n>=0)
            n = ev[n];
      }
      if (ord==2) {
         long n1=hv[i], n2=ev[n1];
         long i1=n1/2, i2=n2/2;
         long j1=n1%2, j2=n2%2; // vertex in the edge
         double angle1 = j1 ? OppositeAngle(eangle[i1]) : eangle[i1];
         double angle2 = !j2 ? OppositeAngle(eangle[i2]) : eangle[i2];
         double da12 = Abs(angle2-angle1);
         if (verbosity>9)
            cout << "     check angle " << i << " " << i1 << " " << i2  << " " << 180*(da12)/Pi 
                 << " " << 180*MaximalAngleOfCorner/Pi << vertices[i] << endl;
         if (  (da12 >= MaximalAngleOfCorner) 
            && (da12 <= 2*Pi-MaximalAngleOfCorner)) {
            vertices[i].SetCorner() ; 
            if (verbosity>7)
               cout << "     The vertex " << i << " is a corner (angle) " 
                    << 180*(da12)/ Pi << " " << 180*MaximalAngleOfCorner/Pi << endl;
         }
//       if the ref a changing then is SetRequired();
         if (edges[i1].flag != edges[i2].flag) {
            vertices[i].SetRequired();
            if (verbosity>7)
               cout << "     The vertex " << i << " is Required the flag change (crack or equi)" << endl;
         }
         if (edges[i1].ref != edges[i2].ref) {
            vertices[i].SetRequired();
            if (verbosity>7)
               cout << "     The vertex " << i << " is Required ref" << endl;
         }
      }

      if (ord != 2) {
         vertices[i].SetCorner();
         if (verbosity>7)
            cout << "     the vertex " << i << " is a corner ordre = " << ord << endl;
      }

//    close the list around the vertex 
      {
         long no = -1, ne = hv[i];
         while (ne>=0)
            ne = ev[no=ne];
         if (no>=0) 
            ev[no] = hv[i];
      } // now the list around the vertex is circular
   }

   k = 0;
   for (i=0; i<nbe; i++) {
      for (jj=0; jj<2; jj++) {
         long n1 = ev[k++]; 
         long i1 = n1/2, j1 = n1%2;
         if (edges[i1].v[j1] != edges[i].v[jj]) {
            cerr << " Bug Adj edge " << i << " " << jj << 
                    " and " << i1 << " " << j1 << " k = " << k;
            cerr << Number(edges[i].v[jj]) <<" <> " << Number(edges[i1].v[j1]) << endl;
            cerr << "edge " << Number(edges[i].v[0]) << " " << Number(edges[i].v[1]) << endl; 
            MeshError(1);
         }
         edges[i1].Adj[j1] = edges + i;
         edges[i1].SensAdj[j1] = jj;
         if (verbosity>10)
	         cout << " edges. Adj " << i1 << " " << j1 << " <--- " << i << " " << jj << endl;
      }
   }

// generation of  all the tangents
   for (i=0; i<nbe; i++) {
      R2 AB = edges[i].v[1]->r - edges[i].v[0]->r;        
      double lAB = Norme2(AB); // length of current edge AB
      double ltg2[2];
      ltg2[0] = 0; ltg2[1] = 0;
      for (jj=0; jj<2; jj++) {
         R2 tg=edges[i].tg[jj];
         double ltg=Norme2(tg); // length of tg
         if (ltg==0) {// no tg
            if (!edges[i].v[jj]->Corner()) {
               tg =  edges[i].v[1-jj]->r 
                   - edges[i].Adj[jj]->v[1-edges[i].SensAdj[jj]]->r;
               ltg =  Norme2(tg);
               tg = tg*(lAB/ltg), ltg = lAB;
            }
         }
         else 
            tg = tg*(lAB/ltg),ltg=lAB;
         ltg2[jj] = ltg;
         if ((tg,AB)<0) 
            tg = -tg;
         edges[i].tg[jj] = tg;
      } 
      if (ltg2[0]!=0)
         edges[i].SetTgA();
      if (ltg2[1]!=0)
         edges[i].SetTgB();
   }

   if (verbosity>7)
      for (i=0; i<nbv; i++)
         if (vertices[i].Required())
            cout << "     The  geo  vertex " << i << " is required " << endl;

   for (int step=0; step<2; step++) {
      for (i=0; i<nbe; i++)
         edges[i].SetUnMark();
      NbOfCurves = 0;
      long nbgem=0;
      for (int level=0; level<2 && nbgem!=nbe; level++) {
         for (i=0; i<nbe; i++) {
            GeometricalEdge &ei = edges[i];   
            for (jj=0; jj<2; jj++) {
               if (!ei.Mark() && (level || ei[jj].Required())) {
                  int k0=jj, k1;
                  GeometricalEdge *e = &ei;
                  GeometricalVertex *a=(*e)(k0);
                  if (curves) {
                     curves[NbOfCurves].be = e;
                     curves[NbOfCurves].kb = k0;
                  }
                  for (;;) { 
                     k1 = 1 - k0; // next vertex of the edge 
                     e->SetMark();
                     nbgem++;
                     e->CurveNumber = NbOfCurves;
                     if (curves) {
                        curves[NbOfCurves].ee = e;
                        curves[NbOfCurves].ke = k1;
                     }
                     GeometricalVertex *b=(*e)(k1);
                     if (a==b || b->Required())
                        break;
                     k0 = e->SensAdj[k1];   // vertex in next edge
                     e = e->Adj[k1];        // next edge 
                  }
                  NbOfCurves++;
                  if (level) {
                     if (verbosity>4)
                        cout << "    Warning: Curve "<< NbOfCurves << " without required vertex " 
                             << "so the vertex " << Number(a) << " becomes required " <<endl;
                     a->SetRequired();
                  }
               }
            }
         }
      }
      assert(nbgem && nbe);
      if (step==0)
         curves = new Curve[NbOfCurves];
   }

   for (int i=0; i<NbOfCurves; i++) {
      GeometricalEdge *be=curves[i].be, *eqbe=be->link;
      curves[i].master = true;
      if (be->Equi() || be->ReverseEqui()) {
         assert(eqbe);
         int nc = eqbe->CurveNumber;
         assert(i!=nc);
         curves[i].next = curves[nc].next;
         curves[i].master = false;
         curves[nc].next = curves + i;
         if (be->ReverseEqui())
           curves[i].Reverse();           
      }
   }
 	 
   if (verbosity>3)
      cout << "    End ReadGeometry: Number of curves in geometry is " << NbOfCurves << endl; 
   if (verbosity>4) {
      for (int i=0; i<NbOfCurves; i++) {
         cout << " Curve " << i << " begin e=" << Number(curves[i].be) << " k=" << curves[i].kb 
              << "  end e= " << Number(curves[i].ee) << " k=" << curves[i].ke << endl;
      }
   }
   delete [] ev;
   delete [] hv;
   delete [] eangle;
}


Geometry::~Geometry() 
{
   assert(NbRef<=0);
   if (vertices)
      delete [] vertices;
   vertices = NULL;
   if (edges)
      delete [] edges;
   edges = NULL;
   if (triangles)
      delete [] triangles;
   triangles = NULL;
   if (_quadtree)
      delete _quadtree;
   _quadtree = NULL;
   if (curves)
      delete [] curves;
   curves = NULL;
   NbOfCurves = 0;
   if (name)
      delete [] name;
   name = NULL;
   if (subdomains)
      delete [] subdomains;
   subdomains = NULL;
   EmptyGeometry();
}


double GeometricalEdge::R1tg(double theta,
                             R2&    t) const // 1/R of radius of curvature
{
   R2 A=v[0]->r, B=v[1]->r;
   double dca, dcb, dcta, dctb, ddca, ddcb, ddcta, ddctb, tt=theta*theta;
   assert(theta>=0);
   assert(theta<=1);
   if (TgA()) { 
      if (TgB()) {// Hermite interpolation
         dcb = 6*theta*(1 - theta);
         ddcb = 6*(1 - 2*theta);
         dca = -dcb;
         ddca = -ddcb;
         dcta = (3*theta - 4)*theta + 1;
         ddcta = 6*theta - 4;
         dctb = 3*tt - 2*theta;
         ddctb = 6*theta-2;
      }
      else { // 1-t*t, t-t*t, t*t
         double t=theta;
         dcb = 2*t;
         ddcb = 2;
         dca = -dcb;
         ddca = -2;
         dcta = 1 - dcb;
         ddcta = -ddcb;
         dctb = 0;    
         ddctb = 0;    
      }
   }
   else {
      if (TgB()) {
         double t=1-theta;
         dca = -2*t;
         ddca = 2;
         dcb = -dca;
         ddcb = -2;
         dctb = 1 + dca;
         ddctb = ddca;
         dcta = 0;
         ddcta = 0;
      }
      else {
         t = B - A;
         return 0;
      }
   } // Lagrange P1
   R2 d = A*dca + B*dcb + tg[0]* dcta + tg[1] * dctb;
   R2 dd = A*ddca + B*ddcb + tg[0]* ddcta + tg[1] * ddctb;
   double d2=(d,d);
   double sd2 = sqrt(d2);
   t = d;
   if (d2>1.0e-20) {
      t /= sd2;
      return Abs(Det(d,dd))/(d2*sd2);
   }
   else return 0;
}


R2 GeometricalEdge::F(double theta) const
{
   R2 A=v[0]->r, B=v[1]->r;
   double ca, cb, cta, ctb;
   assert(theta>=-1e-12);
   assert(theta<=1+1e-12);
   if (TgA()) {
      if (TgB()) {   // Hermite interpolation
         cb = theta*theta*(3-2*theta);
         ca = 1 - cb;
         cta = (1-theta)*(1-theta)*theta;
         ctb = (theta-1)*theta*theta;
      }
      else { // 1-t*t, t-t*t, t*t
         double t=theta;
         cb = t*t;
         ca = 1 - cb;
         cta = t - cb;
         ctb = 0;
      }
   }
   else {
      if (TgB()) {
         double t=1-theta;
         ca = t*t;
         cb = 1 - ca;
         ctb= -t + ca;
         cta = 0;    
      }
      else
         ca = (1-theta), cb = theta, cta = ctb = 0; // Lagrange P1
   }
   return A*ca + B*cb + tg[0]* cta + tg[1] * ctb;
}

}
