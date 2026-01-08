// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY:  Bamg: Bidimensional Anisotrope Mesh Generator
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

#include <stdio.h>
#include "mesh/bamg/Meshio.h"
#include "mesh/bamg/Mesh2.h"

namespace bamg {

inline double det3x3(double A[3],
                     double B[3],
                     double C[3])
{
   return   A[0]*(B[1]*C[2]-B[2]*C[1])
          - A[1]*(B[0]*C[2]-B[2]*C[0])
          + A[2]*(B[0]*C[1]-B[1]*C[0]);
}


SaveMetricInterpole LastMetricInterpole;


void ReductionSimultanee(MetricAnIso M1,
                         MetricAnIso M2,
                         double&     l1,
                         double&     l2,
                         D2xD2&      V) 
{
   double a11=M1.a11, a21=M1.a21, a22=M1.a22;
   double b11=M2.a11, b21=M2.a21, b22=M2.a22;
   const double c21= a21*a21, d21= b21*b21;
   const double a=b11*b22 - d21;
   const double b=-a11*b22-a22*b11+2*a21*b21;
   const double c=-c21+a11*a22;
   const double bb=b*b, ac= a*c;
   const double delta = bb - 4 * ac;
   if (bb + Abs(ac) < 1.0e-20 || (delta< 1.0E-4*bb)) {
      if (Abs(a) < 1.e-30)
         l1 = l2 = 0;
      else 
         l1 = l2 = -b/(2*a); 
      V = D2xD2(1,0,0,1);
   }
   else {
      const double delta2 = sqrt(delta);
      l1 = (-b - delta2)/(2*a);
      l2 = (-b + delta2)/(2*a);
      double v0 = a11-l1*b11, v1 = a21-l1*b21,v2 = a22 - l1*b22;
      double s0 = v0*v0 + v1*v1, s1 = v1*v1 +v2*v2;
      double vp1x, vp1y, vp2x, vp2y;
      if (s1 < s0)
         s0 = sqrt(s0), vp1x = v1/s0, vp1y = -v0/s0;
      else
         s1 = sqrt(s1), vp1x = v2/s1, vp1y = -v1/s1;
      v0 = a11-l2*b11, v1 = a21-l2*b21,v2 = a22 - l2*b22;
      s0 = v0*v0 + v1*v1, s1 = v1*v1 +v2*v2;
      if (s1 < s0)
         s0 = sqrt(s0), vp2x = v1/s0, vp2y = -v0/s0;
      else
         s1 = sqrt(s1), vp2x = v2/s1, vp2y = -v1/s1;
      V = D2xD2(vp1x,vp2x,vp1y,vp2y);
  }
  return;
}


MetricAnIso Intersection(const MetricAnIso M1,
                         const MetricAnIso M2);

MetricAnIso Intersection(const MetricAnIso M1,
                         const MetricAnIso M2) 
{
   D2xD2 M;
   double l1, l2;
   ReductionSimultanee(M1,M2,l1,l2,M);
   R2 v0(M.x.x,M.y.x);
   R2 v1(M.x.y,M.y.y);
   D2xD2 M_1(M.inv());
   D2xD2 D(Max(M1(v0,v0),M2(v0,v0)),0,0,Max(M1(v1,v1),M2(v1,v1)));
   D2xD2 Mi(M_1.t()*D*M_1);
   return MetricAnIso(Mi.x.x,0.5*(Mi.x.y+Mi.y.x),Mi.y.y);
}


MetricAnIso::MetricAnIso(const double      a[3],
                         const MetricAnIso m0,
                         const MetricAnIso m1,
                         const MetricAnIso m2)
{
   MetricAnIso mab(a[0]*m0.a11 + a[1]*m1.a11 + a[2]*m2.a11,
                   a[0]*m0.a21 + a[1]*m1.a21 + a[2]*m2.a21,
                   a[0]*m0.a22 + a[1]*m1.a22 + a[2]*m2.a22);
   MatVVP2x2 vab(mab);
   R2 v1(vab.v.x,vab.v.y);
   R2 v2(-v1.y,v1.x);
   double h1 = a[0] / m0(v1) + a[1] / m1(v1) + a[2] / m2(v1);
   double h2 = a[0] / m0(v2) + a[1] / m1(v2) + a[2] / m2(v2);
   vab.lambda1 =  1 / (h1*h1);
   vab.lambda2 =  1 / (h2*h2);
   *this = vab;
}


MetricAnIso::MetricAnIso(double            a,
                         const MetricAnIso ma,
                         double            b,
                         const MetricAnIso mb)
{
   MetricAnIso mab(a*ma.a11+b*mb.a11,a*ma.a21+b*mb.a21,a*ma.a22+b*mb.a22);
   MatVVP2x2 vab(mab);
   R2 v1(vab.v.x,vab.v.y);
   R2 v2(-v1.y,v1.x);
   double h1=a/ma(v1)+b/mb(v1), h2=a/ma(v2)+b/mb(v2);
   vab.lambda1 = 1/(h1*h1);
   vab.lambda2 = 1/(h2*h2);
   *this = vab;
}


MatVVP2x2::MatVVP2x2(const MetricAnIso M) 
{
   double a11=M.a11, a21=M.a21, a22=M.a22;
   const double eps = 1.e-5;
   double c11=a11*a11, c22=a22*a22, c21=a21*a21;
   double b=-a11-a22,c=-c21+a11*a22;
   double delta = b*b - 4*c;
   double n2=(c11+c22+c21);
   if (n2 < 1e-30)
      lambda1=lambda2=0, v.x=1, v.y=0;
   else if (delta < eps*n2)
      lambda1 = lambda2 = -b/2, v.x = 1, v.y = 0;
   else {
      delta = sqrt(delta);
      lambda1 = (-b-delta)/2.0,lambda2 = (-b+delta)/2.0;
      double v0=a11-lambda1, v1=a21, v2=a22-lambda1;
      double s0=v0*v0+v1*v1, s1=v1*v1+v2*v2;
      if (s1<s0)
         s0 = sqrt(s0), v.x = v1/s0, v.y = -v0/s0;
      else
         s1 = sqrt(s1), v.x = v2/s1, v.y = -v1/s1;
   }
}


int MetricAnIso::IntersectWith(const MetricAnIso M2) 
{
   int r=0;
   MetricAnIso &M1 = *this;
   D2xD2 M;
   double l1, l2;
   ReductionSimultanee(*this,M2,l1,l2,M);
   R2 v1(M.x.x,M.y.x);
   R2 v2(M.x.y,M.y.y);
   double l11=M1(v1,v1), l12=M1(v2,v2), l21=M2(v1,v1), l22=M2(v2,v2);
   if (l11<l21)
      r = 1, l11 = l21;
   if (l12<l22)
      r = 1, l12 = l22; 
   if (r) {
      D2xD2 M_1(M.inv());
      D2xD2 D(l11,0,0,l12); 
      D2xD2 Mi(M_1.t()*D*M_1);
      a11 = Mi.x.x;
      a21 = 0.5*(Mi.x.y+Mi.y.x);
      a22 = Mi.y.y;
   }
   return r;
}


void Triangles::IntersectGeomMetric(const double err=1,
                                    const int    iso=0)
{
   if (verbosity>1)
      cout << " -- IntersectGeomMetric geometric err = " << err << (iso ? " iso " : " aniso ") << endl;
   double ss[2]={0.00001,0.99999};
   double errC=2*sqrt(2*err);
   double hmax = Gh.MaximalHmax(), hmin = Gh.MinimalHmin();
   double maxaniso = 1e6;
   assert(hmax>0);
   SetVertexFieldOn();
   if (errC>1)
      errC = 1;
   for (long i=0; i<nbe; i++) {
      for (int j=0; j<2; j++) {
         Vertex V;
         VertexOnGeom GV;
         Gh.ProjectOnCurve(edges[i],ss[j],V,GV);
         {
            GeometricalEdge *eg = GV;
            double s=GV;
            R2 tg;
            double R1=eg->R1tg(s,tg), ht=hmax;
            if (R1>1.0e-20)
               ht = Min(Max(errC/R1,hmin),hmax);
            double hn = iso? ht : Min(hmax,ht*maxaniso);
            assert(ht>0 && hn>0);
            MatVVP2x2 Vp(1/(ht*ht),1/(hn*hn),tg);
            Metric MVp(Vp);
            edges[i][j].m.IntersectWith(MVp);
         }
      }
   }
}


void Triangles::BoundAnisotropy(double anisomax,
                                double hminaniso)
{
   double lminaniso=1/(Max(hminaniso*hminaniso,1e-100));
   if (verbosity>1)
      cout << " -- BoundAnisotropy by  " << anisomax << endl; 
   double h1=1.e30, h2=1e-30, rx=0;
   double coef=1./(anisomax*anisomax);
   double hn1=1.e30, hn2=1e-30, rnx=1.e-30;  
   for (long i=0; i<nbv; i++) {
      MatVVP2x2 Vp(_vertices[i]);
      double lmax=Vp.lmax();
      h1 = Min(h1,Vp.lmin());
      h2 = Max(h2,Vp.lmax());
      rx = Max(rx,Vp.Aniso2());
      Vp *= Min(lminaniso,lmax)/lmax;
      Vp.BoundAniso2(coef);
      hn1 = Min(hn1,Vp.lmin());
      hn2 = Max(hn2,Vp.lmax());
      rnx = Max(rnx,Vp.Aniso2());
      _vertices[i].m = Vp;
   }

   if (verbosity>2) {
      cout << "     Input:  Hmin = " << sqrt(1/h2)  << " Hmax = " << sqrt(1/h1) 
           << " factor of anisotropy max  = " << sqrt(rx) << endl;
      cout << "     output:  Hmin = " << sqrt(1/hn2) << " Hmax = " << sqrt(1/hn1) 
           << " factor of anisotropy max  = " << sqrt(rnx) << endl;
   }
}


void Triangles::IntersectConsMetric(const double* s,
                                    const long    nbsol,
                                    const int*    typsols,
                                    const double  hmin1,
                                    const double  hmax1,
                                    const double  coef,
                                    const double  anisomax,
                                    const double  CutOff,
                                    const int     NbJacobi,
                                    const int     DoNormalisation,
                                    const double  power,
                                    const int     choice)
{
  //  the array of solution s is store    
  // sol0,sol1,...,soln    on vertex 0
  //  sol0,sol1,...,soln   on vertex 1
  //  etc.
  //  choise = 0 =>  H is computed with green formule
  //   otherwise  => H is computed from P2 on 4T
 
   const int dim = 2;
   int sizeoftype[] = {1, dim, dim*(dim+1)/2, dim*dim};

// computation of the nb of field 
   long ntmp = 0;
   if (typsols) {
      for (long i=0; i<nbsol; i++)
	     ntmp += sizeoftype[typsols[i]];
   }
   else
      ntmp = nbsol;

  const long n = ntmp;
  long i, k, iA, iB, iC, iv;
  R2 O(0,0);
  int RelativeMetric = CutOff>1e-30;
  double hmin=Max(hmin1,MinimalHmin()), hmax=Min(hmax1,MaximalHmax());
  double coef2=1/(coef*coef);
  if (verbosity>1) {
     cout << " -- Construction of Metric: Nb of field. " << n << " nbt = " 
          << nbt << " nbv= " << nbv 
          << " coef = " << coef << endl
          << "     hmin = " << hmin << " hmax=" << hmax 
          << " anisomax = " << anisomax <<  " Nb Jacobi " << NbJacobi << " Power = " << power;
      if (RelativeMetric)
         cout << " RelativeErr with CutOff= "  <<  CutOff << endl;
      else
         cout << " Absolute Err" <<endl;
   }
   double *ss=(double*)s;
   double sA, sB, sC;

   double *detT = new double[nbt];
   double *Mmass = new double[nbv];
   double *Mmassxx = new double[nbv];
   double *dxdx = new double[nbv];
   double *dxdy = new double[nbv];
   double *dydy = new double[nbv];
   double *workT = new double[nbt];
   double *workV = new double[nbv];
   int *OnBoundary = new int[nbv];
   for (iv=0; iv<nbv; iv++) {
      Mmass[iv] = 0;
      OnBoundary[iv] = 0;
      Mmassxx[iv] = 0;
   }
   for (i=0; i<nbt; i++) {
      if (triangles[i].link) {
         const Triangle &t=triangles[i];
         R2 A=t[0], B=t[1], C=t[2];
         iA = Number(t[0]);
         iB = Number(t[1]);
         iC = Number(t[2]);
         double dett = bamg::Area2(A,B,C);
         detT[i] = dett;
         dett /= 6;

//       construction of on boundary 
         int nbb=0;
         for (int j=0; j<3; j++) {
            Triangle *ta=t.Adj(j);
            if (!ta || !ta->link) // no adj triangle => edge on boundary
               OnBoundary[Number(t[VerticesOfTriangularEdge[j][0]])] = 1,
            OnBoundary[Number(t[VerticesOfTriangularEdge[j][1]])] = 1,
            nbb++;
         }
         workT[i] = nbb;
         Mmass[iA] += dett;
         Mmass[iB] += dett;
         Mmass[iC] += dett;
         if (nbb==0|| !choice) {
            Mmassxx[iA] += dett;
            Mmassxx[iB] += dett;
            Mmassxx[iC] += dett;
         }
      }
      else
         workT[i] = -1;
   }

   for (long nusol=0; nusol<nbsol; nusol++) {
      double smin=ss[0], smax=ss[0];
      double h1=1.e30, h2=1e-30, rx=0;
      double coef = 1./(anisomax*anisomax);
      double hn1=1.e30, hn2=1e-30, rnx =1.e-30;  
      int nbfield = typsols? sizeoftype[typsols[nusol]] : 1; 
      if (nbfield==1) {
         for (iv=0, k=0; iv<nbv; iv++, k+=n) {
            dxdx[iv] = dxdy[iv] = dydy[iv] = 0;
            smin = Min(smin,ss[k]);
            smax = Max(smax,ss[k]);
         }
      }
      else {
         for (iv=0, k=0; iv<nbv; iv++, k+=n) {
            double v=0;		     
            for (int i=0; i<nbfield; i++)
               v += ss[k+i]*ss[k+i];
            v = sqrt(v);
            smin = Min(smin,v);
            smax = Max(smax,v);
         }
      }
      double sdelta=smax-smin;
      double absmax=Max(Abs(smin),Abs(smax));
      double cnorm=DoNormalisation ? coef2/sdelta : coef2;
      if (verbosity>2) 
         cout << "    Solution " << nusol <<  " Min = " << smin << " Max = " 
              << smax << " Delta =" << sdelta << " cnorm = " << cnorm
              << " Nb of fields =" << nbfield << endl;
      if (sdelta < 1.0e-10*Max(absmax,1e-20) && (nbfield ==1)) {
         if (verbosity>2)
            cout << "      Solution " << nusol << " is constant. We skip. " 
                 << " Min = " << smin << " Max = " << smax << endl;
         continue;
      }

      double *sf=ss; 
      for (long nufield=0; nufield<nbfield; nufield++, ss++) {
         for (iv=0, k=0; iv<nbv; iv++, k+=n)
            dxdx[iv] = dxdy[iv] = dydy[iv] = 0;
         for (i=0; i<nbt; i++) { 
            if (triangles[i].link) {
               R2 A=triangles[i][0], B=triangles[i][1], C=triangles[i][2];
               R2 nAB=Orthogonal(B-A), nBC=Orthogonal(C-B), nCA=Orthogonal(A-C);
               iA = Number(triangles[i][0]);
               iB = Number(triangles[i][1]);
               iC = Number(triangles[i][2]);
               Triangle *tBC = triangles[i].TriangleAdj(OppositeEdge[0]);
               Triangle *tCA = triangles[i].TriangleAdj(OppositeEdge[1]);
               Triangle *tAB = triangles[i].TriangleAdj(OppositeEdge[2]);
               sA = ss[iA*n], sB = ss[iB*n], sC = ss[iC*n];
               R2 Grads = (nAB*sC + nBC*sA + nCA*sB)/detT[i];
               if (choice) {
                  int nbb=0;
                  double dd=detT[i], taa[3][3], bb[3];
                  for (int j=0; j<3; j++) {
                     int ie = OppositeEdge[j];
                     TriangleAdjacent ta=triangles[i].Adj(ie);
                     Triangle *tt=ta;
                     if (tt && tt->link) {
                        Vertex &v = *ta.OppositeVertex();
                        R2 V = v;
                        long iV = Number(v);
                        double lA = bamg::Area2(V,B,C)/dd;
                        double lB = bamg::Area2(A,V,C)/dd;
                        double lC = bamg::Area2(A,B,V)/dd;
                        taa[0][j] = lB*lC;
                        taa[1][j] = lC*lA;
                        taa[2][j] = lA*lB;
                        bb[j] = ss[iV*n] - (sA*lA + sB*lB + sC*lC);
                     }
                     else {
                        nbb++;
                        taa[0][j] = 0;
                        taa[1][j] = 0;
                        taa[2][j] = 0;
                        taa[j][j] = 1;
                        bb[j] = 0;
                     }
                  }

                  double det33 =  det3x3(taa[0],taa[1],taa[2]);		
                  double cBC   =  det3x3(bb,taa[1],taa[2]);
                  double cCA   =  det3x3(taa[0],bb,taa[2]);
                  double cAB   =  det3x3(taa[0],taa[1],bb);
                  assert(det33);
                  double Hxx = cAB*(nBC.x*nCA.x) + cBC*(nCA.x*nAB.x) + cCA*(nAB.x*nBC.x);
                  double Hyy = cAB*(nBC.y*nCA.y) +  cBC*(nCA.y*nAB.y) + cCA*(nAB.y*nBC.y);
                  double Hxy = cAB*(nBC.y*nCA.x) +  cBC*(nCA.y*nAB.x) + cCA*(nAB.y*nBC.x) 
                             + cAB*(nBC.x*nCA.y) +  cBC*(nCA.x*nAB.y) + cCA*(nAB.x*nBC.y);
                  double coef = 1.0/(3*dd*det33);
                  double coef2 = 2*coef;
                  Hxx *= coef2, Hyy *= coef2, Hxy *= coef2;
                  if (nbb==0) {
                     dxdx[iA] += Hxx;
                     dydy[iA] += Hyy;
                     dxdy[iA] += Hxy;
                     dxdx[iB] += Hxx;
                     dydy[iB] += Hyy;
                     dxdy[iB] += Hxy;
                     dxdx[iC] += Hxx;
                     dydy[iC] += Hyy;
                     dxdy[iC] += Hxy;
                  }
               }
               else {
//                if edge on boundary no contribution  => normal = 0
                  if (!tBC || !tBC->link)
                     nBC = O;
                  if (!tCA || !tCA->link)
                     nCA = O;
                  if (!tAB || !tAB->link)
                     nAB = O;
//                Remark we forgot a 1/2 because
//                \int_{edge} w_i = 1/2 if i is in edge 
//                                  0  if not
//                if we don't take the boundary 
//                dxdx[iA] += ( nCA.x + nAB.x ) *Grads.x;

                  dxdx[iA] += (nCA.x + nAB.x)*Grads.x;
                  dxdx[iB] += (nAB.x + nBC.x)*Grads.x;
                  dxdx[iC] += (nBC.x + nCA.x)*Grads.x;

//                warning optimization (1) the divide by 2 is done on the metrix construction
                  dxdy[iA] += ((nCA.y + nAB.y)*Grads.x + (nCA.x + nAB.x)*Grads.y) ;
                  dxdy[iB] += ((nAB.y + nBC.y)*Grads.x + (nAB.x + nBC.x)*Grads.y) ;
                  dxdy[iC] += ((nBC.y + nCA.y)*Grads.x + (nBC.x + nCA.x)*Grads.y) ; 
                  dydy[iA] += (nCA.y + nAB.y)*Grads.y;
                  dydy[iB] += (nAB.y + nBC.y)*Grads.y;
                  dydy[iC] += (nBC.y + nCA.y)*Grads.y;
               }
            }
         }

         for (iv=0, k=0; iv<nbv; iv++, k+=n) {
            if (Mmassxx[iv]>0) {
               dxdx[iv] /= 2*Mmassxx[iv];
//             warning optimization (1) on term dxdy[iv]*ci/2 
               dxdy[iv] /= 4*Mmassxx[iv];
               dydy[iv] /= 2*Mmassxx[iv];
//             Compute the matrix with abs(eigenvalue)
               Metric M(dxdx[iv], dxdy[iv], dydy[iv]);
               MatVVP2x2 Vp(M);
               Vp.Abs();
               M = Vp;
               dxdx[iv] = M.a11;
               dxdy[iv] = M.a21;
               dydy[iv] = M.a22;
            }
         }

//       Correction of second-order derivative by a laplacian
         double *d2[3] = {dxdx, dxdy, dydy};
         double *dd;
         for (int xy=0; xy<3; xy++) {
            dd = d2[xy];
            for (int ijacobi=0; ijacobi<Max(NbJacobi,2); ijacobi++) {
               for (i=0; i<nbt; i++) {
                  if (triangles[i].link) {
                     iA = Number(triangles[i][0]);
                     iB = Number(triangles[i][1]);
                     iC = Number(triangles[i][2]);
                     double cc=3;
                     if (ijacobi==0)
                        cc = Max((double) ((Mmassxx[iA]>0)+(Mmassxx[iB]>0)+(Mmassxx[iC]>0)),1.);
                     workT[i] = (dd[iA]+dd[iB]+dd[iC])/cc;
                  }
               }
               for (iv=0; iv<nbv; iv++)
                  workV[iv] = 0;
               for (i=0; i<nbt; i++) {
                  if (triangles[i].link) {
                     iA = Number(triangles[i][0]);
                     iB = Number(triangles[i][1]);
                     iC = Number(triangles[i][2]);
                     double cc = workT[i]*detT[i];
                     workV[iA] += cc;
                     workV[iB] += cc;
                     workV[iC] += cc;
                  }
               }
               for (iv=0; iv<nbv; iv++) {
                  if (ijacobi<NbJacobi || OnBoundary[iv])
                     dd[iv] = workV[iv]/(Mmass[iv]*6);
               }
            }
         }

//       Constuction of the matrix from the Hessian dxdx. dxdy,dydy
         double rCutOff=CutOff*absmax;// relative cut off
         for (iv=0,k=0; iv<nbv; iv++, k+=n) { // for all vertices 
            MetricIso Miso;
//          New code to compute ci ---	  
            double ci;
            if (RelativeMetric) { //   compute the norm of the solution
               double xx =0, *sfk=sf+k;
               for (int ifield=0; ifield<nbfield; ifield++ ,sfk++)
                  xx += (*sfk) * (*sfk);	       
               xx = sqrt(xx);
               ci = coef2/Max(xx,rCutOff);
            }
            else
               ci = cnorm;

            Metric Miv(dxdx[iv]*ci,dxdy[iv]*ci,dydy[iv]*ci);
            MatVVP2x2 Vp(Miv);
            Vp.Abs();
            if (power!=1.0) 
               Vp.pow(power);
            h1 = Min(h1,Vp.lmin());
            h2 = Max(h2,Vp.lmax());
            Vp.Maxh(hmin);
            Vp.Minh(hmax);
            rx = Max(rx,Vp.Aniso2());
            Vp.BoundAniso2(coef);
            hn1 = Min(hn1,Vp.lmin());
            hn2 = Max(hn2,Vp.lmax());
            rnx = Max(rnx,Vp.Aniso2());
            Metric MVp(Vp);
            _vertices[iv].m.IntersectWith(MVp);
         }

         if (verbosity>2) {
            cout << "              Field " << nufield << " of solution " << nusol  << endl;
            cout << "              before bounding :  Hmin = " << sqrt(1/h2) << " Hmax = " 
                 << sqrt(1/h1)  << " factor of anisotropy max  = " << sqrt(rx) << endl;
            cout << "              after  bounding :  Hmin = " << sqrt(1/hn2) << " Hmax = " 
                 << sqrt(1/hn1)  << " factor of anisotropy max  = " << sqrt(rnx) << endl;
         }
      }
   }

   delete [] detT;
   delete [] Mmass;
   delete [] dxdx;
   delete [] dxdy;
   delete [] dydy;
   delete [] workT;
   delete [] workV;
   delete [] Mmassxx;
   delete [] OnBoundary;
}


void Triangles::ReadMetric(const char* fmetrix,
                           const double hmin1=1.0e-30,
                           const double hmax1=1.0e30,
                           const double coef=1)
{
   double hmin = Max(hmin1,MinimalHmin());
   double hmax = Min(hmax1,MaximalHmax());
   MeshIstream f_metrix(fmetrix);
   long k, j;
   f_metrix >> k >> j;
   if (verbosity>1)
      cout << " metrix: open " << fmetrix 
           << ", le coef = " << coef
           << ", hmin = " << hmin << ", hmax = " << hmax 
           << ((j==1) ? " Iso " : " AnIso ") << endl;
  
   if (k!=nbv || !(j==1 || j==3)) {
      cerr << " Error Pb metrix " << k << " <> " 
	   <<  nbv << " or  1 or 3 <> " << j << endl;
      MeshError(1002);
   }

   for (long iv=0; iv<nbv; iv++) {
      double h;
      if (j == 1) {
         f_metrix >> h;
         _vertices[iv].m = Metric(Max(hmin,Min(hmax, h*coef)));
      }
      else if (j==3) {
         double a, b, c;	     
         f_metrix >> a >> b >> c;
         MetricAnIso M(a,b,c);
         MatVVP2x2 Vp(M/coef);
         Vp.Maxh(hmin);
         Vp.Minh(hmax);
         _vertices[iv].m = Vp;
      }
   }
}


void Triangles::WriteMetric(ostream& f,
                            int      iso)
{
   if (iso) {
      f << nbv << " " << 1 << endl;
      for (long iv=0; iv<nbv; iv++) {
         MatVVP2x2 V=_vertices[iv].m;
         f << V.hmin() << endl;
      }
   }
   else {
      f << nbv << " " << 3 << endl;
      for (long iv=0; iv<nbv; iv++)
         f << _vertices[iv].m.a11 << " " << _vertices[iv].m.a21 << " " 
           << _vertices[iv].m.a22 << endl;
   }
}


void  Triangles::MaxSubDivision(double maxsubdiv)
{
   const  double maxsubdiv2 = maxsubdiv*maxsubdiv;
   if (verbosity>1)
      cout << " -- Limit the subdivision of a edges in the new mesh by " << maxsubdiv << endl;
   long nbchange=0;    
   double lmax=0;
   for (long it=0; it<nbt; it++) {
      Triangle &t=triangles[it];
      for (int j=0; j<3; j++) {
         Triangle &tt=*t.TriangleAdj(j);
         Triangle *ttt=&tt;
         if ((!ttt || it<Number(tt)) && (tt.link || t.link)) {
            Vertex &v0 = t[VerticesOfTriangularEdge[j][0]];
            Vertex &v1 = t[VerticesOfTriangularEdge[j][1]];
            R2 AB= (R2)v1 - (R2)v0;
            Metric M = v0;
            double l = M(AB,AB);
            lmax = Max(lmax,l);
            if (l> maxsubdiv2) {
               R2 AC = M.Orthogonal(AB);// the ortogonal vector of AB in M
               double lc = M(AC,AC);
               D2xD2 Rt(AB,AC);// Rt.x = AB , Rt.y = AC;
               D2xD2 Rt1(Rt.inv());
               D2xD2 D(maxsubdiv2,0,0,lc);
               D2xD2 MM = Rt1*D*Rt1.t();
               v0.m = M = MetricAnIso(MM.x.x,MM.y.x,MM.y.y);
               nbchange++;
            }
            M = v1;
            l = M(AB,AB);
            lmax = Max(lmax,l);
            if (l> maxsubdiv2) {
               R2 AC = M.Orthogonal(AB);// the orthogonal vector of AB in M
               double lc = M(AC,AC);
               D2xD2 Rt(AB,AC);
               D2xD2 Rt1(Rt.inv());
               D2xD2 D(maxsubdiv2,0,0,lc);
               D2xD2 MM = Rt1*D*Rt1.t();
               v1.m = M = MetricAnIso(MM.x.x,MM.y.x,MM.y.y);
               nbchange++;
            }		
         }
      }
   }
   if (verbosity>3)
      cout << "    Nb of metric change = " << nbchange 
           << " Max of the subdivision of a edges before change = " << sqrt(lmax) << endl;
}


void Triangles::SmoothMetric(double raisonmax) 
{
   if (raisonmax<1.1)
      return;
   if (verbosity > 1)
      cout << " -- Triangles::SmoothMetric raisonmax = " << raisonmax << " " <<nbv <<endl;
   ReMakeTriangleContainingTheVertex();
   long i, j, kch, kk, ip;
   long *first_np_or_next_t0 = new long [nbv];
   long *first_np_or_next_t1 = new long [nbv];
   long Head0 =0, Head1=-1;
   double logseuil=log(raisonmax);

   for (i=0; i<nbv-1; i++)
      first_np_or_next_t0[i] = i+1;
   first_np_or_next_t0[nbv-1] = -1;
   for (i=0; i<nbv; i++)
      first_np_or_next_t1[i] = -1;
   kk = 0;
   while (Head0>=0 && kk++<100) {
      kch = 0;
      for (i=Head0; i>=0; i=first_np_or_next_t0[ip=i], first_np_or_next_t0[ip]=-1) {
//       for all triangles around the vertex s
         Triangle *t=_vertices[i].t;
         assert(t);
         Vertex &vi=_vertices[i];
         TriangleAdjacent ta(t,EdgesVertexTriangle[_vertices[i].vint][0]);
         Vertex *pvj0 = ta.EdgeVertex(0);
         while (1) {
            ta = Previous(Adj(ta));
            assert(_vertices+i == ta.EdgeVertex(1));
            Vertex *vj = &(*ta.EdgeVertex(0));
            if (vj) {
               j = vj - _vertices;
               assert(j>=0 && j<nbv);
               R2 Aij = (R2)(*vj) - (R2)vi;
               double ll=Norme2(Aij);
               if (0) {
                  double hi=ll/vi.m(Aij), hj=ll/vj->m(Aij);
                  if (hi < hj) {
                     double dh=(hj-hi)/ll;
                     if (dh>logseuil) {
                        (*vj).m.IntersectWith(vi.m/(1 +logseuil*ll/hi));
                        if (first_np_or_next_t1[j]<0)
                           kch++, first_np_or_next_t1[j]=Head1, Head1=j;
                     }
                  }
               }
               else {
                  double li = vi.m(Aij);
                  if ((*vj).m.IntersectWith(vi.m/(1 +logseuil*li)))
                     if (first_np_or_next_t1[j]<0) // if the metrix change 
                        kch++, first_np_or_next_t1[j]=Head1, Head1=j;
               }
            }
            if (vj == pvj0)
               break;
         }
      }
      Head0 = Head1;
      Head1 = -1;
      Exchange(first_np_or_next_t0,first_np_or_next_t1);
      if (verbosity>5)
         cout << "     Iteration " << kk << " Number of vertices with change " << kch << endl;
   }
   if (verbosity>2 && verbosity < 5) 
      cout << "    Nb of Loops " << kch << endl;
   delete [] first_np_or_next_t0;
   delete [] first_np_or_next_t1;
}


void Geometry::ReadMetric(const char* fmetrix,
                          double      hmin=1.0e-30,
                          double      hmax=1.0e30,
                          double      coef=1)
{
  hmin = Max(hmin,MinimalHmin());
  MeshIstream f_metrix(fmetrix);
  long k,j;
  f_metrix >>  k >> j ;
  if(verbosity>1)
    cout << " -- ReadMetric  " << fmetrix 
	 << ",  coef = " << coef
	 << ", hmin = " << hmin 
	 << ", hmax = " << hmax 
	 << (  (j == 1)? " Iso " : " AnIso " ) << endl;

  if (k != nbv || !(j==1 || j==3)) {
    cerr << " Error Pb metrix " << k << " <> " 
	 <<  nbv << " or  1 or 3  <> " << j << endl;
   MeshError(1003);}
   for (long iv=0; iv<nbv; iv++) {
      double h;
      if (j == 1) {
         f_metrix >> h;
         _vertices[iv].m = Metric(Max(hmin,Min(hmax, h*coef)));
      }
      else if (j==3) {
         double a, b, c;
         f_metrix >> a >> b >> c;
         MetricAnIso M(a,b,c);
         MatVVP2x2 Vp(M/coef);
         Vp.Maxh(hmin);
         Vp.Minh(hmax);
         _vertices[iv].m = Vp;
      }
   }
}


double LengthInterpole(const MetricAnIso Ma,
                       const MetricAnIso Mb,
                       R2                AB)
{
   double k=1./2.;
   int level=0;
   static int kkk=0;
   static  Metric Ms1[32], Ms2[32];
   static double lMs1[32], lMs2[32];
   static double K[32];
   double l=0,sss=0;
   Ms1[level] = Ma;
   Ms2[level] = Mb;
   double sa = Ma(AB), sb = Mb(AB);
   lMs1[level] = sa;
   lMs2[level] = sb;
   K[level] = k;
   level++;
   int i=0;
   double* L= LastMetricInterpole.L, *S = LastMetricInterpole.S;
   double sstop = 0.1; // Max(0.6,(sa+sb)/5000);
   while (level) {
      level--;
      Metric M1=Ms1[level], M2=Ms2[level];
      k = K[level];
      double s1=lMs1[level], s2=lMs2[level];
      double s= (s1+s2)*k;
      if (s>sstop && level<30 && i<500-level) {
         Metric Mi(0.5,M1,0.5,M2);
         double si = Mi(AB);
         if (Abs((s1+s2)-(si+si)) > s1*0.001) {
            k /= 2;
//          We begin by the end to walk in the correct sens from a to b
//          due to the stack 
            Ms1[level] = Mi;
            Ms2[level] = M2;
            lMs1[level] = si;
            lMs2[level] = s2;
            K[level] = k;
            level++;
            Ms1[level] = M1;
            Ms2[level] = Mi;
            lMs1[level] = s1;
            lMs2[level] = si;
            K[level] = k;
            level++;
         }
         else
            L[i]= l += s, S[i] = sss += k, i++;
      }
      else 
         L[i] = l += s, S[i] = sss += k, i++;
   }

// warning for optimisation S is in [0:0.5] not in [0:1]
   assert(i<512);
   LastMetricInterpole.lab = l;
   LastMetricInterpole.opt = i;
   if (i>200 && kkk++<10)
      cout << "Warning LengthInterpole: ( i = " << i << " l = " << l << " sss "
           << sss << " ) " << sstop << endl;
   return l;
}


double abscisseInterpole(const MetricAnIso Ma,
                        const MetricAnIso  Mb,
                        R2                 AB,
                        double             s,
                        int                optim)
{
   if (!optim)
      LengthInterpole(Ma,Mb,AB);
   double l = s* LastMetricInterpole.lab, r;
   int j=LastMetricInterpole.opt-1, i=0, k;
   double *L= LastMetricInterpole.L, *S = LastMetricInterpole.S;

// warning for optimisation S is the abcisse in [0:0.5]
// and L is the length 
   if (l<=L[0])
      r = 2*S[0]*l/L[0];
   else if (l>=L[j])
      r = 1;
   else {
      while (j-i>1) {
         k = (i+j)/2;
         if (l<=L[k])
            j = k;
         else
            i = k;
      }
      if (i==j)
         r = 2*S[i];
      else
         r =  2*(S[i]*(L[j]-l) + S[j]*(l-L[i]))/(L[j]-L[i]);
   }
   assert(r<=1 && r>=0);
   return r;
}

}   // end of namespace bamg 
