// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY:  Bamg: Bidimensional Anisotrope Mesh Generator
// RELEASE: 0 
// USAGE  : You may copy freely these files and use it for
//          teaching or research. These or part of these may
//          not be sold or used for a commercial purpose without
//          our consent: fax (33) 1 39 63 55 14
//
// AUTHOR:   F. Hecht,    
// ORG    :  INRIA
// E-MAIL :   Frederic.Hecht@Inria.fr   
//
// ORIG-DATE:     Dec 97
// Modif          March 98

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <setjmp.h>
#include <new>
#include <cassert>
#include "mesh/bamg/Meshio.h"
#include <iomanip>
#include "mesh/bamg/Mesh2.h"
#include "mesh/bamg/QuadTree.h"

using namespace bamg;
using namespace std;
#include <fstream>

#define NBVMAX 500000

int initgraph=0;
long verbosity=0;

void MeshErrorIO(ios&)
{
   MeshError(999);
   exit(1);
}


int main_bamg(string input_file,
              string output_file)
{
   static char argv[10][80];
//   2 way for uses ---
//   1 first mesh
//   or adaptation
   MeshIstreamErrorHandler = MeshErrorIO;
   long i;
   hinterpole=1;
   int fileout=0, nbvx=NBVMAX, iso=0, AbsError=0, nbjacoby=1, allquad=0;
   int NoMeshReconstruction=0;
   int Rescaling=1;
   double costheta=2;
   double cutoffradian=-1;
   double anisomax = 1e6;
   double err=0.01, errg=0.1, coef=1, hmin=1.e-100, hmax=1.e17, ratio=0, CutOff=1e-5;
   int KeepBackVertices=1;
   double hminaniso=1e-100; 
   const double boundmaxsubdiv=10;
   double maxsubdiv=boundmaxsubdiv;
   double omega=1.8;
   int NbSmooth=3;
   double *solMbb=0, *solMBB=0;
   int *typesolsBB=0;
   long nbsolbb=0, lsolbb=0, nbsolBB=0, lsolBB=0;
   int SplitEdgeWith2Boundary=0;
   int rbbeqMbb=0, rBBeqMBB=0;
   int ChoiseHessien=0;
   double power=1;
   Triangles *Thr=0, *Thb=0;

   char *fgeom=0, *fmeshback=0, *fmeshout=0, *fmeshr=0, *fmetrix=0, *fmsh=0;
   char *rbb=0, *rBB=0, *wbb=0, *wBB=0;
   char *fMbb=0, *foM=0, *fMBB=0;
   verbosity = 0;
   char *datargv[128];
   int datargc=1;
   datargv[0] = datargv[1] = 0;

   int argc = 7;
   strcpy(argv[1],"-g");
   strcpy(argv[2],input_file.c_str());
   strcpy(argv[3],"-o");
   strcpy(argv[4],output_file.c_str());
   strcpy(argv[5],"-v");
   strcpy(argv[6],"0");

   for (i=1; i<argc; i++)
      if (!strcmp(argv[i],"-v") && ++i<argc)
         verbosity = atoi(argv[i]);
      else if (!strcmp(argv[i],"-NoRescaling")) 
         Rescaling = 0;
      else if (!strcmp(argv[i],"-Rescaling")) 
         Rescaling = 1;
      else if (!strcmp(argv[i],"-H") && ++i<argc) 
         ChoiseHessien = atoi(argv[i]);
      else if (!strcmp(argv[i],"-power") && ++i<argc) 
         power = atof(argv[i]);
      else if (!strcmp(argv[i],"-nbv") && ++i<argc)
         nbvx = atoi(argv[i]);
      else if (!strcmp(argv[i],"-nbs") && ++i<argc)
         nbvx = atoi(argv[i]);
      else if (!strcmp(argv[i],"-g") && ++i<argc)
         fgeom = argv[i];
      else if ((!strcmp(argv[i],"-thetaquad") || (!strcmp(argv[i],"-ThetaQuad"))) && ++i<argc)
         costheta = 1 - Abs(cos(Pi*atof(argv[i])/180.0));
      else if (!strcmp(argv[i],"-b")  && ++i<argc)
         fmeshback = argv[i];
      else if (!strcmp(argv[i],"-r")  && ++i<argc)
         fmeshr = argv[i];
      else if (!strcmp(argv[i],"-M")  && ++i<argc) 
         fmetrix = argv[i];
      else if (!strcmp(argv[i],"-rbb") && ++i<argc) 
         rbb = argv[i];
      else if (!strcmp(argv[i],"-rBB") && ++i<argc) 
         rBB = argv[i];
      else if (!strcmp(argv[i],"-Mbb") && ++i<argc) 
         fMbb = argv[i];
      else if (!strcmp(argv[i],"-MBB") && ++i<argc) 
         fMBB = argv[i];
      else if (!strcmp(argv[i],"-wBB") && ++i<argc) 
         wBB = argv[i];
      else if (!strcmp(argv[i],"-wbb") && ++i<argc) 
         wbb = argv[i];
      else if (!strcmp(argv[i],"-o")  && ++i<argc) 
         fmeshout = argv[i];
      else if (!strcmp(argv[i],"-intm")  && i<argc)
         hinterpole = 0; 
      else if (!strcmp(argv[i],"-omsh") && ++i<argc)
         fmsh = argv[i];
      else if (!strcmp(argv[i],"-oM") && ++i<argc)
         foM = argv[i];
      else if (!strcmp(argv[i],"-coef") && ++i<argc) 
         coef = atof(argv[i]);
      else if (!strcmp(argv[i],"-err") && ++i<argc) 
         err = atof(argv[i]);
      else if (!strcmp(argv[i],"-errg") && ++i<argc) 
         errg = atof(argv[i]);
      else if (!strcmp(argv[i],"-ratio") && ++i<argc) 
         ratio = atof(argv[i]);
      else if (!strcmp(argv[i],"-hmin") && ++i<argc) 
         hmin = atof(argv[i]);
      else if (!strcmp(argv[i],"-hminaniso") && ++i<argc) 
         hminaniso = atof(argv[i]);
      else if (!strcmp(argv[i],"-hmax") && ++i<argc) 
         hmax = atof(argv[i]);
      else if (!strcmp(argv[i],"-thetamax")  && ++i<argc) 
         cutoffradian = atof(argv[i])*Pi/180.0;
      else if (!strcmp(argv[i],"-anisomax") && ++i<argc) 
         anisomax = atof(argv[i]);
      else if (!strcmp(argv[i],"-maxsubdiv") && ++i<argc) 
         maxsubdiv = atof(argv[i]);
      else if (!strcmp(argv[i],"-iso"))
         iso = 1;
      else if (!strcmp(argv[i],"-KeepBackVertices"))
         KeepBackVertices = 1;
      else if (!strcmp(argv[i],"-noKeepBackVertices"))
         KeepBackVertices = 0;
      else if (!strcmp(argv[i],"-KeepBackVertices"))
         KeepBackVertices = 1;
      else if (!strcmp(argv[i],"-noKeepBackVertices"))
         KeepBackVertices = 0;
      else if (!strcmp(argv[i],"-2q") || !strcmp(argv[i],"-2Q"))
         allquad = 1;
      else if (!strcmp(argv[i],"-2") )
         allquad = 2;
      else if (!strcmp(argv[i],"-aniso"))
         iso = 0;
      else if (!strcmp(argv[i],"-RelError"))
         AbsError = 0;
      else if (!strcmp(argv[i],"-AbsError"))
         AbsError = 1;
      else if (!strcmp(argv[i],"-splitpbedge") ) 
         SplitEdgeWith2Boundary=1;
      else if (!strcmp(argv[i],"-nosplitpbedge") ) 
         SplitEdgeWith2Boundary=0;
      else if (!strcmp(argv[i],"-NbJacobi") && ++i<argc) 
         nbjacoby = atoi(argv[i]);
      else if (!strcmp(argv[i],"-CutOff") && ++i<argc) 
         CutOff = atof(argv[i]), AbsError = 0;
      else if (!strcmp(argv[i],"-NbSmooth") && ++i<argc) 
         NbSmooth = atoi(argv[i]);
      else if (!strcmp(argv[i],"-omega") && ++i<argc) 
         omega =  atof(argv[i]);
      else {
         cout << " Usage:" << endl;
         cout << "  Mesh INPUT: The 2 arguments are exclusives" << endl;
         cout << "" << endl;
         cout << "     -g  filename    Set geometry for mesh generation. " << endl;
         cout << "                     DB mesh file." << endl;
         cout << "     -b  filename    Set the background mesh for mesh adaption " << endl;
         cout << "                     (require {\tt -M} or {-Mbb} arguments). " << endl;
         cout << "                     - msh file if filename match  *.msh " << endl;
         cout << "                     - otherwise the file is a BD mesh file " << endl;
         cout << "                     Remark: the geometry is the background geometry." << endl;
         cout << "                     DB mesh file." << endl;
         cout << "     -r  filename    Set the  mesh for modification mesh  " << endl;
         cout << "                     no reconstruction " << endl;
         cout << "                     same as in the case of -b argument" << endl;
         cout << "" << endl;
         cout << "     -thetamax (double)   change the angular limit for a corner in degre " << endl;
         cout << "                       the angle is defined from 2 normals of 2 consecutive edges " << endl;
         cout << "                       if no geometry cf reading an  am_fmt .. file " << endl;

         cout << "" << endl;
         cout << "  METRIC definition or mesh size definition, one of the" << endl;
         cout << "  2 next arguments is need in case of  mesh adaption." << endl;
         cout << "" << endl;
         cout << "     -M filename     Set  the metric which is  defined on the background mesh " << endl;
         cout << "                     or on the geometry. Metric file." << endl;
         cout << "     -Mbb filename   Set the solution  defined on the background mesh for" << endl;
         cout << "                     metric construction, the solutions was FE P1 defined on " << endl;
         cout << "                     the background mesh. bb file." << endl;
         cout << "     -MBB filename   same with -Mbb but with BB file " << endl;
         cout << "" << endl;
         cout << "     -errg (double)  Set the level of error on geometry (0.1) by default" << endl;
         cout << "     -coef (double)  Set the value of mutiplicative " << endl;
         cout << "                     coef on the mesh size (1 by default)." << endl;
         cout << "     -power (double) Set the power parameter of hessien to construct " << endl;
         cout << "                     the metric  (1 by default)" << endl;
         cout << "     -maxsubdiv  (double) Change the metric such that the maximal subdivision  " << endl;
         cout << "                     of a background's edge is bounded by the " <<endl;
         cout << "                     given number (always limited by 10) " << endl;
         cout << "     -ratio (double) Set the ratio for a smoothing of the metric." << endl;
         cout << "                     If ratio is  0 then no smoothing" << endl;
         cout << "	              and if ratio  is in  [1.1,10] then the smoothing " << endl;
         cout << "                     change the metrix such that the greatest geometrical" << endl;
         cout << "                     progression (speed of mesh size variation) " << endl;
         cout << "                     in a mesh is bounded  by ratio (by default no smoothing)." << endl;
         cout << "     -hmin (double)  Set the value of the minimal edge size." << endl;
         cout << "     -hminaniso (double)  Set the value of the minimal edge size and save aniso." << endl;
         cout << "     -hmax (double)  Set the value of the maximal edge size." << endl;
         cout << "     -NbSmooth (int) Number of Smoothing iteration " << endl;
         cout << "                     (3 by default if the metric is set otherwise 0)." << endl;
         cout << "     -omega (double)  relaxation parameter for Smoothing " << endl;
         cout << "     -splitpbedge     split in 2 all internal edges with 2 vertex on boundary" << endl;
         cout << "     -nosplitpbedge   d'ont cut internal edges with 2 vertex on boundary (default)" << endl;
	 cout << "" << endl << "" << endl;
         cout << "        the next arguments are used with the -Mbb argument" << endl;
         cout << "" << endl;
         cout << "     -KeepBackVertices  Try to Keep old vertices (default) " << endl;
         cout << "     -noKeepBackVertices  all vertices are create from scratch  " << endl;
         cout << "     -err   (double)    Set the level of error (default 0.01)"    << endl;
         cout << "     -iso               The constructed metric must be isotropic" << endl;
         cout << "     -aniso             The constructed metric can be anisotropic " << endl;
         cout << "     -anisomax (double) Set  maximum ratio  of anisotropy " << endl;
         cout << "                            1 =>  isotropic (1e6 by default) " << endl;
         cout << "     -RelError          Construct metric with relative  error " << endl;
         cout << "     -AbsError          Construct metric with with abs error " << endl;
         cout << "     -CutOff (double)   Set the limit of the relative error  (1.e-5)" << endl;
         cout << "     -NbJacobi (int)    Set the number of Jacobi iteration for smoothing " << endl;
         cout << "                        the construction of metric  (1 by default)." << endl;
         cout << "     -NoRescaling       Don't rescaling of all solution between [0..1] " 
              << "before metric computation " << endl; 
         cout << "     -Rescaling       Do rescaling of all solution between [0..1] " 
              << "before metric computation ((default case)" << endl; 
         cout << "     -H (int)           choices for computing the hessian (test)" << endl;
         cout << "" << endl << "" << endl;
         cout << "  Definition of some internal variable and limitation." << endl;
         cout << "" << endl;
         cout << "     -v   (int)       Set the level of impression (verbosity) " << endl;
         cout << "                      between 0 and 10 (1 by default)." << endl;
         cout << "     -nbv (int)       Set the maximal of  vertices (" << nbvx  << " by default)." << endl;
         cout << "" << endl << "" << endl;
         cout << "  To interpolate a solution form background mesh to generated " << endl;
         cout << "  mesh (in case of adpatation)" << endl;
         cout << "" << endl;
         cout << "    -rbb filename    Read  solution  file defined on the background " << endl;
         cout << "                      mesh for interpolation on created mesh." << endl;
         cout << "                      bb file. (by default the  -Mbb filename)" << endl;
         cout << "     -wbb filename    Write the file of interpolation of the solutions" << endl;
         cout << "                      read with  -rbb argument. " << endl;
         cout << "                      bb file." << endl;
         cout << "     -rBB filename    same -rbb but with BB file ";
         cout << "     -wBB filename    same -wbb but with BB file ";
         cout << "" << endl;
         cout << "  Output mesh file for adpation or generation." << endl;
         cout << "    Remark: one of output mesh file is require " << endl << endl;
         cout << "     -o       filename Create a DB Mesh file." << endl;
         cout << "     -omsh    filename Create a msh file (freefem3 file)." << endl;
         cout << "     -oM filename      Create a metric file. " << endl;
	 cout << endl << endl;
         cout << "     -thetaquad (double)  minimal angle of a quadrangle " << endl;
         cout << "     -2q                  split triangles in 3 quad and quad in 4 quad  " << endl;
         cout << "     -2                   split triangles in 4 trai and quad in 4 quad  " << endl;
         cout << endl;
         exit(3);
      }

// some verification
   NoMeshReconstruction = fmeshr != 0;
   if (!fmeshback)
      fmeshback = fmeshr;
   fileout = fmeshout || wbb || wBB;
   if (!fileout && !foM ) {
      cerr << "No Output file given" << endl;
      MeshError(1);
   }
   if (maxsubdiv > boundmaxsubdiv || maxsubdiv <= 1.0) {
      cerr << " -maxsubdiv " << maxsubdiv << " is not in ]1,"<< boundmaxsubdiv << "]" << endl;
      exit(3);
   }
   if (iso) 
      anisomax = 1;
   if (!(fmetrix||fMbb))
      NbSmooth = 0; 
   if (!rbb) 
      rbb = fMbb;
   if (!rBB) 
      rBB = fMBB;
   if (fMbb && rbb) 
      rbbeqMbb = !strcmp(rbb,fMbb);
   if (fMBB && rBB) 
      rBBeqMBB = !strcmp(rBB,fMBB);

   if (verbosity) {
      if (fgeom && fileout)
         cout << "Construction of mesh from the geometry file " << fgeom << endl;
      else if (fmeshback && fileout) {
         if (NoMeshReconstruction)
            cout << "Modification of adpated mesh " << fmeshback << endl;
         else
            cout << "Construction of adpated mesh from the background mesh " << fmeshback << endl;
      }
      else if (fmeshback && foM)
         cout << "Construction of the metric file on the background mesh " << fmeshback << endl;
   }

   if (fgeom && fileout) {
      if (verbosity) 
         cout << "Construction of mesh from the geometry file " << fgeom << endl;
      Geometry Gh(fgeom);
      hmin = Max(hmin,Gh.MinimalHmin());
      hmax = Min(hmax,Gh.MaximalHmax());
      if (fmetrix) 
         Gh.ReadMetric(fmetrix,hmin,hmax,coef);
      else {
         for (long iv=0; iv<Gh.nbv; iv++) {
            MetricAnIso M = Gh[iv];
            MatVVP2x2 Vp(M/coef);
            Vp.Maxh(hmin);
            Vp.Minh(hmax);
            Gh.vertices[iv].m = Vp;
         }
      }
      Triangles Th(nbvx,Gh);
      if (costheta<=1)
         Th.MakeQuadrangles(costheta);
      if (allquad)
         Th.SplitElement(allquad==2);
      if (SplitEdgeWith2Boundary)
         Th.SplitInternalEdgeWithBorderVertices();
      Th.ReNumberingTheTriangleBySubDomain();
      if (verbosity>3)
         Th.ShowHistogram();
      if (NbSmooth>0)
         Th.SmoothingVertex(NbSmooth,omega);
      if (verbosity>3 && NbSmooth>0)
         Th.ShowHistogram();
      if (fmeshout)
         Th.Write(fmeshout,Triangles::BDMesh);
      if (fmsh)
         Th.Write(fmsh,Triangles::mshMesh);
      if (verbosity>0) {
         cout << "Number of vertices: " << Th.nbv << endl;
         if (Th.nbt-Th.NbOutT-Th.NbOfQuad*2)
            cout << "Number of triangles: " << (Th.nbt-Th.NbOutT-Th.NbOfQuad*2) << endl;
         if (Th.NbOfQuad)
            cout  << "Number of quadrilaterals: " << Th.NbOfQuad << endl;
      }
   }
   else if (fmeshback && (fmetrix||fMbb||fMBB||NoMeshReconstruction) && 
            (fileout || foM || (rbb && wbb) || (rBB && wBB))) {
      Triangles BTh(fmeshback,cutoffradian); 
      hmin = Max(hmin,BTh.MinimalHmin());
      hmax = Min(hmax,BTh.MaximalHmax());
      BTh.MakeQuadTree();
      if (fmetrix) 
         BTh.ReadMetric(fmetrix,hmin,hmax,coef);
      else {
         Metric Mhmax(hmax);
         for (long iv=0; iv<BTh.nbv; iv++)
            BTh[iv].m = Mhmax;
      }
      if (fMbb) {
         solMbb = ReadbbFile(fMbb,nbsolbb,lsolbb,2,2);
         if (lsolbb != BTh.nbv) {
            cerr << "Fatal error  nbsol " << nbsolbb << " " << lsolbb<< " =! " << BTh.nbv << endl;
            cerr << "Size of sol incorrect " << endl;
            MeshError(99);
         }
         assert(lsolbb==BTh.nbv);
         BTh.IntersectConsMetric(solMbb,nbsolbb,0,hmin,hmax,sqrt(err)*coef,1e30,
                                 AbsError?0.0:CutOff,nbjacoby,Rescaling,power,ChoiseHessien);
         if (!rbbeqMbb)
            delete [] solMbb, solMbb = NULL;
      }
      if (fMBB) {
         solMBB = ReadBBFile(fMBB,nbsolBB,lsolBB,typesolsBB,2,2);
         if (lsolBB != BTh.nbv) {
            cerr << "Fatal error  nbsol " << nbsolBB << " " << lsolBB << " =! " << BTh.nbv << endl;
            cerr << "Size of sol incorrect " << endl;
            MeshError(99);
         }
         assert(lsolBB==BTh.nbv);
         BTh.IntersectConsMetric(solMBB,nbsolBB,0,hmin,hmax,sqrt(err)*coef,1e30,
                                 AbsError?0.0:CutOff,nbjacoby,Rescaling,ChoiseHessien);
         if (!rBBeqMBB)
            delete [] solMBB, solMBB = NULL;
      }
      BTh.IntersectGeomMetric(errg,iso);
      if (ratio) 
         BTh.SmoothMetric(ratio);
      BTh.MaxSubDivision(maxsubdiv);

      if (iso) 
         anisomax = 1;
      BTh.BoundAnisotropy(anisomax,hminaniso);
      if (foM) {
         if (verbosity >2)
            cout << " -- write Metric  file " << foM << endl;
         ofstream f(foM);
         if (f)
            BTh.WriteMetric(f,iso);
      }

      if (fileout) {
         if (NoMeshReconstruction) {
            if ((fmeshback==fmeshr) || (!strcmp(fmeshback,fmeshr)))
               Thr = &BTh, Thb = 0; // back and r mesh are the same 
            else
               Thr = new Triangles(fmeshr,cutoffradian), Thb=&BTh;
         }
         Triangles & Th( *(NoMeshReconstruction 
                         ? new Triangles(*Thr,&Thr->Gh,Thb,nbvx)
                         : new Triangles(nbvx,BTh,KeepBackVertices)));
         if (Thr != &BTh)
            delete Thr;
         if (costheta<=1.0)
            Th.MakeQuadrangles(costheta);
         if (allquad)
            Th.SplitElement(allquad==2);
         if (SplitEdgeWith2Boundary)
            Th.SplitInternalEdgeWithBorderVertices();
         Th.ReNumberingTheTriangleBySubDomain();
         if (verbosity>3)
            Th.ShowHistogram();
         if (NbSmooth>0)
            Th.SmoothingVertex(NbSmooth,omega);
         if (verbosity>3 && NbSmooth>0)
            Th.ShowHistogram();
         Th.BTh.ReMakeTriangleContainingTheVertex();
         if (fmeshout)
            Th.Write(fmeshout,Triangles::BDMesh);
         if (fmsh)
            Th.Write(fmsh,Triangles::mshMesh);

         if ((rbb && wbb) ||(rBB && wBB)) {
            if (verbosity>1) {
               if (rbb) 
                  cout << "Interpolation P1 files: " << rbb << " -> " << wbb << endl;
	       if (rBB)
                  cout << "Interpolation P1 files: " << rBB<< " -> " << wBB << endl;
            }
            const int dim = 2;
            double *solbb=0, *solBB=0;
            if (rbb)
               solbb = rbbeqMbb? solMbb : ReadbbFile(rbb,nbsolbb,lsolbb,2,2);
            if (rBB)
               solBB = rBBeqMBB? solMBB : ReadBBFile(rBB,nbsolBB,lsolBB,typesolsBB,2,2);

            if (!solBB && !solbb) {
               cerr << "Fatal Error " << "solBB = " << solBB << " solbb = " << solbb << endl;
               exit(2);
            }
            ofstream *fbb = wbb ? new ofstream(wbb) : 0;
            ofstream *fBB = wBB ? new ofstream(wBB) : 0;
            long nbfieldBB = 0, nbfieldbb = nbsolbb;
            if (fbb)
               *fbb << dim << " " << nbsolbb << " " << Th.nbv << " " << 2 << endl; 
            if (fBB) {
               int i;
               *fBB << dim << " " << nbsolBB;
               for (i=0; i<nbsolBB; i++)
                  *fBB << " " << (typesolsBB ?typesolsBB[i]+1:1) ;
               *fBB << " " << Th.nbv <<" " << 2 << endl; 
               assert(dim==2);
               for (i=0; i<nbsolBB; i++)
                  nbfieldBB += typesolsBB ? typesolsBB[i]+1 : 1;
            }
            cout << "nb of field BB " << nbfieldBB << endl;
            for (i=0; i<Th.nbv; i++) {
               long i0, i1, i2;
               double a[3];
               Icoor2 dete[3];
               I2 I = Th.BTh.toI2(Th.vertices[i].r);
               Triangle & tb = *Th.BTh.FindTriangleContaining(I,dete);
               if (tb.det>0) {
                  a[0] = (dete[0]) / tb.det;
                  a[1] = double(dete[1]) / tb.det;
                  a[2] = double(dete[2]) / tb.det;
                  i0 = Th.BTh.Number(tb[0]);
                  i1 = Th.BTh.Number(tb[1]);
                  i2 = Th.BTh.Number(tb[2]);
               }
               else {
                  double aa, bb;
                  TriangleAdjacent ta = CloseBoundaryEdge(I,&tb,aa,bb).Adj();
                  int k = ta;
                  Triangle &tc = *(Triangle *)ta;
                  i0 = Th.BTh.Number(tc[0]);
                  i1 = Th.BTh.Number(tc[1]);
                  i2 = Th.BTh.Number(tc[2]);
                  a[VerticesOfTriangularEdge[k][1]] = aa;
                  a[VerticesOfTriangularEdge[k][0]] = bb;
                  a[OppositeVertex[k]] = 1 - aa - bb;
               }
               long ibb0=nbfieldbb*i0, ibb1=nbfieldbb*i1, ibb2=nbfieldbb*i2;
               long iBB0=nbfieldBB*i0, iBB1=nbfieldBB*i1, iBB2=nbfieldBB*i2;
               long j=0;
               for (j=0; j<nbfieldbb; j++) 
                  *fbb << " " << (a[0]*solbb[ibb0++] + a[1]*solbb[ibb1++] + a[2]*solbb[ibb2++]);
               for (j=0; j<nbfieldBB; j++)
                  *fBB << " " << (a[0]*solBB[iBB0++] + a[1]*solBB[iBB1++] + a[2]*solBB[iBB2++]);

               if (fbb)
                  *fbb << endl;
               if (fBB)
                  *fBB << endl;
            }
            if (fbb)
               delete fbb;
            if (fBB)
               delete fBB;
            if (solbb)
               delete [] solbb;
            if (solBB) 
               delete [] solBB;
         }
         if (verbosity>0) {
            if (Th.nbt-Th.NbOutT-Th.NbOfQuad*2) 
               cout << " Number of triangles: " << (Th.nbt-Th.NbOutT-Th.NbOfQuad*2);
            if (Th.NbOfQuad) 
               cout  << "Number of quadrilaterals: " << Th.NbOfQuad << endl;
         }
         delete &Th;
      }
   }

   for (i=1; i<datargc; i++)
      delete [] datargv[i] ;
   cout << flush;
   return 0;
}
