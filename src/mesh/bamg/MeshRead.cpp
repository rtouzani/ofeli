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

void Triangles::Read(MeshIstream& f_in,
                     int          Version,
                     double       cutoffradian)
{
   double hmin=HUGE_VAL;
   long i, dim=0, hvertices=0, ifgeom=0;
   Metric M1(1);
   if (verbosity>1)
      cout << " -- ReadMesh " << f_in.CurrentFile << " Version = " << Version << endl;
   int field=0, showfield=0;
   while (f_in.cm()) {  
      field = 0;
      char fieldname[256];
      if (f_in.eof())
         break;
      f_in.cm() >> fieldname;
      if (f_in.eof())
         break;
      f_in.err() ;
      if (verbosity>9)
         cout <<  "    " << fieldname << endl;
      if (!strcmp(fieldname,"MeshVersionFormatted"))
         f_in >> Version;
      else if (!strcmp(fieldname,"End"))
         break;
      else if (!strcmp(fieldname,"Dimension")) {
         f_in >> dim;
         assert(dim ==2);
      }
      else if (!strcmp(fieldname,"Geometry")) {
         char *fgeom ;
         f_in >> fgeom;
         if (strlen(fgeom))
            Gh.ReadGeometry(fgeom);
         else {
            f_in.cm();
            Gh.ReadGeometry(f_in,fgeom);
         }
         Gh.AfterRead();
         ifgeom = 1;
         delete [] fgeom;
      }
      else if (!strcmp(fieldname,"Identifier")) {
         if (identity)
            delete [] identity;
         f_in >> identity;
      }
      else if (!strcmp(fieldname,"hVertices")) {
         hvertices = 1;
         double h;
         for (i=0; i<nbv; i++) {
            f_in >> h;
            _vertices[i].m = Metric(h);
         }
      }
      else if (!strcmp(fieldname,"Vertices")) { 
         assert(dim ==2);
         f_in >> nbv;
         if (verbosity>3)
            cout << "   Nb of Vertices = " << nbv << endl;
         nbvx = nbv;
         _vertices = new Vertex[nbvx];
         assert(_vertices);
         ordre = new Vertex* [nbvx];
         assert(ordre);
         nbiv = nbv;
         for (i=0; i<nbv; i++) {
            f_in >> _vertices[i].r.x >> _vertices[i].r.y >> _vertices[i].ReferenceNumber;
            _vertices[i].DirOfSearch = NoDirOfSearch;
            _vertices[i].m = M1;
            _vertices[i].color = 0;
         }
         nbtx =  2*nbv - 2; // for filling The Holes and quadrilaterals 
         triangles = new Triangle [nbtx];
         assert(triangles);
         nbt = 0;
      }
      else if (!strcmp(fieldname,"Triangles")) { 
         if (dim !=2)
            cerr << "ReadMesh: Dimension <> 2" <<endl, f_in.ShowIoErr(0);
         if (!_vertices || !triangles || !nbv)
            cerr << "ReadMesh:Triangles before Vertices" <<endl,
         f_in.ShowIoErr(0);
         int NbOfTria;
         f_in >> NbOfTria;
         if (verbosity>3)
            cout << "   NbOfTria = " << NbOfTria << endl;
         if (nbt+NbOfTria >= nbtx)
            cerr << "ReadMesh: We must have 2*NbOfQuad + NbOfTria  = "
                 << nbt + NbOfTria<<" < 2*nbv-2 = "  << nbtx << endl,
         f_in.ShowIoErr(0);
         for (i=0; i<NbOfTria; i++) {
            long i1, i2, i3, iref;
            Triangle& t = triangles[nbt++];
            f_in >> i1 >> i2 >> i3 >> iref;
            t = Triangle(this,i1-1,i2-1,i3-1);
            t.color = iref;
         }
      }
      else if (!strcmp(fieldname,"Quadrilaterals")) {
         if (dim !=2)
            cerr << "ReadMesh: Dimension <> 2" <<endl , f_in.ShowIoErr(0);
         if (!_vertices || !triangles || !nbv )
            cerr << "ReadMesh:Quadrilaterals before Vertices" <<endl,
         f_in.ShowIoErr(0);
         f_in >> NbOfQuad;
         if (verbosity>3)
            cout << "   NbOfQuad = " << NbOfQuad << endl;
         if (nbt+2*NbOfQuad >= nbtx)
            cerr << "ReadMesh: We must have 2*NbOfQuad + NbOfTria  = "
                 << nbt + 2*NbOfQuad <<"  < 2*nbv-2 ="  << nbtx << endl,
         f_in.ShowIoErr(0);
         for (i=0; i<NbOfQuad; i++) {
            long i1, i2, i3, i4, iref;
            Triangle &t1 = triangles[nbt++], &t2 = triangles[nbt++];
            f_in >> i1 >> i2 >> i3 >> i4 >> iref;
            t1 = Triangle(this,i1-1,i2-1,i3-1);
            t1.color = iref;
            t2 = Triangle(this,i3-1,i4-1,i1-1);
            t2.color = iref;
            t1.SetHidden(OppositeEdge[1]); // twice because the adj was not created 
            t2.SetHidden(OppositeEdge[1]); 
         }
      }
      else if (!strcmp(fieldname,"VertexOnGeometricVertex")) {
         f_in >> NbVerticesOnGeomVertex;
         if (verbosity>5)
            cout << "   NbVerticesOnGeomVertex = " << NbVerticesOnGeomVertex << endl
                 << " Gh.vertices " << Gh._vertices << endl;
         if (NbVerticesOnGeomVertex) {
            VerticesOnGeomVertex = new VertexOnGeom[NbVerticesOnGeomVertex];
            if (verbosity>7)
               cout << "   VerticesOnGeomVertex = " << VerticesOnGeomVertex << endl
                    << "   Gh.vertices " << Gh._vertices << endl;
            assert(VerticesOnGeomVertex);
	    for (long i0=0; i0<NbVerticesOnGeomVertex; i0++) { 
               long i1, i2;
               f_in >> i1 >> i2;
               VerticesOnGeomVertex[i0] = VertexOnGeom(_vertices[i1-1],Gh._vertices[i2-1]);
            }
         }
      }
      else if (!strcmp(fieldname,"VertexOnGeometricEdge")) { 
         f_in >> NbVerticesOnGeomEdge;
 	 if (verbosity>3)
            cout << "   NbVerticesOnGeomEdge = " << NbVerticesOnGeomEdge << endl;
         if (NbVerticesOnGeomEdge) {
            VerticesOnGeomEdge = new VertexOnGeom[NbVerticesOnGeomEdge] ;
            assert(VerticesOnGeomEdge);
            for (long i0=0; i0<NbVerticesOnGeomEdge; i0++) { 
               long i1, i2;
               double s=0;
               f_in >> i1 >> i2 >> s;
               VerticesOnGeomEdge[i0] = VertexOnGeom(_vertices[i1-1],Gh.edges[i2-1],s);
            }
         }
      }
      else if (!strcmp(fieldname,"Edges")) { 
         long i, j, i1, i2;
         f_in >> nbe;
         edges = new Edge[nbe];
         if (verbosity>5)
            cout << "     Record Edges: Nb of Edge " << nbe << " edges " <<  edges << endl;
         assert(edges);
         double *len = 0;
         if (!hvertices) {
            len = new double [nbv];
            for (i=0; i<nbv; i++)
               len[i] = 0;
         }
         for (i=0; i<nbe; i++) {
            f_in >> i1 >> i2 >> edges[i].ref;
            assert(i1>0 && i2>0);
            assert(i1<=nbv && i2<=nbv);
            i1--; i2--;
            edges[i].v[0] = _vertices + i1;
            edges[i].v[1] = _vertices + i2;
            edges[i].adj[0] = 0;
            edges[i].adj[1] = 0;
            R2 x12 = _vertices[i2].r-_vertices[i1].r;
            double l12=sqrt((x12,x12));        
            if (!hvertices) {
               _vertices[i1].color++;
               _vertices[i2].color++;
               len[i1] += l12;
               len[i2] += l12;
            }
            hmin = Min(hmin,l12);
         }

//       Definition of the default of the given mesh size 
         if (!hvertices) {
            for (i=0; i<nbv; i++) 
               if (_vertices[i].color > 0) 
                  _vertices[i].m = Metric(len[i]/_vertices[i].color);
               else
                  _vertices[i].m = Metric(hmin);
            delete [] len;
         }
         if (verbosity>5)
            cout << "     hmin " << hmin << endl;

//       construction of edges[].adj
         for (i=0; i<nbv; i++) 
            _vertices[i].color = (_vertices[i].color ==2) ? -1 : -2;
         for (i=0; i<nbe; i++) {
            for (j=0; j<2; j++) {
               Vertex *v=edges[i].v[j];
               long i0=v->color, j0;
               if (i0==-1)
                  v->color = i*2 + j;
               else if (i0>=0) {// i and i0 edge are adjacent by the vertex v
                  j0 = i0%2, i0 = i0/2;
                  assert(v==edges[i0 ].v[j0]);
                  edges[i].adj[j] = edges + i0;
                  edges[i0].adj[j0] = edges + i;
                  assert(edges[i0].v[j0]==v);
                  v->color = -3;
               }
            }
         }
      }
      else if (!strcmp(fieldname,"EdgeOnGeometricEdge")) {
         assert(edges);
         int i1, i2, i, j;
         f_in >> i2;
         if (verbosity>3)
            cout << "     Record EdgeOnGeometricEdge: Nb " << i2 <<endl;
         for (i1=0; i1<i2; i1++) {
            f_in >> i >> j;
            if (!(i>0 && j >0 && i <= nbe && j <= Gh.nbe)) {
               cerr << "      Record EdgeOnGeometricEdge i=" << i << " j = " << j;
               cerr << " nbe = " << nbe << " Gh.nbe = " <<  Gh.nbe << endl;
               cerr << " We must have: (i>0 && j >0 && i <= nbe && j <= Gh.nbe) ";
               cerr << " Fatal error in file " << name << " line " << f_in.LineNumber << endl;
               MeshError(1);
            }
            edges[i-1].on = Gh.edges + j - 1;
         }
      }
      else if (!strcmp(fieldname,"SubDomain") || !strcmp(fieldname,"SubDomainFromMesh")) {
         f_in >> NbSubDomains;
         subdomains = new SubDomain [NbSubDomains];
         for (i=0; i<NbSubDomains; i++) {
            long i3, head, sens;
            f_in >> i3 >> head >> sens >> subdomains[i].ref;
            assert(i3==3);
            head--;
            assert(head<nbt && head>=0);
            subdomains[i].head = triangles + head;
         }
      }
      else { // unknown field
         field = ++showfield;
         if (showfield==1) // just to show one time 
            if (verbosity>5)
               cout << "     Warning we skip the field " << fieldname << " at line "
                    << f_in.LineNumber << endl;
      }
      showfield = field; // just to show one time
   }
   if (ifgeom==0) {
      if (verbosity)
         cout << " ## Warning: Mesh without geometry we construct a geometry (theta =" 
              << cutoffradian*180/Pi << " degres )" << endl;
      ConsGeometry(cutoffradian);	
      Gh.AfterRead();
   }
}


void Triangles::Read_msh(MeshIstream& f_in)
{
   if (verbosity>1)
      cout << " -- ReadMesh .msh file " << f_in.CurrentFile << endl;	
   long i;
   f_in.cm() >> nbv >> nbt;
   while (f_in.in.peek()==' ')
      f_in.in.get();
   if (isdigit(f_in.in.peek())) 
      f_in >> nbe;
   if (verbosity>3)
      cout << "    nbv = " << nbv  << " nbt = " << nbt << " nbe = " << nbe << endl;
   nbvx = nbv;
   nbtx =  2*nbv - 2; // for filling The Holes and quadrilaterals 
   triangles = new Triangle [nbtx];
   assert(triangles);
   _vertices = new Vertex[nbvx];
   ordre = new Vertex* [nbvx];
   edges = new Edge[nbe];
   for (i=0; i<nbv; i++)
      f_in >> _vertices[i].r.x >> _vertices[i].r.y >> _vertices[i].ReferenceNumber;
   for (i=0; i<nbt; i++) {
      long i1, i2, i3, r;
      f_in >> i1 >> i2 >> i3 >> r;
      triangles[i] = Triangle(this,i1-1,i2-1,i3-1);
      triangles[i].color = r;
   }
   for (i=0; i<nbe; i++) {
      long i1, i2, r;
      f_in >> i1 >> i2 >> r;
      edges[i].v[0] = _vertices + i1 - 1;
      edges[i].v[1] = _vertices + i2 - 1;
      edges[i].adj[0] = 0;
      edges[i].adj[1] = 0;
      edges[i].ref = r;
   }
}


Triangles::Triangles(const char*  filename,
                           double cutoffradian) 
          : Gh(*(new Geometry())), BTh(*this)
{
   int lll=strlen(filename);
   int msh=!strcmp(filename+lll-4,".msh");

   char *cname = new char[lll+1];
   strcpy(cname,filename);
   long inbvx =0;
   PreInit(inbvx,cname);
   OnDisk = 1;

   MeshIstream f_in(filename);
   if (f_in.IsString("MeshVersionFormatted")) {
      int version;
      f_in >> version;
      Read(f_in,version,cutoffradian);
   }
   else {
      if (msh)
         Read_msh(f_in);
      else {
         cerr << " Unknown mesh type " << filename << endl;
         MeshError(2);
      }
      ConsGeometry(cutoffradian);
      Gh.AfterRead();    
   }
   SetIntCoor();
   FillHoleInMesh();
}


void Geometry::ReadGeometry(const char* filename)
{
   OnDisk = 1;
   if (verbosity>1)
      cout << " -- ReadGeometry " << filename << endl;
   MeshIstream f_in(filename);
   ReadGeometry(f_in,filename);
}


void Geometry::ReadGeometry(MeshIstream& f_in,
                            const char*  filename)
{
   if (verbosity>1)
      cout << " -- ReadGeometry " << filename << endl;
   assert(empty());
   nbiv = nbv = nbvx = 0;
   nbe = nbt = nbtx = 0;
   NbOfCurves = 0;
   name = new char [strlen(filename)+1];
   strcpy(name,filename);
   double Hmin=HUGE_VAL;
   long hvertices=0, i, Version, dim=0;
   int field=0, showfield=0, NbErr=0;

   while (f_in.cm()) { 
      field = 0;
//    Warning: A trick for on allocate fiedname at each time 
      char fieldname[256];
      f_in.cm() >> fieldname;
      f_in.err();
      if (f_in.eof())
         break;
      if (!strcmp(fieldname,"MeshVersionFormatted") )
         f_in >> Version;
      else if (!strcmp(fieldname,"End"))
         break;
      else if (!strcmp(fieldname,"end"))
         break;
      else if (!strcmp(fieldname,"Dimension")) {
         f_in >> dim;
         if (verbosity>5) 
            cout << "     Geom Record Dimension dim = " << dim << endl;        
         assert(dim==2);
      }
      else if (!strcmp(fieldname,"hVertices")) { 
         if (nbv <=0) {
            cerr << "Error: the field Vertex is not found before hVertices " << filename << endl;
            NbErr++;
         }       
         if (verbosity>5) 
            cout << "     Geom Record hVertices nbv=" << nbv << endl;
         hvertices = 1;
         for (i=0; i<nbv; i++) {
            double h;
            f_in >> h;
            _vertices[i].m = Metric(h);
         }
      }
      else if (!strcmp(fieldname,"MetricVertices")) {
         hvertices = 1;
         if (nbv<=0) {
            cerr << "Error: the field Vertex is not found before MetricVertices " << filename << endl;
            NbErr++;
         }
         if (verbosity>5) 
	    cout << "     Geom Record MetricVertices nbv =" << nbv <<  endl;
         for (i=0; i<nbv; i++) {
            double a11, a21, a22;
            f_in >> a11 >> a21 >> a22; 
            _vertices[i].m = Metric(a11,a21,a22);
         }
      }
      else if (!strcmp(fieldname,"h1h2VpVertices")) {
         hvertices = 1;
	 if (nbv<=0) {
	    cerr << "Error: the field Vertex is not found before h1h2VpVertices " << filename << endl;
            NbErr++;
         }
         if (verbosity>5) 
            cout << "     Geom Record h1h2VpVertices nbv=" << nbv << endl;
         for (i=0; i<nbv; i++) {
            double h1, h2, v1, v2;
            f_in >> h1 >> h2 >> v1 >> v2; 
            _vertices[i].m = Metric(MatVVP2x2(1/(h1*h1),1/(h2*h2),D2(v1,v2)));
         }
      }
      else if (!strcmp(fieldname,"Vertices")) {
         assert(dim==2);
         f_in >> nbv;
         nbvx = nbv;
         _vertices = new GeometricalVertex[nbvx];
         if (verbosity>5) 
            cout << "     Geom Record Vertices nbv = " << nbv << "vertices = " 
                 << _vertices << endl;
         assert(nbvx >= nbv);
         nbiv = nbv;
         for (i=0; i<nbv; i++) {
            f_in >> _vertices[i].r.x;
            f_in >> _vertices[i].r.y;
            f_in >> _vertices[i].ReferenceNumber;
            _vertices[i].DirOfSearch = NoDirOfSearch;
            _vertices[i].color = 0;
            _vertices[i].Set();
         }
         pmin = pmax = _vertices[0].r;
         for (i=0; i<nbv; i++) {
            pmin.x = Min(pmin.x,_vertices[i].r.x);
            pmin.y = Min(pmin.y,_vertices[i].r.y);
            pmax.x = Max(pmax.x,_vertices[i].r.x);
            pmax.y = Max(pmax.y,_vertices[i].r.y);
         }
         R2 DD05 = (pmax-pmin)*0.05;
         pmin -= DD05;
         pmax += DD05;
         coefIcoor = (MaxICoor)/(Max(pmax.x-pmin.x,pmax.y-pmin.y));
         assert(coefIcoor >0);
         if (verbosity>5)
            cout << "     Geom: min=" << pmin << "max =" << pmax
                 << " hmin = " << MinimalHmin() <<  endl;
      }
      else if (!strcmp(fieldname,"MaximalAngleOfCorner") ||
               !strcmp(fieldname,"AngleOfCornerBound")) {
         f_in >> MaximalAngleOfCorner;
         if (verbosity>5) 
            cout << "     Geom Record MaximalAngleOfCorner " << MaximalAngleOfCorner << endl;
         MaximalAngleOfCorner *= Pi/180;
      }
      else if (!strcmp(fieldname,"Edges")) {
         if (nbv <=0) {
            cerr << "Error: the field edges is not found before MetricVertices " << filename << endl;
	    NbErr++;
         }
         else {
            int i1, i2;
            R2 zero2(0,0);
            f_in >> nbe;
            edges = new GeometricalEdge [nbe];
            if (verbosity>5) 
               cout << "     Record Edges: Nb of Edge " << nbe <<endl;
            assert(edges);
            assert (nbv >0); 
            double *len = 0;
            if (!hvertices) {
               len = new double [nbv];
               for (i=0; i<nbv; i++)
                  len[i] = 0;
            }
            for (i=0; i<nbe; i++) {
               f_in >> i1 >> i2 >> edges[i].ref;
               i1--; i2--;
               edges[i].v[0] = _vertices + i1;
               edges[i].v[1] = _vertices + i2;
               R2 x12 = _vertices[i2].r-_vertices[i1].r;
               double l12=sqrt((x12,x12));
               edges[i].tg[0] = zero2;
               edges[i].tg[1] = zero2;
               edges[i].SensAdj[0] = edges[i].SensAdj[1] = -1;
               edges[i].Adj[0] = edges[i].Adj[1] = 0;
               edges[i].flag = 0;
               if (!hvertices) {
                  _vertices[i1].color++;
                  _vertices[i2].color++;
                  len[i1] += l12;
                  len[i2] += l12;
               }
               Hmin = Min(Hmin,l12);
            }
//          definition the default of the given mesh size 
            if (!hvertices) {
               for (i=0; i<nbv; i++) 
                  if (_vertices[i].color > 0) 
                     _vertices[i].m = Metric(len[i]/_vertices[i].color);
                  else
                     _vertices[i].m = Metric(Hmin);
               delete [] len;
               if (verbosity>3) 
                  cout << "     Geom Hmin " << Hmin << endl;
            }
         }
      }
      else if (!strcmp(fieldname,"EdgesTangence") || !strcmp(fieldname,"TangentAtEdges")) { 
         int n, i, j, k;
         R2 tg;
         f_in >> n;
         if (verbosity>5) 
            cout << "     Record TangentAtEdges: Nb of Edge " << n <<endl;   
         for (k=0; k<n; k++) {
            f_in >> i >> j;
            f_in >> tg.x  >> tg.y;
            assert(i<=nbe);
            assert(i>0);
            assert(j==1 || j==2);
            i--; j--;
            edges[i].tg[j] = tg;
         }
      }
      else if (!strcmp(fieldname,"Corners")) { 
         int i, j, n;
         f_in >> n;
         if (verbosity>5) 
            cout << "     Record Corner: Nb of Corner " << n <<endl;
         for (i=0; i<n; i++) {
            f_in >> j;
            assert(j<=nbv);
            assert(j>0);
            j--;
            _vertices[j].SetCorner();
            _vertices[j].SetRequired();
         }
      }
      else if (!strcmp(fieldname,"RequiredVertices")) { 
         int i, j, n;
         f_in >> n;
         for (i=0;i<n;i++) {     
            f_in >> j;
            assert(j<=nbv);
            assert(j>0);
            j--;
            _vertices[j].SetRequired();
         }
      }
      else if (!strcmp(fieldname,"RequiredEdges")) { 
         int i, j, n;
         f_in >> n;
         for (i=0; i<n; i++) {     
            f_in >> j;
            assert(j<=nbe);
            assert(j>0);
            j--;
            edges[j].SetRequired(); 
         }
      }
      else if (!strcmp(fieldname,"SubDomain") || !strcmp(fieldname,"SubDomainFromGeom")) { 
         f_in >> NbSubDomains;
         if (NbSubDomains>0) {
            subdomains = new GeometricalSubDomain[NbSubDomains];
            long i0, i1, i2, i3;
            for (i=0; i<NbSubDomains; i++) {
               f_in >> i0 >> i1 >> i2 >> i3;
               assert(i0 == 2);
               assert(i1<=nbe && i1>0);
               subdomains[i].edge=edges + (i1-1);
               subdomains[i].sens = (int) i2;
               subdomains[i].ref = i3;
            }
         }
      }
      else { // unknown field
         field = ++showfield;
         if (showfield==1) 
            if (verbosity>3)
               cout << "    Warning we skip the field " << fieldname << " at line "
                    << f_in.LineNumber << endl;
      }
      showfield = field;
   }
   if (nbv<=0) {
      cerr << "Error: the field Vertex is not found in " << filename << endl;
      NbErr++;
   }
   if (nbe<=0) {
      cerr << "Error: the field Edges is not found in " << filename << endl;
      NbErr++;
   }
   if (NbErr)
      MeshError(1);
}

}  // end of namespace bamg 
