#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>

#include "mesh/bamg/Meshio.h"
#include "mesh/bamg/Mesh2.h"
#include "mesh/bamg/QuadTree.h"
#include "mesh/bamg/SetOfE4.h"

namespace bamg {

void Triangles::Write(const char*        filename,
                      const TypeFileMesh typein)
{
   TypeFileMesh type = typein;
   const char *gsuffix=".gmsh";
   int ls=0, lll = strlen(filename);
   if (type==AutoMesh) {
      type = BDMesh;
      if (!strcmp(filename + lll - (ls=4),".bamg"))
         type = mshMesh;
      else if (!strcmp(filename + lll - (ls=4),".BAMG"))
         type = mshMesh;
      else
         ls = 0;
   }
   if (verbosity>1) {
      cout << " -- Writing the file " << filename << " of type " ;

      switch (type) {

         case BDMesh:
              cout << " BD Mesh ";
              break;

         case mshMesh:
              cout << " msh ";
              break;

         default: 
              cerr << endl << " Unknown type mesh file " << int(type)
                   << " for Writing "<< filename << endl;
              MeshError(1);
      }     
      long NbOfTria = nbt - 2*NbOfQuad - NbOutT;
      if (NbOfTria)
         cout << " NbOfTria = " << NbOfTria;
      if (NbOfQuad)
         cout << " NbOfQuad = " << NbOfQuad;
      if (nbe)
         cout << " NbOfRefEdge = " << nbe;
      cout << endl;
   }
   ofstream f(filename /*,ios::trunc*/);
   f.precision(12);
   if (f)
      switch (type) {

         case BDMesh:
            {
               if (!Gh.OnDisk) {
                  delete [] Gh.name;
                  Gh.name = new char[lll+1+strlen(gsuffix)];
                  strcpy(Gh.name,filename);
                  if (Gh.name[lll-ls-1]=='.')
                     strcpy(Gh.name+lll-ls, gsuffix+1);
                  else
                     strcpy(Gh.name+lll-ls,gsuffix);
                  Gh.OnDisk = true;
               }
               f << *this;
               break;
            }

         case mshMesh:
            Write_msh(f);
            break;

	 default: 
            cerr << " Unknown mesh type file " << int(type) << " for writing " << filename <<endl;
	    MeshError(1);
      }
      else {
         cerr << " Error opening file " << filename << endl;
         MeshError(1);
      }
}


void Triangles::Write_msh(ostream& f) const 
{
   long i;
   long *reft = new long[nbt];
   long nbInT = ConsRefTriangle(reft);
   f.precision(12);
   f << nbv << " " << nbInT << " " << nbe << endl;
   for (i=0; i<nbv; i++)
      f << _vertices[i].r.x << " " << _vertices[i].r.y << " " 
        << _vertices[i].ref() << endl;
   for (i=0; i<nbt; i++)
      if (reft[i]>=0)
         f << Number(triangles[i][0])+1 << " " 
           << Number(triangles[i][1])+1 << " " 
           << Number(triangles[i][2])+1 << " " 
           << subdomains[reft[i]].ref << endl;
   for (i=0; i<nbe; i++)
      f << Number(edges[i][0]) +1 << " "  << Number(edges[i][1])+1
        << " " << edges[i].ref << endl;      
   delete [] reft;
}


void Triangles::Write(const char* filename)
{
   ofstream f(filename);
   if (f) {
      if (name)
         delete name;
      name = new char[strlen(filename)+1];
      strcpy(name,filename);
      OnDisk = 1;
      f << *this;
   }
}


void Triangles::WriteElements(ostream& f,
                              long*    reft ,
                              long     nbInT) const
{
   const Triangles& Th=*this;
   if (verbosity>9) 
      cout  << " In Triangles::WriteElements " << endl
            << "   Nb of In triangles " << nbInT-Th.NbOfQuad*2 << endl
            << "   Nb of Quadrilaterals " <<  Th.NbOfQuad << endl
            << "   Nb of in+out+quad  triangles " << Th.nbt << " " << nbInT << endl;	 
   long k=nbInT-Th.NbOfQuad*2, num =0;
   if (k>0) {
      f << "\nTriangles\n"<< k << endl;
      for (long i=0; i<Th.nbt; i++) {
         Triangle &t = Th.triangles[i];
         if (reft[i]>=0 && !(t.Hidden(0) || t.Hidden(1) || t.Hidden(2))) {
            k--;
            f << Th.Number(t[0])+1 << " " << Th.Number(t[1])+1 
              << " "  << Th.Number(t[2])+1  << " " << Th.subdomains[reft[i]].ref << endl;
            reft[i] = ++num;
         }
      }
   } 
   if (Th.NbOfQuad>0) {
      f << "\nQuadrilaterals\n"<<Th.NbOfQuad << endl;
      k = Th.NbOfQuad;
      for (long i=0; i<Th.nbt; i++) {
         Triangle & t = Th.triangles[i];
         Triangle *ta;
         Vertex *v0, *v1, *v2, *v3;
         if (reft[i]<0)
            continue;
         if ((ta=t.Quadrangle(v0,v1,v2,v3)) !=0 && &t<ta) {
            k--;
            f << Th.Number(v0)+1 << " " << Th.Number(v1)+1  << " "  
              << Th.Number(v2)+1 << " "  << Th.Number(v3)+1 << " "  
              << Th.subdomains[reft[i]].ref << endl;
            reft[i] = ++num;
            reft[Number(ta)] = num;
         }
      }
      assert(k==0);
   }
// warning reft is now the element number 
}


ostream& operator <<(ostream&         f,
                     const Triangles& Th) 
{
   long* reft = new long[Th.nbt];
   long nbInT = Th.ConsRefTriangle(reft);
   {
      f << "MeshVersionFormatted 0" <<endl;
      f << "\nDimension\n"  << 2 << endl;
      f << "\nIdentifier\n" ;
      WriteStr(f,Th.identity);
      f << "\n\nGeometry\n" ;
      if (Th.Gh.OnDisk)
         WriteStr(f,Th.Gh.name), f << endl;
      else { // empty file name -> geom in same file
         f << "\"\"" << endl << endl;
         f << "# BEGIN of the include geometry file because geometry is not on the disk"
           << Th.Gh << endl;
         f << "End" << endl
           << "# END of the include geometry file because geometry is not on the disk" << endl;
      }
   }
   { 
      f.precision(12);
      f << "\nVertices\n" << Th.nbv << endl;
      for (long i=0; i<Th.nbv; i++) {
         Vertex& v = Th._vertices[i];
         f << v.r.x << " " << v.r.y << " " << v.ref() << endl;
      }
   }
   long ie; 
   {
      f << "\nEdges\n"<< Th.nbe << endl;
      for (ie=0; ie<Th.nbe; ie++) {
	 Edge& e = Th.edges[ie];
	 f << Th.Number(e[0])+1 << " " << Th.Number(e[1])+1 << " " << e.ref << endl;
      }
      if (Th.NbCrackedEdges) {
         f << "\nCrackedEdges\n"<< Th.NbCrackedEdges << endl;
         for (ie=0; ie<Th.NbCrackedEdges; ie++) {
            Edge& e1 = *Th.CrackedEdges[ie].a.edge;
            Edge& e2 = *Th.CrackedEdges[ie].b.edge;
            f << Th.Number(e1)+1 << " " << Th.Number(e2)+1 << endl;
         }
      }
   }

   Th.WriteElements(f,reft,nbInT);
   {
      f << "\nSubDomainFromMesh\n" << Th.NbSubDomains << endl;
      for (long i=0; i<Th.NbSubDomains; i++)
         f << 3 << " " << reft[Th.Number(Th.subdomains[i].head)] << " " << 1 << " "
           <<  Th.subdomains[i].ref << endl;
   }
   if (Th.Gh.NbSubDomains) {
      f << "\nSubDomainFromGeom\n" << Th.Gh.NbSubDomains << endl;
      for (long i=0; i<Th.NbSubDomains; i++) {  
         f << 2 << " " << Th.Number(Th.subdomains[i].edge)+1 << " " 
           << Th.subdomains[i].sens << " " <<  Th.Gh.subdomains[i].ref << endl;
      }
   }
   {
      f << "\nVertexOnGeometricVertex\n"<< Th.NbVerticesOnGeomVertex << endl;
      for (long i0=0; i0<Th.NbVerticesOnGeomVertex; i0++) {
         VertexOnGeom& v =Th.VerticesOnGeomVertex[i0];
         assert(v.OnGeomVertex());
         f << " " << Th.Number((Vertex *)v)+1  
           << " " << Th.Gh.Number(( GeometricalVertex * )v)+1 << endl;
      }
   }
   { 
      f << "\nVertexOnGeometricEdge\n"<< Th.NbVerticesOnGeomEdge << endl;
      for (long i0=0; i0<Th.NbVerticesOnGeomEdge; i0++) {
         const VertexOnGeom& v = Th.VerticesOnGeomEdge[i0];
         assert(v.OnGeomEdge());
         f << " " << Th.Number((Vertex * )v)+1;
         f << " " << Th.Gh.Number((const GeometricalEdge * )v)+1;
         f << " " << (double)v << endl;
      }
   }
   {
      long i0, k=0;
      for (i0=0; i0<Th.nbe; i0++)
         if (Th.edges[i0].on)
            k++;
      f << "\nEdgeOnGeometricEdge\n"<< k << endl;
      for (i0=0; i0<Th.nbe; i0++)
         if (Th.edges[i0].on) 
            f << (i0+1) << " " << (1+Th.Gh.Number(Th.edges[i0].on)) << endl;
      if (Th.NbCrackedEdges) {
         f << "\nCrackedEdges\n"<< Th.NbCrackedEdges << endl;	  
         for (i0=0; i0<Th.NbCrackedEdges; i0++) {
            f << Th.Number(Th.CrackedEdges[i0].a.edge) << " " ;
            f << Th.Number(Th.CrackedEdges[i0].b.edge) << endl;
         }
      }
   }
   if (&Th.BTh != &Th && Th.BTh.OnDisk && Th.BTh.name) {
      int *mark=new int[Th.nbv];
      long i;
      for (i=0; i<Th.nbv; i++)
         mark[i] = -1;
      f << "\nMeshSupportOfVertices\n" << endl;
      WriteStr(f,Th.BTh.name);
      f << endl << "\nIdentityOfMeshSupport" << endl;
      WriteStr(f,Th.BTh.identity);
      f << endl;
      f << "\nVertexOnSupportVertex" << endl;
      f << Th.NbVertexOnBThVertex << endl;
      for (i=0; i<Th.NbVertexOnBThVertex; i++) {
         const VertexOnVertex & vov = Th.VertexOnBThVertex[i];
         long iv = Th.Number(vov.v);
         mark[iv] = 0;
         f << iv+1 << " " << Th.BTh.Number(vov.bv)+1 << endl;
      }
      f << "\nVertexOnSupportEdge" << endl;
      f << Th.NbVertexOnBThEdge << endl;
      for (i=0; i<Th.NbVertexOnBThEdge; i++) {
         const VertexOnEdge& voe = Th.VertexOnBThEdge[i];
         long iv = Th.Number(voe.v);
         mark[iv] = 1;
         f << iv+1 << " " << Th.BTh.Number(voe.be)+1 << " " << voe.abcisse << endl;
      }
      f << "\nVertexOnSupportTriangle" << endl;
      long k = Th.nbv -  Th.NbVertexOnBThEdge - Th.NbVertexOnBThVertex;
      f << k << endl;
      CurrentTh = &Th.BTh;
      for (i=0; i<Th.nbv; i++)
         if (mark[i] == -1) {
            k--;
            Icoor2 dete[3];
            I2 I = Th.BTh.toI2(Th._vertices[i].r);
            Triangle* tb = Th.BTh.FindTriangleContaining(I,dete);
            if (tb->link) {
               double aa=(double)dete[1]/tb->det, bb=(double)dete[2]/tb->det;
               f << i+1 << " " << Th.BTh.Number(tb)+1 << " " << aa << " " << bb << endl;
            }
            else {
               double aa, bb, det[3];
               TriangleAdjacent ta=CloseBoundaryEdgeV2(I,tb,aa,bb);
               int k = ta;
               det[VerticesOfTriangularEdge[k][1]] = aa;
               det[VerticesOfTriangularEdge[k][0]] = bb;
               det[OppositeVertex[k]] = 1 - aa - bb;
               Triangle* tb = ta;
               f << i+1 << Th.BTh.Number(tb)+1 << " " << det[1] << " " << det[2] << endl;
         }
      }
      assert(!k);
      delete [] mark;
   }
   f << "\nEnd" << endl;
   delete [] reft;
   return f;
}


void Geometry::Write(const char* filename)
{
   ofstream f(filename);
   if (f) {
      if (verbosity>1)
         cout << " -- write geometry in file " << filename << endl;
      if (name)
         delete name;
      name = new char[strlen(filename)+1];
      strcpy(name,filename);
      OnDisk =1;
      f << *this;
   }
}


ostream& operator <<(ostream&        f,
                     const Geometry& Gh) 
{
   long NbCorner=0;
   {
      f << "MeshVersionFormatted 0" << endl;
      f << "\nDimension\n"  << 2 << endl;
   }
   int nbreqv=0;
   { 
      f.precision(12);
      f << "\nVertices\n" << Gh.nbv << endl;
      for (long i=0; i<Gh.nbv; i++) {
         GeometricalVertex& v = Gh._vertices[i];
         if (v.Required())
            nbreqv++;
         f << v.r.x << " " << v.r.y << " " << v.ref() << endl;
         if (v.Corner())
            NbCorner++;
      }
   }
   int nbcracked=0;
   {
      int nbreq=0;
      f << "\nEdges\n"<< Gh.nbe << endl;
      for (long ie=0; ie<Gh.nbe; ie++) {
         GeometricalEdge & e = Gh.edges[ie];
         if (e.Required())
            nbreq++;
         if (e.Cracked()) { 
            long ie1 = Gh.Number(e.link);
            if (ie <= ie1)
               ++nbcracked;
         }
         f << Gh.Number(e[0])+1 << " " << Gh.Number(e[1])+1 << " " << e.ref << endl;
      }
     
      if (nbcracked) {
	 f << "\nCrackedEdges\n"<< nbcracked << endl;
	 for (long ie=0; ie<Gh.nbe; ie++) {
            GeometricalEdge & e = Gh.edges[ie];
            if (e.Cracked()) { 
               long ie1 = Gh.Number(e.link);
               if (ie <= ie1)
                  f << ie+1 << " " << ie1+1 << endl;
            }
         }
      }
      if (nbreq) {
         f << "\nRequiredEdges\n" << nbreq << endl;
         for (long ie=0; ie<Gh.nbe; ie++) {
            GeometricalEdge& e = Gh.edges[ie];
            if (e.Required()) 
	       f << ie+1 << endl;
         }
      }
   }
   f << "\nAngleOfCornerBound\n" << Gh.MaximalAngleOfCorner*180/Pi << endl;
   if (NbCorner) {
      f << "\nCorners\n" << NbCorner << endl;
      for (long i=0,j=0; i<Gh.nbv; i++) {
         GeometricalVertex& v = Gh._vertices[i];
         if (v.Corner()) 
            j++, f << Gh.Number(v)+1 << (j % 5 ? ' ' : '\n');
      }
   }

   if (nbreqv) {
      f << "\nRequiredVertices\n"<< nbreqv << endl;
      for (long j=0,i=0; i<Gh.nbv; i++) {
         GeometricalVertex& v = Gh._vertices[i];
         if (v.Required()) 
         j++,f << i+1 << (j % 5 ? ' ' : '\n');
      }
      f << endl;
   }
    
   { 
      long i;
      f << "\nSubDomainFromGeom\n" << Gh.NbSubDomains << endl;
      for (i=0; i<Gh.NbSubDomains; i++)
         f << "2 " << Gh.Number(Gh.subdomains[i].edge)+1 << " " << Gh.subdomains[i].sens
           << " " << Gh.subdomains[i].ref << endl;
   }
   {
      long n=0, i;
      for (i=0; i<Gh.nbe; i++) {
         if (Gh.edges[i].TgA() && Gh.edges[i][0].Corner())
            n++;
         if (Gh.edges[i].TgB() && Gh.edges[i][1].Corner())
            n++;
      }
      if (n) {
	 f << "TangentAtEdges " << n << endl;
	 for (i=0; i<Gh.nbe; i++) {
            if (Gh.edges[i].TgA() && Gh.edges[i][0].Corner())
               f << i+1 << " 1 " << Gh.edges[i].tg[0].x 
                 << " " << Gh.edges[i].tg[0].y << endl;
            if (Gh.edges[i].TgB() && Gh.edges[i][1].Corner()) 
               f << i+1 << " 2 " << Gh.edges[i].tg[1].x
                 << " " << Gh.edges[i].tg[1].y << endl;
         }
      }
   }
   return f;
}

} // end of namespace bamg 
