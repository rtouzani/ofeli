// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY:  Bamg: Bidimensional Anisotrope Mesh Generator
// RELEASE: 0 
// USAGE  : You may copy freely these files and use it for    
//          teaching or research. These or part of these may
//          not be sold or used for a commercial purpose with- 
//          out our consent: fax (33) 1 39 63 55 14       
//
// AUTHOR:   F. Hecht,
// ORG    :  INRIA
// E-MAIL :   Frederic.Hecht@Inria.fr
//
// ORIG-DATE:     Dec 97

#ifndef __MESH2_H
#define __MESH2_H

#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <assert.h>
#if  (defined(unix) || defined(__unix)) && !defined(__AIX)
#define SYSTIMES
#include <sys/times.h>
#include <unistd.h>
#endif

extern long verbosity;
extern int SHOW;
#include "meshtype.h"
#include "R2.h"

namespace bamg {

const double Pi = 3.14159265358979323846264338328;
const float fPi = (float)3.14159265358979323846264338328;

class MeshIstream;
class OFortranUnFormattedFile;
class IFortranUnFormattedFile;

extern int hinterpole;
typedef P2<Icoor1,Icoor2> I2;

inline int BinaryRand() {
#ifdef RAND_MAX
 const long HalfRandMax = RAND_MAX/2;
 return rand() <HalfRandMax;
#else
 return rand() & 16384; // 2^14 (for sun because RAND_MAX is not def in stdlib.h)
#endif
}

typedef P2<double,double> R2;
typedef P2xP2<short,long> I2xI2;
typedef P2<double,double> R2xR2;

}

#include "Metric.h"

namespace bamg {

inline double OppositeAngle(double a) { return a<0 ? Pi + a : a - Pi; }

Icoor2 inline det(const I2& a, const I2& b, const I2& c)
{
   Icoor2 bax=b.x-a.x, bay=b.y-a.y; 
   Icoor2 cax=c.x-a.x, cay=c.y-a.y;
   return bax*cay - bay*cax;
}


// Numbering in a triangle 
static const short VerticesOfTriangularEdge[3][2] = {{1,2},{2,0},{0,1}};
static const short EdgesVertexTriangle[3][2] = {{1,2},{2,0},{0,1}};
static const short OppositeVertex[3] = {0,1,2};
static const short OppositeEdge[3] =  {0,1,2};
static const short NextEdge[3] = {1,2,0};
static const short PreviousEdge[3] = {2,0,1};
static const short NextVertex[3] = {1,2,0};
static const short PreviousVertex[3] = {2,0,1};

long AGoodNumberPrimeWith(long n);

// All angles are in radian beetwen [-Pi,Pi]


class Geometry;
class Triangles;
class Triangle;
class QuadTree;
class GeometricalEdge;
class VertexOnGeom;
class VertexOnEdge;

const int IsVertexOnGeom = 8;
const int IsVertexOnVertex = 16;
const int IsVertexOnEdge = 32;

class Direction {
  private:
  Icoor1 dir;
  public:
  Direction(): dir(MaxICoor){}; //  no direction set
  Direction(Icoor1 i,Icoor1 j) { Icoor2 n2 = 2*(Abs(i)+Abs(j));  
                                 Icoor2 r = MaxICoor* (Icoor2) i;
                                 Icoor1 r1 = (Icoor1) (2*(r/ n2)); // odd number 
                                 dir = (j>0) ? r1 : r1+1; //  odd -> j>0 even -> j<0
                               }
  int sens(Icoor1 i,
           Icoor1 j)
  {
     int r=1;
     if (dir!=MaxICoor) {
        Icoor2 x(dir/2),y1(MaxICoor/2-Abs(x)),y(dir%2?-y1:y1);
        r = (x*i + y*j) >=0;
     }
     return r;
  }                                   
};


class Vertex
{
 public:
  I2 i;  // allow to use i.x, and i.y in long int (beware of scale and centering)
  R2 r;  // allow to use r.x, and r.y in double
  Metric m;
  long ReferenceNumber;
  Direction DirOfSearch;
  union {
    Triangle *t; // one triangle which contains the vertex
    long color;
    Vertex *to;  // use in geometry Vertex to now the mesh Vertex associated 
    VertexOnGeom *on;     // if vint 8; // set with Triangles::SetVertexFieldOn()
    Vertex *onbv; // if vint == 16 on Background vertex Triangles::SetVertexFieldOnBTh()
    VertexOnEdge *onbe;   // if vint == 32 on Background edge
  };
  short vint;  // the vertex number in triangle; varies between 0 and 2 in t
  operator I2() const { return i; } // operator de cast 
  operator const R2 &() const { return r; } 
  double operator()(R2 x) const { return m(x); }
  operator Metric() const { return m; } 
  long Optim(int =1, int =0);
  double Smoothing(Triangles &, const Triangles &, Triangle *&, double=1);
  int ref() const { return ReferenceNumber; }

  friend ostream& operator<<(ostream& s, const Vertex& v)
  { s << "(" << v.i << "," << v.r << MatVVP2x2(v.m) << ")"; return s; }
  inline void Set(const Vertex& rec, const Triangles &, Triangles &);
};

 
double QuadQuality(const Vertex &, const Vertex &, const Vertex &, const Vertex &);


class TriangleAdjacent
{
   friend ostream& operator <<(ostream& f, const TriangleAdjacent &ta)
   {
      f << "{" << ta.t << "," << int(ta.a) << "}";
      return f;
   }

 public:
   Triangle *t;
   int a;
   TriangleAdjacent(Triangle *tt, int aa): t(tt),a(aa &3) { };
   TriangleAdjacent() { };

   operator Triangle* () const { return t; }
   operator Triangle& () const { return *t; }
   operator int() const { return a; }
   TriangleAdjacent& operator++() { a = NextEdge[a]; return *this; }
   TriangleAdjacent operator--() { a = PreviousEdge[a]; return *this; }

  inline TriangleAdjacent Adj() const;
  int swap();
  inline void SetAdj2(const TriangleAdjacent&, int =0);
  inline Vertex* EdgeVertex(const int &) const;
  inline Vertex* OppositeVertex() const;
  inline Icoor2& det() const;
  inline int Locked() const;
  inline int GetAllFlag_UnSwap() const;
  inline void SetLock();
  inline int MarkUnSwap() const;
  inline void SetMarkUnSwap();
  inline void SetCracked();
  inline int Cracked() const;
};


class Edge
{

 public:
   Vertex *v[2];
   long ref;
   GeometricalEdge *on;
   Vertex &operator[](int i) { return *v[i]; }
   Vertex *operator()(int i) { return v[i]; }

   void ReNumbering(Vertex* vb, Vertex* ve, long* renu) 
   {
      if (v[0]>=vb && v[0]<ve)
         v[0] = vb + renu[v[0] - vb];
      if (v[1]>=vb && v[1]<ve)
         v[1] = vb + renu[v[1] - vb];
   }

   const Vertex& operator[](int i) const { return *v[i]; }

// return the point on the curve edge a t in [0:1]
   R2 operator()(double t) const;

// the 2 adj edges if on the same curve
   Edge *adj[2];

   int Intersection(const Edge& e) const
   { 
      if (!(adj[0]==&e || adj[1]==&e)) 
         cerr << "Bug: Intersection " << (void*) &e <<  "  " 
	      << adj[0] << " " <<  adj[1] << endl;
      assert(adj[0]==&e || adj[1]==&e);
      return adj[0]==&e ? 0 : 1;
   }

   double MetricLength() const;

   inline void Set(const Triangles&, long, Triangles&);
};


class GeometricalVertex : public Vertex
{
   int cas;
   friend class Geometry;
   GeometricalVertex *link;  //  link all the same GeometricalVertex circular (Crack) 

 public:

   int Corner() const { return cas&4; }
   int Required() const { return cas&6; }  // a corner is required
   void SetCorner() { cas |= 4; }
   void SetRequired() { cas |= 2; }
   void Set() { cas = 0; }
   GeometricalVertex() : cas(0), link(this) { }
   GeometricalVertex *The() { assert(link); return link; }  // return a unique vertices
   int IsThe() const { return link==this; }

   inline void Set(const GeometricalVertex& rec, const Geometry &Gh, const Geometry &GhNew);

   inline friend ostream& operator <<(ostream& f, const GeometricalVertex &s)
   { f << s.r << "," << s.cas << "."; return f; }
};


class GeometricalEdge
{

 public:

   GeometricalVertex *v[2];
   long ref, CurveNumber;
   R2 tg[2]; // the 2 tangents if tg[0] =0 => no continuity 
   GeometricalEdge *Adj[2]; 
   int SensAdj[2];
   int flag;
   GeometricalEdge *link;
  
   GeometricalVertex & operator[](int i){ return *v[i]; };
   const GeometricalVertex & operator[](int i) const { return *v[i]; };
   GeometricalVertex * operator()(int i){ return v[i]; };
 
  R2 F(double theta) const; // parametrization of the curve edge
  double R1tg(double theta, R2& t) const; // 1/radius of curvature + tangente
  int Cracked() const { return flag & 1; }
  int Dup() const { return flag & 32; }
  int Equi()const { return flag & 2; }
  int ReverseEqui() const { return flag & 128; }
  int TgA()const { return flag & 4; }
  int TgB()const { return flag & 8; }
  int Tg(int i) const { return i==0 ? TgA() : TgB(); }
  int Mark() const { return flag & 16; }
  int Required() { return flag & 64; }
  void SetCracked() { flag |= 1; }
  void SetDup() { flag |= 32; } 
  void SetEqui() { flag |= 2; }
  void SetTgA() { flag|=4; }
  void SetTgB() { flag|=8; }
  void SetMark() { flag|=16; }
  void SetUnMark() { flag &= 1007 /* 1023-16*/; }
  void SetRequired() { flag |= 64; }
  void SetReverseEqui() {flag |= 128; }
  inline void Set(const GeometricalEdge& rec, const Geometry& Th, Geometry& ThNew);
};


class Curve
{
 public:
   GeometricalEdge *be, *ee; // begin et end edge
   int kb, ke;  //  begin vetex and end vertex
   Curve *next; // next curve equi to this
   bool master; // true => of equi curve point on this curve
   inline void Set(const Curve &rec, const Geometry &Th, Geometry &ThNew);
   Curve() : be(0), ee(0), kb(0), ke(0), next(0), master(true) { }
   void Reverse() { Exchange(be,ee); Exchange(kb,ke); } //  revese the sense of the curve
};


class Triangle
{
   friend class TriangleAdjacent;
   friend ostream& operator<<(ostream& f, const Triangle& ta);

 private: // les arete sont opposes a un sommet
   Vertex *ns[3]; // 3 vertices if t is triangle, t[i] allowed by access function, (*t)[i] if pointer
   Triangle *at[3]; // nu triangle adjacent  
   short aa[3];  // les nu des arete dans le triangles (mod 4)

 public:
   Icoor2 det; // determinant du triangle (2 fois l aire des vertices entieres)
   union {
      Triangle *link;
      long color;
   };

   void SetDet() {
     if (ns[0] && ns[1] && ns[2])
        det = bamg::det(*ns[0],*ns[1],*ns[2]);
     else
        det = -1;
   }
 
   Triangle() { }
   Triangle(Triangles *Th, long i, long j, long k);
   Triangle(Vertex *v0, Vertex *v1, Vertex *v2);
   void Set(const Triangle &, const Triangles &, Triangles &);
   int In(Vertex *v) const { return ns[0]==v || ns[1]==v || ns[2]==v; }
   TriangleAdjacent FindBoundaryEdge(int ) const;

   void ReNumbering(Triangle* tb, Triangle* te, long* renu) 
   {
      if (link>=tb && link<te)
         link = tb + renu[link-tb];
      if (at[0]>=tb && at[0]<te)
         at[0] = tb + renu[at[0]-tb];
      if (at[1]>=tb && at[1]<te)
         at[1] = tb + renu[at[1]-tb];
      if (at[2]>=tb && at[2]<te)
         at[2] = tb + renu[at[2]-tb];
   }

   void ReNumbering(Vertex* vb, Vertex* ve, long* renu) 
   {
      if (ns[0]>=vb && ns[0]<ve)
         ns[0] = vb + renu[ns[0]-vb];
      if (ns[1]>=vb && ns[1]<ve)
         ns[1] = vb + renu[ns[1]-vb];
      if (ns[2]>=vb && ns[2]<ve)
         ns[2] = vb + renu[ns[2]-vb];
   }

   const Vertex &operator[](int i) const { return *ns[i]; }
   Vertex & operator[](int i) { return *ns[i]; }
   const Vertex* operator()(int i) const { return ns[i]; }
   Vertex* & operator()(int i) { return ns[i]; }
   TriangleAdjacent Adj(int i) const { return TriangleAdjacent(at[i],aa[i]&3); }
   Triangle* TriangleAdj(int i) const { return at[i&3]; }
   short NuEdgeTriangleAdj(int i) const { return aa[i&3]&3; }
   inline double quality();

   void SetAdjAdj(short a)
   {
      a &= 3;
      Triangle *tt=at[a];
      aa[a] &= 55;
      short aatt=aa[a] & 3;
      if (tt) {
         tt->at[aatt] = this;
         tt->aa[aatt]=a + (aa[a] & 60);
      }
   }

   void SetAdj2(short     a,
                Triangle* t,
                short     aat)
   {
      at[a] = t;
      aa[a] = aat;
      if (t) {
         t->at[aat] = this;
         t->aa[aat] = a;
      }
   }
 
   void SetTriangleContainingTheVertex()
   { 
      if (ns[0])
         ns[0]->t = this, ns[0]->vint = 0;
      if (ns[1])
         ns[1]->t = this, ns[1]->vint = 1;
      if (ns[2])
         ns[2]->t = this, ns[2]->vint = 2;
   }

   int swap(short a1, int=0);
   long Optim(short a, int=0);

   int  Locked(int a) const { return aa[a]&4; } 
   int  Hidden(int a) const { return aa[a]&16; } 
   int  Cracked(int a) const { return aa[a]&32; }

// for optimization
   int  GetAllflag(int a) { return aa[a] & 1020; }
   void SetAllFlag(int a,int f) { aa[a] = (aa[a] &3) + (1020 & f); }

   void SetHidden(int a)
   {
      Triangle *t=at[a];
      if (t)
         t->aa[aa[a] & 3] |= 16;
      aa[a] |= 16;
   }

   void SetCracked(int a)
   {
      Triangle *t=at[a];
      if (t)
         t->aa[aa[a] & 3] |= 32;
      aa[a] |= 32;
   }

   double QualityQuad(int a,int option=1) const;
   Triangle *Quadrangle(Vertex* &v0, Vertex* &v1, Vertex* &v2, Vertex* &v3) const;

   void SetLocked(int a)
   {
      Triangle *t = at[a];
      t->aa[aa[a] & 3] |= 4;
      aa[a] |= 4;
   }

   void SetMarkUnSwap(int a)
   {
      Triangle *t = at[a];
      t->aa[aa[a] & 3] |= 8;
      aa[a] |= 8;
   }

   void SetUnMarkUnSwap(int a)
   {
      Triangle *t = at[a];
      t->aa[aa[a]&3] &= 55; // 23 + 32 
      aa[a] &= 55;
   }
};


class ListofIntersectionTriangles {

  class IntersectionTriangles {

   public:
      Triangle *t;
      double bary[3];  // use if t != 0
      R2 x;
      Metric m;
      double s, p, sn;
  };


  class SegInterpolation {

   public:
      GeometricalEdge *e;
      double sBegin, sEnd; // abscisse of the seg on edge parameter
      double lBegin, lEnd; // length abscisse  set in ListofIntersectionTriangles::Length
      int last;// last index in ListofIntersectionTriangles for this Sub seg of edge

      R2 F(double s) {
         double c01=lEnd-lBegin, c0=(lEnd-s)/c01, c1=(s-lBegin)/c01;
         assert(lBegin<=s && s<=lEnd);
         return e->F(sBegin*c0+sEnd*c1);
      }
  };

  int MaxSize, Size, state, MaxNbSeg, NbSeg;
  double len;
  IntersectionTriangles *lIntTria;
  SegInterpolation *lSegsI;

 public:
  IntersectionTriangles& operator[](int i) { return lIntTria[i]; }
  operator int&() { return Size; }
  ListofIntersectionTriangles(int n=256,
                              int m=16)
  {
     MaxNbSeg = m;
     lIntTria = new IntersectionTriangles [n];
     lSegsI = new SegInterpolation[m];
     state = -1;
     MaxSize = n;
     Size = 0;
     NbSeg = 0; 
     if (verbosity>9) 
        cout << "      construct ListofIntersectionTriangles"
             << MaxSize << " " <<  MaxNbSeg<< endl;};

  ~ListofIntersectionTriangles()
  {
     if (lIntTria)
        delete [] lIntTria, lIntTria=0;
     if (lSegsI)
        delete [] lSegsI, lSegsI=0;
  } 
  
  void init() { state=0; len=0; Size=0; }

  int NewItem(Triangle* tt,
              double    d0,
              double    d1,
              double    d2);

  int NewItem(R2,
              const Metric&);

  void NewSubSeg(GeometricalEdge* e,
                 double           s0,
                 double           s1) 
  {
     if (NbSeg>=MaxNbSeg) {
        int mneo = MaxNbSeg;
        MaxNbSeg *= 2;
       if (verbosity>3) 
          cout <<" reshape lSegsI from " << mneo << " to " 
               << MaxNbSeg <<endl;
       SegInterpolation *lEn = new SegInterpolation[MaxNbSeg];
       assert(lSegsI && NbSeg<MaxNbSeg);
       for (int i=0; i<NbSeg; i++) 
          lEn[i] = lSegsI[MaxNbSeg];
       delete [] lSegsI;
       lSegsI = lEn;        
     }
     if (NbSeg)
        lSegsI[NbSeg-1].last = Size;
     lSegsI[NbSeg].e = e;
     lSegsI[NbSeg].sBegin = s0;
     lSegsI[NbSeg].sEnd = s1;
     NbSeg++;
  }

  void ReShape()
  {
     int newsize = MaxSize*2;
     IntersectionTriangles *nw = new IntersectionTriangles[newsize];
     assert(nw);
     for (int i=0; i<MaxSize; i++)
        nw[i] = lIntTria[i];       
     if (verbosity>3)
        cout << " ListofIntersectionTriangles  ReShape MaxSize "
             << MaxSize << " -> " <<  newsize << endl;
     MaxSize = newsize;
     delete [] lIntTria;
     lIntTria = nw;
  }

  void SplitEdge(const Triangles&,
                 const R2&,
                 const R2&,
                       int         nbegin=0);

  double Length();

  long NewPoints(Vertex* ,
                 long&   nbv,
                 long    nbvx);
};



class GeometricalSubDomain {

 public:
  GeometricalEdge *edge;
  int sens; // -1 or 1
  long ref;
  inline void Set(const GeometricalSubDomain&,
                  const Geometry&            ,
                  const Geometry&             ); 

};


class SubDomain
{
 public:
  Triangle* head;
  long ref;  
  int sens; // -1 or 1
  Edge* edge; 	
  inline void Set(const Triangles &,
                        long,
                        Triangles &);
};


class VertexOnGeom
{
 public:
   double abscisse;  
   Vertex* mv;
   union { 
     GeometricalVertex* gv; // if abscisse <0; 
     GeometricalEdge*   ge;  // if abscisse in [0..1]
   };
   inline void Set(const VertexOnGeom&,const Triangles &, Triangles &);  
   int OnGeomVertex() const { return abscisse<0; }
   int OnGeomEdge() const { return abscisse>=0; }
   VertexOnGeom() { mv = 0; abscisse = 0; gv=0; }
   VertexOnGeom(Vertex& m, GeometricalVertex& g) { mv = &m; abscisse = -1; gv=&g; }
   VertexOnGeom(Vertex& m, GeometricalEdge& g, double s) { mv = &m, abscisse = s; ge = &g; }
   operator Vertex * () const { return mv; }
   operator GeometricalVertex* () const { return gv; }
   operator GeometricalEdge* () const { return ge; }
   operator const double & () const { return abscisse; }
   int IsRequiredVertex() { return ((abscisse<0 ? (gv?gv->Required():0):0)); }
   void SetOn() { mv->on = this; mv->vint = IsVertexOnGeom; }
   friend ostream& operator<<(ostream& f, const VertexOnGeom& vog)
   {
      f << vog.abscisse << " " << vog.mv << " " << vog.gv << " ; ";
      if (vog.abscisse < 0)
         f << *vog.gv << " ;; ";
      return f;
   }

   void Set(const Triangles &,
            long,
            Triangles &);
};


class VertexOnVertex
{
 public:
   Vertex *v, *bv;
   VertexOnVertex(Vertex* w, Vertex* bw) : v(w), bv(bw) {}
   VertexOnVertex() {};
   inline void Set(const Triangles &, long, Triangles &);
   void SetOnBTh() { v->onbv=bv; v->vint = IsVertexOnVertex; }
};


class VertexOnEdge
{
 public:
  Vertex *v;
  Edge *be;
  double abcisse;
  VertexOnEdge(Vertex* w, Edge* bw, double s) : v(w), be(bw), abcisse(s) {}
  VertexOnEdge() { }
  inline void Set(const Triangles &, long, Triangles &);
  void SetOnBTh() { v->onbe = this; v->vint = IsVertexOnEdge; }  
  Vertex & operator[](int i) const { return (*be)[i]; }
  operator double () const { return abcisse; }
  operator  Vertex * () const { return v; } 
  operator  Edge * () const { return be; }
};

inline TriangleAdjacent FindTriangleAdjacent(Edge &E);
inline Vertex* TheVertex(Vertex* a); // remove crak in mesh 


class CrackedEdge { // a small class to store on crack an uncrack information 
  friend class Triangles;
  friend ostream& operator <<(ostream& f, const Triangles& Th) ;  
  class CrackedTriangle {
     friend class Triangles;
     friend class CrackedEdge;
     friend ostream& operator <<(ostream& f, const Triangles& Th) ;  
     Triangle *t; // edge of triangle t
     int i; //  edge number of in triangle
     Edge *edge; // the  2 edge 
     Vertex *New[2]; // new vertex number 
     CrackedTriangle() : t(0),i(0),edge(0) { New[0] = New[1] = 0; }
     CrackedTriangle(Edge * a) : t(0),i(0),edge(a) { New[0] = New[1] = 0; }
     void Crack() {
        Triangle & T(*t); 
        int i0 = VerticesOfTriangularEdge[i][0];
        int i1 = VerticesOfTriangularEdge[i][0];
        assert(New[0] && New[1]);
        T(i0) = New[0];
        T(i1) = New[1];
     }    
     void UnCrack() { 
        Triangle & T(*t); 
       int i0 = VerticesOfTriangularEdge[i][0];
       int i1 = VerticesOfTriangularEdge[i][0];
       assert(New[0] && New[1]);
       T(i0) = TheVertex(T(i0));
       T(i1) = TheVertex(T(i1));}
       void Set() {
          TriangleAdjacent ta(FindTriangleAdjacent(*edge));
          t = ta;
          i = ta;
          New[0] = ta.EdgeVertex(0);
          New[1] = ta.EdgeVertex(1);
       }
  };

 public:
  CrackedTriangle a, b;
  CrackedEdge() : a(), b() {}
  CrackedEdge(Edge *start, long i, long j) : a(start+i), b(start+j) {};
  CrackedEdge(Edge *e0, Edge *e1 ) : a(e0), b(e1) {};

  void Crack() { a.Crack(); b.Crack();}
  void UnCrack() { a.UnCrack(); b.UnCrack();}
  void Set() { a.Set(), b.Set();}
};


class Triangles { 

 public:

  enum TypeFileMesh {
    AutoMesh     = 0,
    BDMesh       = 1,
    mshMesh      = 2
  };

  int static counter; // to kown the number of mesh in memory 
  int OnDisk;       // true if on disk 
  Geometry &Gh;   // Geometry
  Triangles &BTh; // Background Mesh Bth==*this => no background 

  long NbRef; // counter of ref on the this class if 0 we can delete
  long nbvx, nbtx;  // nombre max de sommets, de triangles

  long nt, nbv, nbt, nbiv, nbe; // nb of legal triangles, nb of vertices, of triangles, 
  // of initial vertices, of edges with reference,
  long NbOfQuad; // nb of quadrilaterals 
  long NbSubDomains; // 
  long NbOutT; // Nb of oudeside triangle
  long NbOfTriangleSearchFind;
  long NbOfSwapTriangle;
  char *name, *identity;
  Vertex *_vertices;

  long NbVerticesOnGeomVertex;
  VertexOnGeom* VerticesOnGeomVertex;
  
  long NbVerticesOnGeomEdge;
  VertexOnGeom* VerticesOnGeomEdge;

  long NbVertexOnBThVertex;
  VertexOnVertex *VertexOnBThVertex;

  VertexOnEdge *VertexOnBThEdge;
  long NbVertexOnBThEdge, NbCrackedVertices, NbCrackedEdges;
  CrackedEdge *CrackedEdges;
  R2 pmin, pmax; // extrema
  double coefIcoor;  // coef to integer Icoor1;
  Triangle *triangles;
  Edge *edges; 
  QuadTree *_quadtree;
  Vertex **ordre;
  SubDomain *subdomains;
  ListofIntersectionTriangles lIntTria;
  Triangles(long i);
  ~Triangles(); 
  Triangles(const char *, double=-1) ;

  Triangles(long nbvx,Triangles & BT,int keepBackVertices=1)
         :Gh(BT.Gh),BTh(BT) { GeomToTriangles1(nbvx,keepBackVertices); }
  Triangles(long nbvx,Geometry & G)
         :Gh(G),BTh(*this){GeomToTriangles0(nbvx);}
  Triangles(Triangles &,Geometry * pGh=0,Triangles* pBTh=0,long nbvxx=0 ); // COPY OPERATEUR
  Triangles(const Triangles &, const int *flag, const int *bb); // truncature
  void SetIntCoor(const char *from=0);

  double MinimalHmin() { return 2.0/coefIcoor; }
  double MaximalHmax() { return Max(pmax.x-pmin.x,pmax.y-pmin.y); }
  const Vertex & operator[] (long i) const { return _vertices[i]; }
  Vertex & operator[](long i) { return _vertices[i]; }
  const Triangle & operator()  (long i) const { return triangles[i]; }
  Triangle & operator()(long i) { return triangles[i]; }
  I2 toI2(const R2 & P) const {
          return I2((Icoor1) (coefIcoor*(P.x-pmin.x)),(Icoor1) (coefIcoor*(P.y-pmin.y))); }
  R2 toR2(const I2 & P) const { return R2((double)P.x/coefIcoor+pmin.x, (double)P.y/coefIcoor+pmin.y); }
  void Add(Vertex &s, Triangle *t, Icoor2 *  =0);
  void Insert();
  void ForceBoundary();
  void Heap();
  void FindSubDomain(int );
  long  ConsRefTriangle(long *) const;
  void ShowHistogram() const;

  void ReMakeTriangleContainingTheVertex();
  void UnMarkUnSwapTriangle();
  void SmoothMetric(double raisonmax) ;
  void BoundAnisotropy(double anisomax, double hminaniso=1e-100);
  void MaxSubDivision(double maxsubdiv);
  void WriteMetric(ostream &,int iso) ;
  Edge** MakeGeometricalEdgeToEdge();
  void SetVertexFieldOn();
  void SetVertexFieldOnBTh();
  long SplitInternalEdgeWithBorderVertices();
  void MakeQuadrangles(double costheta);
  int SplitElement(int choice);
  void MakeQuadTree();
  void NewPoints( Triangles &,int KeepBackVertex =1 );
  long InsertNewPoints(long nbvold,long & NbTSwap) ; 
  void NewPointsOld( Triangles & );
  void NewPoints(int KeepBackVertex=1){ NewPoints(*this,KeepBackVertex);}
  void ReNumberingTheTriangleBySubDomain(bool justcompress=false);
  void ReNumberingVertex(long * renu);
  void SmoothingVertex(int =3,double=0.3);
  Metric MetricAt (const R2 &) const;
  GeometricalEdge *ProjectOnCurve(Edge &AB, Vertex &A, Vertex &B, double theta,
                                  Vertex &R, VertexOnEdge &BR, VertexOnGeom &GR);
  
  void WriteElements(ostream &f, long *reft, long nbInT) const;
  
  long Number(const Triangle &t) const { return long(&t - triangles); }
  long Number(const Triangle *t) const { return long(t - triangles); }
  long Number(const Vertex &t) const { return long(&t - _vertices); }
  long Number(const Vertex *t) const { return long(t - _vertices); }
  long Number(const Edge &t) const { return long(&t - edges); }
  long Number(const Edge *t) const { return long(t - edges); }
  long Number2(const Triangle *t) const { return long(t - triangles); }

  Vertex *NearestVertex(Icoor1 i, Icoor1 j);
  Triangle *FindTriangleContaining(const I2 &, Icoor2[3], Triangle *tstart=0) const;
  void Write(const char* filename, const TypeFileMesh type=AutoMesh);
  void Write_msh(ostream &) const;

  void Read(MeshIstream &, int version, double cutoffradian);
  void Read_msh(MeshIstream &);

  void ReadMetric(const char *fmetrix, const double hmin, const double hmax,
                  const double coef);
  void IntersectConsMetric(const double *s, const long nbsol, const int *typsols,
                           const double hmin, const double hmax, const double coef,
                           const double anisomax, const double CutOff=1.e-4,
                           const int NbJacobi=1, const int DoNormalisation=1,
                           const double power=1.0, const int choise=0);
  void IntersectGeomMetric(const double err, const int iso);

  int isCracked() const { return NbCrackedVertices!=0; }
  int Crack();
  int UnCrack();

  friend ostream& operator<<(ostream &f, const Triangles &Th); 
  void Write(const char *filename);
  void ConsGeometry(double =-1.0,int *equiedges=0); // construct a geometry if no geo 
  void FillHoleInMesh() ;
  int CrackMesh();

 private:
  void GeomToTriangles1(long nbvx, int KeepBackVertices=1);// the  real constructor mesh adaption
  void GeomToTriangles0(long nbvx);// the  real constructor mesh generator
  void PreInit(long,char * =0);
};


class Geometry
{ 
 public:
  int OnDisk;
  long NbRef; // counter of ref on the this class if 0 we can delete
  char *name;
  long nbvx,nbtx; // nombre max  de sommets , de  Geometry
  long nbv,nbt,nbiv,nbe; // nombre de sommets, de Geometry, de sommets initiaux,
  long NbSubDomains, NbEquiEdges, NbCrackedEdges, NbOfCurves;
  int empty(){return (nbv ==0) && (nbt==0) && (nbe==0) && (NbSubDomains==0); }
  GeometricalVertex *_vertices;   // data of vertices
  Triangle *triangles; 
  GeometricalEdge *edges;
  QuadTree *_quadtree;
  GeometricalSubDomain *subdomains;
  Curve *curves;
  ~Geometry();
  Geometry(const Geometry &Gh); // Copy  Operator
  Geometry(int nbg, const Geometry **ag); // intersection operator 

  R2 pmin,pmax; // extrema
  double coefIcoor;  // coef to integer Icoor1;
  double MaximalAngleOfCorner;
  
  I2 toI2(const R2 & P) const {
          return I2( (Icoor1)(coefIcoor*(P.x-pmin.x))
                    ,(Icoor1) (coefIcoor*(P.y-pmin.y)) );
  }

  double MinimalHmin() { return 2.0/coefIcoor; }
  double MaximalHmax() { return Max(pmax.x-pmin.x,pmax.y-pmin.y); }
  void ReadGeometry(const char *);
  void ReadGeometry(MeshIstream &, const char *);

  void EmptyGeometry();
  Geometry() {EmptyGeometry();}// empty Geometry
  void AfterRead();
  Geometry(const char * filename) {
      EmptyGeometry();
      OnDisk = 1;
      ReadGeometry(filename);
      AfterRead();
  }

  void ReadMetric(const char *, double hmin, double hmax, double coef);
  const GeometricalVertex & operator[](long i) const { return _vertices[i]; }
  GeometricalVertex & operator[](long i) { return _vertices[i]; }
  const  GeometricalEdge & operator()(long i) const { return edges[i]; }
  GeometricalEdge & operator()(long i) { return edges[i]; }
  long Number(const GeometricalVertex & t) const { return long(&t - _vertices); }
  long Number(const GeometricalVertex * t) const { return long(t - _vertices); }
  long Number(const GeometricalEdge & t) const { return long(&t - edges); }
  long Number(const GeometricalEdge * t) const { return long(t - edges); }
  long Number(const Curve * c) const { return long(c - curves); }

  void UnMarkEdges() {
     for (long i=0;i<nbe;i++)
        edges[i].SetUnMark();
  }

  GeometricalEdge *ProjectOnCurve(const Edge &, double, Vertex &, VertexOnGeom &) const ;
  GeometricalEdge *Containing(const R2 P, GeometricalEdge* start) const;
  friend ostream& operator<<(ostream& f, const Geometry& Gh); 
  void Write(const char *filename);
};

inline Triangles::Triangles(long i) : Gh(*new Geometry()), BTh(*this)
{
   PreInit(i);
}

extern Triangles *CurrentTh;

TriangleAdjacent CloseBoundaryEdge(I2, Triangle *, double &, double &) ;
TriangleAdjacent CloseBoundaryEdgeV2(I2 A, Triangle *t, double &a, double &b);

long FindTriangle(Triangles &Th, double x, double y, double *a, int &inside);

inline Triangle *Triangle::Quadrangle(Vertex* & v0,
                                      Vertex* & v1,
                                      Vertex* & v2,
                                      Vertex* & v3) const
{
// return the other triangle of the quad if a quad or 0 if not a quat
   Triangle* t =0;
   if (link) {
      int a=-1;
      if (aa[0] & 16)
         a = 0;
      if (aa[1] & 16)
         a = 1;
      if (aa[2] & 16)
         a = 2;
      if (a>=0) {
         t = at[a];
         v2 = ns[VerticesOfTriangularEdge[a][0]];
         v0 = ns[VerticesOfTriangularEdge[a][1]];
         v1 = ns[OppositeEdge[a]];
         v3 = t->ns[OppositeEdge[aa[a]&3]];
      }
   }
   return t;
}


inline double Triangle::QualityQuad(int a,
                                    int option) const
{
   double q;
   if (!link || aa[a] &4)
      q = -1;
   else {
      Triangle *t = at[a];
      if (t-this<0)
         q = -1;
      else if (!t->link)
         q = -1;
      else if (aa[0]&16 || aa[1]&16  || aa[2] & 16 || t->aa[0] & 16 || 
               t->aa[1] & 16 || t->aa[2] & 16)
         q = -1;
      else if(option) {
         const Vertex &v2 = *ns[VerticesOfTriangularEdge[a][0]];
         const Vertex &v0 = *ns[VerticesOfTriangularEdge[a][1]];
         const Vertex &v1 = *ns[OppositeEdge[a]];
         const Vertex &v3 = * t->ns[OppositeEdge[aa[a]&3]];
         q = QuadQuality(v0,v1,v2,v3);
      }
      else
         q = 1;
   }
   return q;
}


inline void Vertex::Set(const Vertex& rec,
                        const Triangles&,
                        Triangles&)
{ 
   *this  = rec;
}


inline void GeometricalVertex::Set(const GeometricalVertex& rec,
                                   const Geometry&,
                                   const Geometry&)
{
   *this  = rec;
}


inline void Edge::Set(const Triangles& Th,
                      long             i,
                      Triangles&       ThNew)
{
   *this = Th.edges[i];
   v[0] = ThNew._vertices + Th.Number(v[0]);    
   v[1] = ThNew._vertices + Th.Number(v[1]);
   if (on) 
      on = ThNew.Gh.edges + Th.Gh.Number(on);
   if (adj[0])
      adj[0] = ThNew.edges + Th.Number(adj[0]);
   if (adj[1])
      adj[1] = ThNew.edges + Th.Number(adj[1]);
}


inline void GeometricalEdge::Set(const GeometricalEdge& rec,
                                 const Geometry&        Gh,
                                 Geometry&              GhNew)
{
   *this = rec;
   v[0] = GhNew._vertices + Gh.Number(v[0]);    
   v[1] = GhNew._vertices + Gh.Number(v[1]); 
   if (Adj[0])
      Adj[0] = GhNew.edges + Gh.Number(Adj[0]);     
   if (Adj[1])
      Adj[1] = GhNew.edges + Gh.Number(Adj[1]);     
}

 
inline void Curve::Set(const Curve&    rec,
                       const Geometry& Gh,
                       Geometry&       GhNew)
{
   *this = rec;
   be = GhNew.edges + Gh.Number(be);    
   ee = GhNew.edges + Gh.Number(ee); 
   if (next)
      next = GhNew.curves + Gh.Number(next); 
}


inline void Triangle::Set(const Triangle&  rec,
                          const Triangles& Th,
                          Triangles&       ThNew)
{
   *this = rec;
   if (ns[0])
      ns[0] = ThNew._vertices + Th.Number(ns[0]);
   if (ns[1])
      ns[1] = ThNew._vertices + Th.Number(ns[1]);
   if (ns[2])
      ns[2] = ThNew._vertices + Th.Number(ns[2]);
   if (at[0])
      at[0] = ThNew.triangles + Th.Number(at[0]);
   if (at[1])
      at[1] = ThNew.triangles + Th.Number(at[1]);
   if (at[2])
      at[2] = ThNew.triangles + Th.Number(at[2]);
   if (link>=Th.triangles && link<Th.triangles+Th.nbt)
      link = ThNew.triangles + Th.Number(link);
}


inline void VertexOnVertex::Set(const Triangles& Th,
                                long             i,
                                Triangles&       ThNew)
{
   *this = Th.VertexOnBThVertex[i];  
   v = ThNew._vertices + Th.Number(v);
}


inline void SubDomain::Set(const Triangles& Th,
                           long             i,
                           Triangles&       ThNew)
{
   *this = Th.subdomains[i];
   assert(head-Th.triangles>=0 && head-Th.triangles<Th.nbt);
   head = ThNew.triangles + Th.Number(head);
   assert(edge-Th.edges>=0 && edge-Th.edges<Th.nbe); 
   edge = ThNew.edges + Th.Number(edge);
}


inline void GeometricalSubDomain::Set(const GeometricalSubDomain& rec,
                                      const Geometry&             Gh,
                                      const Geometry&             GhNew)
{
   *this = rec;
   edge = Gh.Number(edge) + GhNew.edges;
}


inline void VertexOnEdge::Set(const Triangles& Th,
                              long             i,
                              Triangles&       ThNew)
{
  *this = Th.VertexOnBThEdge[i];
  v = ThNew._vertices + Th.Number(v);
}


inline void VertexOnGeom::Set(const VertexOnGeom& rec,
                              const Triangles&    Th,
                              Triangles&          ThNew)
{
   *this = rec;
   mv = ThNew._vertices + Th.Number(mv);
   if (gv) {
      if (abscisse<0)
         gv = ThNew.Gh._vertices + Th.Gh.Number(gv);
      else
         ge = ThNew.Gh.edges + Th.Gh.Number(ge);
   }
}


inline double Edge::MetricLength() const
{ 
   return LengthInterpole(v[0]->m,v[1]->m,v[1]->r-v[0]->r);
}


inline void Triangles::ReMakeTriangleContainingTheVertex()
{
   for (long i=0; i<nbv; i++) {
      _vertices[i].vint = 0;
      _vertices[i].t = 0;
   }
   for (long i=0; i<nbt; i++) 
      triangles[i].SetTriangleContainingTheVertex();
}


inline void Triangles::UnMarkUnSwapTriangle()
{
   for (long i=0; i<nbt; i++) 
      for (int j=0; j<3; j++)
         triangles[i].SetUnMarkUnSwap(j);
}


inline void Triangles::SetVertexFieldOn()
{
   for (long i=0; i<nbv; i++) 
      _vertices[i].on = 0;
   for (long j=0; j<NbVerticesOnGeomVertex; j++) 
      VerticesOnGeomVertex[j].SetOn();
   for (long k=0; k<NbVerticesOnGeomEdge; k++)
      VerticesOnGeomEdge[k].SetOn();
}


inline void Triangles::SetVertexFieldOnBTh()
{
   for (long i=0; i<nbv; i++)
      _vertices[i].on = 0;
   for (long j=0; j<NbVertexOnBThVertex; j++ ) 
      VertexOnBThVertex[j].SetOnBTh();
   for (long k=0; k<NbVertexOnBThEdge; k++) 
      VertexOnBThEdge[k].SetOnBTh();      
}


inline void TriangleAdjacent::SetAdj2(const TriangleAdjacent& ta,
                                      int                     l)
{
   if (t) {
      t->at[a] = ta.t;
      t->aa[a]=ta.a|l;
   }
   if (ta.t) {
      ta.t->at[ta.a] = t;
      ta.t->aa[ta.a] = a|l;
   }
}


inline int TriangleAdjacent::Locked() const
{
   return t->aa[a] &4;
}


inline int TriangleAdjacent::Cracked() const
{
   return t->aa[a] &32;
}


inline int TriangleAdjacent::GetAllFlag_UnSwap() const
{
   return t->aa[a] & 1012;
}


inline int TriangleAdjacent::MarkUnSwap() const { return t->aa[a] & 8; }


inline void TriangleAdjacent::SetLock() { t->SetLocked(a); }


inline void TriangleAdjacent::SetCracked() { t->SetCracked(a); }


inline TriangleAdjacent TriangleAdjacent::Adj() const { return t->Adj(a); }


inline Vertex *TriangleAdjacent::EdgeVertex(const int& i) const { return t->ns[VerticesOfTriangularEdge[a][i]]; }


inline Vertex *TriangleAdjacent::OppositeVertex() const { return t->ns[bamg::OppositeVertex[a]]; }


inline Icoor2 &TriangleAdjacent::det() const { return t->det; }


inline TriangleAdjacent Adj(const TriangleAdjacent& a) { return a.Adj(); }


inline TriangleAdjacent Next(const TriangleAdjacent& ta) { return TriangleAdjacent(ta.t,NextEdge[ta.a]); }


inline TriangleAdjacent Previous(const TriangleAdjacent& ta) { return TriangleAdjacent(ta.t,PreviousEdge[ta.a]); }


inline void Adj(GeometricalEdge* & on,
                int&               i)
{
   int j=i;
   i = on->SensAdj[i];
   on = on->Adj[j];
}


inline double quality(const Vertex& va,
                      const Vertex& vb,
                      const Vertex& vc)
{
   double ret;
   I2 ia=va, ib=vb, ic=vc;
   I2 ab=ib-ia, bc=ic-ib, ac=ic-ia;
   Icoor2 deta=Det(ab,ac);
   if (deta<=0)
      ret = -1;
   else {
      double a = sqrt((ac,ac)),
             b = sqrt((bc,bc)),
             c = sqrt((ab,ab));
      double h=Max(Max(a,b),c), ro=deta/(a+b+c);
      ret = ro/h;
   }
   return ret;
}


inline Triangle::Triangle(Triangles* Th,
                          long       i,
                          long       j,
                          long       k)
{
   Vertex *v = Th->_vertices;
   assert(i>=0 && j>=0 && k>=0);
   ns[0] = v + i;
   ns[1] = v + j;
   ns[2] = v + k;
   at[0] = at[1] = at[2] = 0;
   aa[0] = aa[1] = aa[2] = 0;
   det = 0;
}


inline Triangle::Triangle(Vertex* v0,
                          Vertex* v1,
                          Vertex* v2)
{
   ns[0] = v0;
   ns[1] = v1;
   ns[2] = v2;
   at[0] = at[1] = at[2] = 0;
   aa[0] = aa[1] = aa[2] = 0;
   if (v0)
      det = 0;
   else {
      det = -1;
      link = NULL;
   }
}


inline double Triangle::quality()
{
   return det < 0 ? -1 : bamg::quality(*ns[0],*ns[1],*ns[2]);
}


long inline Vertex::Optim(int i,
                          int koption)
{
   long ret=0;
   if (t && vint>=0 && vint<3) {
      ret = t->Optim(vint,koption);
      if (!i) {
         t = 0;
         vint = 0;
      }
   }
   return ret;
}


Icoor2 inline det(const Vertex& a,
                  const Vertex& b,
                  const Vertex& c)
{
   Icoor2 bax=b.i.x-a.i.x, bay=b.i.y-a.i.y; 
   Icoor2 cax=c.i.x-a.i.x, cay=c.i.y-a.i.y; 
   return bax*cay - bay*cax;
}


void swap(Triangle* t1,
          short     a1,
          Triangle* t2,
          short     a2,
          Vertex*   s1,
          Vertex*   s2,
          Icoor2    det1,
          Icoor2    det2);


int inline TriangleAdjacent::swap()
{
   return t->swap(a);
}


int SwapForForcingEdge(Vertex* &         pva,
                       Vertex* &         pvb,
                       TriangleAdjacent& tt1,
                       Icoor2&           dets1,
                       Icoor2&           detsa,
                       Icoor2&           detsb,
                       int&              nbswap);

int ForceEdge(Vertex&           a,
              Vertex&           b,
              TriangleAdjacent& taret);


inline TriangleAdjacent FindTriangleAdjacent(Edge &E)
{
   Vertex *a=E.v[0], *b=E.v[1];
   Triangle* t=a->t;
   int i=a->vint;
   TriangleAdjacent ta(t,EdgesVertexTriangle[i][0]); // Previous edge
   assert(t && i>=0 && i<3);
   assert(a==(*t)(i));
   do { // turn around vertex in direct sens (trigo)
      //  in no crack => ta.EdgeVertex(1) == a otherwise ??? 
      if (ta.EdgeVertex(1)==a && ta.EdgeVertex(0)==b)
         return ta; 
      ta = ta.Adj();
      if (ta.EdgeVertex(0)==a && ta.EdgeVertex(1)==b)
         return ta; 
      --ta;
   } while (t != (Triangle *)ta);
   assert(0);
   return TriangleAdjacent(0,0);// error 
}
  

inline Vertex *TheVertex(Vertex* a) // give a unique vertex with smallest number
{ // in case on crack in mesh 
    Vertex *r(a), *rr;
    Triangle *t=a->t;
    int i=a->vint;
    TriangleAdjacent ta(t,EdgesVertexTriangle[i][0]); // Previous edge
    assert(t && i>=0 && i<3);
    assert(a==(*t)(i));
    do { // turn around vertex in direct sens (trigo)
//     if no crack => ta.EdgeVertex(1) == a
       if ((rr=ta.EdgeVertex(0)) < r)
          r = rr;
       ta = ta.Adj();
       if ((rr=ta.EdgeVertex(1)) < r)
          r = rr;
       --ta;
     } while (t != (Triangle*) ta);  
     return r;
}


inline double CPUtime()
{
#ifdef SYSTIMES
   struct tms buf;
   if (times(&buf)!=-1)
      return (double(buf.tms_utime)+double(buf.tms_stime))/long(sysconf(_SC_CLK_TCK));
   else
#endif
      return (double(clock()))/CLOCKS_PER_SEC;
}

}

#endif
