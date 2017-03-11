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

#include <limits.h>
#include <string.h>
#include <stdlib.h>

#include "mesh/bamg/Meshio.h"
#include "mesh/bamg/Mesh2.h"
#include "mesh/bamg/QuadTree.h"

namespace bamg {

#define INTER_SEG(a,b,x,y) (((y) > (a)) && ((x) <(b)))
#define ABS(i) ((i)<0 ?-(i) :(i))
#define MAX1(i,j) ((i)>(j) ?(i) :(j))
#define NORM(i1,j1,i2,j2) MAX1(ABS((i1)-(j1)),ABS((i2)-(j2)))

#define IJ(i,j,l) ( (j&l) ? ((i&l) ? 3 : 2) :( (i&l)? 1 : 0 ))
#define I_IJ(k,l)  ((k&1) ? l : 0)
#define J_IJ(k,l)  ((k&2) ? l : 0)


Vertex *QuadTree::NearestVertex(Icoor1 i,
                                Icoor1 j)
{
   QuadTreeBox* pb[MaxDeep];
   int pi[MaxDeep];
   Icoor1 ii[MaxDeep], jj[MaxDeep];
   int l = 0;
   IntQuad h = MaxISize, h0, hb = MaxISize;
   Icoor1 i0 = 0, j0 = 0;
   Icoor1 iplus(i<MaxISize?(i<0?0:i):MaxISize-1);
   Icoor1 jplus(j<MaxISize?(j<0?0:j):MaxISize-1);
   Vertex *vn = NULL;

// init for optimization
   QuadTreeBox *b = _root;
   long n0;
   if (!_root->n)
      return vn;
   while ((n0=b->n)<0) {
      Icoor1 hb2 = hb>>1;
      int k = IJ(iplus,jplus,hb2);  // QuadTreeBox number of size hb2 containing i and j
      QuadTreeBox *b0 = b->b[k];
      if (b0==0 || b0->n==0)
         break;
      NbQuadTreeBoxSearch++;
      b = b0;
      i0 += I_IJ(k,hb2); // i origin of QuadTreeBox
      j0 += J_IJ(k,hb2); // j origin of QuadTreeBox
      hb = hb2;
   }

   for (int k=0; k<n0; k++) {
      I2 i2 = b->v[k]->i;
      h0 = NORM(iplus,i2.x,jplus,i2.y);
      if (h0<h) {
         h = h0;
         vn = b->v[k];
      }
      NbVerticesSearch++;
      return vn;
   }

   pb[0] = b;
   pi[0] = b->n>0 ? int(b->n) : 4;
   ii[0] = i0, jj[0] = j0;
   h = hb;
   do {
      b = pb[l];
      while (pi[l]--) {
         int k = pi[l];
         if (b->n>0) { 
            NbVerticesSearch++;
            I2 i2 = b->v[k]->i;
            h0 = NORM(iplus,i2.x,jplus,i2.y);
            if (h0<h) {
               h = h0;
               vn = b->v[k];
            }
         }
         else {
            QuadTreeBox *b0 = b;
            NbQuadTreeBoxSearch++;
            if ((b=b->b[k])) {
               hb >>= 1;
               Icoor1 iii = ii[l]+I_IJ(k,hb), jjj = jj[l]+J_IJ(k,hb);
               if (INTER_SEG(iii,iii+hb,iplus-h,iplus+h) && INTER_SEG(jjj,jjj+hb,jplus-h,jplus+h)) {
                  pb[++l] = b;
                  pi[l] = b->n>0 ? int(b->n) : 4;
                  ii[l] = iii, jj[l] = jjj;
               }
               else
                  b = b0, hb <<= 1;
            }
            else
               b = b0;
         }
      }
      hb <<= 1;
   } while (l--);
   return vn;
}


Vertex* QuadTree::ToClose(Vertex& v,
                          double  seuil,
                          Icoor1  hx,
                          Icoor1  hy)
{
   const Icoor1 i=v.i.x, j=v.i.y;
   const R2 X(v.r);
   const Metric Mx(v.m);
   QuadTreeBox* pb[MaxDeep];
   int pi[MaxDeep];
   Icoor1 ii[MaxDeep], jj[MaxDeep];
   int l = 0;
   QuadTreeBox *b;
   Icoor1 hb = MaxISize, i0 = 0, j0 = 0;
   if (!_root->n)
      return 0;

// general case -----
   pb[0] = _root;
   pi[0] = _root->n>0 ? int(_root->n) : 4;
   ii[0] = i0, jj[0] = j0;
   do {
      b = pb[l];
      while (pi[l]--) {
         int k = pi[l];
         if (b->n>0) {
            NbVerticesSearch++;
            I2 i2 = b->v[k]->i;
            if (ABS(i-i2.x)<hx && ABS(j-i2.y)<hy) {
               R2 XY(X,b->v[k]->r);
               double dd;
               if ((dd=LengthInterpole(Mx(XY),b->v[k]->m(XY)))<seuil)
                  return b->v[k];
            }
         }
         else {
            QuadTreeBox *b0=b;
            NbQuadTreeBoxSearch++;
            if ((b=b->b[k])) {
               hb >>= 1;
               long iii=ii[l]+I_IJ(k,hb), jjj=jj[l]+J_IJ(k,hb);
               if (INTER_SEG(iii,iii+hb,i-hx,i+hx) && INTER_SEG(jjj,jjj+hb,j-hy,j+hy)) {
                  pb[++l] = b;
                  pi[l] = b->n>0 ? int(b->n) : 4;
                  ii[l] = iii, jj[l] = jjj;
               }
               else
                  b = b0, hb <<= 1;
            }
            else
               b = b0;
         }
      }
      hb <<= 1;
   } while (l--);
   return 0;
}


void QuadTree::Add(Vertex& w)
{
   QuadTreeBox **pb, *b;
   long i = w.i.x, j = w.i.y, l = MaxISize;
   pb = &_root;
   while ((b=*pb) && b->n<0) {
      b->n--;
      l >>= 1;
      pb = &b->b[IJ(i,j,l)];
   }
   if (b) {
      if (b->n>3 && b->v[3]==&w)
         return;
      if (b->n>2 && b->v[2]==&w)
         return;
      if (b->n>1 && b->v[1]==&w)
         return;
      if (b->n>0 && b->v[0]==&w)
         return;
   }
   assert(l);
   while ((b=*pb) && (b->n==4)) {  // the QuadTreeBox is full
      Vertex *v4[4]; // copy of the QuadTreeBox vertices
      v4[0] = b->v[0]; v4[1] = b->v[1];
      v4[2] = b->v[2]; v4[3] = b->v[3];
      b->n = -b->n; // mark is pointer QuadTreeBox
      b->b[0] = b->b[1] = b->b[2] = b->b[3] = 0; // set empty QuadTreeBox ptr
      l >>= 1;    // div the size by 2
      for (int k=0; k<4; k++) { // for the 4 vertices find the sub QuadTreeBox ij
         int ij;
         QuadTreeBox *bb = b->b[ij=IJ(v4[k]->i.x,v4[k]->i.y,l)];
         if (!bb)
            bb = b->b[ij] = NewQuadTreeBox(); // alloc the QuadTreeBox 
         bb->v[bb->n++] = v4[k];
      }
      pb = &b->b[IJ(i,j,l)];
   }
   if (!(b=*pb))
      b = *pb = NewQuadTreeBox(); //  alloc the QuadTreeBox 
   b->v[b->n++] = &w; // we add the vertex 
   NbVertices++;
}


QuadTree::QuadTree(Triangles* t,
                   long       nbv)
         : lenStorageQuadTreeBox(t->nbvx/8+10)
{
   _th = t;
   NbVertices = NbQuadTreeBox = NbQuadTreeBoxSearch = NbVerticesSearch = 0;
   if (nbv==-1)
      nbv = t->nbv;
   _sb = new StorageQuadTreeBox(lenStorageQuadTreeBox);
   _root = NewQuadTreeBox();
   assert(MaxISize>MaxICoor);
   for (long i=0; i<nbv; i++) 
      Add(t->_vertices[i]);
}


QuadTree::QuadTree()
         : lenStorageQuadTreeBox(100),
           NbQuadTreeBoxSearch(0),
           NbVerticesSearch(0)
{
   _th = NULL;
   NbQuadTreeBox = 0;
   NbVertices = 0;
   _sb = new StorageQuadTreeBox(lenStorageQuadTreeBox);
   _root = NewQuadTreeBox();
}


QuadTree::StorageQuadTreeBox::StorageQuadTreeBox(long                ll,
                                                 StorageQuadTreeBox* nn)
{
   len = ll;
   n = nn;
   b = new QuadTreeBox[ll];
   for (int i=0; i<ll; i++)
      b[i].n = 0, b[i].b[0] = b[i].b[1] = b[i].b[2] = b[i].b[3] = 0;
   bc = b;
   be = b + ll;
   assert(b);
}


QuadTree::~QuadTree()
{
   delete _sb;
   _root = NULL;
}


Vertex* QuadTree::NearestVertexWithNormal(Icoor1 i,
                                          Icoor1 j)
{
   QuadTreeBox* b;
   IntQuad h = MaxISize, h0, hb = MaxISize;
   Icoor1 i0 = 0, j0 = 0;
   Vertex *vn = NULL;

// init for optimization
   b = _root;
   if (!_root->n)
      return vn;
   Icoor1 iplus(i<MaxISize?(i<0?0:i):MaxISize-1);
   Icoor1 jplus(j<MaxISize?(j<0?0:j):MaxISize-1);
   long n0;
   while((n0=b->n)<0) {
      Icoor1 hb2 = hb >> 1;
      int k = IJ(iplus,jplus,hb2);// QuadTreeBox number of size hb2 containing i;j
      QuadTreeBox *b0=b->b[k];
      if ((b0==0) || (b0->n==0)) 
         break; // null box or empty   => break 	    
      NbQuadTreeBoxSearch++;
      b = b0;	
      i0 += I_IJ(k,hb2); // i orign of QuadTreeBox
      j0 += J_IJ(k,hb2); // j orign of QuadTreeBox 
      hb = hb2; 
   }

   if (n0>0) {
      for (int k=0; k<n0; k++) {
         I2 i2 = b->v[k]->i;
//       try if it is in the right sense 
         h0 = NORM(iplus,i2.x,jplus,i2.y);
         if (h0<h) {
            h = h0;
            vn = b->v[k];
          }
          NbVerticesSearch++;
      }
      if (vn)
         return vn; 
   }

// general case -----
// INITIALIZATION OF THE HEAP 
   QuadTreeBox *pb[MaxDeep];
   pb[0] = b;
   int pi[MaxDeep];
   pi[0] = b->n>0 ? (int)b->n : 4;
   Icoor1 ii[MaxDeep], jj[MaxDeep];
   ii[0] = i0, jj[0] = j0;
   h = hb;
   int l = 0;
   do {   // walk on the tree  
      b = pb[l];
      while (pi[l]--) { // loop on 4 element of the box
         int k = pi[l];
         if (b->n>0) { // Vertex QuadTreeBox none empty
            NbVerticesSearch++;
            I2 i2=b->v[k]->i;
            h0 = NORM(iplus,i2.x,jplus,i2.y);
            if (h0<h) {
               h = h0;
               vn = b->v[k];
            }
         }
         else {
            QuadTreeBox *b0=b;
            NbQuadTreeBoxSearch++;
            if ((b=b->b[k])) {
               hb >>= 1;
               Icoor1 iii=ii[l]+I_IJ(k,hb), jjj=jj[l]+J_IJ(k,hb);
               if (INTER_SEG(iii,iii+hb,iplus-h,iplus+h) && INTER_SEG(jjj,jjj+hb,jplus-h,jplus+h)) {
                  pb[++l] = b;
                  pi[l] = b->n>0 ? (int)b->n : 4;
                  ii[l] = iii, jj[l] = jjj;
               }
               else
                  b = b0, hb <<= 1;
            }
            else
               b = b0;
         }
      }
      hb <<= 1;
   } while (l--);
   return vn;
}


ostream& operator <<(ostream&        f,
                     const QuadTree& qt)
{
   f << " the quadtree "  << endl;
   f << " NbQuadTreeBox = " << qt.NbQuadTreeBox << " Nb Vertices = " << qt.NbVertices << endl;
   f << " NbQuadTreeBoxSearch " << qt.NbQuadTreeBoxSearch
     << " NbVerticesSearch " << qt.NbVerticesSearch << endl;
   f << " SizeOf QuadTree" << qt.SizeOf() << endl;
   return f;
}

}  // end of namespace bamg
