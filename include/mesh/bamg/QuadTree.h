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

#ifndef __QUADTREE_H
#define __QUADTREE_H

namespace bamg {

const int MaxDeep = 30;
typedef  long IntQuad;
const IntQuad MaxISize = (1L << MaxDeep);

class Triangles;
class Vertex;

class QuadTree
{

 public:

   class QuadTreeBox {

     public:

        long n; // if n < 4 => Vertex else =>  QuadTreeBox;
        union {
           QuadTreeBox *b[4];
           Vertex *v[4];
        };
   };

   class StorageQuadTreeBox {

     public:
        QuadTreeBox *b, *bc, *be;
        long len;
        StorageQuadTreeBox *n; // next StorageQuadTreeBox
        StorageQuadTreeBox(long, StorageQuadTreeBox* =0);
        ~StorageQuadTreeBox()
        {
           if (n)
              delete n;
           delete [] b;
        }
        long SizeOf() const { return len*sizeof(QuadTreeBox)+sizeof(StorageQuadTreeBox)+ (n?n->SizeOf():0); }
   };

   long lenStorageQuadTreeBox, NbQuadTreeBoxSearch, NbVerticesSearch, NbQuadTreeBox, NbVertices;
   Vertex* NearestVertex(Icoor1 i, Icoor1 j);
   Vertex* NearestVertexWithNormal(Icoor1 i, Icoor1 j); 
   Vertex* ToClose(Vertex &, double, Icoor1, Icoor1);
   long SizeOf() const { return sizeof(QuadTree)+_sb->SizeOf(); }
   void Add(Vertex& w);

   QuadTreeBox* NewQuadTreeBox()
   {
      if (!(_sb->bc<_sb->be))
         _sb = new StorageQuadTreeBox(lenStorageQuadTreeBox,_sb);
      assert(_sb && (_sb->bc->n == 0));
      NbQuadTreeBox++;
      return _sb->bc++;
   }
   ~QuadTree();
   QuadTree(Triangles* t, long nbv=-1);
   QuadTree();
   friend ostream& operator <<(ostream& f, const QuadTree& qt);

 private:
   StorageQuadTreeBox* _sb;
   QuadTreeBox *_root;
   Triangles *_th;
};

}

#endif
