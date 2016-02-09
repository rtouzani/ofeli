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

#ifndef _SetOfEdge4_h
#define _SetOfEdge4_h

namespace bamg {

class SetOfEdges4;

class longEdge {

  friend class SetOfEdges4;

  public:
    long i, j, next; 

};


class SetOfEdges4 {

 private:
   long nx, nbax, NbOfEdges, *_head;
   longEdge *_edges;

 public:

  SetOfEdges4(long,
              long);

  ~SetOfEdges4() {
     delete [] _head;
     delete [] _edges;
  }

  long add (long ii,
            long jj);

  long addtrie (long ii,
                long jj)
  { return ii <=jj ? add(ii,jj) : add(jj,ii); }

  long nb() const { return NbOfEdges; }

  long find(long ii,
            long jj) const;

  long findtrie(long ii,
                long jj)
  { return ii <=jj ? find(ii,jj) : find(jj,ii); }

  long i(long k) const { return _edges[k].i; }
  long j(long k) const { return _edges[k].j; }
  long newarete(long k) const { return NbOfEdges == k+1; }
  longEdge& operator[](long k) const { return _edges[k]; }
};

}

#endif 
