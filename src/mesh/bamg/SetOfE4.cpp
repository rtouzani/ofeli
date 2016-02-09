#include <iostream> 
using namespace std;

#include "mesh/bamg/meshtype.h"
#include "mesh/bamg/SetOfE4.h"

namespace bamg {

SetOfEdges4::SetOfEdges4(long mmx,
                         long nnx)
{
   nx = nnx;
   nbax = mmx;
   NbOfEdges = 0;
   _head = new long [nx];
   long i=nx;
   while (i--)
      _head[i] = -1;
   _edges = new longEdge[nbax];
}

    
long SetOfEdges4::find(long ii,
                       long jj) const
{
   if (_head==NULL) {
      cerr <<"SetOfEdges4::find \nNo more heap head\n";
      MeshError(888);
   }
   long n = _head[Abs(ii) % nx];

   while (n >= 0) 
      if (ii==_edges[n].i && jj==_edges[n].j)
         return n;
      else
         n = _edges[n].next;
   return -1;
}


long SetOfEdges4::add(long ii,
                      long jj)
{
   if (_head==NULL) {
      cerr << "SetOfEdges4::add\nNo more heap head" << endl;
      MeshError(888);
   }

   long h, n=_head[h=Abs(ii) % nx];
   while (n >= 0) { 
      if (ii==_edges[n].i && jj==_edges[n].j)
         return n;
      else
         n = _edges[n].next;
   }
   if (nbax <=NbOfEdges ) {
      cerr << " SetOfEdges4::add\nHeap overflow " << nbax << " " << NbOfEdges << endl;
      MeshError(888);
   }
   _edges[NbOfEdges].i = ii;
   _edges[NbOfEdges].j = jj;
   _edges[NbOfEdges].next = _head[h];
   _head[h] = NbOfEdges;
   return NbOfEdges++;
}

}  // end of namespace bamg 
