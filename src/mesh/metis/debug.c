/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * debug.c
 *
 * This file contains code that performs self debuging
 *
 * Started 7/24/97
 * George
 *
 * $Id: debug.c,v 1.1 1998/11/27 17:59:13 karypis Exp $
 *
 */

#include "mesh/metis/metis.h"

/*************************************************************************
* This function computes the cut given the graph and a where vector
**************************************************************************/
int ComputeCut(GraphType *graph, idxtype *where)
{
   int i, j, cut;
   if (graph->adjwgt == NULL) {
      for (cut=0, i=0; i<graph->nvtxs; i++) {
         for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++)
            if (where[i] != where[graph->adjncy[j]])
               cut++;
      }
   }
   else {
      for (cut=0, i=0; i<graph->nvtxs; i++) {
         for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++)
            if (where[i] != where[graph->adjncy[j]])
               cut += graph->adjwgt[j];
      }
   }
   return cut/2;
}


/*************************************************************************
* This function checks whether or not the boundary information is correct
**************************************************************************/
int CheckBnd(GraphType *graph) 
{
   int i, j, nvtxs;
   idxtype *xadj, *adjncy, *where, *bndptr, *bndind=0;
   nvtxs = graph->nvtxs;
   xadj = graph->xadj;
   adjncy = graph->adjncy;
   where = graph->where;
   bndptr = graph->bndptr;
   bndind = graph->bndind;
   for (i=0; i<nvtxs; i++) {
      if (xadj[i+1]-xadj[i] == 0)

      for (j=xadj[i]; j<xadj[i+1]; j++) {
         if (where[i] != where[adjncy[j]]) {
            if (bndptr[i]==-1)
               printf("Warning in ChekBnd: bndptr[%d] is equal to -1\n",i);
            if (bndind[bndptr[i]]!=i)
               printf("Warning in ChekBnd: bndind[bndptr[%d]] is different from %d\n",i,i);
            break;
         }
      }
   }
   return 1;
}


/*************************************************************************
* This function checks whether or not the boundary information is correct
**************************************************************************/
int CheckBnd2(GraphType *graph) 
{
   int j, nvtxs, id, ed;
   idxtype *xadj, *adjncy, *where, *bndptr, *bndind;

   nvtxs = graph->nvtxs;
   xadj = graph->xadj;
   adjncy = graph->adjncy;
   where = graph->where;
   bndptr = graph->bndptr;
   bndind = graph->bndind;

   for (int i=0; i<nvtxs; i++) {
      id = ed = 0;
      for (j=xadj[i]; j<xadj[i+1]; j++) {
         if (where[i] != where[adjncy[j]]) 
            ed += graph->adjwgt[j];
         else
            id += graph->adjwgt[j];
      }
      if (ed - id >= 0 && xadj[i] < xadj[i+1]) {
         if (bndptr[i]==-1)
            printf("Warning in ChekBnd2: bndptr[%d] is equal to -1\n",i);
         if (bndind[bndptr[i]]!=i)
            printf("Warning in ChekBnd2: bndind[bndptr[%d]] is different from %d\n",i,i);
      }
   }
   ASSERTP(nbnd == graph->nbnd, ("%d %d\n", nbnd, graph->nbnd));
   return 1;
}

/*************************************************************************
* This function checks whether or not the boundary information is correct
**************************************************************************/
int CheckNodeBnd(GraphType *graph, int onbnd) 
{
  int i, nvtxs, nbnd;
  idxtype *where, *bndptr, *bndind;

  nvtxs = graph->nvtxs;
  where = graph->where;
  bndptr = graph->bndptr;
  bndind = graph->bndind;
  for (nbnd=0, i=0; i<nvtxs; i++) {
    if (where[i] == 2) 
      nbnd++;
  }

  if (nbnd != onbnd)
      printf("Warning: nbnd %d is different from onbnd\n",i);
  for (i=0; i<nvtxs; i++) {
    if (where[i] != 2) {
       if (bndptr[i]!=-1)
          printf("Warning: bndptr[i] %d is equal -1\n",i);
    }
    else {
       if (bndind[bndptr[i]]!=i)
          printf("Warning: bndind array is erroneous\n");
    }
  }
  return 1;
}



/*************************************************************************
* This function checks whether or not the rinfo of a vertex is consistent
**************************************************************************/
int CheckRInfo(RInfoType *rinfo)
{
  int i, j;

  for (i=0; i<rinfo->ndegrees; i++) {
     for (j=i+1; j<rinfo->ndegrees; j++) {
       if (rinfo->edegrees[i].pid == rinfo->edegrees[j].pid)
          printf("Warning: rinfo->edegrees[i].pid == rinfo->edegrees[j].pid, "
                 "%d %d %d %d\n",i,j,rinfo->edegrees[i].pid, rinfo->edegrees[j].pid);
     }
  }

  return 1;
}



/*************************************************************************
* This function checks the correctness of the NodeFM data structures
**************************************************************************/
int CheckNodePartitionParams(GraphType *graph)
{
   int i, j, nvtxs, me, other;
   idxtype *xadj, *adjncy, *vwgt, *where;
   idxtype edegrees[2], pwgts[3];

   nvtxs = graph->nvtxs;
   xadj = graph->xadj;
   vwgt = graph->vwgt;
   adjncy = graph->adjncy;
   where = graph->where;

/*------------------------------------------------------------
/ Compute now the separator external degrees
/------------------------------------------------------------*/
   pwgts[0] = pwgts[1] = pwgts[2] = 0;
   for (i=0; i<nvtxs; i++) {
      me = where[i];
      pwgts[me] += vwgt[i];
      if (me == 2) { /* If it is on the separator do some computations */
         edegrees[0] = edegrees[1] = 0;
         for (j=xadj[i]; j<xadj[i+1]; j++) {
            other = where[adjncy[j]];
            if (other != 2)
               edegrees[other] += vwgt[adjncy[j]];
         }
         if (edegrees[0] != graph->nrinfo[i].edegrees[0] ||
            edegrees[1] != graph->nrinfo[i].edegrees[1]) {
            printf("Something wrong with edegrees: %d %d %d %d %d\n", i,
                   edegrees[0], edegrees[1], graph->nrinfo[i].edegrees[0],
                   graph->nrinfo[i].edegrees[1]);
            return 0;
         }
      }
   }
   if (pwgts[0] != graph->pwgts[0] || pwgts[1] != graph->pwgts[1] ||
       pwgts[2] != graph->pwgts[2])
      printf("Something wrong with part-weights: %d %d %d %d %d %d\n",
             pwgts[0], pwgts[1], pwgts[2], graph->pwgts[0],
            graph->pwgts[1], graph->pwgts[2]);
   return 1;
}


/*************************************************************************
* This function checks if the separator is indeed a separator
**************************************************************************/
int IsSeparable(GraphType *graph)
{
   int i, j, nvtxs, other;
   idxtype *xadj, *adjncy, *where;
   nvtxs = graph->nvtxs;
   xadj = graph->xadj;
   adjncy = graph->adjncy;
   where = graph->where;

   for (i=0; i<nvtxs; i++) {
      if (where[i] == 2)
         continue;
      other = (where[i]+1)%2;
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (where[adjncy[j]] == other)
          printf("Warning: %d %d %d %d %d %d\n", i, where[i], adjncy[j],
                 where[adjncy[j]], xadj[i+1]-xadj[i], xadj[adjncy[j]+1]-xadj[adjncy[j]]);
      }
   }
   return 1;
}
