/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2025 Rachid Touzani

   This file is part of OFELI.

   OFELI is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   OFELI is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with OFELI. If not, see <http://www.gnu.org/licenses/>.

  ==============================================================================

                   Various functions to determine matrix graph

  ==============================================================================*/


#include "linear_algebra/GraphOfMatrix.h"
#include "mesh/Mesh.h"
#include <algorithm>

using std::min;
using std::max;

namespace OFELI {


size_t getMatrixLength(const vector<size_t>& ch)
{
   size_t l=0;
   for (size_t i=0; i<ch.size(); i++)
      l += ch[i];
   return l;
}


size_t FinalizeGraph(vector<size_t>& row_ptr,
                     vector<size_t>& col_ind,
                     vector<RC>&     IJ,
                     vector<size_t>& nbc)
{
   sort(IJ.begin(),IJ.end());
   vector<RC>::iterator new_end=unique(IJ.begin(),IJ.end());
   IJ.erase(new_end,IJ.end());
   StoreGraph(IJ,row_ptr,col_ind);
#ifdef USE_EIGEN
   nbc.clear();
   for (size_t i=0; i<row_ptr.size()-1; i++)
      nbc.push_back(row_ptr[i+1]-row_ptr[i]);
#endif
   return IJ.size();
}


void FinalizeGraph(vector<RC>&     IJ,
                   vector<size_t>& col_ht)
{
   sort(IJ.begin(),IJ.end());
   vector<RC>::iterator new_end=unique(IJ.begin(),IJ.end());
   IJ.erase(new_end,IJ.end());
   for (const auto &a: IJ) {
      size_t i=a.first;
      col_ht[i-1] = std::max(col_ht[i-1],i-a.second+1);
   }
}


size_t SimpleSkyline(const Mesh&     m,
                     vector<size_t>& ch)
{
   int l;
   size_t n, size=m.getNbNodes();
   ch.resize(size,0);

   element_loop(&m) {
      int mini = int(size);
      for (n=1; n<=The_element.getNbNodes(); n++)
         mini = min(int(The_element(n)->n()),mini);
      for (n=1; n<=The_element.getNbNodes(); n++)
         if ((l=int(The_element(n)->n())) != 0)
            ch[l-1] = max(l-mini+1,int(ch[l-1]));
   }
   return getMatrixLength(ch);
}


size_t NodeSkyline(const Mesh&     m,
                   vector<size_t>& ch)
{
   int label_eq;
   size_t size=m.getNbEq();
   ch.resize(size,0);

   element_loop(&m) {
      int mini=int(size);
      for (size_t n=1; n<=The_element.getNbNodes(); n++) {
         the_node = The_element(n);
         for (size_t k=1; k<=The_node.getNbDOF(); k++)
            if (The_node.getDOF(k) != 0)
               mini = min(int(The_node.getDOF(k)),mini);
      }
      for (size_t n=1; n<=The_element.getNbNodes(); n++) {
         Node *nd=The_element(n);
         for (size_t l=1; l<=nd->getNbDOF(); l++)
            if ((label_eq=int(nd->getDOF(l))) != 0)
               ch[label_eq-1] = max(label_eq-mini+1,int(ch[label_eq-1]));
      }
   }
   return getMatrixLength(ch);
}


size_t NodeSkyline(const Mesh&     m,
                   vector<size_t>& ch,
                   size_t          dof1,
                   size_t          dof2)
{
   int i;
   size_t size=m.getNbNodes()*(dof2-dof1+1);
   ch.resize(size,0);

   element_loop(&m) {
      int mini = int(size);
      for (size_t n=1; n<=The_element.getNbNodes(); n++) {
          the_node = The_element(n);
         for (size_t k=dof1; k<=dof2; k++)
            if (The_node.getDOF(k) != 0)
               mini = min(int(The_node.getDOF(k)),mini);
      }
      for (size_t n=1; n<=The_element.getNbNodes(); n++) {
         the_node = The_element(n);
         for (size_t k=dof1; k<=dof2; k++)
            if ((i=int(The_node.getDOF(k))) != 0)
               ch[i-1] = max(i-mini+1,int(ch[i-1]));
      }
   }
   return getMatrixLength(ch);
}


size_t NodeSkyline(const Mesh&     m,
                   vector<size_t>& ch,
                   size_t          dof)
{
   int i;
   size_t size=m.getNbEq();
   ch.resize(size,0);
   element_loop(&m) {
      int mini=int(size);
      for (size_t n=1; n<=The_element.getNbNodes(); n++) {
         the_node = The_element(n);
         if (The_node.getDOF(dof) != 0)
            mini = min(int(node_label),mini);
      }
      for (size_t n=1; n<=The_element.getNbNodes(); n++)
         if ((i=int(The_element(n)->n())) != 0)
             ch[i-1] = max(i-mini+1,int(ch[i-1]));
   }
   return getMatrixLength(ch);
}


size_t SideSkyline(const Mesh&     m,
                   vector<size_t>& ch,
                   size_t          dof1,
                   size_t          dof2)
{
   int label_eq;
   size_t size=m.getNbEq();
   ch.resize(size,0);

   element_loop(&m) {
      int mini=int(size);
      for (size_t n=1; n<=The_element.getNbSides(); n++) {
         the_side = The_element.getPtrSide(n);
         for (size_t k=dof1; k<=dof2; k++)
            if (The_side.getDOF(k) != 0)
               mini = min(int(The_side.getDOF(k)),mini);
      }
      for (size_t n=1; n<=The_element.getNbSides(); n++) {
         the_side = The_element.getPtrSide(n);
         for (size_t l=dof1; l<=dof2; l++)
            if ((label_eq=int(The_side.getDOF(l))) != 0)
               ch[label_eq-1] = max(int(ch[label_eq-1]),label_eq-mini+1);
      }
   }
   return getMatrixLength(ch);
}


size_t SideSkyline(const Mesh&     m,
                   vector<size_t>& ch)
{
   int label_eq;
   size_t size=m.getNbEq();
   ch.resize(size,0);

   element_loop(&m) {
      int mini=int(size);
      for (size_t n=1; n<=The_element.getNbSides(); n++) {
         the_side = The_element.getPtrSide(n);
         for (size_t k=1; k<=The_side.getNbDOF(); k++)
            if (The_side.getDOF(k) != 0)
               mini = min(int(The_side.getDOF(k)),mini);
      }
      for (size_t n=1; n<=The_element.getNbSides(); n++) {
         the_side = The_element.getPtrSide(n);
         for (size_t l=1; l<=The_side.getNbDOF(); l++)
            if ((label_eq=int(The_side.getDOF(l))) != 0)
               ch[label_eq-1] = max(int(ch[label_eq-1]),label_eq-mini+1);
      }
   }
   return getMatrixLength(ch);
}


size_t SideSkyline(const Mesh&     m,
                   vector<size_t>& ch,
                   size_t          dof)
{
   dof = 0;
   int label_eq;
   size_t size=m.getNbEq();
   ch.resize(size,0);

   element_loop(&m) {
      int mini=int(size);
      for (size_t n=1; n<=The_element.getNbSides(); n++)
         mini = min(int(The_element.getSideLabel(n)),mini);
      for (size_t n=1; n<=The_element.getNbSides(); n++) {
         the_side = The_element.getPtrSide(n);
         label_eq = int(side_label);
         ch[label_eq-1] = max(int(ch[label_eq-1]),label_eq-mini+1);
      }
   }
   return getMatrixLength(ch);
}


size_t ElementSkyline(const Mesh&     m,
                      vector<size_t>& ch)
{
   const Element *el1, *el2;
   size_t size=m.getNbElements();
   ch.resize(size,0);

   side_loop(&m) {
      int mini = int(size);
      el1 = The_side.getNeighborElement(1);
      el2 = The_side.getNeighborElement(2);
      if (el1 && el2) {
         int i=int(el1->n()), j=int(el2->n());
         mini = min(mini,min(i,j));
         ch[i-1] = max(int(ch[i-1]),i+1-mini);
         ch[j-1] = max(int(ch[j-1]),j+1-mini);
      }
   }
   return getMatrixLength(ch);
}


size_t ElementSkyline(const Mesh&     m,
                      vector<size_t>& ch,
                      size_t          dof)
{
   dof = 0;
   const Element *el1, *el2;
   size_t size=m.getNbEq();
   ch.resize(size,0);

   side_loop(&m) {
      int mini = int(size);
      if ((el1=The_side.getNeighborElement(1)) &&
          (el2=The_side.getNeighborElement(2))) {
         mini = min(mini,min(int(el1->n()),int(el2->n())));
         int label_eq = int(el1->n());
         ch[label_eq-1] = max(int(ch[label_eq-1]),label_eq-mini+1);
         label_eq = int(el2->n());
         ch[label_eq-1] = max(int(ch[label_eq-1]),label_eq-mini+1);
      }
   }
   return getMatrixLength(ch);
}


size_t SimpleGraph(const Mesh&     m,
                   vector<long>&   xadj,
                   vector<size_t>& adjncy)
{
   size_t SimpleSkyline(const Mesh &m, vector<size_t> &ch);
   size_t k, l, length=0, size=m.getNbNodes();
   vector<size_t> ch(size,0);
   SimpleSkyline(m,ch);

   vector<RC> pp;
   element_loop(&m) {
      for (size_t in=1; in<=The_element.getNbNodes(); in++) {
         if ((k=The_element(in)->n())!=0)
            for (size_t jn=1; jn<=The_element.getNbNodes(); jn++) {
               if ((l=The_element(jn)->n())!=0) {
                  pp.push_back(RC(k-1,l-1));
                  if (k!=l)
                     pp.push_back(RC(l-1,k-1));
               }
            }
      }
   }
   sort(pp.begin(),pp.end());
   vector<RC>::iterator new_end=unique(pp.begin(),pp.end());
   pp.erase(new_end, pp.end());
   length = pp.size();
   xadj.resize(size+1);
   adjncy.resize(length);
   for (size_t i=0; i<size; i++)
      xadj[i] = 0;
   k = 0;
   for (size_t j=0; j<length; j++) {
      if (pp[j].first==pp[j].second)
         k++;
      xadj[k-1]++;
      adjncy[j] = pp[j].second + 1;
   }
   for (k=1; k<size; k++)
      xadj[k] += xadj[k-1];
   for (int ii=int(size); ii>0; ii--)
      xadj[ii] = xadj[ii-1] + 1;
   xadj[0] = 1;
   return length;
}


size_t NodeGraph(const Mesh&     mesh,
                 vector<size_t>& row_ptr,
                 vector<size_t>& col_ind,
                 vector<RC>&     IJ,
                 vector<size_t>& nbc)
{
   element_loop(&mesh) {
      for (size_t in=1; in<=The_element.getNbNodes(); in++) {
         Node *nd1=The_element(in);
         for (size_t k=1; k<=nd1->getNbDOF(); k++) {
            size_t nk=nd1->getDOF(k);
            if (nk) {
               for (size_t jn=1; jn<=The_element.getNbNodes(); jn++) {
                  Node *nd2=The_element(jn);
                  for (size_t l=1; l<=nd2->getNbDOF(); l++) {
                     size_t nl=nd2->getDOF(l);
                     if (nl) {
                        IJ.push_back(RC(nk-1,nl-1));
                        if (nk!=nl)
                           IJ.push_back(RC(nl-1,nk-1));
                     }
                  }
               }
            }
         }
      }
   }
   return FinalizeGraph(row_ptr,col_ind,IJ,nbc);
}


size_t NodeGraph(const Mesh&     mesh,
                 size_t          dof1,
                 size_t          dof2,
                 vector<size_t>& row_ptr,
                 vector<size_t>& col_ind,
                 vector<RC>&     IJ,
                 vector<size_t>& nbc)
{
   element_loop(&mesh) {
      for (size_t in=1; in<=The_element.getNbNodes(); in++) {
         Node *nd1=The_element(in);
         for (size_t k=dof1; k<=dof2; k++) {
            if (nd1->getDOF(k)!=0) {
               for (size_t jn=1; jn<=The_element.getNbNodes(); jn++) {
                  Node *nd2=The_element(jn);
                  for (size_t l=dof1; l<=dof2; l++) {
                     if (nd2->getDOF(l)!=0) {
                        IJ.push_back(RC(nd1->getDOF(k)-1,nd2->getDOF(l)-1));
                        if (nd1->getDOF(k)!=nd2->getDOF(l))
                           IJ.push_back(RC(nd2->getDOF(l)-1,nd1->getDOF(k)-1));
                     }
                  }
               }
            }
         }
      }
   }
   return FinalizeGraph(row_ptr,col_ind,IJ,nbc);
}


size_t SideGraph(const Mesh&     mesh,
                 vector<size_t>& row_ptr,
                 vector<size_t>& col_ind,
                 vector<RC>&     IJ,
                 vector<size_t>& nbc)
{
   element_loop(&mesh) {
      for (size_t in=1; in<=The_element.getNbSides(); in++) {
         Side *sd1=The_element.getPtrSide(in);
         for (size_t k=1; k<=sd1->getNbDOF(); k++) {
            if (sd1->getDOF(k)!=0) {
               for (size_t jn=1; jn<=The_element.getNbSides(); jn++) {
                  Side *sd2=The_element.getPtrSide(jn);
                  for (size_t l=1; l<=sd2->getNbDOF(); l++) {
                     if (sd2->getDOF(l)!=0) {
                        IJ.push_back(RC(sd1->getDOF(k)-1,sd2->getDOF(l)-1));
                        if (sd1->getDOF(k)!=sd2->getDOF(l))
                           IJ.push_back(RC(sd2->getDOF(l)-1,sd1->getDOF(k)-1));
                     }
                  }
               }
            }
         }
      }
   }
   return FinalizeGraph(row_ptr,col_ind,IJ,nbc);
}


size_t SideGraph(const Mesh&     mesh,
                 size_t          dof1,
                 size_t          dof2,
                 vector<size_t>& row_ptr,
                 vector<size_t>& col_ind,
                 vector<RC>&     IJ,
                 vector<size_t>& nbc)
{
   element_loop(&mesh) {
      for (size_t in=1; in<=The_element.getNbSides(); in++) {
         Side *sd1=The_element.getPtrSide(in);
         for (size_t k=dof1; k<=dof2; k++) {
            if (sd1->getDOF(k)!=0) {
               for (size_t jn=1; jn<=The_element.getNbSides(); jn++) {
                  Side *sd2=The_element.getPtrSide(jn);
                  for (size_t l=dof1; l<=dof2; l++) {
                     if (sd2->getDOF(l)!=0) {
                        IJ.push_back(RC(sd1->getDOF(k)-1,sd2->getDOF(l)-1));
                        if (sd1->getDOF(k)!=sd2->getDOF(l))
                           IJ.push_back(RC(sd2->getDOF(l)-1,sd1->getDOF(k)-1));
                     }
                  }
               }
            }
         }
      }
   }
   return FinalizeGraph(row_ptr,col_ind,IJ,nbc);
}


size_t ElementGraph(const Mesh&     mesh,
                    vector<size_t>& row_ptr,
                    vector<size_t>& col_ind,
                    vector<RC>&     IJ,
                    vector<size_t>& nbc)
{
   Element *el1, *el2;
   side_loop(&mesh) {
      if ((el1=The_side.getNeighborElement(1)) && (el2=The_side.getNeighborElement(2))) {
         IJ.push_back(RC(el1->n()-1,el1->n()-1));
         IJ.push_back(RC(el2->n()-1,el2->n()-1));
         IJ.push_back(RC(el1->n()-1,el2->n()-1));
         IJ.push_back(RC(el2->n()-1,el1->n()-1));
      }
   }
   return FinalizeGraph(row_ptr,col_ind,IJ,nbc);
}


size_t SideNodeGraph(const Mesh&     mesh,
                     vector<size_t>& row_ptr,
                     vector<size_t>& col_ind,
                     vector<RC>&     IJ,
                     vector<size_t>& nbc)
{
   side_loop(&mesh) {
      for (size_t k=1; k<=The_side.getNbDOF(); k++) {
         if (The_side.getCode(k)!=0) {
            for (size_t in=1; in<=The_side.getNbNodes(); in++) {
               the_node = The_side(in);
               for (size_t l=1; l<=The_node.getNbDOF(); l++)
                  if (The_node.getDOF(l)!=0)
                     IJ.push_back(RC(The_node.getDOF(k)-1,The_node.getDOF(l)-1));
            }
         }
      }
   }
   return FinalizeGraph(row_ptr,col_ind,IJ,nbc);
}


size_t NodeSideGraph(const Mesh&     mesh,
                     vector<size_t>& row_ptr, 
                     vector<size_t>& col_ind,
                     vector<RC>&     IJ,
                     vector<size_t>& nbc)
{
   side_loop(&mesh) {
      for (size_t k=1; k<=The_side.getNbDOF(); k++) {
         if (theSide->getCode(k)!=0) {
            for (size_t in=1; in<=The_side.getNbNodes(); in++) {
               the_node = The_side(in);
               for (size_t l=1; l<=The_node.getNbDOF(); l++)
                  if (The_node.getDOF(l)!=0)
                     IJ.push_back(RC(The_node.getDOF(l)-1,The_node.getDOF(k)-1));
            }
         }
      }
   }
   return FinalizeGraph(row_ptr,col_ind,IJ,nbc);
}


size_t NodeGraphScal(const Mesh&     mesh,
                     vector<size_t>& row_ptr,
                     vector<size_t>& col_ind,
                     vector<RC>&     IJ,
                     vector<size_t>& nbc)
{
   element_loop(&mesh) {
      for (size_t in=1; in<=The_element.getNbNodes(); in++) {
         for (size_t jn=1; jn<=The_element.getNbNodes(); jn++) {
            IJ.push_back(RC(The_element(in)->n()-1,The_element(jn)->n()-1));
            if (The_element(in)->n() != The_element(jn)->n())
               IJ.push_back(RC(The_element(jn)->n()-1,The_element(in)->n()-1));
         }
      }
   }
   return FinalizeGraph(row_ptr,col_ind,IJ,nbc);
}


size_t NodeGraphScal(const Mesh&     mesh,
                     size_t          dof,
                     size_t          nb_eq,
                     vector<size_t>& row_ptr,
                     vector<size_t>& col_ind,
                     vector<RC>&     IJ,
                     vector<size_t>& nbc)
{
   element_loop(&mesh) {
      for (size_t in=1; in<=The_element.getNbNodes(); in++) {
         size_t ii=The_element(in)->getDOF(dof);
         if (ii) {
            for (size_t jn=1; jn<=The_element.getNbNodes(); jn++) {
                size_t jj=The_element(jn)->getDOF(dof);
                if (jj) {
                   IJ.push_back(RC(ii-1,jj-1));
                   if (ii != jj)
                      IJ.push_back(RC(jj-1,ii-1));
                }
            }
         }
      }
   }
   return FinalizeGraph(row_ptr,col_ind,IJ,nbc);
}


size_t XGraph(const Mesh&     mesh,
              vector<size_t>& row_ptr,
              vector<size_t>& col_ind,
              vector<RC>&     IJ,
              vector<size_t>& nbc)
{
   element_loop(&mesh) {
      for (size_t in=1; in<=The_element.getNbNodes(); in++) {
         for (size_t jn=1; jn<=The_element.getNbNodes(); jn++) {
            IJ.push_back(RC(The_element(in)->n()-1,The_element(jn)->n()-1));
            if (The_element(in)->n() != The_element(jn)->n())
               IJ.push_back(RC(The_element(jn)->n()-1,The_element(in)->n()-1));
         }
      }
   }
   FinalizeGraph(row_ptr,col_ind,IJ,nbc);
   IJ.clear();

// Extended graph
   size_t k, l, size=mesh.getNbNodes();
   for (size_t i=0; i<size; i++) {
      IJ.push_back(RC(i,i));
      for (size_t j=row_ptr[i]-1; j<row_ptr[i+1]-1; j++) {
         if ((k=col_ind[j]) != i+1) {
            IJ.push_back(RC(k-1,i));
            if (k-1 != i)
               IJ.push_back(RC(i,k-1));
            for (size_t m=row_ptr[k-1]-1; m<row_ptr[k]-1; m++) {
               if ((l=col_ind[m])!=i+1 && l!=k) {
                  IJ.push_back(RC(i,l-1));
                  if (l-1 != i)
                     IJ.push_back(RC(l-1,i));
               }
            }
         }
      }
   }
   return FinalizeGraph(row_ptr,col_ind,IJ,nbc);
}


void StoreGraph(const vector<RC>& IJ,
                vector<size_t>&   row_ptr,
                vector<size_t>&   col_ind)
{
   size_t length=IJ.size();
   col_ind.clear();
   row_ptr.clear();
   col_ind.push_back(IJ[0].second+1);
   row_ptr.push_back(0);
   for (size_t k=1; k<length; k++) {
      col_ind.push_back(IJ[k].second+1);
      if (IJ[k].first>IJ[k-1].first)
         row_ptr.push_back(k);
   }
   row_ptr.push_back(length);
}


size_t DG0Graph(const Mesh&     mesh,
                vector<RC>&     I,
                vector<size_t>& nbc)
{
   I.clear();
   Element *el1, *el2, *el;
   vector<RC> IJ;
   element_loop(&mesh) {
      for (size_t i=1; i<=The_element.getNbSides(); i++) {
         if ((el=The_element.getNeighborElement(i))) {
            IJ.push_back(RC(el->n(),element_label));
            for (size_t j=1; j<=el->getNbSides(); j++) {
               el1 = el->getPtrSide(j)->getNeighborElement(1);
               el2 = el->getPtrSide(j)->getNeighborElement(2);
               if (el1 && el1!=el)
                  IJ.push_back(RC(el1->n(),element_label));
               if (el2 && el2!=el)
                  IJ.push_back(RC(el2->n(),element_label));
            }
         }
      }
   }
   for (vector<RC>::iterator it=IJ.begin(); it!=IJ.end(); it++)
      I.push_back(*it);
   return IJ.size();
}


size_t DGGraph(const Mesh&     mesh,
               vector<RC>&     I,
               vector<size_t>& nbc)
{
   I.clear();
   size_t fd=1;
   element_loop(&mesh) {
      The_element.setFirstDOF(fd);
      for (size_t i=1; i<=The_element.getNbDOF(); i++)
        for (size_t j=1; j<=The_element.getNbDOF(); j++)
           I.push_back(RC(fd+i-1,fd+j-1));
      fd += The_element.getNbDOF();
   }

   side_loop(&mesh) {
      Element *el1=The_side.getNeighborElement(1),
              *el2=The_side.getNeighborElement(2);
      if (el1 && el2) {
         for (size_t i=1; i<=el1->getNbDOF(); i++) {
            for (size_t j=1; j<=el2->getNbDOF(); j++) {
               I.push_back(RC(el1->getFirstDOF()+i-1,el2->getFirstDOF()+j-1));
               I.push_back(RC(el2->getFirstDOF()+j-1,el1->getFirstDOF()+i-1));
            }
         }
      }
   }
   return I.size();
}

} /* namespace OFELI */
