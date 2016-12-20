/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2017 Rachid Touzani

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

                     Implementation of utility mesh functions

  ==============================================================================*/

#include <algorithm>
using std::min;
using std::max;
using std::pair;

#include <vector>
using std::vector;

#include <string>
using std::string;

#include "mesh/MeshUtil.h"
#include "linear_algebra/Vect.h"
#include "mesh/Mesh.h"
#include "mesh/Grid.h"
#include "shape_functions/Line2.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Quad4.h"
#include "shape_functions/Tetra4.h"
#include "shape_functions/Hexa8.h"
#include "shape_functions/Penta6.h"
#include "linear_algebra/LocalVect.h"
#include "post/Reconstruction.h"

namespace OFELI {


bool operator==(const Element& el1,
                const Element& el2)
{
#ifdef _OFELI_RANGE_CHECK
   assert(el1.getNbNodes() == el2.getNbNodes());
#endif
   size_t nn=el1.getNbNodes();
   vector<size_t> en1(nn), en2(nn);
   for (size_t i=0; i<nn; i++) {
      en1[i] = el1(i+1)->n();
      en2[i] = el2(i+1)->n();
   }
   std::sort(en1.begin(),en1.end());
   std::sort(en2.begin(),en2.end());
   if (en1==en2)
      return true;
   else
      return false;
}


bool operator==(const Side& sd1,
                const Side& sd2)
{
#ifdef _OFELI_RANGE_CHECK
   assert(sd1.getNbNodes() == sd2.getNbNodes());
#endif
   size_t nn=sd1.getNbNodes();
   vector<size_t> sn1(nn), sn2(nn);
   for (size_t i=0; i<nn; i++) {
      sn1[i] = sd1(i+1)->n();
      sn2[i] = sd2(i+1)->n();
   }
   std::sort(sn1.begin(),sn1.end());
   std::sort(sn2.begin(),sn2.end());
   if (sn1==sn2)
      return true;
   else
      return false;
}


void QuadrilateralsToTriangles(Mesh& m1,
                               Mesh& m2)
{
   if (m1.getDim()!=2) {
      cerr << "Valid for 2-D meshes only !!" << endl;
      exit(0);
   }
   m2.setDim(2);
   mesh_nodes(m1)
      m2.Add(the_node);
   int code;
   size_t nn[4], n=1;
   mesh_elements(m1) {
      if (the_element->getShape()==QUADRILATERAL) {
         for (size_t k=0; k<4; k++)
            nn[k] = The_element(k+1)->n();
         code = The_element.getCode();
         Element *el = new Element(n++,TRIANGLE,code);
         el->Add(m2[nn[0]]);
         el->Add(m2[nn[1]]);
         el->Add(m2[nn[2]]);
         m2.Add(el);
         el = new Element(n++,TRIANGLE,code);
         el->Add(m2[nn[2]]);
         el->Add(m2[nn[3]]);
         el->Add(m2[nn[0]]);
         m2.Add(el);
      }
      else
         m2.Add(the_element);
   }
}


void HexahedraToTetrahedra(Mesh& m1,
                           Mesh& m2)
{
   if (m1.getDim()!=3) {
      cerr << "Valid for 3-D meshes only !!" << endl;
      exit(0);
   }
   m2.setDim(3);
   mesh_nodes(m1)
      m2.Add(the_node);

   int code;
   size_t nn[8], n=1;
   mesh_elements(m1) {
      if (the_element->getShape()==HEXAHEDRON) {
         for (size_t k=0; k<8; k++)
            nn[k] = The_element(k+1)->n();
         code = the_element->getCode();
         Element *el = new Element(n++,TETRAHEDRON,code);
         el->Add(m1[nn[4]]);
         el->Add(m1[nn[1]]);
         el->Add(m1[nn[3]]);
         el->Add(m1[nn[0]]);
         m2.Add(el);
         el = new Element(n++,TETRAHEDRON,code);
         el->Add(m1[nn[7]]);
         el->Add(m1[nn[6]]);
         el->Add(m1[nn[5]]);
         el->Add(m1[nn[2]]);
         m2.Add(el);
         el = new Element(n++,TETRAHEDRON,code);
         el->Add(m1[nn[7]]);
         el->Add(m1[nn[4]]);
         el->Add(m1[nn[5]]);
         el->Add(m1[nn[3]]);
         m2.Add(el);
         el = new Element(n++,TETRAHEDRON,code);
         el->Add(m1[nn[7]]);
         el->Add(m1[nn[5]]);
         el->Add(m1[nn[2]]);
         el->Add(m1[nn[3]]);
         m2.Add(el);
         el = new Element(n++,TETRAHEDRON,code);
         el->Add(m1[nn[5]]);
         el->Add(m1[nn[1]]);
         el->Add(m1[nn[2]]);
         el->Add(m1[nn[3]]);
         m2.Add(el);
         el = new Element(n++,TETRAHEDRON,code);
         el->Add(m1[nn[5]]);
         el->Add(m1[nn[4]]);
         el->Add(m1[nn[1]]);
         el->Add(m1[nn[3]]);
         m2.Add(el);
      }
      else
         m2.Add(the_element);
   }
}


void DeformMesh(      Mesh&         mesh,
                const Vect<real_t>& u,
                      real_t        a)
{
   mesh_nodes(mesh) {
      size_t n=node_label;
      for (size_t m=1; m<=mesh.getDim(); m++)
         The_node.setCoord(m,The_node.getCoord(m)+a*u(n,m));
   }
}


#ifdef USE_PETSC
void DeformMesh(      Mesh&              mesh,
                const PETScVect<real_t>& u,
                      real_t             a)
{
   mesh_nodes(mesh) {
      size_t n=node_label, m=The_node.getNbDOF();
      The_node.setCoord(1,The_node.getX()+a*u(n,1));
      The_node.setCoord(2,The_node.getY()+a*u(n,2));
      if (m>2)
         The_node.setCoord(3,The_node.getZ()+a*u(n,3));
   }
}
#endif

void Refine(Mesh& in_mesh,
            Mesh& out_mesh)
{
   int code[MAX_NBDOF_NODE];
   size_t i;
   in_mesh.getAllSides();
   out_mesh._first_dof = 1;
   out_mesh._nb_dof = 0;
   out_mesh._nb_nodes = out_mesh._nb_elements = out_mesh._nb_sides = out_mesh._nb_vertices = 0;
   out_mesh._dim = in_mesh._dim;
   out_mesh._verb = in_mesh._verb;

// Copy nodes
   size_t first_dof = 1;
   mesh_nodes(in_mesh) {
      Node *nd = new Node(the_node->n(),the_node->getCoord());
      size_t nb_dof = the_node->getNbDOF();
      nd->setNbDOF(nb_dof);
      for (i=0; i<nb_dof; i++)
         code[i] = the_node->getCode(i+1);
      nd->setDOF(first_dof,nb_dof);
      nd->setCode(code);
      out_mesh.Add(nd);
   }

// Create midside nodes
   size_t nd_label = in_mesh.getNbNodes();
   mesh_sides(in_mesh) {
      Node *nd1 = The_side(1), *nd2 = The_side(2);
      Point<real_t> x = 0.5*(nd1->getCoord() + nd2->getCoord());
      the_node = new Node(++nd_label,x);
      size_t nb_dof = nd1->getNbDOF();
      the_node->setNbDOF(nb_dof);
      for (i=0; i<nb_dof; i++) {
         int c1=nd1->getCode(i+1), c2=nd2->getCode(i+1);
         code[i] = c1;
         if (c2==0)
            code[i] = c2;
      }
      the_node->setDOF(in_mesh._first_dof,nb_dof);
      the_node->setCode(code);
      out_mesh.Add(the_node);
      The_side.Add(the_node);
   }

// Divide elements
   size_t el_label=0;
   mesh_elements(in_mesh) {
      try {
         if (The_element.getShape() != TRIANGLE)
            THROW_RT(" Element " + itos(element_label) + " is not a triangle.");
      }
      CATCH("Refine(Mesh,Mesh):");
      Node *nd[6]={NULL,NULL,NULL,NULL,NULL,NULL};
      nd[0] = The_element(1); 
      nd[1] = The_element(2); 
      nd[2] = The_element(3);
      for (i=1; i<=3; i++) {
         Side *sd=The_element.getPtrSide(i);
         size_t n1=(*sd)(1)->n(), n2=(*sd)(2)->n();
         if ((n1==nd[0]->n() && n2==nd[1]->n()) || (n1==nd[1]->n() && n2==nd[0]->n()))
            nd[3] = (*sd)(3);
         else if ((n1==nd[1]->n() && n2==nd[2]->n()) || (n1==nd[2]->n() && n2==nd[1]->n()))
            nd[4] = (*sd)(3);
         else if ((n1==nd[2]->n() && n2==nd[0]->n()) || (n1==nd[0]->n() && n2==nd[2]->n()))
            nd[5] = (*sd)(3);
      }
      Element *el = new Element(++el_label,The_element.getShape(),The_element.getCode());
      Element *e = in_mesh(element_label);
      int code = The_element.getCode();
      el->Add(nd[0]); el->Add(nd[3]); el->Add(nd[5]);
      el->setCode(code);
      out_mesh.Add(el);
      e->setChild(el);
      el = new Element(++el_label,The_element.getShape(),The_element.getCode());
      el->Add(nd[3]); el->Add(nd[1]); el->Add(nd[4]);
      el->setCode(code);
      out_mesh.Add(el);
      e->setChild(el);
      el = new Element(++el_label,The_element.getShape(),The_element.getCode());
      el->Add(nd[4]); el->Add(nd[2]); el->Add(nd[5]);
      el->setCode(code);
      out_mesh.Add(el);
      e->setChild(el);
      el = new Element(++el_label,The_element.getShape(),The_element.getCode());
      el->Add(nd[3]); el->Add(nd[4]); el->Add(nd[5]);
      el->setCode(code);
      out_mesh.Add(el);
      e->setChild(el);
   }
}


void Refine(Mesh& in_mesh,
            Mesh& out_mesh,
            int   n)
{
   Mesh ms(in_mesh);
   ms.getAllSides();
   for (int i=1; i<=n; i++) {
      Refine(ms,out_mesh);
      ms = out_mesh;
      ms.getAllSides();
   }
}


size_t init_side_node_numbering(int                      shape,
                                vector<vector<size_t> >& nsd,
                                int&                     sh)
{
   size_t ns=2;
   sh = LINE;
   if (shape == TRIANGLE) {
      nsd.resize(3);
      nsd[0].push_back(1); nsd[0].push_back(2);
      nsd[1].push_back(2); nsd[1].push_back(3);
      nsd[2].push_back(3); nsd[2].push_back(1);
   }
   else if (shape == QUADRILATERAL) {
      nsd.resize(4);
      nsd[0].push_back(1); nsd[0].push_back(2);
      nsd[1].push_back(2); nsd[1].push_back(3);
      nsd[2].push_back(3); nsd[2].push_back(4);
      nsd[3].push_back(4); nsd[3].push_back(1);
   }
   if (shape == TETRAHEDRON) {
      ns = 3;
      sh = TRIANGLE;
      nsd.resize(4);
      nsd[0].push_back(1); nsd[0].push_back(2); nsd[0].push_back(3);
      nsd[1].push_back(2); nsd[1].push_back(4); nsd[1].push_back(3);
      nsd[2].push_back(4); nsd[2].push_back(1); nsd[2].push_back(3);
      nsd[3].push_back(1); nsd[3].push_back(4); nsd[3].push_back(2);
   }
   return ns;
}


void getWindowField(const Vect<real_t>& u,
                          Vect<real_t>& v)
{
   const Mesh *ms = &(v.getMesh());
   for (size_t n=1; n<=ms->getNbNodes(); n++)
      for (size_t i=1; i<=v.getNbDOF(); i++)
         v(n,i) = u(ms->getNodeNewLabel(n),i);
}


int equal_sides(const Side* sd1,
                const Side* sd2)
{
   size_t n = sd1->getNbNodes();
   vector<size_t> s1(n), s2(n);
   for (size_t i=1; i<=n; i++) {
      s1[i-1] = (*sd1)(i)->n();
      s2[i-1] = (*sd2)(i)->n();
   }
   std::sort(s1.begin(),s1.end());
   std::sort(s2.begin(),s2.end());
   if (s1==s2)
      return 1;
   else
      return 0;
}


int equal_sides(const Side*           sd,
                      vector<size_t>& s)
{
   size_t n = sd->getNbNodes();
   vector<size_t> s1(n);
   for (size_t i=1; i<=n; i++)
      s1[i-1] = (*sd)(i)->n();
   std::sort(s1.begin(),s1.end());
   std::sort(s.begin(),s.end());
   if (s1==s)
      return 1;
   else
      return 0;
}


void order_side_nodes(size_t ns,
                      ND&    s)
{
   size_t n=0, t = s.nd[0];
   s.n = ns;
   ND ss(s);
   if (t > s.nd[1]) {
      t = s.nd[1];
      n = 1;
   }
   if (ns>2 && t>s.nd[2])
      t = s.nd[2], n = 2;
   if (ns>3 && t>s.nd[3])
      t = s.nd[3], n = 3;
   s.nd[0] = ss.nd[n];
   s.nd[1] = ss.nd[(n+1)%ns];
   if (ns > 2)
      s.nd[2] = ss.nd[(n+2)%ns];
   ss = s;

   if (ns > 2) {
      n = 1;
      t = s.nd[1];
      if (ns > 2 && t > s.nd[2]) {
         s.nd[1] = ss.nd[2];
         s.nd[2] = ss.nd[1];
      }
   }
}


void order_edge_nodes(ND& s)
{
   size_t n=0, t = s.nd[0];
   s.n = 2;
   ND ss(s);
   if (t > s.nd[1])
      t = s.nd[1], n = 1;
   s.nd[0] = ss.nd[n];
   s.nd[1] = ss.nd[(n+1)%2];
}


bool compare_sides(const ND& s1,
                   const ND& s2)
{
   if (s1.nd[0] > s2.nd[0])
      return false;
   else if (s1.nd[0] < s2.nd[0])
      return true;
   else {
      if (s1.nd[1] > s2.nd[1])
         return false;
      else if (s1.nd[1] < s2.nd[1])
         return true;
      else {
         if (s1.nd[2] > s2.nd[2])
            return false;
         else if (s1.nd[2] < s2.nd[2])
            return true;
         else {
            if (s1.nd[0] > s2.nd[3])
               return false;
            else
               return true;
         }
      }
   }
}


bool compare_edges(const ND& s1,
                   const ND& s2)
{
   if (s1.nd[0] > s2.nd[0])
      return false;
   else if (s1.nd[0] < s2.nd[0])
      return true;
   else {
      if (s1.nd[1] > s2.nd[1])
         return false;
      else if (s1.nd[1] < s2.nd[1])
         return true;
      else ;
   }
   return false;
}


size_t remove_internal_sides(vector<ND>& p)
{
   size_t k=0, n=p.size();
   bool f;
   for (size_t i=0; i<n; i++) {
      size_t j = i+1;
      f = true;
      while (j<n && p[j]==p[i]) {
         p[i].e2 = p[j].e1;
         j++;
         f = false;
      }
      if (f == true)
         p[k++] = p[i];
      i = j - 1;
   }
   p.erase(p.begin()+k,p.end());
   return k;
}


void complete_sides(vector<ND>& p)
{
   size_t n = p.size();
   for (size_t i=0; i<n; i++) {
      size_t j = i+1;
      while (j<n && p[j]==p[i]) {
         p[i].e2 = p[j].e1;
         j++;
      }
      i = j - 1;
   }
}


void order_edge_nodes(pair<size_t,size_t>& ed)
{
   size_t t = ed.first;
   if (t > ed.second) {
      t = ed.second;
      ed.second = ed.first;
      ed.first = t;
   }
}


size_t clean_edges(      vector<pair<size_t,size_t> >& p,
                         vector<pair<size_t,size_t> >& q,
                   const size_t&                       n)
{
   size_t k=0;
   for (size_t i=0; i<n; i++) {
      size_t j = i+1;
      q[k++] = p[i];
      while (j<n && p[j]==p[i]) {
//        q[k-1].e2 = p[j].e1;
        j++;
      }
      i = j - 1;
   }
   return k;
}


void FindRoot(size_t&         root,
              vector<long>&   xadj,
              vector<size_t>& adjncy,
              vector<size_t>& mask,
              long&           nlvl,
              vector<size_t>& xls,
              size_t*         ls)
{
   size_t node, mindeg, ndeg;
   long nabor, kstop, jstrt, kstrt, nunlvl, ccsize, k;
   real_t csize=0.;

// Determine the level structure rooted at root
   RootLs(root, xadj, adjncy, mask, nlvl, xls, ls);
   ccsize = xls[nlvl] - 1;
   if (nlvl==1 || nlvl==ccsize)
      return;

// Pick a node with minimum degree from the last level

   do {
      jstrt = xls[nlvl-1];
      mindeg = ccsize;
      root = ls[jstrt-1];
      if (ccsize != jstrt) {
         for (long jj=jstrt; jj<=ccsize; ++jj) {
            node = ls[jj-1];
            ndeg = 0;
            kstrt = xadj[node-1];
            kstop = xadj[node] - 1;
            for (k=kstrt; k<=kstop; ++k) {
               nabor = adjncy[k-1];
               if (mask[nabor-1] > 0)
                  ++ndeg;
            }
            if (ndeg < mindeg) {
               root = node;
               mindeg = ndeg;
            }
         }
      }

//    And generate its rooted level structure
      RootLs(root, xadj, adjncy, mask, nunlvl, xls, ls);
      if (nunlvl <= nlvl)
         return;
      nlvl = nunlvl;
   } while (real_t(nlvl) < csize);
}


void RootLs(size_t&         root,
            vector<long>&   xadj,
            vector<size_t>& adjncy,
            vector<size_t>& mask,
            long&           nlvl,
            vector<size_t>& xls,
            size_t*         ls)
{
    size_t nbr;
    long node, lbegin, ccsize, lvlend, lvsize, i;
    long unsigned jstrt, jstop;

    mask[root-1] = 0;
    ls[0] = root;
    nlvl = 0;
    lvlend = 0;
    ccsize = 1;

//  lbegin is the pointer to the beginning of the current level,
//  and lvlend points to the end of this level

    do {
       lbegin = lvlend + 1;
       lvlend = ccsize;
       xls[nlvl++] = lbegin;

//     Generate the next level by finding all the masked neighbors
//     of nodes in the current level

       for (i=lbegin-1; i<lvlend; ++i) {
          node = ls[i];
          jstrt = xadj[node-1];
          jstop = xadj[node] - 1;
          if (jstop >= jstrt) {
             for (unsigned long jj=jstrt-1; jj<jstop; ++jj) {
                nbr = adjncy[jj];
                if (mask[nbr-1] != 0) {
                   ls[ccsize++] = nbr;
                   mask[nbr-1] = 0;
                }
             }
          }
       }

//     Compute the current level width. if it is nonzero, Generate the next level
       lvsize = ccsize - lvlend;
    } while (lvsize>0);

//  Reset mask to one for the nodes in the level structure
    xls[nlvl] = lvlend + 1;
    for (long j=0; j<ccsize; j++)
       mask[ls[j]-1] = 1;
}


void RCM(size_t& root,
         vector<long>&   xadj,
         vector<size_t>& adjncy,
         vector<size_t>& mask,
         size_t*         perm,
         size_t&         ccsize,
         vector<size_t>& deg)
{
    unsigned long jstrt, jstop;
    long l;
    size_t node, lbegin, lvlend, nbr, lnbr, fnbr, lperm, k;

    Degree(root, xadj, adjncy, mask, deg, ccsize, perm);
    mask[root-1] = 0;
    if (ccsize <= 1)
       return;
    lvlend = 0;
    lnbr = 1;

    do {
       lbegin = lvlend + 1;
       lvlend = lnbr;
       for (size_t i=lbegin; i<=lvlend; ++i) {
          node = perm[i-1];
          jstrt = xadj[node-1];
          jstop = xadj[node] - 1;
          fnbr = lnbr + 1;
          for (unsigned long jj=jstrt; jj<=jstop; ++jj) {
             nbr = adjncy[jj-1];
             if (mask[nbr-1] != 0) {
                mask[nbr-1] = 0;
                perm[lnbr++] = nbr;
             }
          }
          if (fnbr < lnbr) {
             k = fnbr;
             while (k >= lnbr) {
                l = k;
                nbr = perm[k++];

L400:
                if (l >= long(fnbr)) {
                   lperm = perm[l-1];
                   if (deg[lperm-1] > deg[nbr-1]) {
                      perm[l--] = lperm;
                      goto L400;
                   }
                }
                perm[l] = nbr;
             }
          }
       }
    } while (lnbr > lvlend);

//  Reverse the Cuthill McKee ordering
    k = ccsize/2; l = ccsize - 1;
    for (size_t i=1; i<=k; ++i) {
       lperm = perm[l];
       perm[l--] = perm[i-1];
       perm[i-1] = lperm;
    }
}


void Degree(size_t&         root,
            vector<long>&   xadj,
            vector<size_t>& adjncy,
            vector<size_t>& mask,
            vector<size_t>& deg,
            size_t&         ccsize,
            size_t*         ls)
{
   size_t lbegin, lvlend, lvsize, node, nbr;
   long jstrt, jstop, j, ideg;

   ls[0] = root;
   xadj[root-1] = -xadj[root-1];
   lvlend = 0; ccsize = 1;

   lvsize = 1;
   while (lvsize > 0) {
     lbegin = lvlend + 1;
     lvlend = ccsize;

//   Find the degrees of nodes in the current level and generate the next level
     for (size_t i=lbegin; i<=lvlend; ++i) {
        node = ls[i-1];
        jstrt = -xadj[node-1];
        long xa = xadj[node];
        jstop = ((xa > 0) ? xa : -xa) - 1;
        ideg = 0;
        if (jstop >= jstrt) {
          for (j=jstrt; j<=jstop; ++j) {
             nbr = adjncy[j-1];
             if (mask[nbr-1] != 0) {
               ++ideg;
               if (xadj[nbr-1] >= 0) {
                  xadj[nbr-1] = -xadj[nbr-1];
                  ls[ccsize++] = nbr;
               }
             }
          }
        }
        deg[node-1] = ideg;
     }
     lvsize = ccsize - lvlend;
   }

   for (size_t ii=0; ii<ccsize; ++ii)
      xadj[ls[ii]-1] = -xadj[ls[ii]-1];
}


size_t compare_node_list(size_t* nod1,
                         size_t* nod2,
                         size_t  n)
{
   size_t ret=0;
   for (size_t i=0; i<n; i++)
      for (size_t j=i+1; j<n; j++)
         if (nod1[i] == nod2[j])
           return 1;
   return ret;
}


int BoundaryConditionCode(vector<string>& str,
                          int*            code,
                          string          s)
{
   if (str[0] == s)
      return code[0];
   else if (str[1] == s)
      return code[1];
   else if (str[2] == s)
      return code[2];
   else if (str[3] == s)
      return code[2];
   else if (str[4] == s)
      return code[3];
   else if (str[5] == s)
      return code[4];
   else return atoi(s.c_str());
}


void DofCode(int    mark,
             size_t nb_dof,
             int*   code)
{
   int m;
   int j = mark;
   for (size_t k=0; k<nb_dof; k++) {
      int kk = int(pow(10.,real_t(nb_dof-k-1)));
      code[k] = m = j/kk;
      j -= m*kk;
   }
}


void MeshToGrid(      Mesh&         m,
                      Grid&         g,
                const Vect<real_t>& u,
                      Vect<real_t>& ug,
                      size_t        dof)
{
   Point<real_t> xm, xM;
   int nx=g.getNx(), ny=g.getNy(), nz=g.getNz();
   Point<real_t> h = Point<real_t>(g.getHx(),g.getHy(),g.getHz());
   size_t dim = m.getDim();

   mesh_elements(m) {
      try {
         if (dim==1 && The_element.getShape() != LINE)
            THROW_RT(" Element " + itos(element_label) + " is not a line.");
      }
      CATCH("MeshToGrid(Mesh,Grid,Vect<real_t>,Vect<real_t>):");
      try {
         if (dim==2 && The_element.getShape() != TRIANGLE)
            THROW_RT(" Element "+itos(element_label)+" is not a triangle.");
      }
      CATCH("MeshToGrid(Mesh,Grid,Vect<real_t>,Vect<real_t>)");
      try {
         if (dim==3 && The_element.getShape() != TETRAHEDRON)
            THROW_RT(" Element "+itos(element_label)+" is not a tetrahedron.");
      }
      CATCH("MeshToGrid(Mesh,Grid,Vect<real_t>,Vect<real_t>):");

      xm = xM = The_element(1)->getCoord();
      for (size_t l=2; l<=The_element.getNbNodes(); l++) {
         xm.x = std::min(xm.x,The_element(l)->getX());
         xM.x = std::max(xM.x,The_element(l)->getX());
         if (dim>1) {
            xm.y = std::min(xm.y,The_element(l)->getY());
            xM.y = std::max(xM.y,The_element(l)->getY());
         }
         if (dim>2) {
            xm.z = std::min(xm.z,The_element(l)->getZ());
            xM.z = std::max(xM.z,The_element(l)->getZ());
         }
      }
      int i1 = std::max(int(1+(xm.x-g.getXMin().x)/h.x),1);
      int i2 = std::min(int(2+(xM.x-g.getXMin().x)/h.x),nx+1);
      int j1=1, j2=1, k1=1, k2=1;
      if (dim>1) {
         j1 = std::max(int(1+(xm.y-g.getXMin().y)/h.y),1);
         j2 = std::min(int(2+(xM.y-g.getXMin().y)/h.y),ny+1);
      }
      if (dim>2) {
         k1 = std::max(int(1+(xm.z-g.getXMin().z)/h.z),1);
         k2 = std::min(int(2+(xM.z-g.getXMin().z)/h.z),nz+1);
      }
      Point<real_t> x = 0;
      for (int ii=i1; ii<=i2; ii++) {
         x.x = g.getXMin().x + (ii-1)*h.x;
         for (int jj=j1; jj<=j2; jj++) {
            x.y = g.getXMin().y + (jj-1)*h.y;
            for (int kk=k1; kk<=k2; kk++) {
               x.z = g.getXMin().z + (kk-1)*h.z;
               switch (The_element.getShape()) {

                  case LINE:
                     {
                        Line2 ln(the_element);
                        if (ln.isIn(x)) {
                           Point<real_t> s = ln.getRefCoord(x);
                           ug(ii) = ln.Sh(1,s)*u(The_element(1)->n(),dof) + 
                                    ln.Sh(2,s)*u(The_element(2)->n(),dof);
                        }
                     }
                     break;

                  case TRIANGLE:
                     {
                        Triang3 tr(the_element);
                        if (tr.isIn(x)) {
                           Point<real_t> s = tr.getRefCoord(x);
                           ug(ii,jj) = tr.Sh(1,s)*u(The_element(1)->n(),dof) +
                                       tr.Sh(2,s)*u(The_element(2)->n(),dof) +
                                       tr.Sh(3,s)*u(The_element(3)->n(),dof);
                        }
                        break;
                     }

                  case TETRAHEDRON:
                     {
                        Tetra4 te(the_element);
                        if (te.isIn(x)) {
                           Point<real_t> s = te.getRefCoord(x);
                           ug(ii,jj,kk) = te.Sh(1,s)*u(The_element(1)->n(),dof) +
                                          te.Sh(2,s)*u(The_element(2)->n(),dof) +
                                          te.Sh(3,s)*u(The_element(3)->n(),dof) +
                                          te.Sh(4,s)*u(The_element(4)->n(),dof);
                        }
                        break;
                     }
               }
            }
         }
      }
   }
}


void GridToMesh(      Grid&         g,
                      Mesh&         m,
                const Vect<real_t>& ug,
                      Vect<real_t>& u,
                      size_t        dof)
{
   Point<real_t> x, xi, s;
   size_t dim = m.getDim();
   int nx=g.getNx(), ny=g.getNy(), nz=g.getNz();
   int i, j, k;
   Point<real_t> h = Point<real_t>(g.getHx(),g.getHy(),g.getHz());

   mesh_nodes(m) {
      size_t n = node_label;
      Point<real_t> x = The_node.getCoord();
      i = std::min(1+int((x.x-g.getXMin().x)/h.x),nx+1);
      xi.x = g.getXMin().x + (i-1)*h.x;
      s.x = (x.x-xi.x)/h.x;
      if (dim==1) {
         if (u.getDOFType()==NODE_DOF) {
            u(n,dof) = (1-s.x)*ug(i);
            if (i<=nx)
               u(n,dof) += s.x*ug(i+1);
         }
         else if (u.getDOFType()==ELEMENT_DOF) {
            u(n,dof) = (1-s.x)*ug(i);
            if (i<=nx)
               u(n,dof) += s.x*ug(i+1);
         }
      }
      else {
         j = std::min(1+int((x.y-g.getXMin().y)/h.y),ny+1);
         xi.y = g.getXMin().y + (j-1)*h.y;
         s.y = (x.y-xi.y)/h.y;
         if (dim==2) {
            u(n,dof) = (1-s.x)*(1-s.y)*ug(i,j);
            if (i<=nx)
               u(n,dof) += s.x*(1-s.y)*ug(i+1,j);
            if (j<=ny)
               u(n,dof) += (1-s.x)*s.y*ug(i,j+1);
            if (i<=nx && j<=ny)
               u(n,dof) += s.x*s.y*ug(i+1,j+1);
	 }
      }
      if (dim==3) {
         k = std::min(1+int((x.z-g.getXMin().z)/h.z),nz+1);
         xi.z = g.getXMin().z + (k-1)*h.z;
         s.z = (x.z-xi.z)/h.z;
         u(n,dof) = (1-s.x)*(1-s.y)*(1-s.z)*ug(i,j,k);
         if (i<=nx)
            u(n,dof) += s.x*(1-s.y)*(1-s.z)*ug(i+1,j,k);
         if (j<=ny)
            u(n,dof) += (1-s.x)*s.y*(1-s.z)*ug(i,j+1,k);
         if (k<=nz)
            u(n,dof) += (1-s.x)*(1-s.y)*s.z*ug(i,j,k+1);
         if (i<=nx && j<=ny)
            u(n,dof) += s.x*s.y*(1-s.z)*ug(i+1,j+1,k);
         if (i<=nx && k<=nz)
            u(n,dof) += s.x*(1-s.y)*s.z*ug(i+1,j,k+1);
         if (j<=ny && k<=nz)
            u(n,dof) += (1-s.x)*s.y*s.z*ug(i,j+1,k+1);
         if (i<=nx && j<=ny && k<=nz)
            u(n,dof) += s.x*s.y*s.z*ug(i+1,j+1,k+1);
      }
   }
}


void MeshToMesh(      Mesh&         m1,
                      Mesh&         m2,
                const Vect<real_t>& u1,
                      Vect<real_t>& u2, 
                      size_t        nx,
                      size_t        ny,
                      size_t        nz,
                      size_t        dof)
{
   MeshToMesh(m1,m2,u1,u2,m1.getMinCoord(),m1.getMaxCoord(),nx,ny,nz);
}


void MeshToMesh(      Mesh&          m1,
                      Mesh&          m2,
                const Vect<real_t>&  u1,
                      Vect<real_t>&  u2, 
                const Point<real_t>& xmin,
                const Point<real_t>& xmax,
                      size_t         nx,
                      size_t         ny,
                      size_t         nz,
                      size_t         dof)
{
   Grid g;
   g.setDomain(xmin,xmax);
   g.setN(nx,ny,nz);
      //   Vect<real_t> ug((nx+1)*(ny+1)*(nz+1)*u1.getNbDOF());
   Vect<real_t> ug(nx+1,ny+1,nz+1);

   try {
      if (u1.getDOFType()==NODE_FIELD) {
         MeshToGrid(m1,g,u1,ug);
         GridToMesh(g,m2,ug,u2);
      }
      else if (u1.getDOFType()==ELEMENT_FIELD) {
         int nb_dof = u1.getNbDOF();
         Vect<real_t> nu1(m1,nb_dof,NODE_FIELD), nu2(m2,nb_dof,NODE_FIELD);
         Reconstruction r(m1);
         r.P0toP1(u1,nu1);
         MeshToGrid(m1,g,nu1,ug);
         GridToMesh(g,m2,ug,nu2);
         mesh_elements(m2) {
            real_t w = 0;
            for (size_t i=1; i<=The_element.getNbNodes(); i++)
               w += nu2(The_element(i)->n(),dof);
            u2(element_label,dof) = w/The_element.getNbNodes();
         }
      }
      else
         THROW_RT(" DOF data type not implemented.");
   }
   CATCH("MeshToMesh(...):");
}


int NodeInElement(const Node*    nd,
                  const Element* el)
{
   for (size_t i=1; i<=el->getNbNodes(); i++)
      if (el->getPtrNode(i)==nd)
         return int(i);
   return 0;
}


int NodeInSide(const Node* nd,
               const Side* sd)
{
   for (size_t i=1; i<=sd->getNbNodes(); i++)
      if (sd->getPtrNode(i)==nd)
         return int(i);
   return 0;
}


int SideInElement(const Side*    sd,
                  const Element* el)
{
   for (size_t i=1; i<=el->getNbSides(); i++)
      if (el->getPtrSide(i)==sd)
         return int(i);
   return 0;
}


real_t getMinElementMeasure(const Mesh& m)
{
   real_t a = 1./OFELI_EPSMCH;
   mesh_elements(m)
      a = std::min(a,The_element.getMeasure());
   return a;
}


real_t getMaxElementMeasure(const Mesh& m)
{
   real_t a = 0;
   mesh_elements(m)
      a = std::max(a,The_element.getMeasure());
   return a;
}


real_t getMinSideMeasure(const Mesh& m)
{
   real_t a = 1./OFELI_EPSMCH;
   mesh_sides(m)
      a = std::min(a,The_side.getMeasure());
   return a;
}


real_t getMaxSize(const Mesh& m)
{
   real_t a=0, b;
   mesh_elements(m) {
      if (The_element.getShape()==LINE)
         b = Line2(the_element).getLength();
      else if (The_element.getShape()==TRIANGLE)
         b = Triang3(the_element).getMaxEdgeLength();
      else if (The_element.getShape()==QUADRILATERAL)
         b = Quad4(the_element).getMaxEdgeLength();
      else if (The_element.getShape()==TETRAHEDRON)
         b = Tetra4(the_element).getMaxEdgeLength();
      else if (The_element.getShape()==HEXAHEDRON)
         b = Hexa8(the_element).getMaxEdgeLength();
      else if (The_element.getShape()==PENTAHEDRON)
         b = Penta6(the_element).getMaxEdgeLength();
      else
         ;
      a = std::max(a,b);
   }
   return a;
}


real_t getMinSize(const Mesh& m)
{
   real_t a=1./OFELI_EPSMCH, b;
   mesh_elements(m) {
      if (The_element.getShape()==LINE)
         b = Line2(the_element).getLength();
      else if (The_element.getShape()==TRIANGLE)
         b = Triang3(the_element).getMaxEdgeLength();
      else if (The_element.getShape()==QUADRILATERAL)
         b = Quad4(the_element).getMaxEdgeLength();
      else if (The_element.getShape()==TETRAHEDRON)
         b = Tetra4(the_element).getMaxEdgeLength();
      else if (The_element.getShape()==HEXAHEDRON)
         b = Hexa8(the_element).getMaxEdgeLength();
      else if (The_element.getShape()==PENTAHEDRON)
         b = Penta6(the_element).getMaxEdgeLength();
      else
         ;
      a = std::min(a,b);
   }
   return a;
}


real_t getMaxSideMeasure(const Mesh& m)
{
   real_t a = 0;
   mesh_sides(m)
      a = max(a,The_side.getMeasure());
   return a;
}


real_t getMeanElementMeasure(const Mesh& m)
{
   real_t a = 0;
   mesh_elements(m)
      a += The_element.getMeasure();
   return a/m.getNbElements();
}


real_t getMeanSideMeasure(const Mesh& m)
{
   real_t a = 0;
   mesh_sides(m)
      a += The_side.getMeasure();
   return a/m.getNbSides();
}


void setNodeCodes(     Mesh&   m,
                 const string& exp,
                       int     code,
                       size_t  dof)
{
   mesh_nodes(m)
      The_node.setCode(exp,code,dof);
}


void setBoundaryNodeCodes(      Mesh&   m,
                          const string& exp,
                                int     code,
                                size_t  dof)
{
   MeshBoundaryNodes(m)
      theNode->setCode(exp,code,dof);
}


void setSideCodes(      Mesh&   m,
                  const string& exp,
                        int     code,
                        size_t  dof)
{
   mesh_sides(m)
      The_side.setCode(exp,code,dof);
}


void setBoundarySideCodes(      Mesh&   m,
                          const string& exp,
                                int     code,
                                size_t  dof)
{
   MeshBoundarySides(m)
      theSide->setCode(exp,code,dof);
}


void setElementCodes(      Mesh&   m,
                     const string& exp,
                           int     code)
{
   mesh_elements(m)
      The_element.setCode(exp,code);
}

} /* namespace OFELI */
