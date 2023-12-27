/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2024 Rachid Touzani

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

                         Implementation of class 'Partition'

  ==============================================================================*/


#include "mesh/Partition.h"
#include "mesh/MeshUtil.h"
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>

using std::cout;
using std::ifstream;
using std::ostream;
using std::endl;

#include <iomanip>
using std::setw;

#include <algorithm>
#include <vector>
#include "util/util.h"

namespace OFELI {


Partition::Partition(Mesh&        mesh,
                     int          n,
                     vector<int>& epart)
          : _theMesh(&mesh), _nb_submesh(n)
{
   _theMesh->getAllSides();
   Init();
   NodeNeighborList();
}


Partition::Partition(Mesh&  mesh,
                     size_t n)
          : _theMesh(&mesh), _nb_submesh(n)
{
   Prepare();
   NodeNeighborList();
}


Partition::~Partition()
{
   if (_nb_submesh>0) {
      for (int i=0; i<_nb_submesh; i++) {
         delete _theSubMesh[i];
         delete [] _sm2m_node[i];
         delete [] _m2sm_node[i];
         delete [] _sm2m_element[i];
         delete [] _m2sm_element[i];
         delete _interface_side[i];
      }
      delete [] _sm2m_element;
      delete [] _m2sm_element;
      delete [] _sm2m_node;
      delete [] _m2sm_node;
   }
}


void Partition::set(Mesh&  mesh,
                    size_t n)
{
   _theMesh = &mesh;
   _nb_submesh = n;
   Prepare();
}


void Partition::Prepare()
{
   _theMesh->getAllSides();
   int ne=_theMesh->getNbElements(), nn=_theMesh->getNbNodes();
   vector<size_t> type(ne);

// Prepare data for Metis
   if (Verbosity)
      cout << "Preparing data for Metis ..." << endl;
   int i=0, esize=0;
   MESH_EL {
      type[i] = 0;
      int sh=the_element->getShape();
      if (sh==TRIANGLE)
         type[i] = 1;
      else if (sh==TETRAHEDRON)
         type[i] = 2;
      else if (sh==HEXAHEDRON)
         type[i] = 3;
      else if (sh==QUADRILATERAL)
         type[i] = 4;
      if (i==0) {
         switch (type[0]) {

            case 0:
               cout << "Unknown element type" << endl;
               exit(0);
               break;

            case 1:
               esize = 3;
               break;

            case 2:
               esize = 4;
               break;

            case 3:
               esize = 8;
               break;

            case 4:
               esize = 4;
               break;
         }
         if (int(the_element->getNbNodes())!=esize) {
            cout << "Error: Illegal number of element nodes" << endl;
            exit(0);
         }
      }

      else {
         if (type[i]!=type[0]) {
            cout << "Only one element shape is allowed." << endl;
            exit(0);
         }
      }
      i++;
   }

   vector<int> elmnts(esize*ne);
   MESH_EL {
      for (int j=0; j<esize; j++)
         elmnts[(element_label-1)*esize+j] = The_element(j+1)->n() - 1;
   }

// Mesh partitioning by Metis
   if (Verbosity)
      cout << "Mesh partitioning by Metis ..." << endl;
   _npart.resize(nn), _epart.resize(ne);
   int etype=type[0], numflag=0, edgecut;
   METIS_PartMeshDual(&ne,&nn,&elmnts[0],&etype,&numflag,&_nb_submesh,&edgecut,&_epart[0],&_npart[0]);
   Init();
}


void Partition::Init()
{
   size_t ne=_theMesh->getNbElements(), nn=_theMesh->getNbNodes();
   if (Verbosity)
      cout << "Creating submesh data ..." << endl;
   _sm2m_element = new size_t *[_nb_submesh];
   _m2sm_element = new size_t *[_nb_submesh];
   _sm2m_node = new size_t *[_nb_submesh];
   _m2sm_node = new size_t *[_nb_submesh];
   _nb_interface_sides.resize(_nb_submesh);
   for (int i=0; i<_nb_submesh; i++) {
      _nb_interface_sides[i] = 0;
      _m2sm_element[i] = new size_t [ne];
      _sm2m_element[i] = new size_t [ne];
      for (size_t j=0; j<ne; j++)
         _m2sm_element[i][j] = _sm2m_element[i][j] = 0;
      _m2sm_node[i] = new size_t [nn];
      _sm2m_node[i] = new size_t [nn];
      for (size_t j=0; j<nn; j++)
         _m2sm_node[i][j] = _sm2m_node[i][j] = 0;
   }

// Number locally nodes and elements in submeshes
   _theSubMesh.resize(_nb_submesh);
   for (int sd=0; sd<_nb_submesh; sd++) {
      size_t el_label=1, nd_label=1;
      Mesh ms;
      _theSubMesh[sd] = new Mesh;
      _theSubMesh[sd]->setDim(_theMesh->getDim());
      MESH_EL {
         size_t n=element_label;
         if (_epart[n-1]==sd) {
            Element *el=new Element(The_element);
            _m2sm_element[sd][n-1] = el_label;
            _sm2m_element[sd][el_label-1] = n;
            el->setLabel(el_label++);
            _theSubMesh[sd]->Add(el);
            for (size_t i=1; i<=the_element->getNbNodes(); i++) {
               Node *nd=new Node(*The_element(i));
               size_t m=The_element(i)->n();
               if (_m2sm_node[sd][m-1]==0) {
                  _m2sm_node[sd][m-1] = nd_label;
                  _sm2m_node[sd][nd_label-1] = m;
                  nd->setLabel(nd_label++);
                  _theSubMesh[sd]->Add(nd);
               }
               else
                  nd->setLabel(_m2sm_node[sd][m-1]);
               el->Replace(i,nd);
            }
         }
      }
      if (_theMesh->NodesAreDOF()==true)
         _theSubMesh[sd]->setNodesForDOF();
      else if (_theMesh->SidesAreDOF()==true)
         _theSubMesh[sd]->setSidesForDOF();
      else if (_theMesh->ElementsAreDOF()==true)
         _theSubMesh[sd]->setElementsForDOF();
      _theSubMesh[sd]->NumberEquations();
      _theSubMesh[sd]->getBoundarySides();
   }

// Create Interface Information
   if (Verbosity)
      cout << "Creating interface information ..." << endl;
   _interface_side = new Interface *[_nb_submesh];
   for (int i=0; i<_nb_submesh; i++)
      _interface_side[i] = new Interface [_theSubMesh[i]->getNbSides()];
   for (int sd1=0; sd1<_nb_submesh; sd1++) {
      for (int sd2=sd1+1; sd2<_nb_submesh; sd2++) {
         for (size_t i=1; i<=_theSubMesh[sd1]->getNbSides(); i++) {
            Side *ss1=_theSubMesh[sd1]->getPtrSide(i);
            for (size_t j=1; j<=_theSubMesh[sd2]->getNbSides(); j++) {
               Side *ss2=_theSubMesh[sd2]->getPtrSide(j);
               if (*ss1==*ss2) {
                  size_t k=_nb_interface_sides[sd1];
                  _interface_side[sd1][k].side1 = ss1->n();
                  _interface_side[sd1][k].sd = sd2;
                  _interface_side[sd1][k].side = ss2->n();
                  _nb_interface_sides[sd1]++;
                  k = _nb_interface_sides[sd2];
                  _interface_side[sd2][k].side1 = ss2->n();
                  _interface_side[sd2][k].sd = sd1;
                  _interface_side[sd2][k].side = ss1->n();
                  _nb_interface_sides[sd2]++;
               }
            }
         }
      }
   }
}


void Partition::NodeNeighborList()
{
   _node_neig.resize(_theMesh->getNbNodes());
   _nnz.resize(_theMesh->getNbNodes());
   MESH_ND {
      _nnz[node_label-1].resize(_nb_submesh+1);
      clear(_nnz[node_label-1]);
   }
   MESH_SD {
      for (size_t i=1; i<=The_side.getNbNodes(); i++) {
         the_node = The_side(i);
         size_t n=The_node.n()-1;
         for (size_t j=1; j<=The_side.getNbNodes(); j++) {
            if (i!=j) {
               size_t ns=The_side(j)->n();
               _node_neig[n].push_back(ns);
               _nnz[n][_npart[ns-1]]++;
            }
         }
      }
   }
   MESH_ND {
      size_t n=node_label-1;
      for (int s=0; s<_nb_submesh; s++)
         if (s!=_npart[n])
            _nnz[n][_nb_submesh] += _nnz[n][s];
   }
}


int Partition::getNbConnectInSubMesh(int n,
                                     int s) const
{
   return _nnz[n-1][s];
}


int Partition::getNbConnectOutSubMesh(int n,
                                      int s) const
{
   int m=0;
   for (int i=0; i<_nb_submesh; i++)
      if (i!=s)
         m += _nnz[n-1][i];
   return m;
}



ostream& operator<<(ostream&         s,
                    const Partition& p)
{
   s << "\nM E S H   P A R T I T I O N   I N F O R M A T I O N" << endl << endl;
   s << "Number of subdomains:   " << setw(6) << p._nb_submesh << endl << endl;
   for (int i=0; i<p._nb_submesh; i++) {
      s << endl << "Sub-domain:             " << setw(6) << i << endl;
      s << "Nb. of nodes:           " << setw(6) << p.getNbNodes(i) << endl;
      s << "Nb. of elements:        " << setw(6) << p.getNbElements(i) << endl;
      s << "Nb. of interface sides: " << setw(6) << p.getNbInterfaceSides(i) << endl;
      if (Verbosity > 1) {
        s << "Side label   Opposite submesh   Opposite Side" << endl;
        for (size_t j=1; j<=p._nb_interface_sides[i]; j++)
           s << setw(8)  << p._interface_side[i][j-1].side1
             << setw(14) << p._interface_side[i][j-1].sd
             << setw(18) << p._interface_side[i][j-1].side
             << endl;
      }
      s << endl;
   }
   return s;
}

}  /* namespace OFELI */
