/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2019 Rachid Touzani

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

                            Implementation of class DG

  ==============================================================================*/


#include "equations/DG.h"

namespace OFELI {

DG::DG(Mesh&  ms,
       size_t degree)
{
   _degree = degree;
   _nb_el_dof = _nb_sd_dof = 0;
   _theMesh = &ms;
   switch (_theMesh->getDim()) {

      case 1: _nb_el_dof = _degree + 1;
              _nb_sd_dof = 2;
              break;

      case 2: if (_degree<10)
                 _nb_el_dof = (_degree+1)*(_degree+2)/2;
              else
                 _nb_el_dof = (_degree+1)*(_degree+1);
              _nb_sd_dof = _degree + 1;
              break;

      case 3: if (_degree<10) {
                 _nb_el_dof = (_degree+1)*(_degree+2)/2;
                 _nb_sd_dof = (_degree+1)*(_degree+2)/2;
              }
              else if (_degree<20) {
                 _nb_el_dof = (_degree+1)*(_degree+1)*(_degree+1);
                 _nb_sd_dof = (_degree+1)*(_degree+1);
              }
              else {
                 _nb_el_dof = (_degree+1)*(_degree+1)*(_degree+2)/2;
                 _nb_sd_dof = (_degree+1)*(_degree+1);
              }
              break;
   }
   _theMesh->getAllSides();
   _nb_eq = _nb_el_dof*_theMesh->getNbElements();
   setDGLabel();
   setGraph();
   _b.setSize(_nb_eq);
}


DG::~DG()
{
}


int DG::setGraph()
{
   Vect<std::pair<size_t,size_t> > ij;
   MESH_EL {
      _ne = theElementLabel;
      for (size_t i=1; i<=_nb_el_dof; i++) {
         for (size_t j=1; j<=_nb_el_dof; j++)
            ij.push_back(RC(II(i),II(j)));
      }
      for (size_t k=1; k<=TheElement.getNbSides(); k++) {
         theSide = TheElement.getPtrSide(k);
         if (TheSide.isOnBoundary()==false) {
            Element *el = TheSide.getOtherNeighborElement(theElement);
            _nf = el->n();
            for (size_t i=1; i<=_nb_el_dof; i++)
               for (size_t j=1; j<=_nb_el_dof; j++)
                  ij.push_back(RC(II(i),IJ(j)));
         }
      }
   }
   _A.setGraph(ij,0);
   return ij.size();
}


void DG::setDGLabel()
{
   _g2l.resize(_theMesh->getNbElements());
   MESH_EL {
      for (size_t i=1; i<=TheElement.getNbNodes(); i++)
         _g2l[theElementLabel-1][TheElement(i)->n()-1] = i;
   }
}

} /* namespace OFELI */
